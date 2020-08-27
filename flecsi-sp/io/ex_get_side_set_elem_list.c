/*
 * Copyright (c) 2005 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.  
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
/*****************************************************************************
*
* exgssn - ex_get_side_set_node_list
*
* entry conditions - 
*   input parameters:
*       int     exoid                   exodus file id
*       int     side_set_id             side set id
*
* exit conditions - 
*       int     *side_set_node_cnt_list returned array of number of nodes for
*                                       side or face
*       int     *side_set_node_list     array of nodes
*
* revision history - 
*
*
*****************************************************************************/

#include <ctype.h>                      // for toupper
#include <inttypes.h>                   // for PRId64
#include <stddef.h>                     // for size_t
#include <stdio.h>                      // for sprintf
#include <stdlib.h>                     // for malloc, NULL, free
#include <string.h>                     // for strncmp, strlen
#include <sys/types.h>                  // for int64_t
#include "exodusII.h"                   // for ex_err, exerrval, ex_block, etc
#include "exodusII_int.h"               // for elem_blk_parm, EX_FATAL, etc

int ex_get_set_range (int   exoid,
		ex_entity_type set_type,
		ex_entity_id   set_id,
		void_int  *set_entry_list, 
		void_int  *set_extra_list, /* NULL if dont want to retrieve data */
    size_t begin,
    size_t end);

int ex_get_block_parm(int exoid,
    struct elem_blk_parm * elem_blk_parms)
{
  int ids_size;
  size_t m;
  size_t i;
  size_t elem_ctr;
  int64_t num_elem_blks, ndim;
  char errmsg[MAX_ERR_LENGTH];
  void_int *elem_blk_ids = NULL;
  int err_stat = EX_NOERR;
  
  ids_size = sizeof(int);
  if (ex_int64_status(exoid) & EX_IDS_INT64_API) {
    ids_size = sizeof(int64_t);
  }

  
  num_elem_blks = ex_inquire_int(exoid, EX_INQ_ELEM_BLK);
  if (num_elem_blks < 0) {
    sprintf(errmsg,
	    "Error: failed to get number of element blocks in file id %d",exoid);
    ex_err("ex_get_side_set_node_list_range",errmsg,exerrval);
    return(EX_FATAL);
  }
  
  /* Allocate space for the element block ids */
  if (!(elem_blk_ids= malloc(num_elem_blks*ids_size))) {
    exerrval = EX_MEMFAIL;
    sprintf(errmsg,
	    "Error: failed to allocate space for element block ids for file id %d",
	    exoid);
    ex_err("ex_get_side_set_node_list_range",errmsg,exerrval);
    err_stat = EX_FATAL;
    goto cleanup;
  }
  
  if (ex_get_elem_blk_ids(exoid, elem_blk_ids) == -1) {
    sprintf(errmsg,
	    "Error: failed to get element block ids in file id %d",
	    exoid);
    ex_err("ex_get_side_set_node_list_range",errmsg,EX_MSG);
    err_stat = EX_FATAL;
    goto cleanup;
  }

  elem_ctr = 0;
  for (i=0; i<num_elem_blks; i++) {
    ex_block block;

    if (ids_size == sizeof(int64_t)) {
      block.id = ((int64_t*)elem_blk_ids)[i];
    } else {
      block.id = ((int*)elem_blk_ids)[i];
    }
    block.type = EX_ELEM_BLOCK;

    /* read in an element block parameter */
    if ((ex_get_block_param (exoid, &block)) == -1) {
      sprintf(errmsg,
	      "Error: failed to get element block %"PRId64" parameters in file id %d",
              block.id, exoid);
      ex_err("ex_get_side_set_node_list_range",errmsg,EX_MSG);
      err_stat = EX_FATAL;
      goto cleanup;
    }

    elem_blk_parms[i].num_elem_in_blk = block.num_entry;
    elem_blk_parms[i].num_nodes_per_elem = block.num_nodes_per_entry;
    elem_blk_parms[i].num_attr = block.num_attribute;
    elem_blk_parms[i].elem_blk_id = block.id;

    for (m=0; m < strlen(block.topology); m++) {
      elem_blk_parms[i].elem_type[m] = toupper(block.topology[m]);
    }
    elem_blk_parms[i].elem_type[m] = '\0';

    if (strncmp(elem_blk_parms[i].elem_type,"CIRCLE",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_CIRCLE;
	/* set side set node stride */
        elem_blk_parms[i].num_nodes_per_side[0] = 1;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"SPHERE",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_SPHERE;
	/* set side set node stride */
        elem_blk_parms[i].num_nodes_per_side[0] = 1;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"QUAD",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_QUAD;
	/* determine side set node stride */
	if (elem_blk_parms[i].num_nodes_per_elem == 4)
	  elem_blk_parms[i].num_nodes_per_side[0] = 2;
	else if (elem_blk_parms[i].num_nodes_per_elem == 5)
	  elem_blk_parms[i].num_nodes_per_side[0] = 2;
	else
	  elem_blk_parms[i].num_nodes_per_side[0] = 3;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"TRIANGLE",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_TRIANGLE;
	/* set default side set node stride */
	if (ndim == 2)  /* 2d TRIs */
	  {
	    if (elem_blk_parms[i].num_nodes_per_elem == 3)
	      elem_blk_parms[i].num_nodes_per_side[0] = 2;
	    else 
	      elem_blk_parms[i].num_nodes_per_side[0] = 3;
	  }
	else if (ndim == 3)  /* 3d TRIs */
	  {
	    if (elem_blk_parms[i].num_nodes_per_elem == 3)
	      elem_blk_parms[i].num_nodes_per_side[0] = 3;
	    else 
	      elem_blk_parms[i].num_nodes_per_side[0] = 6;
	  }
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"SHELL",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_SHELL;
	/* determine side set node stride */
	if (elem_blk_parms[i].num_nodes_per_elem == 2) /* KLUDGE for 2D Shells*/
	  elem_blk_parms[i].num_nodes_per_side[0] = 2;
	else if (elem_blk_parms[i].num_nodes_per_elem == 4)
	  elem_blk_parms[i].num_nodes_per_side[0] = 4;
	else
	  elem_blk_parms[i].num_nodes_per_side[0] = 8;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"HEX",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_HEX;
	/* determine side set node stride */
	if (elem_blk_parms[i].num_nodes_per_elem == 8)  /* 8-node bricks */
	  elem_blk_parms[i].num_nodes_per_side[0] = 4;
	else if (elem_blk_parms[i].num_nodes_per_elem == 9)  /* 9-node bricks */
	  elem_blk_parms[i].num_nodes_per_side[0] = 4;
	else if (elem_blk_parms[i].num_nodes_per_elem == 12)  /* HEXSHELLS */
	  elem_blk_parms[i].num_nodes_per_side[0] = 4;
	else if (elem_blk_parms[i].num_nodes_per_elem == 27)  /* 27-node bricks */
	  elem_blk_parms[i].num_nodes_per_side[0] = 9;
	else 
	  elem_blk_parms[i].num_nodes_per_side[0] = 8;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"TETRA",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_TETRA;
	/* determine side set node stride */
	if (elem_blk_parms[i].num_nodes_per_elem == 4)
	  elem_blk_parms[i].num_nodes_per_side[0] = 3;
	else if (elem_blk_parms[i].num_nodes_per_elem == 8)
	  elem_blk_parms[i].num_nodes_per_side[0] = 4;
	else
	  elem_blk_parms[i].num_nodes_per_side[0] = 6;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"WEDGE",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_WEDGE;
	/* determine side set node stride */
	if (elem_blk_parms[i].num_nodes_per_elem == 6)
	  elem_blk_parms[i].num_nodes_per_side[0] = 4;
	else
	  elem_blk_parms[i].num_nodes_per_side[0] = 8;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"PYRAMID",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_PYRAMID;
	/* determine side set node stride */
	if (elem_blk_parms[i].num_nodes_per_elem == 5)
	  elem_blk_parms[i].num_nodes_per_side[0] = 4;
	else
	  elem_blk_parms[i].num_nodes_per_side[0] = 8;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"BEAM",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_BEAM;
	/* determine side set node stride */
	if (elem_blk_parms[i].num_nodes_per_elem == 2)
	  elem_blk_parms[i].num_nodes_per_side[0] = 2;
	else 
	  elem_blk_parms[i].num_nodes_per_side[0] = 3;
      }
    else if ( (strncmp(elem_blk_parms[i].elem_type,"TRUSS",3) == 0) ||
              (strncmp(elem_blk_parms[i].elem_type,"BAR",3) == 0) ||
              (strncmp(elem_blk_parms[i].elem_type,"EDGE",3) == 0) )
      {
	elem_blk_parms[i].elem_type_val = EX_EL_TRUSS;
	/* determine side set node stride */
	if (elem_blk_parms[i].num_nodes_per_elem == 2)
	  elem_blk_parms[i].num_nodes_per_side[0] = 2;
	else 
	  elem_blk_parms[i].num_nodes_per_side[0] = 3;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"NULL",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_NULL_ELEMENT;
	elem_blk_parms[i].num_nodes_per_side[0] = 0;
	elem_blk_parms[i].num_elem_in_blk = 0;
      }
    else
      { /* unsupported element type; no problem if no sides specified for
	   this element block */
	elem_blk_parms[i].elem_type_val = EX_EL_UNK;
	elem_blk_parms[i].num_nodes_per_side[0] = 0;
      }
    elem_blk_parms[i].elem_blk_id = block.id;    /* save id */
    elem_ctr += elem_blk_parms[i].num_elem_in_blk;
    elem_blk_parms[i].elem_ctr = elem_ctr;      /* save elem number max */
  }

 cleanup:
  ex_safe_free(elem_blk_ids);
}

int ex_get_side_set_elem_list_range(int exoid,
			      ex_entity_id side_set_id,
            const struct elem_blk_parm * elem_blk_parms,
			      void_int *side_set_elem_list,
			      void_int *side_set_side_list,
            void_int *ss_elem_ndx,
            void_int* ss_elem_node_ndx,
            void_int* ss_parm_ndx,
            size_t begin,
            size_t end)
{
  size_t m;
  size_t i, j; 
  int64_t elem, side;
  int64_t num_side_sets, num_elem_blks, num_df, ndim;
  int64_t tot_num_elem = 0, tot_num_ss_elem = 0, elem_num = 0;
  size_t connect_offset, side_num, node_pos;
  void_int *elem_blk_ids = NULL;
  void_int *connect = NULL; 
  //void_int *ss_elem_ndx = NULL;
  //void_int *ss_elem_node_ndx = NULL;
  //void_int *ss_parm_ndx = NULL;
  //void_int *side_set_elem_list = NULL;
  //void_int *side_set_side_list = NULL;
  size_t elem_ctr, node_ctr, elem_num_pos;
  size_t num_nodes_per_elem;
  int int_size, ids_size;

  int err_stat = EX_NOERR;
  int status;

  //struct elem_blk_parm  *elem_blk_parms = NULL;

  size_t num_ss_elem = end - begin;

  char errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  /* first check if any side sets are specified */
  /* inquire how many side sets have been stored */

  num_side_sets = ex_inquire_int(exoid, EX_INQ_SIDE_SETS);
  if (num_side_sets < 0) {
    sprintf(errmsg,
	    "Error: failed to get number of side sets in file id %d",exoid);
    ex_err("ex_get_side_set_node_list_range",errmsg,exerrval);
    return(EX_FATAL);
  }

  if (num_side_sets == 0) {
    sprintf(errmsg,
	    "Warning: no side sets defined in file id %d",exoid);
    ex_err("ex_get_side_set_node_list_range",errmsg,EX_WARN);
    return(EX_WARN);
  }

  /* Lookup index of side set id in VAR_SS_IDS array */
  ex_id_lkup(exoid,EX_SIDE_SET,side_set_id);
  if (exerrval != 0)  {
    if (exerrval == EX_NULLENTITY) {
      sprintf(errmsg,
              "Warning: side set %"PRId64" is NULL in file id %d",
	      side_set_id,exoid);
      ex_err("ex_get_side_set_node_list_range",errmsg,EX_NULLENTITY);
      return (EX_WARN);
    }
    else {

      sprintf(errmsg,
	      "Error: failed to locate side set %"PRId64" in VAR_SS_IDS array in file id %d",
	      side_set_id,exoid);
      ex_err("ex_get_side_set_node_list_range",errmsg,exerrval);
      return (EX_FATAL);
    }
  }

  num_elem_blks = ex_inquire_int(exoid, EX_INQ_ELEM_BLK);
  if (num_elem_blks < 0) {
    sprintf(errmsg,
	    "Error: failed to get number of element blocks in file id %d",exoid);
    ex_err("ex_get_side_set_node_list_range",errmsg,exerrval);
    return(EX_FATAL);
  }

  tot_num_elem = ex_inquire_int(exoid, EX_INQ_ELEM);
  if (tot_num_elem < 0) {
    sprintf(errmsg,
	    "Error: failed to get total number of elements in file id %d",exoid);
    ex_err("ex_get_side_set_node_list_range",errmsg,exerrval);
    return(EX_FATAL);
  }

  /* get the dimensionality of the coordinates;  this is necessary to
     distinguish between 2d TRIs and 3d TRIs */
  ndim = ex_inquire_int(exoid, EX_INQ_DIM);
  if (ndim < 0) {
    sprintf(errmsg,
	    "Error: failed to get dimensionality in file id %d",exoid);
    ex_err("ex_cvt_nodes_to_sides",errmsg,exerrval);
    return(EX_FATAL);
  }

  int_size = sizeof(int);
  if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
    int_size = sizeof(int64_t);
  }

  ids_size = sizeof(int);
  if (ex_int64_status(exoid) & EX_IDS_INT64_API) {
    ids_size = sizeof(int64_t);
  }

  /* First determine the  # of elements in the side set*/
  if (int_size == sizeof(int64_t)) {
    status = ex_get_set_param(exoid,EX_SIDE_SET, side_set_id,&tot_num_ss_elem,&num_df);
  } else {
    int tot, df;
    status = ex_get_set_param(exoid,EX_SIDE_SET, side_set_id,&tot,&df);
    tot_num_ss_elem = tot;
    num_df = df;
  }
  
  if (end>tot_num_ss_elem || begin<0)
  {
      exerrval = 1;    
      sprintf(errmsg,
              "Error: request is outside acceptable domain for file id %d",
              exoid);
      ex_err("ex_get_side_set_node_list_range",errmsg,exerrval);
      return (EX_FATAL);
   }

  if (status != EX_NOERR) {
    sprintf(errmsg,
	    "Error: failed to get number of elements in side set %"PRId64" in file id %d",
            side_set_id, exoid);
    ex_err("ex_get_side_set_node_list_range",errmsg,exerrval);
    return(EX_FATAL);
  }

#if 0
  /* Allocate space for the side set element list */
  if (!(side_set_elem_list=malloc(num_ss_elem*int_size))) {
    exerrval = EX_MEMFAIL;
    sprintf(errmsg,
	    "Error: failed to allocate space for side set element list for file id %d",
	    exoid);
    ex_err("ex_get_side_set_node_list_range",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* Allocate space for the side set side list */
  if (!(side_set_side_list=malloc(num_ss_elem*int_size))) {
    exerrval = EX_MEMFAIL;
    sprintf(errmsg,
	    "Error: failed to allocate space for side set side list for file id %d",
	    exoid);
    ex_err("ex_get_side_set_node_list_range",errmsg,exerrval);
    err_stat = EX_FATAL;
    goto cleanup;
  }
#endif

  if (ex_get_set_range(exoid, EX_SIDE_SET, side_set_id, 
		 side_set_elem_list, side_set_side_list, begin, end) == -1) {
    sprintf(errmsg,
	    "Error: failed to get side set %"PRId64" in file id %d",
	    side_set_id, exoid);
    ex_err("ex_get_side_set_node_list_range",errmsg,exerrval);
    err_stat = EX_FATAL;
    goto cleanup;
  }

#if 0
  /* Allocate space for the ss element index array */
  if (!(ss_elem_ndx= malloc(num_ss_elem*int_size))) {
    exerrval = EX_MEMFAIL;
    sprintf(errmsg,
	    "Error: failed to allocate space for side set elem sort array for file id %d",
	    exoid);
    ex_err("ex_get_side_set_node_list_range",errmsg,exerrval);
    err_stat = EX_FATAL;
    goto cleanup;
  }
#endif

  /* Sort side set element list into index array  - non-destructive */
  if (int_size == sizeof(int64_t)) {
    /* Sort side set element list into index array  - non-destructive */
    int64_t *elems = (int64_t*)ss_elem_ndx;
    for (i=0;i<num_ss_elem;i++) {
      elems[i] = i; /* init index array to current position */
    }
    ex_iqsort64(side_set_elem_list, ss_elem_ndx,num_ss_elem);
  } else {
    /* Sort side set element list into index array  - non-destructive */
    int *elems = (int*)ss_elem_ndx;
    for (i=0;i<num_ss_elem;i++) {
      elems[i] = i; /* init index array to current position */
    }
    ex_iqsort(side_set_elem_list, ss_elem_ndx,num_ss_elem);
  }

#if 0
  /* Allocate space for the element block ids */
  if (!(elem_blk_ids= malloc(num_elem_blks*ids_size))) {
    exerrval = EX_MEMFAIL;
    sprintf(errmsg,
	    "Error: failed to allocate space for element block ids for file id %d",
	    exoid);
    ex_err("ex_get_side_set_node_list_range",errmsg,exerrval);
    err_stat = EX_FATAL;
    goto cleanup;
  }
  
  if (ex_get_elem_blk_ids(exoid, elem_blk_ids) == -1) {
    sprintf(errmsg,
	    "Error: failed to get element block ids in file id %d",
	    exoid);
    ex_err("ex_get_side_set_node_list_range",errmsg,EX_MSG);
    err_stat = EX_FATAL;
    goto cleanup;
  }

  /* Allocate space for the element block params */
  if (!(elem_blk_parms= malloc(num_elem_blks*sizeof(struct elem_blk_parm)))) {
    exerrval = EX_MEMFAIL;
    sprintf(errmsg,
	    "Error: failed to allocate space for element block params for file id %d",
            exoid);
    ex_err("ex_get_side_set_node_list_range",errmsg,exerrval);
    err_stat = EX_FATAL;
    goto cleanup;
  }

  elem_ctr = 0;
  for (i=0; i<num_elem_blks; i++) {
    ex_block block;

    if (ids_size == sizeof(int64_t)) {
      block.id = ((int64_t*)elem_blk_ids)[i];
    } else {
      block.id = ((int*)elem_blk_ids)[i];
    }
    block.type = EX_ELEM_BLOCK;

    /* read in an element block parameter */
    if ((ex_get_block_param (exoid, &block)) == -1) {
      sprintf(errmsg,
	      "Error: failed to get element block %"PRId64" parameters in file id %d",
              block.id, exoid);
      ex_err("ex_get_side_set_node_list_range",errmsg,EX_MSG);
      err_stat = EX_FATAL;
      goto cleanup;
    }

    elem_blk_parms[i].num_elem_in_blk = block.num_entry;
    elem_blk_parms[i].num_nodes_per_elem = block.num_nodes_per_entry;
    elem_blk_parms[i].num_attr = block.num_attribute;
    elem_blk_parms[i].elem_blk_id = block.id;

    for (m=0; m < strlen(block.topology); m++) {
      elem_blk_parms[i].elem_type[m] = toupper(block.topology[m]);
    }
    elem_blk_parms[i].elem_type[m] = '\0';

    if (strncmp(elem_blk_parms[i].elem_type,"CIRCLE",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_CIRCLE;
	/* set side set node stride */
        elem_blk_parms[i].num_nodes_per_side[0] = 1;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"SPHERE",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_SPHERE;
	/* set side set node stride */
        elem_blk_parms[i].num_nodes_per_side[0] = 1;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"QUAD",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_QUAD;
	/* determine side set node stride */
	if (elem_blk_parms[i].num_nodes_per_elem == 4)
	  elem_blk_parms[i].num_nodes_per_side[0] = 2;
	else if (elem_blk_parms[i].num_nodes_per_elem == 5)
	  elem_blk_parms[i].num_nodes_per_side[0] = 2;
	else
	  elem_blk_parms[i].num_nodes_per_side[0] = 3;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"TRIANGLE",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_TRIANGLE;
	/* set default side set node stride */
	if (ndim == 2)  /* 2d TRIs */
	  {
	    if (elem_blk_parms[i].num_nodes_per_elem == 3)
	      elem_blk_parms[i].num_nodes_per_side[0] = 2;
	    else 
	      elem_blk_parms[i].num_nodes_per_side[0] = 3;
	  }
	else if (ndim == 3)  /* 3d TRIs */
	  {
	    if (elem_blk_parms[i].num_nodes_per_elem == 3)
	      elem_blk_parms[i].num_nodes_per_side[0] = 3;
	    else 
	      elem_blk_parms[i].num_nodes_per_side[0] = 6;
	  }
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"SHELL",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_SHELL;
	/* determine side set node stride */
	if (elem_blk_parms[i].num_nodes_per_elem == 2) /* KLUDGE for 2D Shells*/
	  elem_blk_parms[i].num_nodes_per_side[0] = 2;
	else if (elem_blk_parms[i].num_nodes_per_elem == 4)
	  elem_blk_parms[i].num_nodes_per_side[0] = 4;
	else
	  elem_blk_parms[i].num_nodes_per_side[0] = 8;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"HEX",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_HEX;
	/* determine side set node stride */
	if (elem_blk_parms[i].num_nodes_per_elem == 8)  /* 8-node bricks */
	  elem_blk_parms[i].num_nodes_per_side[0] = 4;
	else if (elem_blk_parms[i].num_nodes_per_elem == 9)  /* 9-node bricks */
	  elem_blk_parms[i].num_nodes_per_side[0] = 4;
	else if (elem_blk_parms[i].num_nodes_per_elem == 12)  /* HEXSHELLS */
	  elem_blk_parms[i].num_nodes_per_side[0] = 4;
	else if (elem_blk_parms[i].num_nodes_per_elem == 27)  /* 27-node bricks */
	  elem_blk_parms[i].num_nodes_per_side[0] = 9;
	else 
	  elem_blk_parms[i].num_nodes_per_side[0] = 8;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"TETRA",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_TETRA;
	/* determine side set node stride */
	if (elem_blk_parms[i].num_nodes_per_elem == 4)
	  elem_blk_parms[i].num_nodes_per_side[0] = 3;
	else if (elem_blk_parms[i].num_nodes_per_elem == 8)
	  elem_blk_parms[i].num_nodes_per_side[0] = 4;
	else
	  elem_blk_parms[i].num_nodes_per_side[0] = 6;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"WEDGE",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_WEDGE;
	/* determine side set node stride */
	if (elem_blk_parms[i].num_nodes_per_elem == 6)
	  elem_blk_parms[i].num_nodes_per_side[0] = 4;
	else
	  elem_blk_parms[i].num_nodes_per_side[0] = 8;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"PYRAMID",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_PYRAMID;
	/* determine side set node stride */
	if (elem_blk_parms[i].num_nodes_per_elem == 5)
	  elem_blk_parms[i].num_nodes_per_side[0] = 4;
	else
	  elem_blk_parms[i].num_nodes_per_side[0] = 8;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"BEAM",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_BEAM;
	/* determine side set node stride */
	if (elem_blk_parms[i].num_nodes_per_elem == 2)
	  elem_blk_parms[i].num_nodes_per_side[0] = 2;
	else 
	  elem_blk_parms[i].num_nodes_per_side[0] = 3;
      }
    else if ( (strncmp(elem_blk_parms[i].elem_type,"TRUSS",3) == 0) ||
              (strncmp(elem_blk_parms[i].elem_type,"BAR",3) == 0) ||
              (strncmp(elem_blk_parms[i].elem_type,"EDGE",3) == 0) )
      {
	elem_blk_parms[i].elem_type_val = EX_EL_TRUSS;
	/* determine side set node stride */
	if (elem_blk_parms[i].num_nodes_per_elem == 2)
	  elem_blk_parms[i].num_nodes_per_side[0] = 2;
	else 
	  elem_blk_parms[i].num_nodes_per_side[0] = 3;
      }
    else if (strncmp(elem_blk_parms[i].elem_type,"NULL",3) == 0)
      {
	elem_blk_parms[i].elem_type_val = EX_EL_NULL_ELEMENT;
	elem_blk_parms[i].num_nodes_per_side[0] = 0;
	elem_blk_parms[i].num_elem_in_blk = 0;
      }
    else
      { /* unsupported element type; no problem if no sides specified for
	   this element block */
	elem_blk_parms[i].elem_type_val = EX_EL_UNK;
	elem_blk_parms[i].num_nodes_per_side[0] = 0;
      }
    elem_blk_parms[i].elem_blk_id = block.id;    /* save id */
    elem_ctr += elem_blk_parms[i].num_elem_in_blk;
    elem_blk_parms[i].elem_ctr = elem_ctr;      /* save elem number max */
  }


  /* Allocate space for the ss element to element block parameter index array */
  if (!(ss_parm_ndx=malloc(num_ss_elem*int_size)))
    {
      exerrval = EX_MEMFAIL;
      sprintf(errmsg,
	      "Error: failed to allocate space for side set elem parms index for file id %d",
	      exoid);
      ex_err("ex_get_side_set_node_list_range",errmsg,exerrval);
      err_stat = EX_FATAL;
      goto cleanup;
    }


  /* Allocate space for the ss element to node list index array */
  if (!(ss_elem_node_ndx=malloc(num_ss_elem*int_size)))
    {
      exerrval = EX_MEMFAIL;
      sprintf(errmsg,
	      "Error: failed to allocate space for side set elem to node index for file id %d",
	      exoid);
      ex_err("ex_get_side_set_node_list_range",errmsg,exerrval);
           
      err_stat = EX_FATAL;
      goto cleanup;
    }
#endif

  /* Build side set element to node list index and side set element 
     parameter index.
  */
  for (i=0;i<num_ss_elem;i++) {
    node_ctr = 0;

    if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
      elem = ((int64_t*)side_set_elem_list)[i];
      side = ((int64_t*)side_set_side_list)[i];
    } else {
      elem = ((int*)side_set_elem_list)[i];
      side = ((int*)side_set_side_list)[i];
    }

    for (j=0; j<num_elem_blks; j++) {
      if (elem_blk_parms[j].elem_type_val != EX_EL_NULL_ELEMENT)
	if (elem <= elem_blk_parms[j].elem_ctr)
	  break;
    }

    if (j >= num_elem_blks) {
      exerrval = EX_BADPARAM;
      sprintf(errmsg,
	      "Error: Invalid element number %"PRId64" found in side set %"PRId64" in file %d",
              elem, side_set_id, exoid);
      ex_err("ex_get_side_set_node_list_range",errmsg,EX_MSG);
      err_stat = EX_FATAL;
      goto cleanup;
    }

    /* Update node_ctr (which points to next node in chain */

    /* WEDGEs with 3 node sides (side 4 or 5) are special cases */
    if (elem_blk_parms[j].elem_type_val == EX_EL_WEDGE &&
        (side == 4 || side == 5))
      {
	if (elem_blk_parms[j].num_nodes_per_elem == 6)
	  node_ctr += 3;  /* 3 node side */
	else
	  node_ctr += 6;  /* 6 node side */
      }
    /* PYRAMIDSs with 3 node sides (sides 1,2,3,4) are also special */
    else if (elem_blk_parms[j].elem_type_val == EX_EL_PYRAMID &&
             (side < 5))
      {
	if (elem_blk_parms[j].num_nodes_per_elem == 5)
	  node_ctr += 3;  /* 3 node side */
	else
	  node_ctr += 6;  /* 6 node side */
      }
    /* side numbers 3,4,5,6 for SHELLs are also special */
    else if (elem_blk_parms[j].elem_type_val == EX_EL_SHELL &&
	     (side > 2 ))
      {
	if (elem_blk_parms[j].num_nodes_per_elem == 4)
	  node_ctr += 2;  /* 2 node side */
	else
	  node_ctr += 3;  /* 3 node side */
      }
    /* side numbers 3,4,5 for 3d TRIs are also special */
    else if (elem_blk_parms[j].elem_type_val == EX_EL_TRIANGLE &&
             ndim == 3 &&
             side > 2 )
      {
	if (elem_blk_parms[j].num_nodes_per_elem == 3)  /* 3-node TRI */
	  node_ctr += 2;  /* 2 node side */
	else   /* 6-node TRI */
	  node_ctr += 3;  /* 3 node side */
      }
    else /* all other element types */
      node_ctr += elem_blk_parms[j].num_nodes_per_side[0];
    
    if (int_size == sizeof(int64_t)) {
      ((int64_t*)ss_parm_ndx)[i] = j; /* assign parameter block index */
      ((int64_t*)ss_elem_node_ndx)[i] = node_ctr;     /* assign node list index */
    } else {
      ((int*)ss_parm_ndx)[i] = j; /* assign parameter block index */
      ((int*)ss_elem_node_ndx)[i] = node_ctr;     /* assign node list index */
    }

  }


  /* All done: release connectivity array space, element block ids array,
     element block parameters array, and side set element index array */
 cleanup:
  //ex_safe_free(connect);
  //ex_safe_free(ss_parm_ndx);
  //ex_safe_free(elem_blk_ids);
  //ex_safe_free(elem_blk_parms);
  //ex_safe_free(ss_elem_ndx);
  //ex_safe_free(ss_elem_node_ndx);
  //ex_safe_free(side_set_side_list);
  //ex_safe_free(side_set_elem_list);

  return(err_stat);
}
