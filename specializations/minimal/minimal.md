<!-- CINCHDOC DOCUMENT(Developer Guide) SECTION(Minimal Mesh) -->

# Minimal Mesh

This example defines the minimum set of types to implement a FleCSI mesh
specialization. This section describes the steps taken in creating it.
**This description is not a substitute for actual review of the files
that define the specialization, and the reader is encouraged to
carefully review the actual source code before embarking on their own
development.**

## Files

The following files are included for this specialization:

* **minimal_mesh.h** This file defines the mesh type *minimal_mesh_t*
  with a simple mesh interface, including a constructor, mesh
  initialization, and methods to add vertices and cells. Additionally,
  this type defines the required *indeces* method that associates a
  generic integer value to a particular index space. This type derives
  from topology::mesh_topology_t.

* **minimal_entity_types.h** This file defines the mesh entity types
  *vertex_t*, *edge_t*, *face_t*, and *cell_t*. Actual specializations
  would define interfaces for these types that reflect the needs of the
  sovlers for which they are designed.

* **minimal_types.h** This file defines the mesh type *MT*
  parameterization used by topology::mesh_topology_t. The requirements
  for this type are defined below.

* **minimal_config.h** This file defines the mesh properties
  *num_dimensions*, and *num_domains*, and some basic types.

## Mesh Adjacency Information

The FleCSI topology type defines two methods for associating mesh
entities to create adjacency information for the user: *connectivities*
and *bindings*. In order to understand the distinction between these
methods, we must first consider how FleCSI thinks about entities.

### Mesh Domains

A FleCSI mesh domain is allowed to have one type associated with each
topological dimension. Mesh specializations that have a single mesh
domain can be termed *simple*, and provide adjacency information between
distinct dimensional types, e.g., vertex (topological dimension 0) ->
cell (topological dimension 2) for a two-dimensional mesh. Adjacencies
between entities within a single mesh domain are defined in the
*connectivies* type discussed below.

A mesh specialization that requires more than one mesh domain is termed
*complex*, and provides adjacency information between distinct
dimensional types *and* between types in different mesh domains. The need
for mesh domains arises when a mesh specialization exposes more than one
type for a given toplogical dimension, e.g., edges and corners. Edges
have toplogical dimension 1 because they are defined by two vertices of
topological dimension 0. However, corners also have topological
dimension 1 because they are defined by a vertex and a cell center.
Defining two entity types with the same topological dimension would
break the underlying logic used by the mesh topology to track adjacency
information. To fix this problem, FleCSI defines mesh domains. Adjacency
information between entities in different mesh domains is defined in the
*bindings* type discussed below.

## Mesh Type

The FleCSI mesh topology is parameterized on a mesh type *MT* that defines
several subtypes:

* **num_dimensions** This must be a static constexpr size_t that
  defines the number of topological dimensions of the mesh
  specialization.

      static constexpr size_t num_dimensions = 3;

* **num_domains** This must be a static constexpr size_t that defines
  the number of mesh domains included in the specialization.

      static constexpr size_t num_domains = 1;

* **entity_types** This must be a std::tuple of std::pair, each of which
  associate a mesh domain with a user-defined entity type for a
  particular topological dimension. As an example, if the user defines
  entity types like:

      struct vertex_t
        : public topology::mesh_entity_t<0, num_domains>
      {
      }; // struct vertex_t

      struct edge_t
        : public topology::mesh_entity_t<1, num_domains>
      {
      }; // struct edge_t

      struct cell_t
        : public topology::mesh_entity_t<2, num_domains>
      {
      }; // struct cell_t

      struct corner_t
        : public topology::mesh_entity_t<1, num_domains>
      {
      }; // struct corner_t

  These can be used to construct the *entity_types* tuple like:

      using entity_types =
        std::tuple<
          std::pair<domain_<0>, vertex_t>,
          std::pair<domain_<0>, edge_t>,
          std::pair<domain_<0>, cell_t>,
          std::pair<domain_<1>, corner_t>
        >;

  The *domain_<N>* construct is simply a way to typeify an integer
  value, so *domain_<0>* simply identifies domain 0 to the compiler.

* **connectivities** This must be a std::tuple of std::tuple, each of
  which defines an adjacency relationship that should be computed and
  stored by the unerlying mesh topology. As an example, using the entity
  types defined above, the user can define connectivities like:

      using connectivities =
        std::tuple<
          std::tuple<domain_<0>, vertex_t, cell_t>,
          std::tuple<domain_<0>, cell_t, vertex_t>
        >;

  This tells the topology that it should compute and store adjacency
  information from a vertex to a cell, and from a cell to a vertex. This
  information is available statically (at compile time), meaning that
  the mesh topology data structure only needs to create storage for
  these specific adjacency tables. It is still possible for a mesh
  specialization to provide other adjacency information by following the
  information stored in the mesh. This makes it possible to employ
  different strategies for computed vs. stored information so that a
  mesh specialization can be tuned to make the best choices for the
  tradeoffs associated with memory usage or compute time.

* **bindings** This must be a std::tuple of std::tuple, each of which
  defines an adjacency relationship between entities in different mesh
  domains that should be computed and stored by the underlying mesh
  topology. As an example, using the entity types defined above, the
  user can define bindings like:

      using bindings =
        std::tuple<
          std::tuple<domain_<1>, domain_<0>, corner_t, cell_t>,
          std::tuple<domain_<0>, domain_<1>, cell_t, corner_t>
        >;

  Note that bindings can be empty for *simple* mesh specializations:

      using bindings = std::tuple<>;

* **create_entity** An implementation of this method must be provided.
  The method creates the appropriate entity type given the domain *M*,
  dimension *D*, and the number of vertices *num_vertices*. The method
  has the form:

      template<
        size_t M,
        size_t D
      >
      static
      topology::mesh_entity_base_t<num_domains> *
      create_entity(
        topology::mesh_topology_base_t * mesh,
        size_t num_vertices
      )
      {
      } // create_entity

  **Explanation:** Initialization of the specialized mesh involves
  direct creation of entities for topological dimension 0 (vertices) and
  dimension N, where N is the highest dimension, e.g., dimension 3
  (cells) for a three-dimensional mesh. Entities of interior topological
  dimension, e.g., of dimension 1 and dimension 2 in this example, are
  created indirectly by the mesh topology itself. To fascilitate this,
  the specialization provides the create_entity method so that the
  correct entity type is created for a given domain, dimension and
  number of vertices. This allows the selection of different types for
  complex, heterogeneous meshes.

<!-- vim: set tabstop=2 shiftwidth=2 expandtab fo=cqt tw=72 : -->
