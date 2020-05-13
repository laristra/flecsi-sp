# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack import *


class FlecsiSpDeps(BundlePackage):
    '''Flecsi-SP contains various specializations for use with the FleCSI core programming system.
    '''
    git = 'https://github.com/laristra/flecsi-sp.git'

    version('1.4', branch='1.4', submodules=False, preferred=True)

    variant('build_type', default='Release', values=('Debug', 'Release', 'RelWithDebInfo', 'MinSizeRel'),
            description='The build type to build', multi=False)
    variant('backend', default='mpi', values=('serial', 'mpi', 'legion', 'charmpp', 'hpx'),
            description='Backend to use for distributed memory', multi=False)
    variant('cinch', default=True,
            description='Enable External Cinch')
    variant('shared', default=True,
            description='Build shared libraries')
    variant('hdf5', default=True,
            description='Enable HDF5 Support')
    variant('caliper', default=False,
            description='Enable Caliper Support')
    variant('graphviz', default=False,
            description='Enable GraphViz Support')
    variant('tutorial', default=False,
            description='Build FleCSI Tutorials')
    variant('portage', default=False,
            description='Enable Portage Support')

    for b in ['serial', 'mpi', 'legion', 'charm++', 'hpx', 'trilinos']:
        depends_on("flecsi-deps@1.4 backend=%s" % b,
            when="backend=%s" % b)
    for b in ['Debug', 'Release', 'RelWithDebInfo', 'MinSizeRel']:
        depends_on("flecsi-deps@1.4 build_type=%s" % b,
            when="build_type=%s" % b)
    for v in ['shared', 'hdf5', 'caliper', 'graphviz', 'tutorial']:
        depends_on("flecsi-deps@1.4 +%s" % v, when="+%s" % v)
        depends_on("flecsi-deps@1.4 ~%s" % v, when="~%s" % v)

    for b in ['Debug', 'Release', 'RelWithDebInfo', 'MinSizeRel']:
        depends_on("libristra build_type=%s" % b,
            when="build_type=%s" % b)

    depends_on('cmake@3.12:')
    # Requires cinch > 1.0 due to cinchlog installation issue
    depends_on('cinch@1.01:', type='build', when='+cinch')
    depends_on('exodusii')
    #depends_on('mpi', when='backend=mpi')
    #depends_on('mpi', when='backend=legion')
    #depends_on('mpi', when='backend=hpx')
    #depends_on('mpi', when='backend=charmpp')
    #depends_on('legion@ctrl-rep-5 +shared +mpi +hdf5', when='backend=legion')
    #depends_on('hpx@1.3.0 cxxstd=14 build_type=Release', when='backend=hpx')
    #depends_on('boost@1.70.0: cxxstd=14 +program_options')
    #depends_on('metis@5.1.0:')
    #depends_on('parmetis@4.0.3:')
    #depends_on('hdf5+hl+mpi')
    #depends_on('lua@5.3.5')
    
    #portage requires LAPACKE
    depends_on('netlib-lapack lapacke=true', when='+portage')
    