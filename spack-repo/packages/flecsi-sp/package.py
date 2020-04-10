# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack import *


class FlecsiSp(CMakePackage):
    '''Flecsi-SP contains various specializations for use with the FleCSI core programming system.
    '''
    git      = 'https://https://github.com/laristra/flecsi-sp.git'

    version('1.4', branch='1.4', submodules=False, preferred=True)

    variant('build_type', default='Release', values=('Debug', 'Release'),
            description='The build type to build', multi=False)
    variant('backend', default='mpi', values=('serial', 'mpi', 'legion', 'charmpp', 'hpx'),
            description='Backend to use for distributed memory', multi=False)
    variant('cinch', default=True,
            description='Enable External Cinch')

    depends_on('cmake@3.12:')
    # Requires cinch > 1.0 due to cinchlog installation issue
    depends_on('cinch@1.01:', type='build', when='+cinch')
    depends_on('mpi', when='backend=mpi')
    depends_on('mpi', when='backend=legion')
    depends_on('mpi', when='backend=hpx')
    depends_on('mpi', when='backend=charmpp')
    depends_on('legion@ctrl-rep-5 +shared +mpi +hdf5', when='backend=legion')
    depends_on('hpx@1.3.0 cxxstd=14 build_type=Release', when='backend=hpx')
    depends_on('boost@1.70.0: cxxstd=14 +program_options')
    depends_on('metis@5.1.0:')
    depends_on('parmetis@4.0.3:')
    depends_on('hdf5+hl+mpi')
    depends_on('lua@5.3.5')
    depends_on('exodusii')
    #depends_on('libristra +cinch')
    #depends_on('flecsi +cinch backend=mpi', when='backend=mpi')
    #depends_on('flecsi +cinch +hdf5 backend=legion', when='backend=legion')

    def cmake_args(self):
        spec = self.spec
        options = []

        if '+cinch' in spec:
            options.append('-DCINCH_SOURCE_DIR=' + spec['cinch'].prefix)

        if self.run_tests:
            options.append('-DENABLE_UNIT_TESTS=ON')
        else:
            options.append('-DENABLE_UNIT_TESTS=OFF')

        if spec.variants['backend'].value == 'legion':
            options.append('-DFLECSI_RUNTIME_MODEL=legion')
            options.append('-DENABLE_MPI=ON')
        elif spec.variants['backend'].value == 'mpi':
            options.append('-DFLECSI_RUNTIME_MODEL=mpi')
            options.append('-DENABLE_MPI=ON')
        elif spec.variants['backend'].value == 'hpx':
            options.append('-DFLECSI_RUNTIME_MODEL=hpx')
            options.append('-DENABLE_MPI=ON')
        elif spec.variants['backend'].value == 'charmpp':
            options.append('-DFLECSI_RUNTIME_MODEL=charmpp')
            options.append('-DENABLE_MPI=ON')
        else:
            options.append('-DFLECSI_RUNTIME_MODEL=serial')
            options.append('-DENABLE_MPI=OFF')

        return options
