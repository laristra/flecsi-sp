# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack import *


class FlecsiSp(CMakePackage):
    '''Flecsi-SP contains various specializations for use with the FleCSI core programming system.
    '''
    git      = 'https://https://github.com/laristra/flecsi-sp.git'

    version('master', branch='master', submodules=False, preferred=True)

    variant('build_type', default='Release', values=('Debug', 'Release'),
            description='The build type to build', multi=False)
    variant('backend', default='mpi', values=('serial', 'mpi', 'legion', 'charm++', 'hpx'),
            description='Backend to use for distributed memory', multi=False)
    variant('cinch', default=False,
            description='Enable External Cinch')

    depends_on('cmake@3.12:',  type='build')
    # Requires cinch > 1.0 due to cinchlog installation issue
    depends_on('cinch@1.01:', type='build', when='+cinch')
    #depends_on('mpi', when='backend=mpi')
    #depends_on('mpi', when='backend=legion')
    #depends_on('mpi', when='backend=hpx')
    #depends_on('legion@ctrl-rep +shared +mpi +hdf5', when='backend=legion')
    #depends_on('boost@1.70.0: cxxstd=14 +program_options')
    #depends_on('metis@5.1.0:')
    #depends_on('parmetis@4.0.3:')
    depends_on('libristra +cinch')
    depends_on('flecsi +cinch backend=mpi', when='backend=mpi')
    depends_on('flecsi +cinch +hdf5 backend=legion', when='backend=legion')
    depends_on('exodusii +mpi')

    def cmake_args(self):
        spec = self.spec
        options = []

        if '+cinch' in spec:
            options.append('-DCINCH_SOURCE_DIR=' + spec['cinch'].prefix)

        if self.run_tests:
            options.append('-DENABLE_UNIT_TESTS=ON')
        else:
            options.append('-DENABLE_UNIT_TESTS=OFF')

        return options
