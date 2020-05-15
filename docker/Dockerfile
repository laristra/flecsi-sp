FROM laristra/libristra:fedora30

SHELL ["/bin/bash", "-c"]	
		
ENV DEBIAN_FRONTEND=noninteractive \
    HOME=/home/laristra \
    BASH_ENV=~/.bashrc

ARG CC
ARG CXX
ARG CXXFLAGS
ARG BUILD_TYPE
ARG MPI_PACKAGE
ARG RUNTIME

#for coverage
ARG CI
ARG TRAVIS
ARG TRAVIS_BRANCH
ARG TRAVIS_JOB_NUMBER
ARG TRAVIS_PULL_REQUEST 
ARG TRAVIS_JOB_ID
ARG TRAVIS_TAG
ARG TRAVIS_REPO_SLUG
ARG TRAVIS_COMMIT
ARG TRAVIS_OS_NAME

COPY flecsi-sp $HOME/flecsi-sp
USER root
RUN chown -R laristra:laristra $HOME/flecsi-sp
USER laristra

WORKDIR $HOME
RUN module list && \
    echo $CMAKE_PREFIX_PATH && \
    git clone https://github.com/laristra/ristra_spackages.git

RUN export SPACK_SPEC="flecsi-sp-deps%gcc backend=$RUNTIME ^$MPI_PACKAGE" && \
    spack repo add --scope site $HOME/ristra_spackages/spack-repo && \
    spack repo add --scope site $HOME/flecsi-sp/spack-repo && \
    spack find && \
    echo "spack install --show-log-on-error $SPACK_SPEC" && \
    spack install --show-log-on-error $SPACK_SPEC && \
    rm deps.sh && \
    rm $HOME/load_deps.sh && \
    spack module tcl loads python | tee -a deps.sh && \
    spack module tcl loads parmetis ^$MPI_PACKAGE | tee -a deps.sh && \
    spack module tcl loads --dependencies $SPACK_SPEC | tee -a deps.sh
RUN echo ". $HOME/spack/share/spack/setup-env.sh && . $HOME/deps.sh" >> $HOME/load_deps.sh

RUN source $HOME/load_deps.sh && module list

WORKDIR $HOME/flecsi-sp
RUN mkdir build
WORKDIR build

RUN source $HOME/load_deps.sh && \
    cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
          -DFLECSI_RUNTIME_MODEL=$RUNTIME \
          -DENABLE_UNIT_TESTS=ON \
          -DENABLE_FLECSIT=OFF \
          -DENABLE_FLECSI_TUTORIAL=OFF \
          .. && \
    make -j2 || make VERBOSE=1 && \
    ctest --output-on-failure && \
    make install DESTDIR=${PWD}/install && rm -rf ${PWD}/install

USER root
RUN make install

USER laristra
WORKDIR $HOME
