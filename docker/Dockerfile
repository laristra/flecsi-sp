FROM laristra/libristra:fedora_serial

ARG CC
ARG CXX
ARG CXXFLAGS
ARG BUILD_TYPE

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

COPY flecsi-sp /home/flecsi/flecsi-sp
USER root
RUN chown -R flecsi:flecsi /home/flecsi/flecsi-sp
USER flecsi

WORKDIR /home/flecsi/flecsi-sp
RUN mkdir build

WORKDIR build

RUN cmake -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          -DENABLE_UNIT_TESTS=ON \
          ..
RUN make -j2 || make VERBOSE=1
#RUN ctest --output-on-failure
RUN make install DESTDIR=${PWD}

WORKDIR /home/flecsi
