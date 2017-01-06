FROM laristra/flecsi:fedora_serial

ARG MPI

RUN rm -rf /home/flecsi/.ccache
COPY flecsi-sp /home/flecsi/flecsi-sp
COPY ccache/ /home/flecsi/.ccache
USER root
RUN chown -R flecsi:flecsi /home/flecsi/flecsi-sp /home/flecsi/.ccache
USER flecsi

WORKDIR /home/flecsi/flecsi-sp
RUN mkdir build

WORKDIR build
RUN cmake -DENABLE_MPI=$MPI \
          -DENABLE_UNIT_TESTS=ON ..

RUN make VERBOSE=1 -j2
RUN make test
RUN make install DESTDIR=${PWD}
USER root
RUN make install
USER flecsi
WORKDIR /home/flecsi
