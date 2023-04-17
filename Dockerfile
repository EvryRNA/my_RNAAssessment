FROM gcc:9.5.0
WORKDIR src
COPY Makefile .
COPY src src
RUN make
COPY . .
CMD /bin/bash
