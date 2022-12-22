FROM ubuntu:22.04

RUN apt-get update
RUN apt-get install ffmpeg libsm6 libxext6 -y
RUN apt-get install freeglut3

RUN mkdir -p /AMTOPTW/instances
COPY /build/AMTOPTW /AMTOPTW/AMTOPTW
COPY /build/instances /AMTOPTW/instances

WORKDIR /AMTOPTW
ENTRYPOINT [ "AMTOPTW" ]