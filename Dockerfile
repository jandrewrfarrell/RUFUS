#
# A Dockerfile to get RUFUS running
#

# GCC 4.9 only available up to 16.04
FROM ubuntu:16.04

ARG DEBIAN_FRONTEND=noninteractive 

COPY . /RUFUS

RUN set -ex; \
# Dependencies
	BUILD_DEPS="cmake build-essential g++-4.9 zlib1g-dev libbz2-dev libbz2-dev liblzma-dev libncurses5-dev"; \
	apt-get update; \
	apt-get install -y software-properties-common; \
	add-apt-repository ppa:ubuntu-toolchain-r/test; \
	apt-get install -y python wget git bc $BUILD_DEPS; \
# Build
	mkdir -p /RUFUS/bin; \
	cd /RUFUS/bin; \
	cmake ../ -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++); \
	make; \
# Cleanup
	apt-get purge -y --auto-remove $BUILD_DEPS; \
	apt-get clean; \
	echo done
	
# Runtime tools
RUN set -ex; \
	apt install samtools; \
	echo done


