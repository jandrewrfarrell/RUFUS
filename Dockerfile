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
	apt-get install -y python wget git bc libgomp1 $BUILD_DEPS; \
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

# Setup samtools 1.11
RUN set -ex; \
	BUILD_DEPS="cmake build-essential libncurses5-dev zlib1g-dev libbz2-dev libbz2-dev liblzma-dev"; \
	apt-get install -y $BUILD_DEPS; \
# Get samtools and build it
	cd /; \
	wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2; \
	tar -xjf samtools-1.11.tar.bz2; \
	cd samtools*; \
	./configure; \
	make; \
	make install; \
	cd /; \
	rm -rf samtools*; \
# Cleanup
	apt-get purge -y --auto-remove $BUILD_DEPS; \
	apt-get clean; \
	echo done

