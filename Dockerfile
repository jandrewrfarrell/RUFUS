FROM ubuntu:16.04

# A quick and dirty Dockerfile to get RUFUS running

ENV DEBIAN_FRONTEND=noninteractive 

# Dependencies
RUN set -ex; \
	apt-get update; \
	apt-get install -y python cmake wget build-essential software-properties-common; \
	add-apt-repository ppa:ubuntu-toolchain-r/test; \
	apt-get install -y g++-4.9 zlib1g-dev libbz2-dev libbz2-dev liblzma-dev bc libncurses5-dev git; \
	echo done

# Building
# Using kohrar/RUFUS because of a single line preventing compilation
RUN set -ex; \
	git clone https://github.com/kohrar/RUFUS.git; \
	mkdir -p /RUFUS/bin; \
	cd /RUFUS/bin; \
	cmake ../ -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++); \
	make ; \
	echo done
	
# Runtime tools
RUN set -ex; \
	apt install samtools; \
	echo done


