.PHONY = all clean

CC = g++

INCLUDES = ./include/

SRCS := $(wildcard *.cpp)

BINS := $(SRCS:%.cpp=%)

all: ${BINS}

%: %.cpp
	${CC} -I ${INCLUDES} $< -o $@

clean:
	rm -rvf ${BINS}
