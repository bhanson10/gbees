# makefile, https://github.com/bhanson10/gbees/tree/main/examples/PCR3BP
# Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

CC=gcc
CFLAGS=-pedantic -Wall -std=c99 -D_GNU_SOURCE -fPIC <path_to_lapack_include>
LDFLAGS=-shared <path_to_lapack_lib> <path_to_lapack_openblas>
LDLIBS=-lm -llapack -lopenblas

SRCS = $(wildcard *.c) $(wildcard ../../*.c)
OBJS = $(SRCS:%.c=%.o)

TARGET = PCR3BP

MODE ?= c

.DEFAULT: all

all: $(MODE)

c: link folders_c
	./${TARGET}

python: folders_py
	python3 -u ${TARGET}.py

folders_c:
	mkdir -p results/c/P0 results/c/P1 results/c/P2 results/c/P3 results/c/C0 results/c/C1 results/c/C2 results/c/C3

folders_py:
	mkdir -p results/python/P0 results/python/P1 results/python/P2 results/python/P3 results/python/C0 results/python/C1 results/python/C2 results/python/C3

link: ${OBJS}
	${CC} ${LDFLAGS} -o ${TARGET} ${OBJS} ${LDLIBS}

%.c.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

clean:
	find ../.. -type f -name '*.o' -delete
	rm -f ${TARGET}