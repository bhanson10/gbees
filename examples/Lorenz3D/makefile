# makefile, https://github.com/bhanson10/gbees/tree/main/examples/Lorenz3D
# Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

CC=gcc

# DEBUG
#CFLAGS=-pedantic -Wall -std=c99 -D_GNU_SOURCE -static -g

# RELEASE
CFLAGS=-pedantic -Wall -std=c99 -D_GNU_SOURCE -fPIC -I<path_to_lapack_include> 

LDFLAGS=-L<path_to_lapack_lib>  -L<path_to_openblas_include> 
LDLIBS=-lm -llapack -lopenblas

SRCS = $(wildcard *.c) $(wildcard ../../*.c)
OBJS = $(SRCS:%.c=%.o)

TARGET = Lorenz3D

MODE ?= c

.DEFAULT: all

all: $(MODE)

c: link folders_c 
	./${TARGET}

python: folders_py
	python3 -u ${TARGET}.py

folders_c:
	mkdir -p results/c/P0 results/c/P1 results/c/C0 results/c/C1 

folders_py:
	mkdir -p results/python/P0 results/python/P1 results/python/C0 results/python/C1 

link: ${OBJS}
	${CC} ${LDFLAGS} -o ${TARGET} ${OBJS} ${LDLIBS}

%.c.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

clean:
	find ../.. -type f -name '*.o' -delete
	rm -f ${TARGET}