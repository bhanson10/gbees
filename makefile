# makefile, 
# Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

CC=gcc
CFLAGS=-pedantic -Wall -std=c99 -D_GNU_SOURCE -fPIC <path_to_lapack_include>
LDFLAGS=-shared <path_to_lapack_lib> <path_to_lapack_openblas>
LDLIBS=-lm -llapack -lopenblas

SRCS = $(wildcard *.c) $(wildcard ../../*.c)
OBJS = $(SRCS:%.c=%.o)

TARGET = gbees.so

.DEFAULT: all

all: link

link: ${OBJS}
	${CC} ${LDFLAGS} -o ${TARGET} ${OBJS} ${LDLIBS}

%.c.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

clean:
	find ./ -type f -name '*.o' -delete