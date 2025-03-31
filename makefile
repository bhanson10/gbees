# makefile, 
# Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

CC=gcc

# DEBUG
#CFLAGS=-pedantic -Wall -std=c99 -D_GNU_SOURCE -static -g

# RELEASE
CFLAGS=-pedantic -Wall -std=c99 -D_GNU_SOURCE -fPIC -I<path_to_lapack_include> 

LDFLAGS=-shared -L<path_to_lapack_lib>  -L<path_to_openblas_lib> 
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
