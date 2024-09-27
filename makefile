# makefile, 
# Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

CC=gcc

CFLAGS=-pedantic -Wall -std=c99 -D_GNU_SOURCE -fPIC

# remove -g for release version, and add -DRELEASE
LDFLAGS=-shared
LDLIBS=-lm

SRCS = $(wildcard *.c) $(wildcard ../../*.c)
OBJS = $(SRCS:%.c=%.o)

TARGET = gbees-hash.so

.DEFAULT: all

all: link

link: ${OBJS}
	${CC} ${LDFLAGS} -o ${TARGET} ${OBJS} ${LDLIBS}

%.c.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

clean:
	find ./ -type f -name '*.o' -delete