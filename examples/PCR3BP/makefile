# makefile, https://github.com/bhanson10/gbees/tree/main/examples/PCR3BP
# Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

CC=gcc

CFLAGS=-pedantic -Wall -std=c99 -D_GNU_SOURCE

# remove -g for release version, and add -DRELEASE -fPIC
LDFLAGS=
LDLIBS=-lm

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
	mkdir -p results/c/P0 results/c/P1 results/c/P2 results/c/P3

folders_py:
	mkdir -p results/python/P0 results/python/P1 results/python/P2 results/python/P3

link: ${OBJS}
	${CC} ${LDFLAGS} -o ${TARGET} ${OBJS} ${LDLIBS}

%.c.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

clean:
	find ../.. -type f -name '*.o' -delete
	rm -f ${TARGET}