# makefile, https://github.com/bhanson10/gbees/tree/main/examples/Aizawa
# Copyright 2025 by Benjamin Hanson, published under BSD 3-Clause License.

CC=gcc
HOME_DIR=/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/Utilities/Functions/c

# DEBUG
#CFLAGS=-pedantic -Wall -std=c99 -D_GNU_SOURCE -static -g

# RELEASE
CFLAGS=-pedantic -Wall -std=c99 -D_GNU_SOURCE -fPIC

LDFLAGS=
LDLIBS=-lm

SRCS = $(wildcard *.c) $(wildcard ../../*.c)
OBJS = $(SRCS:%.c=%.o)

TARGET = Aizawa

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

mc_Aizawa: frame_conversion.so
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I$(HOME_DIR)/ -Wl,-rpath,./ frame_conversion.so mc_Aizawa.c -o mc_Aizawa
	@echo ""
	@echo "Problem file compiled successfully."
	./mc_Aizawa

frame_conversion.so: 
	@-rm -f frame_conversion.so
	$(CC) -pedantic -Wall -std=c99 -D_GNU_SOURCE -fPIC -shared -o frame_conversion.so $(HOME_DIR)/frame_conversion.c -lm

clean:
	find ../.. -type f -name '*.o' -delete
	rm -f ${TARGET}
	@-rm -vf mc_Aizawa
	@-rm -vf frame_conversion.so