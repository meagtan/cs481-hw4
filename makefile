all: main.c nw.c upgma.c
	gcc -O0 -o main main.c nw.c upgma.c
