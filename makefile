all: main.c nw.c upgma.c
	gcc -O3 -o main main.c nw.c upgma.c
