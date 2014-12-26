# debug is default because that will usually be what I want until I have a
# working alpha.
debug:
	gcc -lgmp -g -DDEBUG shlife.c

build:
	gcc -lgmp shlife.c
