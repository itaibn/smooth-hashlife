# debug is default because that will usually be what I want until I have a
# working alpha.
debug:
	gcc -lgmp -g -DDEBUG -Wall shlife.c -o shlife

build:
	gcc -lgmp shlife.c
