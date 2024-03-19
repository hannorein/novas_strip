all: 
	gcc -lm test.c novas.c -o novas.bin
strip:
	gcc -fdata-sections -ffunction-sections -c novas.c -o novas.o 
	gcc -lm -Wl,--gc-sections -Wl,--print-gc-sections -o novas.bin novas.o
