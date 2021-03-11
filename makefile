LD = gcc
CFLAGS = -Werror -lm -std=c99 -g 
VFLAGS =  -O2 -ftree-vectorize -fopt-info-vec-missed  $(CFLAGS)
GPFLAGS = -Wall -lm -std=c99 -pg 
OBJS = strassen_algorithm.c strassen_help.c
PNAME = strassen
LIBS = -lm

$(PNAME): strassen_algorithm.o strassen_help.o strassen_help.h
	$(LD) -o $@ $^ $(LIBS)

strassen_algorithm.o: strassen_algorithm.c strassen_help.h
	$(LD)  $(VFLAGS) -c -o $@  $<

strassen_help.o: strassen_help.c strassen_help.h
	$(LD) $(VFLAGS) -c -o $@  $<

strassen_algorithm.c: strassen_help.h

strassen_help.c: strassen_help.h

test: test_main.o test_help.o test_help.h
	$(LD) -o $@ $^ $(LIBS)

test_main.o: test_main.c test_help.h
	$(LD) -c -o $@  $(CFLAGS) $<

test_help.o: test_help.c test_help.h
	$(LD) -c -o $@  $(CFLAGS) $<

test_main.c: test_help.h

test_help.c: test_help.h

clean_test:
	rm *.o test

clean:
	rm *.o $(PNAME)
