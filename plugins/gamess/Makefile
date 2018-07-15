include ../../Makefile.inc

OBJS=gamess2lus.o mystring.o
EXE1=gamess2lus.exe

all:	$(EXE1) 

$(EXE1):	$(OBJS)
	$(CC) -o $(EXE1) $(OBJS) -lm

gamess2lus.o:	gamess2lus.c mystring.h
	$(CC) -c gamess2lus.c

mystring.o:	mystring.h
	$(CC) -c mystring.c

clean:	
	rm -f *.exe $(OBJS)
