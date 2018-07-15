include ../../Makefile.inc

OBJS=gauss2lus.o mystring.o
EXE1=gauss2lus.exe

all:	$(EXE1) 

$(EXE1):	$(OBJS)
	$(CC) -o $(EXE1) $(OBJS) -lm

gauss2lus.o:	gauss2lus.c mystring.h
	$(CC) -c gauss2lus.c

mystring.o:	mystring.h
	$(CC) -c mystring.c

clean:	
	rm -f *.exe $(OBJS)
