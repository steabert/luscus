include ../../Makefile.inc

OBJS=qe2lus.o mystring.o
EXE1=qe2lus.exe

all:	$(EXE1) 

$(EXE1):	$(OBJS)
	$(CC) -o $(EXE1) $(OBJS) -lm

qe2lus.o:	qe2lus.c mystring.h
	$(CC) -c qe2lus.c

mystring.o:	mystring.c mystring.h
	$(CC) -c mystring.c

clean:	
	rm -f *.exe $(OBJS)
