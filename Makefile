CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
CFLAGS = -Wall -I$(CFITSIO) 
LIBS = -L$(CFITSIO) -lcfitsio -lm
GLIBS = 
GLIBS += 
OBJECTS = subtractOverscan.o 
HEADERS = globalConstants.h

ALL : subtractOverscan.exe
	@echo "Listo!"

subtractOverscan.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o subtractOverscan.exe $(LIBS) $(GLIBS) $(CFLAGS)

subtractOverscan.o : subtractOverscan.cc $(HEADERS)
	$(CPP) -c subtractOverscan.cc -o subtractOverscan.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe
