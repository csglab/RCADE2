TARGET = ../../bin/FASTAtoRF
CPP01 = commandline
CPP02 = find_zf
CPP03 = main
CPP04 = misc
CPP05 = open_read
CPP06 = read_fasta
CPP07 = write

H01 = declarations
H02 = structures
H03 = commandline

OBJS = $(CPP01).o $(CPP02).o $(CPP03).o $(CPP04).o $(CPP05).o $(CPP06).o $(CPP07).o
CC = g++	# change this line if you are using a different compiler

$(TARGET): $(OBJS)
	$(CC) -O3 -o $(TARGET) $(OBJS)

$(CPP01).o: $(CPP01).cpp
	$(CC) -O3 -c $(CPP01).cpp

$(CPP02).o: $(CPP02).cpp
	$(CC) -O3 -c $(CPP02).cpp

$(CPP03).o: $(CPP03).cpp
	$(CC) -O3 -c $(CPP03).cpp

$(CPP04).o: $(CPP04).cpp
	$(CC) -O3 -c $(CPP04).cpp

$(CPP05).o: $(CPP05).cpp
	$(CC) -O3 -c $(CPP05).cpp

$(CPP06).o: $(CPP06).cpp
	$(CC) -O3 -c $(CPP06).cpp

$(CPP07).o: $(CPP07).cpp
	$(CC) -O3 -c $(CPP07).cpp


$(OBJS): $(H01).h $(H02).h $(H03).h
