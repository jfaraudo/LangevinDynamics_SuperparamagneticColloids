CC=g++
CFLAGS=-c -Wall -O2
LDFLAGS=
SOURCES=app.aggregation.cpp aggregate.cpp bead.cpp dump.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXE=aggregation-28Mar2010.ex

	
$(EXE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

clean:
	rm $(EXE) $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
