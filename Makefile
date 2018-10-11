# c++ compiler
CXX=clang++
# arguments to c++ compiler
CXXFLAGS = -Wall -Wextra -O0 -g3 -std=c++11
# the name of the executable
PROGRAM:=project1

#
SOURCES:= main.cc code.cc 
OBJECTS:=$(addsuffix .o, $(basename $(SOURCES)))
LOGFILE:= logfile.txt TotalE.txt
all: $(PROGRAM)

%.o: %.cc
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

$(PROGRAM): $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS)

clean:
	rm -f $(OBJECTS) $(PROGRAM) $(LOGFILE)
clear: clean
