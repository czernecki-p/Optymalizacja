
# Makefile for compiling the optimization project

CXX = g++
CXXFLAGS = -Wall -O2 -std=c++11
TARGET = main_program

# Source files
SOURCES = main.cpp matrix.cpp ode_solver.cpp opt_alg.cpp solution.cpp user_funs.cpp
HEADERS = matrix.h ode_solver.h opt_alg.h solution.h user_funs.h

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: clean

