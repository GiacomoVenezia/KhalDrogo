CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall
SRCDIR = src
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
TARGET = drogo
INCLUDES = -I$(SRCDIR) -I.

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SOURCES) -o $(TARGET)

clean:
	rm -f $(TARGET)

.PHONY: clean
