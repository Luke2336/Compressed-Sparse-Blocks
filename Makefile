CXX := g++
TARGET := csb
OPENMPFLAG := 
CXXFLAGS := -std=c++11 -O3
WARNINGFLAGS := -Wall -Wextra
INCLUDE := src
SRC_DIRS := src
SRCS := $(wildcard $(SRC_DIRS:=/*.cpp))
OBJS := $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)

ifndef BOOST_LIBRARY_PATH
BOOST_LIBRARY_PATH := ""
endif

all: $(TARGET)

$(TARGET): $(OBJS) $(FLUTE_WRAPPER_LIB)
	$(CXX) -o $@ $^ -lpthread $(OPENMPFLAG)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(WARNINGFLAGS) $(OPENMPFLAG) -I $(INCLUDE) -isystem $(BOOST_LIBRARY_PATH) -MMD -c $< -o $@

clean:
	rm -rf $(TARGET) $(OBJS) $(DEPS)

.PHONY: all clean test
-include $(DEPS)