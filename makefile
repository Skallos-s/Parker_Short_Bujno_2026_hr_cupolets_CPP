# File: makefile
# Purpose: compiling
# Author: Daniel Bujno
# Followed tutorial by Ktitagya Agarwal (https://www.geeksforgeeks.org/cpp/makefile-in-c-and-its-applications/)
# Followed tutorial by MestreLionn and Chnossos (https://stackoverflow.com/questions/30573481/how-to-write-a-makefile-with-separate-source-and-header-directories)

# Compiler
CXX := g++

# Location of source, header, and object files
SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := .

# Compilar flags
CXXFLAGS := -Wall -g
CPPFLAGS := -Iinclude -MMD -MP # -I is a preprocessor flag, not a compiler flag
LDFLAGS  := -Llib              # -L is a linker flag

# Target executable
TARGET := $(BIN_DIR)/GENERATE.exe

# For deleting the target
TARGET_DEL := $(BIN_DIR)/GENERATE.exe

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)

# Object files
OBJS := $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

.PHONY: all clean

# Default rule to build and run the executable
all: $(TARGET) run

# Rule to link object files into the target executable
$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CXX) $(LDLIBS) -o $(TARGET) $(OBJS) -std=c++17 -lstdc++fs

# Rule to compile .cpp files into .o files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Make binary directory if it does not exist
$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

# Rule to run the executable
run: $(TARGET)
	./$(TARGET)

# Clean rule to remove generated files
clean:
	@$(RM) -rv $(BIN_DIR) $(OBJ_DIR)

-include $(OBJ:.o=.d)
