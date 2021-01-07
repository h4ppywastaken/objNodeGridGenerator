CC = g++
LD = g++

TARGET = objNodeGridGenerator

CPPFLAGS = -O3 -Wall -std=c++11 -fopenmp
LDLIBS = -fopenmp
INCLUDES = -Isource

SRC_DIR = source
BUILD_DIR = build
VPATH = source

# Rules
all: $(TARGET)

$(TARGET): $(BUILD_DIR)/$(TARGET).o
	$(CC) $(CPPFLAGS) $(INCLUDES) $^ -o $(TARGET)

$(BUILD_DIR)/$(TARGET).o: $(TARGET).cpp
	$(CC) $(CPPFLAGS) $(INCLUDES) -c $^ -o $@

clean:
	rm -f $(BUILD_DIR)/*.o $(TARGET) 

.PHONY: clean

# Dependencies
$(TARGET): $(OBJ)
