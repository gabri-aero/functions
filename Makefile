# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -O2 -std=c++17

# Library definition
LIB = functions
LIB_FILE = lib$(LIB).a

# Directories
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build
FFT_DIR = external/fft/include

# Source and object files
HEADERS = $(wildcard $(INCLUDE_DIR)/*.hpp)
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Default target
all: $(LIB_FILE)

# Create static library
$(LIB_FILE): $(OBJECTS)
	ar rcs $@ $^

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(INCLUDE_DIR)/%.hpp | $(BUILD_DIR)
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE_DIR) -I$(FFT_DIR) -c $< -o $@

# Create build directory if it does not exist
$(BUILD_DIR): 
	@mkdir -p $(BUILD_DIR)


# Clean build files
clean:
	rm -rf $(BUILD_DIR)
