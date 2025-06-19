# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -O3 -ffast-math -march=native -std=c++17

# Library definition
LIB = functions
LIB_FILE = lib$(LIB).a

# Directories
SRC_DIR = src
HEADERS_DIR = functions
BUILD_DIR = build
TEST_DIR = tests
FFT_DIR = external/fft

# Define GTest libraries
GTEST_LIBS = -lgtest -lgtest_main -pthread
GTEST_DIR = /usr/include/gtest

# Source and object files
HEADERS = $(wildcard $(HEADERS_DIR)/*.hpp)
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
TEST_SOURCES = $(HEADERS:$(HEADERS_DIR)/%.hpp=$(TEST_DIR)/Test%.cpp)
TEST_EXES = $(TEST_SOURCES:$(TEST_DIR)/%.cpp=$(BUILD_DIR)/%.exe)

# Default target
all: $(LIB_FILE)

# Create static library
$(LIB_FILE): $(OBJECTS)
	ar rcs $@ $^

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS_DIR)/%.hpp | $(BUILD_DIR)
	$(CXX) $(CXX_FLAGS) -I$(HEADERS_DIR) -I$(FFT_DIR) -c $< -o $@

# Create build directory if it does not exist
$(BUILD_DIR): 
	@mkdir -p $(BUILD_DIR)

# Build tests
test: $(LIB_FILE) $(TEST_EXES)
	echo $(TEST_EXES)
	@for test_exe in $(TEST_EXES); do \
		./$$test_exe || exit 1; \
	done

# Test targets
$(TEST_EXES): $(TEST_SOURCES) $(BUILD_DIR)
	$(CXX) $(CXX_FLAGS) -I. -I$(FFT_DIR) -I$(GTEST_DIR) $< -L. -l$(LIB) -o $@ $(GTEST_LIBS)


# Clean build files
clean:
	rm -rf $(BUILD_DIR)
	rm $(LIB_FILE)
