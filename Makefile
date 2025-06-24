# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -O3 -ffast-math -march=native -std=c++17

# Library definition
LIB = functions
LIB_DIR = lib
LIB_FILE = $(LIB_DIR)/lib$(LIB).a

# Directories
SRC_DIR = src
INCLUDE_DIR = include
HEADERS_DIR = $(INCLUDE_DIR)/functions
BUILD_DIR = build
TEST_DIR = tests
EXTERNAL_DIRS = $(filter %/,$(wildcard external/*/))

# Define include flags
INCLUDE_FLAGS := -I$(HEADERS_DIR) $(addprefix -I, $(EXTERNAL_DIRS))

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
$(LIB_FILE): $(OBJECTS) | $(LIB_DIR)
	ar rcs $@ $^

# Create lib directory if it does not exist
$(LIB_DIR): 
	@mkdir -p $(LIB_DIR)

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS_DIR)/%.hpp | $(BUILD_DIR)
	$(CXX) $(CXX_FLAGS) $(INCLUDE_FLAGS) -c $< -o $@

# Create build directory if it does not exist
$(BUILD_DIR): 
	@mkdir -p $(BUILD_DIR)

# Build tests
test: $(LIB_FILE) $(TEST_EXES)
	echo $(TEST_EXES)
	@for test_exe in $(TEST_EXES); do \
		./$$test_exe; \
	done

# Test targets
$(BUILD_DIR)/%.exe: $(TEST_DIR)/%.cpp $(BUILD_DIR)
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE_DIR) $(INCLUDE_FLAGS) -I$(GTEST_DIR) $< -L$(LIB_DIR) -l$(LIB) -o $@ $(GTEST_LIBS)


# Clean build files
clean:
	rm -rf $(BUILD_DIR)
	rm -rf $(LIB_DIR)
