# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -O3 -ffast-math -march=native -std=c++17

# Directories
INCLUDE_DIR = include
HEADERS_DIR = $(INCLUDE_DIR)/functions
BUILD_DIR = build
TEST_DIR = tests
EXTERNAL_DIRS = $(filter %/,$(wildcard external/*/))

# Define include flags
INCLUDE_FLAGS := -I. $(addprefix -I, $(EXTERNAL_DIRS))

# Define GTest libraries
GTEST_LIBS = -lgtest -lgtest_main -pthread
GTEST_DIR = /usr/include/gtest

# Source and object files
HEADERS = $(wildcard $(HEADERS_DIR)/*.hpp)
TEST_SOURCES = $(HEADERS:$(HEADERS_DIR)/%.hpp=$(TEST_DIR)/Test%.cpp)
TEST_EXES = $(TEST_SOURCES:$(TEST_DIR)/%.cpp=$(BUILD_DIR)/%.exe)

# Create build directory if it does not exist
$(BUILD_DIR): 
	@mkdir -p $(BUILD_DIR)

# Build tests
test: $(TEST_EXES)
	echo $(HEADERS_DIR)
	echo $(TEST_EXES)
	@for test_exe in $(TEST_EXES); do \
		./$$test_exe; \
	done

# Test targets
$(BUILD_DIR)/%.exe: $(TEST_DIR)/%.cpp $(BUILD_DIR)
	$(CXX) $(CXX_FLAGS) $(INCLUDE_FLAGS) -I$(GTEST_DIR) $< -L$(LIB_DIR) -l$(LIB) -o $@ $(GTEST_LIBS)


# Clean build files
clean:
	rm -rf $(BUILD_DIR)
