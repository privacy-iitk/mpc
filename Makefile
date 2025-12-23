CXX = g++
CXXFLAGS = -std=c++20 -march=native -O3 -g -Wall -pedantic -Wno-ignored-attributes \
           -fopenmp -msse2 -mssse3 -maes -mpclmul -DBOOST_ERROR_CODE_HEADER_ONLY
LIBS = -lbsd -fopenmp -lboost_system -pthread

# ---------------------------------------------------------------------
# Runtime parameters (NO compile-time effect)
# ---------------------------------------------------------------------
PREPROCESS ?= 0
ROLE ?= 0

ifeq ($(PREPROCESS),1)
RUNTIME_ARGS = -p --role $(ROLE)
else
RUNTIME_ARGS = --role $(ROLE)
endif

# ---------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------
TARGET = mpc_async
SRC = main.cpp

TEST1 = mpcops_test
TEST1_SRC = mpcops_test.cpp

TEST2 = test_shares
TEST2_SRC = shares_test.cpp

TEST3 = locoram_test
TEST3_SRC = locoram_test.cpp

TEST4 = dpf_test
TEST4_SRC = dpf_test.cpp

# ---------------------------------------------------------------------
# Default rule
# ---------------------------------------------------------------------
all: $(TARGET) $(TEST1) $(TEST2)

# ---------------------------------------------------------------------
# Main binary
# ---------------------------------------------------------------------
$(TARGET): $(SRC) common.hpp shares.hpp mpcops.hpp dpf/dpf.h doram.hpp doram.tpp
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LIBS)

# ---------------------------------------------------------------------
# MPC operations test
# ---------------------------------------------------------------------
$(TEST1): $(TEST1_SRC) shares.hpp mpcops.hpp types.hpp prg.hpp
	$(CXX) $(CXXFLAGS) -DBOOST_ASIO_HAS_CO_AWAIT $(TEST1_SRC) -o $(TEST1) $(LIBS)

# ---------------------------------------------------------------------
# Shares test
# ---------------------------------------------------------------------
$(TEST2): $(TEST2_SRC)
	$(CXX) $(CXXFLAGS) $(TEST2_SRC) -o $(TEST2) $(LIBS)

# ---------------------------------------------------------------------
# Locoram test
# ---------------------------------------------------------------------
$(TEST3): $(TEST3_SRC) locoram.tpp
	$(CXX) $(CXXFLAGS) $(TEST3_SRC) -o $(TEST3) $(LIBS)

# ---------------------------------------------------------------------
# DPF test
# ---------------------------------------------------------------------
$(TEST4): $(TEST4_SRC) dpf.hpp types.hpp
	$(CXX) $(CXXFLAGS) $(TEST4_SRC) -o $(TEST4) $(LIBS)

# ---------------------------------------------------------------------
# Run helpers (THIS IS THE NEW PART)
# ---------------------------------------------------------------------
run: $(TARGET)
	./$(TARGET) $(RUNTIME_ARGS)

run-test1: $(TEST1)
	./$(TEST1) $(RUNTIME_ARGS)

run-test2: $(TEST2)
	./$(TEST2)

# ---------------------------------------------------------------------
# Clean
# ---------------------------------------------------------------------
clean:
	rm -f $(TARGET) $(TEST1) $(TEST2) $(TEST3) $(TEST4)


# CXX = g++
# CXXFLAGS = -std=c++20 -march=native -O3 -g -Wall -pedantic -Wno-ignored-attributes \
#            -fopenmp -msse2 -mssse3 -maes -mpclmul -DBOOST_ERROR_CODE_HEADER_ONLY
# LIBS = -lbsd -fopenmp -lboost_system -pthread

# # Main target
# TARGET = mpc_async
# SRC = main.cpp

# # Test targets
# TEST1 = mpcops_test
# TEST1_SRC = mpcops_test.cpp

# TEST2 = test_shares
# TEST2_SRC = shares_test.cpp

# TEST3 = locoram_test
# TEST3_SRC = locoram_test.cpp

# TEST4 = dpf_test
# TEST4_SRC = dpf_test.cpp

# # ---------------------------------------------------------------------
# # Default rule: build everything
# # ---------------------------------------------------------------------
# all: $(TARGET) $(TEST1) $(TEST2)

# # ---------------------------------------------------------------------
# # Main binary
# # ---------------------------------------------------------------------
# $(TARGET): $(SRC) common.hpp shares.hpp mpcops.hpp dpf/dpf.h doram.hpp doram.tpp
# 	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LIBS)

# # ---------------------------------------------------------------------
# # MPC operations test
# # ---------------------------------------------------------------------
# $(TEST1): $(TEST1_SRC) shares.hpp mpcops.hpp types.hpp prg.hpp
# 	$(CXX) $(CXXFLAGS) -DBOOST_ASIO_HAS_CO_AWAIT $(TEST1_SRC) -o $(TEST1) $(LIBS)

# # ---------------------------------------------------------------------
# # Shares test
# # ---------------------------------------------------------------------
# $(TEST2): $(TEST2_SRC)
# 	$(CXX) $(CXXFLAGS) $(TEST2_SRC) -o $(TEST2) $(LIBS)


# # ---------------------------------------------------------------------
# # Locoram test
# # ---------------------------------------------------------------------
# $(TEST3): $(TEST3_SRC) locoram.tpp
# 	$(CXX) $(CXXFLAGS) $(TEST3_SRC) -o $(TEST3) $(LIBS)


# # ---------------------------------------------------------------------
# # DPF test
# # ---------------------------------------------------------------------
# $(TEST4): $(TEST4_SRC) dpf.hpp types.hpp
# 	$(CXX) $(CXXFLAGS) $(TEST4_SRC) -o $(TEST4) $(LIBS)

# # ---------------------------------------------------------------------
# # Clean rule
# # ---------------------------------------------------------------------
# clean:
# 	rm -f $(TARGET) $(TEST1) $(TEST2) $(TEST3)
