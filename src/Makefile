# Multi-step build: compile .cpp -> .o, then link to produce `fast2snp`

CC := g++
CXXFLAGS := -std=c++17 -fopenmp -O3 -I. -I../vendor/parallel-hashmap -Wall -Wextra -MMD -MP
LDFLAGS :=

# 目录定义
BUILD_DIR := ../build
BIN_DIR := $(BUILD_DIR)/bin
OBJ_DIR := $(BUILD_DIR)/obj

# 确保目录存在
$(shell mkdir -p $(BIN_DIR) $(OBJ_DIR))

SRCS := fast2snp.cpp main.cpp
OBJS := $(addprefix $(OBJ_DIR)/, $(notdir $(SRCS:.cpp=.o)))
DEPS := $(patsubst %.o,%.d, $(OBJS))

TARGET := $(BIN_DIR)/fast2snp
SAM2VCF := $(BIN_DIR)/sam2vcf

.PHONY: all clean bin

all: $(TARGET) $(SAM2VCF)

$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $@ $(CXXFLAGS) $(LDFLAGS)

$(SAM2VCF): sam2vcf.cpp
	$(CC) $(CXXFLAGS) -o $@ $<

$(OBJ_DIR)/%.o: %.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)

# Convenience target kept for compatibility
bin: $(TARGET)

-include $(DEPS)