CC      = gcc
CFLAGS  = -D_GNU_SOURCE -Iinclude -march=native -O3 -std=c2x -Wall -Wextra
LDFLAGS = -lm

SRC_DIR     = src
OBJ_DIR     = obj
BIN_DIR     = bin
COMMON_DIR  = $(SRC_DIR)/common
Q1_DIR      = $(SRC_DIR)/q1
Q2_DIR      = $(SRC_DIR)/q2

COMMON_SRCS = multigrid.c

COMMON_OBJS = $(patsubst %.c,$(OBJ_DIR)/common/%.o,$(COMMON_SRCS))
Q1_OBJ      = $(OBJ_DIR)/q1/main.o
Q2_OBJ      = $(OBJ_DIR)/q2/main.o

TARGETS = $(BIN_DIR)/q1 $(BIN_DIR)/q2

all: directories $(TARGETS)

directories:
	@mkdir -p $(OBJ_DIR)/common $(OBJ_DIR)/q1 $(OBJ_DIR)/q2 $(BIN_DIR)

$(OBJ_DIR)/common/%.o: $(COMMON_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/q1/main.o: $(Q1_DIR)/main.c
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/q2/main.o: $(Q2_DIR)/main.c
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN_DIR)/q1: $(Q1_OBJ) $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(BIN_DIR)/q2: $(Q2_OBJ) $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

run: bin/q1 bin/q2
	@mkdir -p data figs
	./bin/q1
	./bin/q2
	python3 scripts/plot.py

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: all clean directories
