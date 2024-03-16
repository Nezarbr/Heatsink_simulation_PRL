CC=mpicc
CFLAGS=-Wall -O3 -lm
SRC_DIR=src
BIN_DIR=bin
EXE_FILES=$(patsubst $(SRC_DIR)/%.c,$(BIN_DIR)/%.exe,$(wildcard $(SRC_DIR)/*.c))

all: $(EXE_FILES)

$(BIN_DIR)/%.exe: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(BIN_DIR)/*.exe