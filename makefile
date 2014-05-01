CC=gcc

CFLAGS=-O3 -Wall -fopenmp

LIBS=-lm -llapack -lblas

SRC_DIR=./src
INC_DIR=./include
OBJ_DIR=./obj
LIB_DIR=./lib
BIN_DIR=./bin

BIN=${BIN_DIR}/main

SRC=$(wildcard ${SRC_DIR}/*.c)
OBJ=$(patsubst ${SRC_DIR}/%.c,${OBJ_DIR}/%.o,$(SRC))

INC=-I${INC_DIR}
LDFLAGS=-L${LIB_DIR}

all: $(BIN)

$(OBJ): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $< $(LIBS)
	
$(BIN): $(OBJ)
	$(CC) $(CFLAGS) $(INC) -o $@ $^ $(LIBS)
clean: 
	rm ${OBJ_DIR}/*.o
	rm ${BIN_DIR}/*
