CC = clang
CFLAGS = -Wall -Wextra -Werror -Wno-newline-eof -Wno-unused-parameter -pedantic -std=c99 -ferror-limit=1
LDFLAGS = -lm
DIFF = diff
# -y diffs side by side.

PROGRAM = program
SOURCE_FILES = $(shell find . -type f -name '*.c')
HEADER_FILES = $(shell find . -type f -name '*.h')
OUTPUT_FILES = $(shell find . -type f -name '*myoutput.txt')
OBJECT_FILES = $(SOURCE_FILES:.c=.o)

.PHONY: all build clean run test0 test1 test2 diff0 diff1 diff2

all: build

diffall: diff0 diff1 diff2

build: $(PROGRAM)

clean:
	rm -f $(PROGRAM) $(OBJECT_FILES) $(OUTPUT_FILES)

run: test0

diff0: build
	./$(PROGRAM) < test0.txt > test0-myoutput.txt
	$(DIFF) test0-myoutput.txt test0-out.txt

diff1: build
	./$(PROGRAM) < test1.txt > test1-myoutput.txt
	$(DIFF) test1-myoutput.txt test1-out.txt

diff2: build
	./$(PROGRAM) < test2.txt > test2-myoutput.txt
	$(DIFF) test2-myoutput.txt test2-out.txt

test0: build
	./$(PROGRAM) < test0.txt

test1: build
	./$(PROGRAM) < test1.txt

test2: build
	./$(PROGRAM) < test2.txt

$(PROGRAM): $(OBJECT_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)