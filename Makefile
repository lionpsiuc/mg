CC      = gcc
CFLAGS  = -Wall -Wextra -O3 -std=c2x -march=native -D_GNU_SOURCE
LDFLAGS = -lm

TARGETS     = mg1 mg2
CSV_FILES   = summary.csv residuals.csv comparison.csv
PLOT_SCRIPT = plot.py
PLOTS       = 1_perf-vs-lmax.png 1_residuals-*.png 2_comparison.png

all: $(TARGETS)

mg1: mg1.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

mg2: mg2.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

run: mg1 mg2
	./mg1
	./mg2
	python3 $(PLOT_SCRIPT)

clean:
	rm -f $(TARGETS) $(CSV_FILES) $(PLOTS)

.PHONY: all run clean
