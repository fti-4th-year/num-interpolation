.phony: all

all: run

run: main.cpp
	g++ -std=c++11 -o $@ $<
