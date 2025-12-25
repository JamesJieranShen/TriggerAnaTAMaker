# Makefile for the lazy.
BUILD_DIR := build
INSTALL_DIR := install
CMAKE_FLAGS := -DCMAKE_BUILD_TYPE=Release
TARGET ?= MakeTriggerActivity

.PHONY: all configure build install clean run

all: install

build:
	@mkdir -p $(BUILD_DIR)
	cmake -S . -B $(BUILD_DIR) $(CMAKE_FLAGS)
	cmake --build $(BUILD_DIR) -- $(MAKEFLAGS)

install: build
	@mkdir -p $(INSTALL_DIR)
	cmake --install $(BUILD_DIR) --prefix "$(abspath $(INSTALL_DIR))"

clean:
	rm -rf $(BUILD_DIR)
	rm -rf $(INSTALL_DIR)

run: all
	./$(BUILD_DIR)/$(TARGET)

