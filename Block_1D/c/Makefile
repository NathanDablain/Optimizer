CC = gcc
CCFLAGS = -Wall -O2
CCDBGFLAGS = -g -Wall -O0
CPPFLAGS = -MMD -MP

BUILD_DIR_DBG = Build_Debug
BUILD_DIR_REL = Build_Release

CCSOURCES = $(wildcard *.c)
CCOBJECTS = $(CCSOURCES:%.c=$(BUILD_DIR_REL)/%.o)
CCDBGOBJECTS = $(CCSOURCES:%.c=$(BUILD_DIR_DBG)/%.o)
CCDEPENDS = $(CCOBJECTS:.o=.d)
CCDBGDEPENDS = $(CCDBGOBJECTS:.o=.d)

default: all

all: $(BUILD_DIR_REL)/Block_1D.exe $(BUILD_DIR_DBG)/Block_1D_DBG.exe

release: $(BUILD_DIR_REL)/Block_1D.exe

debug: $(BUILD_DIR_DBG)/Block_1D_DBG.exe

$(BUILD_DIR_REL)/Block_1D.exe: $(CCOBJECTS) 
	$(CC) $^ -o $@

$(BUILD_DIR_DBG)/Block_1D_DBG.exe: $(CCDBGOBJECTS) 
	$(CC) $^ -o $@

$(BUILD_DIR_REL)/%.o: %.c
	mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CCFLAGS) -c $< -o $@

$(BUILD_DIR_DBG)/%.o: %.c
	mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CCDBGFLAGS) -c $< -o $@

-include $(CCDEPENDS) $(CCDBGDEPENDS)


.PHONY: clean

clean:
	rm -rf $(BUILD_DIR_REL)/Block_1D.exe $(BUILD_DIR_DBG)/Block_1D_DBG.exe $(CCOBJECTS) $(CCDBGOBJECTS)

.PHONY: clean_dep

clean_dep:
	rm -rf$(CCDEPENDS) $(CCDBGDEPENDS)

.PHONY: pristine

pristine:
	rm -rf $(BUILD_DIR_REL)/Block_1D.exe $(BUILD_DIR_DBG)/Block_1D_DBG.exe $(CCDEPENDS) $(CCOBJECTS) $(CCDBGOBJECTS) $(CCDBGDEPENDS)