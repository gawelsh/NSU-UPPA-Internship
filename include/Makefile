CC = gcc
override CFLAGS += -g -O2 -Wall -pedantic
override LDFLAGS += -lm -lgsl -lblas
HDF5INC := $(shell pkg-config --cflags hdf5)
HDF5LIB := $(shell pkg-config --libs hdf5)

GCCVERSION = $(shell gcc -dumpversion | cut -f1 -d.)
ifeq "$(GCCVERSION)" "4"
    override CFLAGS += -std=c99
endif

DEPDIR := .deps
DEPFLAGS = -E -MMD -MP -MF $@

allsrc := $(wildcard *.c)
allobj := $(patsubst %.c,%.o,$(allsrc))
pattern=" main ?\("
maintp := $(shell grep -lE $(pattern) *c)
src := $(filter-out $(maintp),$(allsrc))
obj := $(patsubst %.c,%.o,$(src))
objtp := $(patsubst %.c,%.o,$(maintp))
exetp := $(patsubst %.c,%.exe,$(maintp))
DEPFILES := $(allsrc:%.c=$(DEPDIR)/%.d)

all : $(DEPFILES) $(exetp)

.SECONDARY: $(allobj)
.PHONY: $(DEPDIR)

%.o : %.c
%.o : %.c $(DEPDIR)/%.d | $(DEPDIR)
	$(CC) -c $(CFLAGS) $(HDF5INC) -o $@ $<

$(exetp): %.exe : %.o
ifeq ($(MAKELEVEL),1)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(HDF5LIB)
else
	@make --no-print-directory $@
endif

clean:
	rm -f $(allobj) $(exetp)
	rm -rf $(DEPDIR)
clean-data:
	rm -f *.dat *.h5 *.xmf

$(DEPDIR): ; @mkdir -p $@

$(DEPFILES): $(DEPDIR)/%.d : %.c| $(DEPDIR)
	@$(CC) $(DEPFLAGS) -c $(CFLAGS) $(HDF5INC) $< > $<.p
	@echo '' >> $@
	@awk '/^$(patsubst $(DEPDIR)/%.d,%.o,$@):/,/^$$/' $@ | sed -E 's/\.o/.exe/g;s/\.(c|h)/.o/g;' >> $@
	@rm -f $<.p


include $(wildcard $(DEPFILES))
