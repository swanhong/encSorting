RM = rm -rf
CC = g++
CFLAGS = -std=c++11 -O3
HEAANDIR = ../HEAAN
OBJS = CyclicArray.o SortingAlgorithm.o TestAlgorithm.o
LOCALLIBS = -lntl -lgmp -lm -lpthread
BOOTLIBS = -I../eigen
LIBS = -L/usr/local/lib $(HEAANDIR)/lib/libHEAAN.a -I$(HEAANDIR)/src  $(LOCALLIBS)

# All Target
all: main run

# Tool invocations
main: $(OBJS)
	@echo 'Building target: $@'
	$(CC) main.cpp -o output $(OBJS) $(LIBS) $(CFLAGS) $(BOOTLIBS)
	@echo 'Finished building target: $@'
	@echo ' '

.cpp.o:
	@echo 'Building file: $<'
	$(CC) -c $< -o $@ $(LIBS) $(CFLAGS) $(BOOTLIBS)
	@echo 'Finished building: $<'
	@echo ' '

# Other Targets
clean:
	$(RM) output $(OBJS)

run:
	./output