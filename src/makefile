RM = rm -rf
CC = g++
CFLAGS = -std=c++11 -O3
HEAANDIR = ../HEAAN
OBJS = CyclicArray.o SortingAlgorithm.o
LOCALLIBS = -lntl -lgmp -lm -lpthread
BOOTLIBS = ../eigen
LIBS = -L/usr/local/lib -I $(BOOTLIBS) $(HEAANDIR)/lib/libHEAAN.a $(LOCALLIBS)

# All Target
all: main run


# g++ main.cpp -o output -L/usr/local/lib -I ../eigen ../HEAAN/lib/libHEAAN.a -lntl -lgmp -lm -std=c++11 -O3 -I/../../src -pthread
# Tool invocations
main: $(OBJS)
	@echo 'Building target: $@'
	$(CC) main.cpp -o output $(OBJS) $(LIBS) $(CFLAGS)	
	@echo 'Finished building target: $@'
	@echo ' '

.cpp.o:
	@echo 'Building file: $<'
	$(CC) -c $< -o $@ $(LIBS) $(CFLAGS)
	@echo 'Finished building: $<'
	@echo ' '

# Other Targets
clean:
	$(RM) output $(OBJS)

run:
	./output