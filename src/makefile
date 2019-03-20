RM = rm -rf
CC = g++
CFLAGS = -std=c++11 -O3

HEAAN_DIR = ../HEAAN
BOOT_DIR = ../bootsrc

BOOT_OBJS = $(BOOT_DIR)/common.o $(BOOT_DIR)/new_bootstrapping.o
OUR_OBJS = CyclicArray.o SortingAlgorithm.o EncContext.o BootContext.o TestAlgorithm.o TestBootstrapping.o
OBJS =  $(BOOT_OBJS) $(OUR_OBJS)

LOCAL_LIBS = -lntl -lgmp -lm -lpthread
BOOT_LIBS = ../eigen
LIBS = -L/usr/local/lib

# All Target
all: main run


# g++ main.cpp -o output -L/usr/local/lib -I ../eigen ../HEAAN/lib/libHEAAN.a -lntl -lgmp -lm -std=c++11 -O3 -I/../../src -pthread
# Tool invocations
main: $(OBJS)
	@echo 'Building target: $@'
	$(CC) main.cpp -o output $(OBJS) $(LIBS) $(CFLAGS) -I $(BOOT_LIBS) $(HEAAN_DIR)/lib/libHEAAN.a $(LOCAL_LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

.cpp.o:
	@echo 'Building file: $<'
	$(CC) -c $< -o $@ $(CFLAGS) -I $(BOOT_LIBS)
	@echo 'Finished building: $<'
	@echo ' '

# Other Targets
clean:
	$(RM) output $(OBJS)

run:
	./output