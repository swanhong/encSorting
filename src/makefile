RM = rm -rf
CC = g++
CFLAGS = -std=c++11 -O3

HEAAN_DIR = ../HEAAN
BOOT_DIR = ../bootsrc
PLAIN_DIR = plain
ENC_DIR = enc
TEST_DIR = test

BOOT_OBJS = $(BOOT_DIR)/common.o $(BOOT_DIR)/new_bootstrapping.o
OUR_OBJS = PrintUtils.o MaskingGenerator.o \
			$(PLAIN_DIR)/CyclicArray.o $(PLAIN_DIR)/PlainSorting.o \
			$(ENC_DIR)/EncAlgorithm.o $(ENC_DIR)/BootScheme.o $(ENC_DIR)/BootAlgorithm.o $(ENC_DIR)/EncSorting.o \
			$(TEST_DIR)/TestPlain.o $(TEST_DIR)/TestEnc.o $(TEST_DIR)/TestBoot.o $(TEST_DIR)/TestSort.o
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