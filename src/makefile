RM = rm -rf
CC = g++
CFLAGS = -std=c++11 -O3
TERGET = main.cpp
OUTPUT = output

HEAAN_DIR = ../HEAAN/

DEFAULT_FILES = PrintUtils MaskingGenerator

BOOT_DIR = ../bootsrc/
BOOT_FILES = common new_bootstrapping

PLAIN_DIR = plain/
PLAIN_FILES = CyclicArray PlainSorting

ENC_DIR = enc/
ENC_FILES = EncAlgorithm BootScheme BootAlgorithm EncSorting

TEST_DIR = test/
TEST_FILES = TestPlain TestEnc TestBoot TestSort

BOOT_OBJS = $(addprefix $(BOOT_DIR), $(addsuffix .o, $(BOOT_FILES)))
OUR_OBJS = $(addsuffix .o, $(DEFAULT_FILES)) \
			$(addprefix $(PLAIN_DIR), $(addsuffix .o, $(PLAIN_FILES))) \
			$(addprefix $(ENC_DIR), $(addsuffix .o, $(ENC_FILES))) \
			$(addprefix $(TEST_DIR), $(addsuffix .o, $(TEST_FILES)))
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
	$(CC) $(TERGET) -o $(OUTPUT) $(OBJS) $(LIBS) $(CFLAGS) -I $(BOOT_LIBS) $(HEAAN_DIR)lib/libHEAAN.a $(LOCAL_LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

.cpp.o:
	@echo 'Building file: $<'
	$(CC) -c $< -o $@ $(CFLAGS) -I $(BOOT_LIBS)
	@echo 'Finished building: $<'
	@echo ' '

# Other Targets
clean:
	$(RM) $(OUTPUT) $(OBJS)

run:
	./$(OUTPUT)