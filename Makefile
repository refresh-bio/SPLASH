all: satc satc_dump satc_merge sig_anch download_kmc

SPLASH_LIBS_DIR = libs
LIBS_DIR = . #/usr/local/lib
INCLUDE_DIR= libs
OPENBLAS_INCLUDE_DIR = libs/openblas/include

SATC_MAIN_DIR=src/satc
SATC_MERGE_MAIN_DIR=src/satc_merge
SATC_DUMP_MAIN_DIR=src/satc_dump
SIG_ANCH_MAIN_DIR=src/sig_anch
COMMON_DIR=src/common

OUT_BIN_DIR=bin

CC 	= g++
CFLAGS	= -fPIC -Wall -O3 -m64 -std=c++17 -mavx -pthread -I $(INCLUDE_DIR) -I $(OPENBLAS_INCLUDE_DIR) -fpermissive
CLINK	= -lm -std=c++17 -lpthread -static-libstdc++ -lgfortran

release: CLINK	= -lm -std=c++17 -static -lgfortran -lquadmath -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
release: all

ifdef MSVC     # Avoid the MingW/Cygwin sections
    uname_S := Windows
else                          # If uname not available => 'not'
    uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
endif
ifeq ($(uname_S),Linux)
	CLINK+=-fabi-version=6
	LIB_ZLIB=cloudflare-zlib/libz.a
	LIB_ZSTD=zstd/linux/libzstd.a
	OPENBLAS_LIB_DIR=libs/openblas/lib/libopenblas.a
endif

ifeq ($(uname_S),Darwin)
	LIB_ZLIB=cloudflare-zlib/libz.mac.a
	LIB_ZSTD=zstd/mac/libzstd.a
	#OPENBLAS_LIB_DIR = #TODO
	CC = g++-11
endif

# default install location (binary placed in the /bin folder)
prefix      = /usr/local

# optional install location
exec_prefix = $(prefix)


%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@


satc: $(OUT_BIN_DIR)/satc

$(OUT_BIN_DIR)/satc: $(SATC_MAIN_DIR)/main.o \
	$(SATC_MAIN_DIR)/kmc_api/kmc_file.o \
	$(SATC_MAIN_DIR)/kmc_api/mmer.o \
	$(SATC_MAIN_DIR)/kmc_api/kmer_api.o
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) -o $@ $^ \
	$(SPLASH_LIBS_DIR)/$(LIB_ZSTD) \
	$(CLINK)

satc_merge: $(OUT_BIN_DIR)/satc_merge

$(OUT_BIN_DIR)/satc_merge: $(SATC_MERGE_MAIN_DIR)/main.o \
	$(SATC_MERGE_MAIN_DIR)/pvals.o \
	$(SATC_MERGE_MAIN_DIR)/extra_stats.o \
	$(SATC_MERGE_MAIN_DIR)/ob_utils.o \
	$(OPENBLAS_LIB_DIR)
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) -o $@ $^ \
	$(SPLASH_LIBS_DIR)/$(LIB_ZSTD) \
	$(CLINK)

satc_dump: $(OUT_BIN_DIR)/satc_dump

$(OUT_BIN_DIR)/satc_dump: $(SATC_DUMP_MAIN_DIR)/main.o
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) -o $@ $^ \
	$(SPLASH_LIBS_DIR)/$(LIB_ZSTD) \
	$(CLINK)

sig_anch: $(OUT_BIN_DIR)/sig_anch

$(OUT_BIN_DIR)/sig_anch: $(SIG_ANCH_MAIN_DIR)/sig_anch.o \
	$(COMMON_DIR)/csv.o
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) -o $@ $^ \
	$(CLINK)

download_kmc:
	-mkdir -p $(OUT_BIN_DIR)
	./download_kmc.sh $(OUT_BIN_DIR)

clean:
	-rm $(SATC_MAIN_DIR)/*.o
	-rm $(SATC_MAIN_DIR)/kmc_api/*.o	
	-rm $(SATC_MERGE_MAIN_DIR)/*.o
	-rm $(SATC_DUMP_MAIN_DIR)/*.o
	-rm $(COMMON_DIR)/*.o
	-rm $(SIG_ANCH_MAIN_DIR)/*.o
	-rm $(OUT_BIN_DIR)/satc	
	-rm $(OUT_BIN_DIR)/satc_merge
	-rm $(OUT_BIN_DIR)/sig_anch
	-rm $(OUT_BIN_DIR)/kmc
	-rm $(OUT_BIN_DIR)/kmc_tools
