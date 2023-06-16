all: satc satc_dump satc_merge sig_anch download_kmc splash

SPLASH_LIBS_DIR = libs
LIBS_DIR = . #/usr/local/lib
INCLUDE_DIR= libs
MIMALLOC_INLUCDE_DIR = libs/mimalloc/include

SATC_MAIN_DIR=src/satc
SATC_MERGE_MAIN_DIR=src/satc_merge
SATC_DUMP_MAIN_DIR=src/satc_dump
SIG_ANCH_MAIN_DIR=src/sig_anch
COMMON_DIR=src/common

OUT_BIN_DIR=bin

CC 	= g++
CFLAGS	= -fPIC -Wall -O3 -m64 -std=c++17 -pthread -I $(INCLUDE_DIR) -I $(MIMALLOC_INLUCDE_DIR) -fpermissive
CLINK	= -lm -std=c++17 -lpthread -static-libstdc++

MIMALLOC_OBJ=libs/mimalloc/mimalloc.o

release: CLINK	= -lm -std=c++17 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
release: CFLAGS	= -fPIC -Wall -O3 -DNDEBUG -m64 -std=c++17 -pthread -I $(INCLUDE_DIR) -I $(MIMALLOC_INLUCDE_DIR) -fpermissive
release: all

debug: CFLAGS	= -fPIC -Wall -O0 -g -m64 -std=c++17 -pthread -I $(INCLUDE_DIR) -I $(MIMALLOC_INLUCDE_DIR) -fpermissive
debug: all

ifdef MSVC     # Avoid the MingW/Cygwin sections
    uname_S := Windows
else                          # If uname not available => 'not'
    uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
endif
ifeq ($(uname_S),Linux)
	CLINK+=-fabi-version=6
	LIB_ZLIB=cloudflare-zlib/libz.a
	LIB_ZSTD=zstd/linux/libzstd.a
endif

ifeq ($(uname_S),Darwin)
	LIB_ZLIB=cloudflare-zlib/libz.mac.a
	LIB_ZSTD=zstd/mac/libzstd.a
	CC = g++-11
endif

# default install location (binary placed in the /bin folder)
prefix      = /usr/local

# optional install location
exec_prefix = $(prefix)

$(MIMALLOC_OBJ):
	$(CC) -DMI_MALLOC_OVERRIDE -O3 -DNDEBUG -fPIC -Wall -Wextra -Wno-unknown-pragmas -fvisibility=hidden -ftls-model=initial-exec -fno-builtin-malloc -c -I libs/mimalloc/include libs/mimalloc/src/static.c -o $(MIMALLOC_OBJ)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

satc: $(OUT_BIN_DIR)/satc

$(OUT_BIN_DIR)/satc: $(SATC_MAIN_DIR)/satc.o \
	$(COMMON_DIR)/kmc_api/kmc_file.o \
	$(COMMON_DIR)/kmc_api/mmer.o \
	$(COMMON_DIR)/kmc_api/kmer_api.o \
	$(COMMON_DIR)/illumina_adapters_static.o
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) -o $@ $^ \
	$(SPLASH_LIBS_DIR)/$(LIB_ZSTD) \
	$(CLINK)

satc_merge: $(OUT_BIN_DIR)/satc_merge

$(OUT_BIN_DIR)/satc_merge: $(SATC_MERGE_MAIN_DIR)/satc_merge.o \
	$(SATC_MERGE_MAIN_DIR)/pvals.o \
	$(SATC_MERGE_MAIN_DIR)/anchor.o \
	$(SATC_MERGE_MAIN_DIR)/extra_stats.o
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) -o $@ $^ \
	$(SPLASH_LIBS_DIR)/$(LIB_ZSTD) \
	$(CLINK)

satc_dump: $(OUT_BIN_DIR)/satc_dump

$(OUT_BIN_DIR)/satc_dump: $(SATC_DUMP_MAIN_DIR)/satc_dump.o
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

splash:
	cp src/splash.py bin/splash

install: all
	install bin/* /usr/local/bin

uninstall:
	-rm /usr/local/bin/satc
	-rm /usr/local/bin/satc_dump
	-rm /usr/local/bin/satc_merge
	-rm /usr/local/bin/sig_anch
	-rm /usr/local/bin/splash
	-rm /usr/local/bin/kmc
	-rm /usr/local/bin/kmc_tools

clean:
	-rm $(SATC_MAIN_DIR)/*.o
	-rm $(COMMON_DIR)/kmc_api/*.o
	-rm $(SATC_MERGE_MAIN_DIR)/*.o
	-rm $(SATC_DUMP_MAIN_DIR)/*.o
	-rm $(COMMON_DIR)/*.o
	-rm $(SIG_ANCH_MAIN_DIR)/*.o
	-rm -rf $(OUT_BIN_DIR)
