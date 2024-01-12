all: satc satc_dump satc_merge sig_anch compactors download_kmc splash supervised_test

SPLASH_LIBS_DIR = libs
LIBS_DIR = . #/usr/local/lib
INCLUDE_DIR= libs
MIMALLOC_INLUCDE_DIR = libs/mimalloc/include

SATC_MAIN_DIR=src/satc
SATC_MERGE_MAIN_DIR=src/satc_merge
SATC_DUMP_MAIN_DIR=src/satc_dump
SIG_ANCH_MAIN_DIR=src/sig_anch
COMMON_DIR=src/common

COMPACTORS_MAIN_DIR=src/compactors
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

ifeq ($(PREFIX),)
	PREFIX ?= /usr/local
endif

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

compactors: $(OUT_BIN_DIR)/compactors

$(OUT_BIN_DIR)/compactors: $(COMPACTORS_MAIN_DIR)/main.o \
	$(COMPACTORS_MAIN_DIR)/console.o \
	$(COMPACTORS_MAIN_DIR)/engine.o \
	$(COMPACTORS_MAIN_DIR)/io.o \
	$(COMPACTORS_MAIN_DIR)/main.o \
	$(COMPACTORS_MAIN_DIR)/read_select.o \
	$(COMMON_DIR)/edit_distance.o \
	$(COMMON_DIR)/illumina_adapters_static.o \
	$(SPLASH_LIBS_DIR)/cdflib/cdflib.o
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) -o $@ $^ \
	$(SPLASH_LIBS_DIR)/$(LIB_ZLIB) \
	$(CLINK) 
download_kmc:
	-mkdir -p $(OUT_BIN_DIR)
	./download_kmc.sh $(OUT_BIN_DIR)

splash:
	-mkdir -p $(OUT_BIN_DIR)
	cp src/splash.py bin/splash

supervised_test:
	cp src/supervised_test/supervised_test.R bin

install: all
	install -d $(PREFIX)/bin
	install bin/* $(PREFIX)/bin

uninstall:
	-rm $(PREFIX)/bin/satc
	-rm $(PREFIX)/bin/compactors
	-rm $(PREFIX)/bin/satc_dump
	-rm $(PREFIX)/bin/satc_merge
	-rm $(PREFIX)/bin/sig_anch
	-rm $(PREFIX)/bin/splash
	-rm $(PREFIX)/bin/kmc
	-rm $(PREFIX)/bin/kmc_tools
	-rm $(PREFIX)/bin/supervised_test.R

clean:
	-rm $(SATC_MAIN_DIR)/*.o
	-rm $(COMMON_DIR)/kmc_api/*.o
	-rm $(SATC_MERGE_MAIN_DIR)/*.o
	-rm $(SATC_DUMP_MAIN_DIR)/*.o
	-rm $(COMPACTORS_MAIN_DIR)/*.o
	-rm $(COMMON_DIR)/*.o
	-rm $(SIG_ANCH_MAIN_DIR)/*.o
	-rm -rf $(OUT_BIN_DIR)
