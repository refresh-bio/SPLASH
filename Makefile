all: satc satc_dump satc_merge sig_anch compactors download_kmc splash supervised_test dsv_manip

dummy := $(shell git submodule update --init --recursive)

SPLASH_LIBS_DIR = libs
LIBS_DIR = . #/usr/local/lib
INCLUDE_DIR= libs
ZLIB_INCLUDE_DIR= libs/cloudflare-zlib

SATC_MAIN_DIR=src/satc
SATC_MERGE_MAIN_DIR=src/satc_merge
SATC_DUMP_MAIN_DIR=src/satc_dump
SIG_ANCH_MAIN_DIR=src/sig_anch
DSV_MANIP_MAIN_DIR=src/dsv_manip
COMMON_DIR=src/common
COMPACTORS_MAIN_DIR=src/compactors

OUT_BIN_DIR=bin



ifdef MSVC     # Avoid the MingW/Cygwin sections
    UNAME_S := Windows
else                          # If uname not available => 'not'
    UNAME_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
    UNAME_M := $(shell uname -m)
endif

D_OS =
D_ARCH =

ifeq ($(UNAME_S),Darwin)
	D_OS=MACOS
	ifeq ($(UNAME_M),arm64)
		D_ARCH=ARM64
	else
		D_ARCH=X64
	endif
else
	D_OS=LINUX
	D_ARCH=X64
	ifeq ($(UNAME_M),arm64)
		D_ARCH=ARM64
	endif
	ifeq ($(UNAME_M),aarch64)
		D_ARCH=ARM64
	endif
endif

CPU_FLAGS =
STATIC_LFLAGS =
PLATFORM_SPECIFIC_FLAGS =

#in some cases we can have different results on ARM
#I guess this is exactly the same as here: https://bugs.mysql.com/bug.php?id=82760
ifeq ($(D_ARCH),ARM64)
	PLATFORM_SPECIFIC_FLAGS = -ffp-contract=off
endif

ifeq ($(D_OS),MACOS)
	CC = g++-11

	ifeq ($(D_ARCH),ARM64)
		CPU_FLAGS = -march=armv8.4-a
	else
		CPU_FLAGS = -m64
	endif
	STATIC_LFLAGS = -static-libgcc -static-libstdc++ -pthread
else
	CC 	= g++

	ifeq ($(D_ARCH),ARM64)
		CPU_FLAGS = -march=armv8-a
		STATIC_LFLAGS = -static-libgcc -static-libstdc++ -lpthread
	else
		CPU_FLAGS = -m64
		STATIC_LFLAGS = -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
	endif
endif


CFLAGS	= -fPIC -Wall -O3 $(PLATFORM_SPECIFIC_FLAGS) $(CPU_FLAGS) -std=c++17 -pthread -I $(INCLUDE_DIR) -I $(ZLIB_INCLUDE_DIR) -fpermissive
CLINK	= -lm -std=c++17 -lpthread

release: CLINK = -lm -std=c++17 $(STATIC_LFLAGS)
release: CLINK = -lm -std=c++17 $(STATIC_LFLAGS)

release: CFLAGS	= -fPIC -Wall -O3 -DNDEBUG $(PLATFORM_SPECIFIC_FLAGS) $(CPU_FLAGS) -std=c++17 -pthread -I $(INCLUDE_DIR) -I $(ZLIB_INCLUDE_DIR) -fpermissive
release: all

debug: CFLAGS	= -fPIC -Wall -O0 -g $(PLATFORM_SPECIFIC_FLAGS) $(CPU_FLAGS) -std=c++17 -pthread -I $(INCLUDE_DIR) -I $(ZLIB_INCLUDE_DIR) -fpermissive
debug: all

ifeq ($(UNAME_S),Linux)
	CLINK+=-fabi-version=6
endif


LIB_ZLIB=$(SPLASH_LIBS_DIR)/cloudflare-zlib/libz.a
LIB_ZSTD=$(SPLASH_LIBS_DIR)/zstd/lib/libzstd.a

# default install location (binary placed in the /bin folder)

ifeq ($(PREFIX),)
	PREFIX ?= /usr/local
endif

$(LIB_ZLIB):
	cd $(SPLASH_LIBS_DIR)/cloudflare-zlib; ./configure; make libz.a

$(LIB_ZSTD):
	cd $(SPLASH_LIBS_DIR)/zstd; make -j

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

satc: $(OUT_BIN_DIR)/satc

$(OUT_BIN_DIR)/satc: $(SATC_MAIN_DIR)/satc.o \
	$(COMMON_DIR)/kmc_api/kmc_file.o \
	$(COMMON_DIR)/kmc_api/mmer.o \
	$(COMMON_DIR)/kmc_api/kmer_api.o \
	$(COMMON_DIR)/illumina_adapters_static.o \
	$(LIB_ZSTD)
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) -o $@ $^ \
	$(LIB_ZSTD) \
	$(CLINK)

satc_merge: $(OUT_BIN_DIR)/satc_merge

$(OUT_BIN_DIR)/satc_merge: $(SATC_MERGE_MAIN_DIR)/satc_merge.o \
	$(SATC_MERGE_MAIN_DIR)/pvals.o \
	$(SATC_MERGE_MAIN_DIR)/anchor.o \
	$(SATC_MERGE_MAIN_DIR)/extra_stats.o \
	$(LIB_ZSTD)
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) -o $@ $^ \
	$(LIB_ZSTD) \
	$(CLINK)

satc_dump: $(OUT_BIN_DIR)/satc_dump

$(OUT_BIN_DIR)/satc_dump: $(SATC_DUMP_MAIN_DIR)/satc_dump.o $(LIB_ZSTD)
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) -o $@ $^ \
	$(LIB_ZSTD) \
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
	$(LIB_ZLIB) \
	$(SPLASH_LIBS_DIR)/cdflib/cdflib.o

	-mkdir -p $(OUT_BIN_DIR)
	$(CC) -o $@ $^ \
	$(LIB_ZLIB) \
	$(CLINK) 

dsv_manip: $(OUT_BIN_DIR)/dsv_manip

$(OUT_BIN_DIR)/dsv_manip: $(DSV_MANIP_MAIN_DIR)/dsv_manip.o \
	$(DSV_MANIP_MAIN_DIR)/dsv_common.o \
	$(DSV_MANIP_MAIN_DIR)/limit_mode.o \
	$(DSV_MANIP_MAIN_DIR)/sort_mode.o
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) -o $@ $^ \
	$(CLINK)
download_kmc:
	-mkdir -p $(OUT_BIN_DIR)
	./download_kmc.py $(OUT_BIN_DIR)

splash:
	-mkdir -p $(OUT_BIN_DIR)
	cp src/splash.py bin/splash

supervised_test:
	-mkdir -p $(OUT_BIN_DIR)
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
	-rm $(PREFIX)/bin/dsv_manip
	-rm $(PREFIX)/bin/supervised_test.R

clean:
	-rm $(SATC_MAIN_DIR)/*.o
	-rm $(COMMON_DIR)/kmc_api/*.o
	-rm $(SATC_MERGE_MAIN_DIR)/*.o
	-rm $(SATC_DUMP_MAIN_DIR)/*.o
	-rm $(COMPACTORS_MAIN_DIR)/*.o
	-rm $(COMMON_DIR)/*.o
	-rm $(SIG_ANCH_MAIN_DIR)/*.o
	-rm $(DSV_MANIP_MAIN_DIR)/*.o
	-rm -rf $(OUT_BIN_DIR)
	-rm $(SPLASH_LIBS_DIR)/cdflib/*.o
