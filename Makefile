all: bkc satc satc_dump satc_undump satc_filter satc_to_fasta satc_merge sig_anch read_selector compactors download_kmc lookup_table splash supervised_test dsv_manip gap_shortener fafq_filter tsv_to_fasta


dummy := $(shell git submodule update --init --recursive)

SPLASH_LIBS_DIR = libs
LIBS_DIR = . #/usr/local/lib
INCLUDE_DIR= libs
ZLIB_INCLUDE_DIR = libs/zlib-ng/build-g++
ZLIB_INCLUDE_DIR_FOR_FILE_WRAPPER = libs/zlib-ng/build-g++/zlib-ng
ZSTD_INCLUDE_DIR = libs/zstd/lib
LIBDEFLATE_INCLUDE_DIR = libs/libdeflate
ISAL_INCLUDE_DIR = libs/isa-l/include
MIMALLOC_INLUCDE_DIR = libs/mimalloc/include

SATC_MAIN_DIR=src/satc
SATC_MERGE_MAIN_DIR=src/satc_merge
SATC_DUMP_MAIN_DIR=src/satc_dump
SATC_UNDUMP_MAIN_DIR=src/satc_undump
SATC_FILTER_MAIN_DIR=src/satc_filter
SATC_TO_FASTA_MAIN_DIR=src/satc_to_fasta
SIG_ANCH_MAIN_DIR=src/sig_anch
LOOKUP_TABLE_MAIN_DIR=src/lookup_table
DSV_MANIP_MAIN_DIR=src/dsv_manip
GAP_SHORTENER_MAIN_DIR=src/gap_shortener
FAFQ_FILTER_MAIN_DIR=src/fafq_filter
TSV_TO_FASTA_MAIN_DIR=src/tsv_to_fasta
COMMON_DIR=src/common

READ_SELECTOR_MAIN_DIR=src/read_selector
COMPACTORS_MAIN_DIR=src/compactors

OUT_BIN_DIR=bin


MIMALLOC_OBJ=libs/mimalloc/mimalloc.o

# Check if the selected C/C++ compilers are GNU gcc/g++
CXX_VERSION := $(shell $(CXX) --version 2>/dev/null | grep 'g++')
CC_VERSION := $(shell $(CC) --version 2>/dev/null | grep 'gcc')

#This is mainly for MAC OS, lets try to get some installeg g++
ifeq ($(CXX_VERSION),)
CXX := $(shell basename $$(which g++-14 || which g++-13 || which g++-12 || which g++-11 || which g++-10 || which g++-9 || echo "g++"))
CXX_VERSION := $(shell $(CXX) --version 2>/dev/null | grep 'g++')
endif
ifeq ($(CC_VERSION),)
CC := $(shell basename $$(which gcc-14 || which gcc-13 || which gcc-12 || which gcc-11 || which gcc-10 || which gcc-9 || echo "gcc"))
CC_VERSION := $(shell $(CC) --version 2>/dev/null | grep 'gcc')
endif

ifeq ($(CXX_VERSION),)
	WRONG_CXX_OR_CC = yes
endif
ifeq ($(CC_VERSION),)
	WRONG_CXX_OR_CC = yes
endif
ifdef WRONG_CXX_OR_CC
$(error The selected C++ compiler ($(CXX)) or C compiler ($(CC)) is not GNU g++/gcc. Please specify a GNU g++/gcc compiler. If you are using MACOS you may get one with brew (https://brew.sh/), and after installing run make CC=gcc-<version> CXX=g++-<version>)
endif

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

ifeq ($(D_ARCH),X64)
	dummy_install_nasm := $(shell \
	if [ ! -f build_tools/nasm/nasm ]; then \
		cd build_tools/nasm && ./autogen.sh && ./configure && make -j; \
	fi)
endif

NASM_V := $(shell build_tools/nasm/nasm --version 2>/dev/null)

CPU_FLAGS =
STATIC_LFLAGS =
PLATFORM_SPECIFIC_FLAGS =

#in some cases we can have different results on ARM
#I guess this is exactly the same as here: https://bugs.mysql.com/bug.php?id=82760
ifeq ($(D_ARCH),ARM64)
	PLATFORM_SPECIFIC_FLAGS = -ffp-contract=off
endif

ifeq ($(D_OS),MACOS)

	ifeq ($(D_ARCH),ARM64)
		CPU_FLAGS = -march=armv8.4-a
	else
		CPU_FLAGS = -m64
	endif
	STATIC_LFLAGS = -static-libgcc -static-libstdc++ -pthread
else

	ifeq ($(D_ARCH),ARM64)
		CPU_FLAGS = -march=armv8-a
		STATIC_LFLAGS = -static-libgcc -static-libstdc++ -lpthread
	else
		CPU_FLAGS = -m64
		STATIC_LFLAGS = -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
	endif
endif




LIB_ZLIB=$(SPLASH_LIBS_DIR)/zlib-ng/build-g++/zlib-ng/libz.a
LIB_ISAL=$(SPLASH_LIBS_DIR)/isa-l/bin/isa-l.a
LIB_ZSTD=$(SPLASH_LIBS_DIR)/zstd/lib/libzstd.a
LIB_LIBDEFLATE=$(SPLASH_LIBS_DIR)/libdeflate/build/libdeflate.a

LIB_GZ=$(LIB_ZLIB)

REFRESH_FLAGS = 
ifeq ($(UNAME_S),Linux)
	ifeq ($(UNAME_M),x86_64)
		REFRESH_FLAGS +=-DARCH_X64
		ifdef NASM_V
			LIB_GZ=$(LIB_ISAL)
			REFRESH_FLAGS += -DREFRESH_USE_IGZIP -I $(ISAL_INCLUDE_DIR)
		else
			REFRESH_FLAGS +=-DREFRESH_USE_ZLIB
		endif
	else
		REFRESH_FLAGS +=-DREFRESH_USE_ZLIB
	endif
else
	REFRESH_FLAGS +=-DREFRESH_USE_ZLIB
endif

CFLAGS	= -fPIC -Wall -O3 $(PLATFORM_SPECIFIC_FLAGS) $(CPU_FLAGS) $(REFRESH_FLAGS) -std=c++17 -pthread -I $(ZLIB_INCLUDE_DIR) -I $(ZLIB_INCLUDE_DIR_FOR_FILE_WRAPPER) -I $(INCLUDE_DIR) -I $(MIMALLOC_INLUCDE_DIR) -I $(ZSTD_INCLUDE_DIR) -I $(LIBDEFLATE_INCLUDE_DIR) -fpermissive
CLINK	= -lm -std=c++17 -lpthread

release: CLINK = -lm -std=c++17 $(STATIC_LFLAGS)

release: CFLAGS	= -fPIC -Wall -O3 -DNDEBUG $(PLATFORM_SPECIFIC_FLAGS) $(CPU_FLAGS) $(REFRESH_FLAGS) -std=c++17 -pthread -I $(ZLIB_INCLUDE_DIR) -I $(ZLIB_INCLUDE_DIR_FOR_FILE_WRAPPER) -I $(INCLUDE_DIR) -I $(MIMALLOC_INLUCDE_DIR) -I $(ZSTD_INCLUDE_DIR) -I $(LIBDEFLATE_INCLUDE_DIR) -fpermissive
#release: satc_undump satc_filter satc_to_fasta fafq_filter
release: all

debug: CFLAGS	= -fPIC -Wall -O0 -g $(PLATFORM_SPECIFIC_FLAGS) $(CPU_FLAGS) $(REFRESH_FLAGS) -std=c++17 -pthread -I $(ZLIB_INCLUDE_DIR) -I $(ZLIB_INCLUDE_DIR_FOR_FILE_WRAPPER) -I $(INCLUDE_DIR) -I $(MIMALLOC_INLUCDE_DIR) -I $(ZSTD_INCLUDE_DIR) -I $(LIBDEFLATE_INCLUDE_DIR) -fpermissive
debug: all

ifeq ($(UNAME_S),Linux)
	CLINK+=-fabi-version=6
endif


# default install location (binary placed in the /bin folder)
prefix      = /usr/local

# optional install location
exec_prefix = $(prefix)

$(LIB_ZLIB):
	cd $(SPLASH_LIBS_DIR)/zlib-ng; cmake -DCMAKE_CXX_COMPILER=$(CXX) -DCMAKE_C_COMPILER=$(CC) -B build-g++/zlib-ng -S . -DZLIB_COMPAT=ON; cmake --build build-g++/zlib-ng --config Release

$(LIB_ISAL):
	cd $(SPLASH_LIBS_DIR)/isa-l && PATH=../../build_tools/nasm:$$PATH make -f Makefile.unx

$(LIB_ZSTD):
	cd $(SPLASH_LIBS_DIR)/zstd; make -j
	
$(LIB_LIBDEFLATE):
	cd $(SPLASH_LIBS_DIR)/libdeflate; cmake -B build && cmake --build build

$(MIMALLOC_OBJ):
	$(CXX) -DMI_MALLOC_OVERRIDE -O3 -DNDEBUG -fPIC -Wall -Wextra -Wno-unknown-pragmas -fvisibility=hidden -ftls-model=initial-exec -fno-builtin-malloc -c -I libs/mimalloc/include libs/mimalloc/src/static.c -o $(MIMALLOC_OBJ)

src/lookup_table/%.o: src/lookup_table/%.cpp $(SPLASH_LIBS_DIR)/SBWT/build/external/sdsl-lite/build/lib/libsdsl.a
	$(CXX) $(CFLAGS) -I libs/SBWT/include -I libs/SBWT/sdsl-lite/include -I libs/SBWT/SeqIO/include -I libs/SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -c $< -o $@

%.o: %.cpp $(LIB_ZLIB) $(LIB_GZ)
	$(CXX) $(CFLAGS) -c $< -o $@

bkc: $(OUT_BIN_DIR)/bkc

$(OUT_BIN_DIR)/bkc:
	(cd $(SPLASH_LIBS_DIR)/bkc; make CC=$(CC) CXX=$(CXX) release SHARED_INCLUDE_DIR=../../src/common BKC_OUT_BIN_DIR=../../bin -j)


satc: $(OUT_BIN_DIR)/satc

$(OUT_BIN_DIR)/satc: $(SATC_MAIN_DIR)/satc.o \
	$(COMMON_DIR)/kmc_api/kmc_file.o \
	$(COMMON_DIR)/kmc_api/mmer.o \
	$(COMMON_DIR)/kmc_api/kmer_api.o \
	$(LIB_ZSTD)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ $^ \
	$(LIB_ZSTD) \
	$(CLINK)

satc_merge: $(OUT_BIN_DIR)/satc_merge

$(OUT_BIN_DIR)/satc_merge: $(SATC_MERGE_MAIN_DIR)/satc_merge.o \
	$(SATC_MERGE_MAIN_DIR)/pvals.o \
	$(SATC_MERGE_MAIN_DIR)/anchor.o \
	$(SATC_MERGE_MAIN_DIR)/extra_stats.o \
	$(LIB_LIBDEFLATE) \
	$(LIB_ZSTD)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ $^ \
	$(LIB_ZSTD) \
	$(LIB_LIBDEFLATE) \
	$(CLINK)

satc_dump: $(OUT_BIN_DIR)/satc_dump

$(OUT_BIN_DIR)/satc_dump: $(SATC_DUMP_MAIN_DIR)/satc_dump.o $(LIB_ZSTD)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ $^ \
	$(LIB_LIBDEFLATE) \
	$(LIB_ZSTD) \
	$(CLINK)

satc_undump: $(OUT_BIN_DIR)/satc_undump

$(OUT_BIN_DIR)/satc_undump: $(SATC_UNDUMP_MAIN_DIR)/satc_undump.o $(LIB_ZSTD)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ $^ \
	$(LIB_ZSTD) \
	$(CLINK)

satc_filter: $(OUT_BIN_DIR)/satc_filter

$(OUT_BIN_DIR)/satc_filter: $(SATC_FILTER_MAIN_DIR)/satc_filter.o $(LIB_ZSTD)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ $^ \
	$(LIB_ZSTD) \
	$(CLINK)

satc_to_fasta: $(OUT_BIN_DIR)/satc_to_fasta

$(OUT_BIN_DIR)/satc_to_fasta: $(SATC_TO_FASTA_MAIN_DIR)/satc_to_fasta.o $(LIB_ZSTD)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ $^ \
	$(LIB_ZSTD) \
	$(CLINK)

tsv_to_fasta: $(OUT_BIN_DIR)/tsv_to_fasta

$(OUT_BIN_DIR)/tsv_to_fasta: $(TSV_TO_FASTA_MAIN_DIR)/tsv_to_fasta.o
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ $^

fafq_filter: $(OUT_BIN_DIR)/fafq_filter

$(OUT_BIN_DIR)/fafq_filter: $(FAFQ_FILTER_MAIN_DIR)/fafq_filter.o \
	$(FAFQ_FILTER_MAIN_DIR)/app.o \
	$(FAFQ_FILTER_MAIN_DIR)/worker.o \
	$(LIB_GZ) \
	$(LIB_ZSTD)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ $^ \
	$(LIB_GZ) \
	$(LIB_ZSTD) \
	$(CLINK)

sig_anch: $(OUT_BIN_DIR)/sig_anch

$(OUT_BIN_DIR)/sig_anch: $(SIG_ANCH_MAIN_DIR)/sig_anch.o \
	$(COMMON_DIR)/csv.o
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ $^ \
	$(CLINK)

read_selector: $(OUT_BIN_DIR)/read_selector

$(OUT_BIN_DIR)/read_selector: $(READ_SELECTOR_MAIN_DIR)/read_selector.o \
	$(LIB_ZLIB)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ $^ \
	$(LIB_ZLIB) \
	$(CLINK) 

compactors: $(OUT_BIN_DIR)/compactors

$(OUT_BIN_DIR)/compactors: $(COMPACTORS_MAIN_DIR)/main.o \
	$(COMPACTORS_MAIN_DIR)/console.o \
	$(COMPACTORS_MAIN_DIR)/engine.o \
	$(COMPACTORS_MAIN_DIR)/io.o \
	$(COMPACTORS_MAIN_DIR)/main.o \
	$(COMPACTORS_MAIN_DIR)/read_select.o \
	$(COMMON_DIR)/edit_distance.o \
	$(LIB_ZLIB) \
	$(SPLASH_LIBS_DIR)/cdflib/cdflib.o

	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ $^ \
	$(LIB_ZLIB) \
	$(CLINK) 

$(SPLASH_LIBS_DIR)/SBWT/build/libsbwt_static.a:
	(cd $(SPLASH_LIBS_DIR)/SBWT/build; cmake -DCMAKE_CXX_COMPILER=$(CXX) -DCMAKE_C_COMPILER=$(CC) .. -DMAX_KMER_LENGTH=32; make -j)

$(SPLASH_LIBS_DIR)/SBWT/build/external/sdsl-lite/build/lib/libsdsl.a: $(SPLASH_LIBS_DIR)/SBWT/build/libsbwt_static.a
$(SPLASH_LIBS_DIR)/SBWT/build/external/KMC/build/libkmc_core.a: $(SPLASH_LIBS_DIR)/SBWT/build/libsbwt_static.a
$(SPLASH_LIBS_DIR)/SBWT/build/external/KMC/build/libkmc_tools.a: $(SPLASH_LIBS_DIR)/SBWT/build/libsbwt_static.a



lookup_table: $(OUT_BIN_DIR)/lookup_table splash
	cp $(LOOKUP_TABLE_MAIN_DIR)/build_lookup_table.py bin/build_lookup_table.py

$(OUT_BIN_DIR)/lookup_table: \
	$(MIMALLOC_OBJ) \
	$(LOOKUP_TABLE_MAIN_DIR)/lookup_table.o \
	$(LOOKUP_TABLE_MAIN_DIR)/lookup.o \
	$(LOOKUP_TABLE_MAIN_DIR)/kmer_index_sbwt.o \
	$(LOOKUP_TABLE_MAIN_DIR)/build.o \
	$(LOOKUP_TABLE_MAIN_DIR)/query.o \
	$(LOOKUP_TABLE_MAIN_DIR)/query_many.o \
	$(LOOKUP_TABLE_MAIN_DIR)/query_common.o \
	$(COMMON_DIR)/kmc_api/kmc_file.o \
	$(COMMON_DIR)/kmc_api/mmer.o \
	$(COMMON_DIR)/kmc_api/kmer_api.o \
	$(SPLASH_LIBS_DIR)/refresh/parallel_queues/lib/parallel-queues-monitor.o \
	$(SPLASH_LIBS_DIR)/SBWT/build/libsbwt_static.a \
	$(SPLASH_LIBS_DIR)/SBWT/build/external/sdsl-lite/build/lib/libsdsl.a \
	$(SPLASH_LIBS_DIR)/SBWT/build/external/KMC/build/libkmc_core.a \
	$(SPLASH_LIBS_DIR)/SBWT/build/external/KMC/build/libkmc_tools.a \
	$(LIB_ZLIB)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ $^ \
	$(CLINK)

dsv_manip: $(OUT_BIN_DIR)/dsv_manip

$(OUT_BIN_DIR)/dsv_manip: $(DSV_MANIP_MAIN_DIR)/dsv_manip.o \
	$(DSV_MANIP_MAIN_DIR)/dsv_common.o \
	$(DSV_MANIP_MAIN_DIR)/limit_mode.o \
	$(DSV_MANIP_MAIN_DIR)/sort_mode.o
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ $^ \
	$(CLINK)

gap_shortener: $(OUT_BIN_DIR)/gap_shortener

$(OUT_BIN_DIR)/gap_shortener: $(GAP_SHORTENER_MAIN_DIR)/gap_shortener.o \
	$(GAP_SHORTENER_MAIN_DIR)/reader.o \
	$(GAP_SHORTENER_MAIN_DIR)/pseudoreads_maker.o \
	$(LIB_ZLIB)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ $^ \
	$(CLINK)

download_kmc:
	-mkdir -p $(OUT_BIN_DIR)
	./download_kmc.py $(OUT_BIN_DIR)

splash:
	mkdir -p bin
	cp src/splash.py bin/splash
	cp -R postprocessing bin/

supervised_test:
	cp src/supervised_test/supervised_test.R bin

install: all
	cp -r bin/* /usr/local/bin

uninstall:
	-rm /usr/local/bin/satc
	-rm /usr/local/bin/bkc
	-rm /usr/local/bin/bkc_dump
	-rm /usr/local/bin/read_selector
	-rm /usr/local/bin/compactors
	-rm /usr/local/bin/satc_dump
	-rm /usr/local/bin/satc_merge
	-rm /usr/local/bin/sig_anch
	-rm /usr/local/bin/splash
	-rm /usr/local/bin/kmc
	-rm /usr/local/bin/kmc_tools
clean:
	-rm $(SPLASH_LIBS_DIR)/refresh/*.o
	-rm $(SATC_MAIN_DIR)/*.o
	-rm $(COMMON_DIR)/kmc_api/*.o
	-rm $(SATC_MERGE_MAIN_DIR)/*.o
	-rm $(SATC_DUMP_MAIN_DIR)/*.o
	-rm $(FAFQ_FILTER_MAIN_DIR)/*.o
	-rm $(SATC_FILTER_MAIN_DIR)/*.o
	-rm $(READ_SELECTOR_MAIN_DIR)/*.o
	-rm $(COMPACTORS_MAIN_DIR)/*.o
	-rm $(COMMON_DIR)/*.o
	-rm $(SIG_ANCH_MAIN_DIR)/*.o
	-rm $(LOOKUP_TABLE_MAIN_DIR)/*.o
	-rm $(DSV_MANIP_MAIN_DIR)/*.o
	-rm $(GAP_SHORTENER_MAIN_DIR)/*.o
	-rm -rf $(OUT_BIN_DIR)
	-rm $(SPLASH_LIBS_DIR)/cdflib/*.o
	(cd $(SPLASH_LIBS_DIR)/bkc; make clean)
