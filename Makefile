all: bkc satc satc_dump satc_undump satc_filter satc_to_fasta satc_merge sig_anch read_selector compactors download_kmc lookup_table splash supervised_test dsv_manip gap_shortener fafq_filter tsv_to_fasta

# splash is quite sensitive to the order of operations, which may be different on different archs due to how SIDM is handled
# to keep results the same across various platforms by default, we set generic platform
PLATFORM?=generic

# *** REFRESH makefile utils
include refresh.mk

$(call INIT_SUBMODULES_FAST)
$(call INIT_GLOBALS)
$(call CHECK_OS_ARCH, $(PLATFORM))

# *** Project directories
$(call SET_SRC_OBJ_BIN,src,obj,bin)
3RD_PARTY_DIR := ./libs
BUILD_TOOLS := ./build_tools

# *** Project configuration
$(call ADD_NASM,$(BUILD_TOOLS)/nasm)
$(call CHECK_NASM)
$(call ADD_MIMALLOC,$(3RD_PARTY_DIR)/mimalloc)
$(call ADD_CDFLIB, $(3RD_PARTY_DIR)/cdflib)
$(call ADD_SBWT, $(3RD_PARTY_DIR)/SBWT)					# For some reason must be before zlib
$(call PROPOSE_ISAL, $(3RD_PARTY_DIR)/isa-l)
$(call ADD_ZLIB_NG, $(3RD_PARTY_DIR)/zlib-ng)
$(call CHOOSE_GZIP_DECOMPRESSION)
$(call ADD_LIBDEFLATE, $(3RD_PARTY_DIR)/libdeflate)
$(call ADD_LIBZSTD, $(3RD_PARTY_DIR)/zstd)
$(call SET_STATIC, $(STATIC_LINK))

$(call SET_C_CPP_STANDARDS, c11, c++17)
$(call ADD_REFRESH_LIB, $(3RD_PARTY_DIR))
$(call ADD_REFRESH_PARALLEL_QUEUES_MONITOR, $(3RD_PARTY_DIR))

$(call SET_GIT_COMMIT)

$(call SET_FLAGS, $(TYPE))

$(call SET_COMPILER_VERSION_ALLOWED, GCC, Linux_x86_64, 10, 20)
$(call SET_COMPILER_VERSION_ALLOWED, GCC, Linux_aarch64, 11, 20)
$(call SET_COMPILER_VERSION_ALLOWED, GCC, Darwin_x86_64, 11, 13)
$(call SET_COMPILER_VERSION_ALLOWED, GCC, Darwin_arm64, 11, 13)

ifneq ($(MAKECMDGOALS),clean)
$(call CHECK_COMPILER_VERSION)
endif

# *** Source files and rules
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,SATC_MAIN,satc))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,SATC_MERGE_MAIN,satc_merge))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,SATC_DUMP_MAIN,satc_dump))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,SATC_UNDUMP_MAIN,satc_undump))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,SATC_FILTER_MAIN,satc_filter))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,SATC_TO_FASTA_MAIN,satc_to_fasta))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,SIG_ANCH_MAIN,sig_anch))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,LOOKUP_TABLE_MAIN,lookup_table))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,DSV_MANIP_MAIN,dsv_manip))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,GAP_SHORTENER_MAIN,gap_shortener))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,FAFQ_FILTER_MAIN,fafq_filter))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,TSV_TO_FASTA_MAIN,tsv_to_fasta))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,COMMON,common))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,READ_SELECTOR_MAIN,read_selector))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,COMPACTORS_MAIN,compactors))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,COMMON_FILTERS,common/filters))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,COMMON_KMC_API,common/kmc_api))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,COMMON_TYPES,common/types))

# *** Targets
splash: $(OUT_BIN_DIR)/splash
$(OUT_BIN_DIR)/splash: $(SRC_DIR)/splash.py
	mkdir -p $(OUT_BIN_DIR)
	cp $(SRC_DIR)/splash.py $(OUT_BIN_DIR)/splash
	cp -R postprocessing $(OUT_BIN_DIR)/

satc: $(OUT_BIN_DIR)/satc
$(OUT_BIN_DIR)/satc: \
	$(OBJ_SATC_MAIN) $(OBJ_COMMON_KMC_API)
	-mkdir -p $(OUT_BIN_DIR)	
	$(CXX) -o $@  \
	$(MIMALLOC_OBJ) \
	$(OBJ_SATC_MAIN) $(OBJ_COMMON_KMC_API) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

satc_merge: $(OUT_BIN_DIR)/satc_merge
$(OUT_BIN_DIR)/satc_merge: \
	$(OBJ_SATC_MERGE_MAIN)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ \
	$(OBJ_SATC_MERGE_MAIN) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

satc_dump: $(OUT_BIN_DIR)/satc_dump
$(OUT_BIN_DIR)/satc_dump: \
	$(OBJ_SATC_DUMP_MAIN)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ \
	$(OBJ_SATC_DUMP_MAIN) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

satc_undump: $(OUT_BIN_DIR)/satc_undump
$(OUT_BIN_DIR)/satc_undump: \
	$(OBJ_SATC_UNDUMP_MAIN)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ \
	$(OBJ_SATC_UNDUMP_MAIN) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

satc_filter: $(OUT_BIN_DIR)/satc_filter
$(OUT_BIN_DIR)/satc_filter: \
	$(OBJ_SATC_FILTER_MAIN)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ \
	$(OBJ_SATC_FILTER_MAIN) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

satc_to_fasta: $(OUT_BIN_DIR)/satc_to_fasta
$(OUT_BIN_DIR)/satc_to_fasta: \
	$(OBJ_SATC_TO_FASTA_MAIN)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ \
	$(OBJ_SATC_TO_FASTA_MAIN) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

tsv_to_fasta: $(OUT_BIN_DIR)/tsv_to_fasta
$(OUT_BIN_DIR)/tsv_to_fasta: \
	$(OBJ_TSV_TO_FASTA_MAIN)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ \
	$(OBJ_TSV_TO_FASTA_MAIN) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

fafq_filter: $(OUT_BIN_DIR)/fafq_filter
$(OUT_BIN_DIR)/fafq_filter: \
	$(OBJ_FAFQ_FILTER_MAIN) 
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ \
	$(OBJ_FAFQ_FILTER_MAIN) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

sig_anch: $(OUT_BIN_DIR)/sig_anch
$(OUT_BIN_DIR)/sig_anch: \
	$(OBJ_SIG_ANCH_MAIN) $(OBJ_COMMON)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ \
	$(OBJ_SIG_ANCH_MAIN) $(OBJ_COMMON) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

read_selector: $(OUT_BIN_DIR)/read_selector
$(OUT_BIN_DIR)/read_selector: \
	$(OBJ_READ_SELECTOR_MAIN)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ \
	$(OBJ_READ_SELECTOR_MAIN) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

compactors: $(OUT_BIN_DIR)/compactors
$(OUT_BIN_DIR)/compactors: \
	$(OBJ_COMPACTORS_MAIN)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ \
	$(OBJ_COMPACTORS_MAIN) $(CDFLIB_OBJ) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

lookup_table: $(OUT_BIN_DIR)/lookup_table $(OUT_BIN_DIR)/build_lookup_table.py
$(OUT_BIN_DIR)/build_lookup_table.py: $(SRC_LOOKUP_TABLE_MAIN_DIR)/build_lookup_table.py
	cp $(SRC_LOOKUP_TABLE_MAIN_DIR)/build_lookup_table.py $(OUT_BIN_DIR)/build_lookup_table.py
$(OUT_BIN_DIR)/lookup_table: \
	$(OBJ_LOOKUP_TABLE_MAIN) $(OBJ_COMMON_KMC_API)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@  \
	$(MIMALLOC_OBJ) $(REFRESH_PARALLEL_QUEUES_MONITOR_OBJ) \
	$(OBJ_LOOKUP_TABLE_MAIN) $(OBJ_COMMON_KMC_API) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

dsv_manip: $(OUT_BIN_DIR)/dsv_manip
$(OUT_BIN_DIR)/dsv_manip: \
	$(OBJ_DSV_MANIP_MAIN)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ \
	$(OBJ_DSV_MANIP_MAIN) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

gap_shortener: $(OUT_BIN_DIR)/gap_shortener
$(OUT_BIN_DIR)/gap_shortener: \
	$(OBJ_GAP_SHORTENER_MAIN)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@ \
	$(OBJ_GAP_SHORTENER_MAIN) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

download_kmc: $(OUT_BIN_DIR)/kmc_marker
$(OUT_BIN_DIR)/kmc_marker: ./download_kmc.py
	mkdir -p $(OUT_BIN_DIR)
	./download_kmc.py $(OUT_BIN_DIR)
	touch $@
$(OUT_BIN_DIR)/kmc $(OUT_BIN_DIR)/kmc_tools: $(OUT_BIN_DIR)/kmc_marker

supervised_test: $(OUT_BIN_DIR)/supervised_test.R
$(OUT_BIN_DIR)/supervised_test.R: $(SRC_DIR)/supervised_test/supervised_test.R
	cp $(SRC_DIR)/supervised_test/supervised_test.R $(OUT_BIN_DIR)

bkc: $(OUT_BIN_DIR)/bkc
$(OUT_BIN_DIR)/bkc:
	(cd $(3RD_PARTY_DIR)/bkc; make CC=$(CC) CXX=$(CXX) release SHARED_INCLUDE_DIR=../../src/common BKC_OUT_BIN_DIR=../../bin -j)

clean-bkc:
	(cd $(3RD_PARTY_DIR)/bkc; make CC=$(CC) CXX=$(CXX) clean)

# *** Cleaning
.PHONY: clean init install uninstall
clean: clean-libzstd clean-zlib-ng clean-isa-l clean-libdeflate clean-mimalloc_obj clean-cdflib_obj clean-sbwt clean-refresh_parallel_queues_monitor_obj clean-bkc
	-rm -r $(OBJ_DIR)
	-rm -r $(OUT_BIN_DIR)

install: all
	$(eval SPLASH_PREFIX?=/usr/local/bin)
	cp -r bin/* $(SPLASH_PREFIX)

uninstall:
	$(eval SPLASH_PREFIX?=/usr/local/bin)
	$(eval FILES_TO_DELETE := bkc bkc_dump build_lookup_table.py compactors dsv_manip fafq_filter gap_shortener kmc kmc_tools lookup_table read_selector satc satc_dump satc_filter satc_merge satc_to_fasta satc_undump sig_anch splash supervised_test.R tsv_to_fasta)
	-rm -f $(addprefix $(SPLASH_PREFIX)/, $(FILES_TO_DELETE))
	-rm -r $(SPLASH_PREFIX)/postprocessing

init:
	$(call INIT_SUBMODULES)
