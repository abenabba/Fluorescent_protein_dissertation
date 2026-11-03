#===============================================================================
# edited to include new motifs for barcode extraction 
SHELL := /bin/bash

# -------------------------
# Directories
# -------------------------
PROJ_DIR := /projects/plesalab/shared/ab_project/illumina_shuffle_analysis/GC3F-AB-11718/makefile_run
CUT3_DIR := $(PROJ_DIR)/cut3
FINAL_BC_3_DIR := $(CUT3_DIR)/Final_BC_3
LOG_3_DIR := $(PROJ_DIR)/logs_3

# Create directories if needed
$(shell mkdir -p $(CUT3_DIR) $(FINAL_BC_3_DIR) $(LOG_3_DIR))

# -------------------------
# Tools
# -------------------------
STARCODE := starcode
SORT := sort -T $(PROJ_DIR) --parallel=8

# -------------------------
# Barcode extraction settings
# Barcode must be exactly BC_LEN nt between:
# LEFT_FLANK  (<= LEFT_MM mismatches)
# RIGHT_FLANK (<= RIGHT_MM mismatches)
# -------------------------
LEFT_FLANK := TGGCTGCGGAAC
RIGHT_FLANK := GCACGACGTCAG
LEFT_MM := 1
RIGHT_MM := 2
BC_LEN := 20

# -------------------------
# Samples
# -------------------------
SAMPLES := \
TTGTGATCCT-ATGCCTTCAT-Medium-Illegal-Char-Slash-brighter-library_S32
#  B12_S11 B13_S12 B14_S13 B15_S14 B16_S15 B17_S16 B18_S17 B19_S18 B20_S19 \
#  Brighter-Library_S20 \
#  MB12_S21 MB13_S22 MB14_S23 MB15_S24 MB16_S25 MB17_S26 MB18_S27 MB19_S28 MB20_S29 MB21_S30 MB22_S31 \
#  Medium-Library_S10 \
#  P12_S1 P13_S2 P14_S3 P15_S4 P16_S5 P17_S6 P18_S7 P19_S8 P20_S9

# -------------------------
# Output file groups
# -------------------------
TRUEBC_FILES       := $(foreach s,$(SAMPLES),$(CUT3_DIR)/$(s)_mrg_flip_trim_SE.truebc.txt)
TRUEBC_COUNT_FILES := $(foreach s,$(SAMPLES),$(CUT3_DIR)/$(s).truebc.count.txt)
TRUEBC_FIX_FILES   := $(foreach s,$(SAMPLES),$(CUT3_DIR)/$(s).truebc.fix.count.txt)
TRUEBC_STAR_FILES  := $(foreach s,$(SAMPLES),$(FINAL_BC_3_DIR)/$(s)_truebc_collapse_d1.tsv)

.SECONDARY: $(TRUEBC_FILES) $(TRUEBC_COUNT_FILES) $(TRUEBC_FIX_FILES)

.PHONY: truebc counts truebc_starcode all

# Default target
all: truebc_starcode

# -------------------------
# Step 1: Extract barcode between two flanks
# Accept only if:
#   left flank matches with <= LEFT_MM mismatches
#   barcode is exactly BC_LEN nt
#   right flank matches with <= RIGHT_MM mismatches
# -------------------------
$(CUT3_DIR)/%_mrg_flip_trim_SE.truebc.txt: $(PROJ_DIR)/%_mrg_flip_trim_SE.fastq.gz
	gzip -cd $< | awk 'NR%4==2 {print toupper($$0)}' | \
	awk 'function hdist(a,b,i,d){d=0; for(i=1;i<=length(a);i++){if(substr(a,i,1)!=substr(b,i,1)) d++} return d} {seq=$$0; left="$(LEFT_FLANK)"; right="$(RIGHT_FLANK)"; llen=length(left); rlen=length(right); for(pos=1; pos<=length(seq)-llen-$(BC_LEN)-rlen+1; pos++){leftwin=substr(seq,pos,llen); if(index(leftwin,"N")>0) continue; if(hdist(leftwin,left) <= $(LEFT_MM)){ bcstart=pos+llen; bc=substr(seq,bcstart,$(BC_LEN)); if(length(bc)!=$(BC_LEN) || bc !~ /^[ACGT]+$$/) continue; rightwin=substr(seq,bcstart+$(BC_LEN),rlen); if(length(rightwin)!=rlen) continue; if(index(rightwin,"N")>0) continue; if(hdist(rightwin,right) <= $(RIGHT_MM)){ print bc; break } } }}' > $@

# -------------------------
# Step 2: Count extracted barcodes
# Output format: barcode <tab> count
# -------------------------
$(CUT3_DIR)/%.truebc.count.txt: $(CUT3_DIR)/%_mrg_flip_trim_SE.truebc.txt
	cat $< | $(SORT) | uniq -c | awk '{print $$2 "\t" $$1}' > $@

# -------------------------
# Step 3: Prepare starcode input
# -------------------------
$(CUT3_DIR)/%.truebc.fix.count.txt: $(CUT3_DIR)/%.truebc.count.txt
	tail -n +1 $< > $@

# -------------------------
# Convenience targets
# -------------------------
truebc: $(TRUEBC_FILES)

counts: $(TRUEBC_COUNT_FILES)

# -------------------------
# Step 4: Collapse corrected barcodes with starcode
# Pattern rule: allows building one specific sample output
# -------------------------
$(FINAL_BC_3_DIR)/%_truebc_collapse_d1.tsv: $(CUT3_DIR)/%.truebc.fix.count.txt
	$(STARCODE) -d1 --sphere \
	  -i $< \
	  --output $@ \
	  --print-clusters \
	  1>$(LOG_3_DIR)/$*_truebc_starcode.stdout.log \
	  2>$(LOG_3_DIR)/$*_truebc_starcode.stderr.log

# Build all corrected starcode outputs
truebc_starcode: $(TRUEBC_STAR_FILES)