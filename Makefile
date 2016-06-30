export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
.SECONDEXPANSION:

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

SAMPLES=1_PGA1_GD_TruSeq \
        1_PGA1_TD_TruSeq \
        2_PGA1_GD_TruSeq \
        2_PGA1_TD_TruSeq \
        3_PGA1_GD_Nextera \
        3_PGA1_TD_Nextera \
        4_PGA1_GD_Nextera \
        4_PGA1_TD_Nextera

all: coverage-read.pdf coverage-fragment.pdf fastqc picard flagstat

#---
#--- READ COVERAGE
#---

coverage-read.pdf: $(foreach S, $(SAMPLES), $S.bedtools-genomecov.reads.txt) /mnt/projects/wgs/scripts/plot-coverage.R
	Rscript /mnt/projects/wgs/scripts/plot-coverage.R \
		--filename-suffix .bedtools-genomecov.reads.txt \
		--title-1 "Coverage" \
		--title-2 "Cumulative coverage" \
		--sort-point 6 \
		--max-coverage 30 \
		--output-pdf coverage-read.pdf.part \
		--output-png coverage-read.png.part
	mv coverage-read.pdf.part coverage-read.pdf
	mv coverage-read.png.part coverage-read.png

%.bedtools-genomecov.reads.txt: /mnt/projects/wgs/data/bam/variant_calling_process_lane_SET_53_C248CACXX_%_recalibrated.bam
	~/tools/samtools-0.1.19/samtools view -b -f 2 -F 1792 -q 1 $< 20 \
		| /data_synology/software/bedtools-2.17.0/bin/bedtools genomecov \
			-ibam stdin \
			-g <(echo -e "20\t59505520\n") \
		| grep ^20 \
		> $@.part
	mv $@.part $@	

#---
#--- FRAGMENT COVERAGE
#---

coverage-fragment.pdf: $(foreach S, $(SAMPLES), $S.bedtools-genomecov.fragments.txt) /mnt/projects/wgs/scripts/plot-coverage.R
	Rscript /mnt/projects/wgs/scripts/plot-coverage.R \
		--filename-suffix .bedtools-genomecov.fragments.txt \
		--title-1 "Physical (fragment) coverage" \
		--title-2 "Cumulative physical (fragment) coverage" \
		--sort-point 18 \
		--max-coverage 30 \
		--output-pdf coverage-fragment.pdf.part \
		--output-png coverage-fragment.png.part
	mv coverage-fragment.pdf.part coverage-fragment.pdf
	mv coverage-fragment.png.part coverage-fragment.png
	
%.bedtools-genomecov.fragments.txt: /mnt/projects/wgs/data/bam/variant_calling_process_lane_SET_53_C248CACXX_%_recalibrated.bam
	~/tools/samtools-0.1.19/samtools sort -@ 10 -no <(~/tools/samtools-0.1.19/samtools view -bh -f 2 -F 1792 -q 1 $< 20) bla \
		| /data_synology/software/bedtools-2.17.0/bin/bamToBed -i stdin -bedpe \
		| cut -f 1,2,6 \
		| sort -k 1,1 \
		| /data_synology/software/bedtools-2.17.0/bin/bedtools genomecov \
			-i stdin \
			-g <(echo -e "20\t59505520\n") \
		| grep ^20 \
		> $@.part
	mv $@.part $@

#---
#--- FASTQC
#---

.PHONY: fastqc
fastqc: $(foreach S, $(SAMPLES), fastqc/$S_fastqc.html)
	
fastqc/%_fastqc.html: /mnt/projects/wgs/data/bam/variant_calling_process_lane_SET_53_C248CACXX_%_recalibrated.bam
	mkdir -p fastqc/$*.part
	/data_synology/software/FastQC-0.11.2/fastqc --outdir fastqc/$*.part --threads 5 $<
	mv fastqc/$*.part/* fastqc
	rmdir fastqc/$*.part

#---
#--- PICARD INSERT SIZE
#---
.PHONY: picard
picard: $(foreach S, $(SAMPLES), picard/$S.chr20.picard.insertsize.out)

picard/%.chr20.picard.insertsize.out: /mnt/projects/wgs/data/bam/variant_calling_process_lane_SET_53_C248CACXX_%_recalibrated.bam ~/tools/picard-tools-1.114/CollectInsertSizeMetrics.jar
	mkdir -p picard
	java -jar ~/tools/picard-tools-1.114/CollectInsertSizeMetrics.jar \
		INPUT=<(~/tools/samtools-0.1.19/samtools view -bh -f 2 -F 1792 -q 1 $< 20) \
		HISTOGRAM_FILE=picard/$*.chr20.picard.insertsize.pdf \
		OUTPUT=$@.part \
		STOP_AFTER=10000000
	mv $@.part $@

#---
#--- SAMTOOLS FLAGSTAT
#---
.PHONY: flagstat
flagstat: $(foreach S, $(SAMPLES), flagstat/$S.chr20.samtools.flagstat)
flagstat/%.chr20.samtools.flagstat: /mnt/projects/wgs/data/bam/variant_calling_process_lane_SET_53_C248CACXX_%_recalibrated.bam
	mkdir -p flagstat
	~/tools/samtools-0.1.19/samtools flagstat <(~/tools/samtools-0.1.19/samtools view -bh -f 2 -F 1792 -q 1 $< 20) 2>&1 1>$@.part | $(LOG)
	mv $@.part $@
