#--- compiler and linker
cc_flag	:= -c -O3 -D_BIOSTAT
l_flag	:= -lz
cc	:= g++
#--- print message
define printMSG
	@echo -e "-----" $(1) "compilation complete -----\\n"
endef

#--- programmes
progMAIN:= biostat
progCC	:= bioCC
progFD	:= fragDist
progRD	:= readDens
progVA	:= vAlign
progFQS	:= fqStatN
#--- sources directories
dirROOT	:= src
dirBAM	:= $(dirROOT)/bam
dirCMN	:= $(dirROOT)/common
dirCC	:= $(dirROOT)/biocc
dirFD	:= $(dirROOT)/fragdist
dirRD	:= $(dirROOT)/readdens
dirVA	:= $(dirROOT)/valign
dirFQS	:= $(dirROOT)/fqstatn
dirALL	:= $(dirCC) $(dirFD) $(dirRD) $(dirVA) $(dirFQS)
#--- sources
srcBAM	:= $(wildcard $(dirBAM)/*.cpp)
srcCMN	:= $(wildcard $(dirCMN)/*.cpp)
srcCMN_	:= $(notdir $(srcCMN))
srcMAIN	:= $(wildcard $(dirROOT)/*.cpp)
srcCC	:= $(addprefix $(dirCC)/, $(srcCMN_) Calc.cpp bioCC.cpp)
srcFD	:= $(addprefix $(dirFD)/, $(srcCMN_) fragDist.cpp)
srcRD	:= $(addprefix $(dirRD)/, $(srcCMN_) readDens.cpp)
srcVA	:= $(addprefix $(dirVA)/, $(srcCMN_) vAlign.cpp)
srcFQS	:= $(addprefix $(dirFQS)/, $(addsuffix .cpp, fqStatN common TxtFile))

#--- objects
objBAM	:= $(srcBAM:.cpp=.o)
objMAIN	:= $(srcMAIN:.cpp=.o)
objCC	:= $(srcCC:.cpp=.o)
objFD	:= $(srcFD:.cpp=.o)
objRD	:= $(srcRD:.cpp=.o)
objVA	:= $(srcVA:.cpp=.o)
objFQS	:= $(srcFQS:.cpp=.o)

.PHONY: all

# all: copySRC $(srcBAM) $(srcMAIN) $(srcCMN) $(progFD) $(progMAIN)
# all: copySRC $(srcBAM) $(srcMAIN) $(srcCMN) $(progCC) $(progFD) $(progRD) $(progVA) $(progFQS) $(progMAIN)
all: copySRC $(progCC) $(progFD) $(progRD) $(progVA) $(progFQS) $(progMAIN)

$(progCC): $(srcCC) $(objCC) $(objBAM)
	$(cc) $(l_flag) $(objCC) $(objBAM) -o $@
	$(call printMSG,$(progCC))

$(progFD): $(srcFD) $(objFD) $(objBAM)
	$(cc) $(l_flag) $(objFD) $(objBAM) -o $@
	$(call printMSG,$(progFD))

$(progRD): $(srcRD) $(objRD) $(objBAM)
	$(cc) $(l_flag) $(objRD) $(objBAM) -o $@
	$(call printMSG,$(progRD))

$(progVA): $(srcVA) $(objVA) $(objBAM)
	$(cc) $(l_flag) $(objVA) $(objBAM) -o $@
	$(call printMSG,$(progVA))

$(progFQS): $(srcFQS) $(objFQS)
	$(cc) $(l_flag) $(objFQS) -o $@
	$(call printMSG,$(progFQS))

$(progMAIN): $(srcMAIN) $(objMAIN)
	$(cc) $(l_flag) $(objMAIN) -o $@
	$(call printMSG,$(progMAIN))

.cpp.o:
	$(cc) $(cc_flag) $< -o $@

copySRC:	# copy newer common files to prog folders
	@for file in $(dirCMN)/*; do \
	if [ ! -s $(dirCC)/$$(basename $$file) -o $$file -nt $(dirCC)/$$(basename $$file) ]; then \
	echo $(dirALL) | xargs -n 1 cp $$file; echo cp $$file; fi; \
	done

clean:
	find $(dirALL) -type f -name '*.o' -exec rm {} +
