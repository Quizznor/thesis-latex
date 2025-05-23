# Makefile for building latex documents via biber

SHELL := bash
TARGET := main
BIBTEX := biber
LATEX := latexmk -pdf

DEPENDENCIES := $(wildcard *.tex)
DEPENDENCIES += $(wildcard abstract/*.tex)
DEPENDENCIES += $(wildcard chap/*.tex)
DEPENDENCIES += $(wildcard app/*.tex)
DEPENDENCIES += $(wildcard include/*.tex)
DEPENDENCIES += $(wildcard *.sty *.bib)

##########################################################

.PHONY: all clean

all:
	scripts/sort_acronyms.py
	scripts/gather_acronyms
	scripts/gather_chapters
	scripts/gather_appendices
	#$(LATEX) $*
clean:
	- $(RM) -f $(TARGET){.fls,.bbl,.synctex.gz,.fdb_latexmk,.pdf,.tex~,.ps,.dvi,.log,.aux,.toc,.bcf,.nav,.out,.snm,.bbl,.blg,-blx.bib,.tex.bbl,.tex.blg,.run.xml,.spl,.lof,.maf,.fff,.lot,.ttt,.mtc*,.ptc,.loc,.tdo,.soc} */*.aux
