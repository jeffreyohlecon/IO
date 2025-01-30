SHELL=bash
PAPER_SECTIONS = $(shell grep -v '%.*input' paper.tex | grep -o 'sections/[A-Za-z0-9_]*\.tex')
PAPER_INPUTS = $(shell grep --no-filename 'input/' $(PAPER_SECTIONS) | grep -v '^%.*\\input' | grep -o 'input/[A-Za-z0-9_\.\-]*\.[a-z]*')

## RECIPES
all: paper.pdf 

paper.pdf: paper.tex preamble.tex $(PAPER_SECTIONS) $(PAPER_INPUTS)
	pdflatex -shell-escape -draftmode $(basename $<).tex
	pdflatex -draftmode $(basename $<).tex
	pdflatex $(basename $<).tex
	rm $(basename $<).log $(basename $<).out
	rm $(basename $<).bbl $(basename $<).blg

input/%.tex: ../tasks/q2/output/%.tex | input/
	ln -sf ../$< $@
