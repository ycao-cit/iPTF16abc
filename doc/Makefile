#
# LaTex_Makefile
# Yi Cao, 2016-12-08 09:45
#

OBJ = ms
LATEX = /usr/bin/pdflatex
BIBTEX = /usr/bin/bibtex

all:
	for item in $(OBJ); \
		do \
		$(LATEX) $$item; \
		$(BIBTEX) $$item; \
		$(LATEX) $$item; \
		$(LATEX) $$item; \
		done;

clean:
	rm -f *.{blg,aux,out,log}

# vim:ft=make
#
