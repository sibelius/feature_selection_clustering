all:
	pdflatex fsc.tex
	rm -rf *.aux *.log *.bak *.out *.bbl *.blg	
bib:
	bibtex fsc
	pdflatex fsc.tex
	pdflatex fsc.tex
	open fsc.pdf
	rm -rf *.aux *.log *.bak *.out *.bbl *.blg	
clean:
	rm -rf *.aux *.log *.bak *.out *.bbl *.blg
