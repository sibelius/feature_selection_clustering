all:
	pdflatex paper.tex
	rm -rf *.aux *.log *.bak *.out *.bbl *.blg	
bib:
	bibtex paper
	pdflatex paper.tex
	pdflatex paper.tex
	open paper.pdf
	rm -rf *.aux *.log *.bak *.out *.bbl *.blg	
clean:
	rm -rf *.aux *.log *.bak *.out *.bbl *.blg
