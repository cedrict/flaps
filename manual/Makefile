paper = manual
all:
	pdflatex  --shell-escape $(paper).tex
	bibtex $(paper)
	pdflatex  --shell-escape $(paper).tex
	pdflatex  --shell-escape $(paper).tex
	rm -f *.out *.log *.aux *.bbl *.blg *.spl
clean:
	rm -f *.out *.log *.aux *.bbl *.blg *.spl *.synctex.* 

cleanall:
	rm -f *.aux *.bbl *.blg *.log *.out *.spl *.synctex.* $(paper).pdf
