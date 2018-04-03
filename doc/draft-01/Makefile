
FNAME=biharm-evals

TAR_DIR=tardir
ARXIV_DIR=arxiv
FIG_DIR=fig
FIG_GEN_DIR=../../fig

pdf $(FNAME).bbl: $(FNAME).tex refs.bib
	pdflatex $(FNAME)
	bibtex8 $(FNAME) || :
	pdflatex $(FNAME)
	pdflatex $(FNAME)

fast: $(FNAME).tex
	pdflatex $(FNAME)

latest-figs:
	@mkdir -p $(FIG_DIR)
	cp $(FIG_GEN_DIR)/*.png $(FIG_DIR)/
	cp $(FIG_GEN_DIR)/*.pdf $(FIG_DIR)/

tar: pdf
	@mkdir -p $(TAR_DIR)
	cd $(TAR_DIR); rm -rf *
	cp *.tex $(TAR_DIR) > /dev/null || :
	cp -r fig/ $(TAR_DIR) > /dev/null || :
	cp *.cls $(TAR_DIR) > /dev/null || :
	cp *.bst $(TAR_DIR) > /dev/null || :
	cp *.sty $(TAR_DIR) > /dev/null || :
	cp *.bib $(TAR_DIR) > /dev/null || :
	cp *.clo $(TAR_DIR) > /dev/null || :
	cp $(FNAME).pdf $(TAR_DIR) > /dev/null || :
	cp Makefile $(TAR_DIR) > /dev/null || :
	tar -cvf $(FNAME).tar -C $(TAR_DIR) .

tar-flat: pdf
	@mkdir -p $(TAR_DIR)
	cd $(TAR_DIR); rm -rf *
	sed 's_fig/__g' $(FNAME).tex > $(TAR_DIR)/$(FNAME).tex
	sed '/\\bibliography{refs}/ r $(FNAME).bbl' $(TAR_DIR)/$(FNAME).tex > $(TAR_DIR)/$(FNAME).tex.tmp
	sed 's_\\bibliography{refs}__g' $(TAR_DIR)/$(FNAME).tex.tmp > $(TAR_DIR)/$(FNAME).tex
	rm -f $(TAR_DIR)/*.tmp
	cp fig/* $(TAR_DIR) > /dev/null || :
	cp *.cls $(TAR_DIR) > /dev/null || :
	cp *.bst $(TAR_DIR) > /dev/null || :
	cp *.sty $(TAR_DIR) > /dev/null || :
	cp *.bib $(TAR_DIR) > /dev/null || :
	cp *.clo $(TAR_DIR) > /dev/null || :
	cp *.bbl $(TAR_DIR) > /dev/null || :
	cp $(FNAME).pdf $(TAR_DIR) > /dev/null || :
	@printf 'all: \n' > $(TAR_DIR)/Makefile
	@printf '\t pdflatex $(FNAME)\n' >> $(TAR_DIR)/Makefile
	@printf '\t pdflatex $(FNAME)\n' >> $(TAR_DIR)/Makefile
	@printf '\t pdflatex $(FNAME)\n' >> $(TAR_DIR)/Makefile
	tar -cvf $(FNAME).tar -C $(TAR_DIR) .

tar-arxiv: $(FNAME).bbl
	@mkdir -p $(ARXIV_DIR)
	cd $(ARXIV_DIR); rm -f *
	sed 's_fig/__g' $(FNAME).tex > $(ARXIV_DIR)/$(FNAME).tex
	cp fig/* $(ARXIV_DIR) || :
	cp *.cls $(ARXIV_DIR) || :
	cp *.bst $(ARXIV_DIR) || :
	cp *.sty $(ARXIV_DIR) || :
	cp *.bib $(ARXIV_DIR) || :
	cp *.bbl $(ARXIV_DIR) || :
	cd $(ARXIV_DIR); tar -cvf $(FNAME).tar *

clean:
	rm -f *.log
	rm -f *.pdf
	rm -f *.spl
	rm -f *.blg
	rm -f *.aux
	rm -f *.bbl
