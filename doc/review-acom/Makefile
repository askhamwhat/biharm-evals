
FNAME=stokes-eig

TAR_DIR=tardir
ARXIV_DIR=arxiv
FIG_DIR=fig

FIGEPS= annulus_discretization.eps \
	conv-studies.eps \
	det_spur.eps \
	ex_barbell_001_bdry.eps \
	ex_barbell_001_coeffs_plot.eps \
	ex_barbell_001_sings_plot.eps \
	ex_many_holes_004_bdry.eps \
	ex_many_holes_004_coeffs_plot.eps \
	ex_many_holes_004_sings_plot.eps \
	ex_many_holes_004_sings_plot_ref.eps \
	ex_many_holes_004_sings_plot_west.eps \
	ext_dom.eps \
	mc_dom.eps \
	multiply_final.eps \
	speed_res.eps


pdf $(FNAME).bbl: $(FNAME).tex refs.bib
	pdflatex $(FNAME)
	bibtex $(FNAME) || :
	pdflatex $(FNAME)
	pdflatex $(FNAME)

fast: $(FNAME).tex
	pdflatex $(FNAME)

fixenviros:
	sed -i 's/{lem}/{lemma}/g' *.tex
	sed -i 's/{thrm}/{theorem}/g' *.tex
	sed -i 's/{cor}/{corollary}/g' *.tex
	sed -i 's%fig/%%g' *.tex
	sed -i 's%fig/%%g' *.tex
	sed -i 's%media/%%g' *.tex

epsfigs: $(FIGEPS)

%.eps:%.pdf
	inkscape $< --export-ps=$@
