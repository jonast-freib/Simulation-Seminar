# Seminar paper

The paper entry point is `main.tex`.

## Data and figures

The figures in `figures/` are generated from CSV data in `data/`. From the
repository root:

```sh
cmake --build build --target sph_experiments
./build/sph_experiments Latex/data        # about 75 s on 12 cores
python3 Latex/plots/make_plots.py         # needs matplotlib and numpy
```

`data/martin_moyce_n2.csv` holds the experimental dam-break front of
Martin & Moyce (1952) for n^2 = 2; `make_plots.py` adds it to
`figures/dam_front.pdf`.

## Build

Build the PDF from this directory with:

```sh
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

The section sources are in `sections/`, and references are stored in
`references.bib`. Replace the placeholder author. The cited bibliography entries were checked
against publisher records; `ritter1892` only against a secondary citation.