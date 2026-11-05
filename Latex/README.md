# Seminar paper

The paper entry point is `main.tex`.

## Data and figures

The figures in `figures/` are generated from CSV data in `data/`. From the
repository root:

```sh
cmake --build build --target sph_experiments
./build/sph_experiments Latex/data        # about 40 s on 12 cores
python3 Latex/plots/make_plots.py         # needs matplotlib and numpy
```

`data/martin_moyce_n2.csv` is a template for the experimental dam-break front
of Martin & Moyce (1952). Fill in the measured points from the original paper;
`make_plots.py` adds them to `figures/dam_front.pdf` automatically.

## Build

Build the PDF from this directory with:

```sh
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

The section sources are in `sections/`, and references are stored in
`references.bib`. Replace the placeholder author. The bibliography entries added
for the methods and results sections were written from memory; check their
volume, page and year details before submission.