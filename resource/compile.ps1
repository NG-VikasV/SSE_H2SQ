$ResourceDir = Split-Path -Parent $MyInvocation.MyCommand.Definition
Write-Host "Compiling main.tex using Docker (texlive/texlive:latest)..."

docker run --rm -v "$($ResourceDir):/workdir" -w /workdir texlive/texlive:latest pdflatex -interaction=nonstopmode main.tex
docker run --rm -v "$($ResourceDir):/workdir" -w /workdir texlive/texlive:latest bibtex main
docker run --rm -v "$($ResourceDir):/workdir" -w /workdir texlive/texlive:latest pdflatex -interaction=nonstopmode main.tex
docker run --rm -v "$($ResourceDir):/workdir" -w /workdir texlive/texlive:latest pdflatex -interaction=nonstopmode main.tex

Write-Host "Compilation finished! Check main.pdf in the resource directory."
