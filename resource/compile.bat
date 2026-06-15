@echo off
echo Compiling main.tex using Docker (texlive/texlive:latest)...

docker run --rm -v "%cd%":/workdir -w /workdir texlive/texlive:latest pdflatex -interaction=nonstopmode main.tex
docker run --rm -v "%cd%":/workdir -w /workdir texlive/texlive:latest bibtex main
docker run --rm -v "%cd%":/workdir -w /workdir texlive/texlive:latest pdflatex -interaction=nonstopmode main.tex
docker run --rm -v "%cd%":/workdir -w /workdir texlive/texlive:latest pdflatex -interaction=nonstopmode main.tex

echo.
echo Compilation finished! Check main.pdf in the resource directory.
pause
