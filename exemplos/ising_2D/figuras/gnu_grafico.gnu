# gera gráficos da magnetizacao e susceptibilidade magnetica

# mude o nome dos arquivos aqui
mag='mag_N400'
chi='chi_N400'

set terminal epslatex standalone color 12
set output(sprintf('ising2D_%s.tex',mag))
set xlabel('$\beta$')
set ylabel('$M/N$')
plot sprintf('../%s.dat',mag) u 1:($2-$3):($2+$3) with filledcurve lc 0 fs transparent solid .25 notitle, '' u 1:2 with linespoints pt 7 pi 1 lc 2 notitle

set output(sprintf('ising2D_%s.tex',chi))
set xlabel('$\beta$')
set ylabel('$\chi/N$')
plot sprintf('../%s.dat',chi) u 1:($2-$3):($2+$3) with filledcurve lc 0 fs transparent solid .25 notitle, '' u 1:2 with linespoints pt 7 pi 1 lc 2 notitle

