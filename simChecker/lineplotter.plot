set macros

linefiles=system('ls *.fitt | tail -n +7|  head -n 1 ')

#set terminal size 1000,1000
#set terminal x11 geometry 1000x1000
#set size 2,1

set key outside

set title "Fittness of the best organism in each jobs"

YLABS="unset ylabel; set ylabel 'Fittness'"
XLABS="unset xlabel; set xlabel 'Iteration (*e7'"

set multiplot layout 3,1 columnsfirst

#plot for [file in linefiles] file w l smooth sbezier title file
@YLABS; @XLABS
plot for [col=1:5] 'fittavgs.fitt' using col with lines title columnheader

YLABS="unset ylabel; set label 'IPR'"
XLABS="unset xlabel; set label 'Checkpoint'"
plot for [col=6:10] 'fittavgs.fitt' using col with lines title columnheader

YLABS="unset ylabel; set label 'IPR'"
XLABS="unset xlabel; set label 'Checkpoint'"
plot for [col=11:15] 'fittavgs.fitt' using col with lines title columnheader
#pause -1
