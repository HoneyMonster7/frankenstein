

linefiles=system('ls *.fitt')

#set terminal size 1000,1000
#set terminal x11 geometry 1000x1000
#set size 2,1

set title "Fittness of the best organism in each jobs"

plot for [file in linefiles] file w l smooth sbezier title file

#pause -1
