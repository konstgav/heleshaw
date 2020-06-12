gnuplot fieldsplot.plt
cd animation
mencoder.exe mf://*.png -mf fps=5:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o motion.avi
mplayer motion.avi