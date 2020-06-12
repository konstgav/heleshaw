reset

L = 20.0
H = 40.0
Nsteps = 40
minVal = -200.0
maxVal = 200.0
numOfContours = 60.0
step = (maxVal -  minVal)/numOfContours

set term png size 300,600 crop
#set xrange [0:L]
#set yrange [0:H]
#set size ratio H/L

do for [ii=0:Nsteps] {
	
	set size ratio H/L
	
	outputFileName = sprintf('animation/%08.0f.png',ii)
	inputFileName  = sprintf('psi/%08.0f.dat',ii)
	inputFileNameT  = sprintf('T/%08.0f.dat',ii)
	set output outputFileName

	set multiplot

	# Surface
	set surface
	set pm3d map interpolate 0,0
	unset colorbox
	unset key
	set palette rgbformulae 33,13,10
	splot inputFileNameT

	# Contour
	set contour base
	set cntrparam level incremental minVal, step, maxVal
	unset surface
	unset key
	unset clabel
	splot inputFileName w l lw 1.5 lt -1

	unset multiplot
	reset
}