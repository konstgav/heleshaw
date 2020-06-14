LDFLAGS =
LDFLAGS +=
FCFLAGS = -ffree-line-length-512 -O3 -Wall -fbounds-check
FCFLAGS +=

hs.run: main.f90 sip.f90
	gfortran $(FCFLAGS) $^ -o $@ $(LDFLAGS)