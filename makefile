#
# makefile to generate TEMSIM multislice package.
#
# Put this file in the same directory as the TEMSIM
# C source files and type "make all" from a command line,
# to compile all of the programs.
#
#  type "make all" to compile everthing
#  type "make remove" to remove all of the compiled files.
#
# Each file name is assumed to be all lower case.
#
# last modified 8-feb-2006 ejk
#

# define compiler with optimize flag
#CC = gcc -ansi -pedantic -O 		# original
CC = condor_compile icc -ansi -O3 -ipo -static		# intel compiler
#DEL = del  # windows/mingw - doesn't work without .exe in file name
DEL = rm  # unix

# define libraries
MYLIBS = fft2dc.o slicelib.o tiffsubs.o writegfx.o 
LIBS = ${MYLIBS}$ -lm

#
#  entry point to build everything
#
all:
	make atompot
	make autoslice
	make display
	make image
	make mulslice
	make probe
	make slicview
	make stemslic
	make sumpix
	make autostem
	make ktiff2dm
	make autowave
	make autopacbed
	make autocbed
#
#  entry point to remove compiled files
#
remove:
	${DEL}$ atompot
	${DEL}$ autoslice
	${DEL}$ display
	${DEL}$ image
	${DEL}$ mulslice
	${DEL}$ probe
	${DEL}$ slicview
	${DEL}$ stemslic
	${DEL}$ sumpix
	${DEL}$ ktiff2dm
	${DEL}$ autowave
	${DEL}$ fft2dc.o
	${DEL}$ slicelib.o
	${DEL}$ tiffsubs.o
	${DEL}$ autopacbed
	${DEL}$ autostem
	${DEL}$ writegfx.o
#
#  main programs
#

atompot: atompot.c  ${MYLIBS}
	${CC} -o atompot atompot.c ${LIBS}

autoslice: autoslice.c  ${MYLIBS}
	${CC} -o autoslice autoslice.c ${LIBS}

display: display.c  ${MYLIBS}
	${CC} -o display display.c ${LIBS}

image: image.c  ${MYLIBS}
	${CC} -o image image.c ${LIBS}

mulslice: mulslice.c  ${MYLIBS}
	${CC} -o mulslice mulslice.c ${LIBS}

probe: probe.c ${MYLIBS}
	${CC} -o probe probe.c ${LIBS}

slicview: slicview.c ${MYLIBS}
	${CC} -o slicview slicview.c ${LIBS}

stemslic: stemslic.c ${MYLIBS}
	${CC} -o stemslic stemslic.c ${LIBS}

sumpix: sumpix.c ${MYLIBS}
	${CC} -o sumpix sumpix.c ${LIBS}

autostem: autostem.c ${MYLIBS}
	${CC} -o autostem autostem.c ${LIBS}

ktiff2dm: ktiff2dm.c ${MYLIBS}
	${CC} -o ktiff2dm ktiff2dm.c ${LIBS}

autowave: autowave.c ${MYLIBS}
	${CC} -o autowave autowave.c ${LIBS}

autopacbed: autopacbed.c ${MYLIBS}
	${CC} -o autopacbed autopacbed.c ${LIBS}

autocbed: autocbed.c ${MYLIBS}
	${CC} -o autocbed autocbed.c ${LIBS}
#
# define subroutine library
#

slicelib.o: slicelib.c
	${CC} -c slicelib.c

tiffsubs.o: tiffsubs.c
	${CC} -c tiffsubs.c

fft2dc.o: fft2dc.c
	${CC} -c fft2dc.c

writegfx.o: writegfx.c
	${CC} -c writegfx.c
