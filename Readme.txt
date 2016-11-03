Multislice package updated 11/3/16, cz

Updates in 11/3/16 version:
-minor bug fixed for using three or more collection angle ranges
-use inner collection anlge to determine the mrad/px ratio on detector sensitivity map
-would ask user whether outer collection anlge is determined by aperture
	- if yes, a cutoff would be applied at outer collection angle
	- if no, the whole detecor would be used regardless of the outer collection angle entered by user




To install the package, type 

	make

to compile everything with icc compiler

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

to compile part of the package, type

	make remove

to remove all compiled files. 

