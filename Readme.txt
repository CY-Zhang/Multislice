Multislice package updated 11/3/16, cz

Updates in 06/09/17 version:
- Major bug fixed in Autocbed which moves the CBED assignment loop outside wobbler loop, CBED pattern should converge with number of phonons now

Updates in 3/17/17 version:
- new function 'autocbed' to simulate bunch of cbed pattern added
- one cbed would be simulated for each probe position and thickness (total nxout*nyout*nThick)
- PACBED pattern and STEM image would be simulated too

Updates in 11/3/16 version:
-minor bug fixed for using three or more collection angle ranges
-use inner collection anlge to determine the mrad/px ratio on detector sensitivity map
-would ask user whether outer collection anlge is determined by aperture
	- if yes, a cutoff would be applied at outer collection angle
	- if no, the whole detecor would be used regardless of the outer collection angle entered by user
-simulation code supports the use of a normalized detector sensitivity map with rotation option availiable




To install the package, type 

	make all

to compile everything with icc compiler on Condor system

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

to compile part of the package, type

	make remove

to remove all compiled files. 

