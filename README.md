# Quasi-Love-Wave-Detections
- Matlab code for the analysis of Quasi-Love waves from seismic surface wave data. 
- Written by Caroline Eakin, ANU (caroline.eakin@anu.edu.au)
- Please acknowledge and refer to the following paper for further details: Eakin (202-) xxx, Communications Earth & Environment, xxx. 
- Pre-print available here: https://www.researchsquare.com/article/rs-121788/v1
- Code reproduces Figure S1 of the above paper. 

**Overview of Quasi-Love Waves**

Quasi-Love waves are thought to be generated by Love-to-Rayleigh scattering due to lateral gradients in seismic anisotropy in the upper mantle. Quasi-Love waves have a similar waveform shape to the fundamental Love wave (G1) but display Rayleigh wave motion, i.e. they appear on the vertical and radial components, and have elliptical particle motion. Quasi-Love waves should arrive after the fundamental Love wave but before the fundamental Rayleigh wave. The delay time between the fundamental Love wave arrival and the Quasi-Love wave can be used to determine the distance to the scattering source. By back-projecting along the great-circle path the location at which the scattering occured (i.e. the location of the lateral gradient in seismic anisotropy) can be determined. 

**Overview of Files**
1) Inspect4QLwave.m 
- Matlab script to plot waveforms and assess the presence of a Quasi-Love wave detection. If detected, the coordinates of the Quasi-Love scatterer, its amplitude, and polarity are calculated and stored. Requires visual inspection and user input. For each event will ask if you wish to "Keep (1), or Discard (0)?", in other words has a clear Quasi-Love wave been detected? 

2) TestData_MORW_FigS1.mat
- Example input data for 1 event in order to run Inspect4QLwave.m 
- See Figure S1 of associated paper for further event details

**Steps to run**
1. Open Matlab window
2. Type the following 2 lines:
   load('TestData_MORW_FigS1.mat');
   out=Inspect4QLwave(E,N,Z,t,sdt,dis,bazi,slat,slong);
