#ifndef CONSTANTS_H
#define CONSTANTS_H

#ifndef __CINT__

#ifndef TRUE
#define TRUE			1
#endif

#ifndef FALSE
#define FALSE			0
#endif

#define WHERE __FILE__ << ":" << __LINE__ << ": "

#endif

// from MWPC
const int nMWPCHitMax = 10000;
const int NMWPCPLANES = 7;

// for the mapping and geometry constants
const int ncards = 2;
const int rows = 8;
const int cols = 8;
const int gains	= 2;
const int FEC[ncards] = {5, 6}; // cards in the crate
// we should have one charge-sensitive pre-amp per module
const int NCSP = rows*cols;

// for the readout
const int NFEC = 2;                             // Max 2 FrontEndCard
const int NCHIP = 4;                            // 4 ALTROs per FEC
const int NCHAN = 16;                           // Channels per ALTRO
const int CHANNELS = NFEC*NCHIP*NCHAN;		// Max uses ALTRO channels 
const int TOTCHAN = CHANNELS;

const int REQSAMPLES = 50;                     // This is what we ask for in the RCU
const int EXTRASAMPLES = 15;                    // We get a few bonus ones
const int REALSAMPLES = REQSAMPLES + EXTRASAMPLES;	// This is the reported wordcount-2
const int SAMPLES = REQSAMPLES + EXTRASAMPLES;		        // Number of readout ALTRO words per channel (last 4 words is trailer)

// general flags
const int debug = 0;            		// Set 0/1/2/3 to print debug messages of different detail levels

// put in nominal errors for pedestal for fits
// low, and high gains are separate
const double  NOMINAL_PEDESTAL_RMS[NFEC][2] = { {0.3, 2.0},   // FEC 1
						{0.6, 4.0} }; // FEC 2 

#endif
