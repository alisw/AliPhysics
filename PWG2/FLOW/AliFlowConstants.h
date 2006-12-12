//////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//
// Description: constants for the flow makers
//
// Original Authors:                 Art Poskanzer & Raimond Snellings
//
//////////////////////////////////////////////////////////////////////

#ifndef AliFlowConstants_h
#define AliFlowConstants_h

#include <TROOT.h>

class Flow {

 public:

  enum {
    nHars        =   2, 					// number of harmonics (>= 2)
    nSels        =   2, 					// number of selections (for each harmonics) 
    nSubs        =   2, 					// number of sub-events
    nPhiBins     = 360, 					// number of phi bins in weighting histograms
    nEtaBins     = 100,						// number of eta bins in FlowAnalysis histograms
    nPtBins      = 100,						// number of pT bins in FlowAnalysis histograms
    nPtBinsPart  = 100,						// number of pT bins in flow histograms
    nCumulIntegOrders =   3, 					// order of the cumulant analysis
    nCumulInteg_qMax  =   8,					// ..
    nCumulDiffOrders  =   2,					// ..
    nCumulDiff_qMax   =   8,					// ..
    nCents       = 9,						// total number of centrality classes
    nPid         = 6						// total number of p.id. hypotesis
  };

  typedef Double_t PhiWgt_t[nSels][nHars][nPhiBins]; 		// intermediate type to import weights from histograms' file 

  static Int_t   fMClabel ;	  				// checking the simulation: pTrack->Label()<fMClabel = primary track 
  static Bool_t  fDebug ;  	  				// for more cout statements (debugging purpose)

 // Histograms limits
  static Float_t fEtaMin ;					// eta lower limit for FlowAnalysis histograms
  static Float_t fEtaMax ;					// eta upper limit for FlowAnalysis histograms
  static Float_t fPtMin ;					// pT lower limit for FlowAnalysis histograms
  static Float_t fPtMax ;					// pT upper limit for yield histograms
  static Float_t fPtMaxPart ;					// pT upper limit for flow histograms
  static Float_t fPtWgtSaturation ;				// flow(pT) saturation value
  static Float_t fEtaMinTpcOnly ;				// eta lower limit of the full TPC acceptance (unused)
  static Float_t fEtaMaxTpcOnly ;				// eta upper limit of the full TPC acceptance (unused)

 // Centrality Measurement
  static Float_t  fEtaMid ;              		    	// Mid-Rapidity interval, used for Centrality measurement
  static Float_t  fEtaGood ;             		    	// Good Rapidity interval (TPC acceptance)
  static Int_t    fCent0[Flow::nCents] ;			// Expected Multiplicity for each Centrality class
  static Double_t fBayesian[Flow::nPid] ;			// Expected particles' abundance
  static Float_t  fCentNorm[Flow::nCents] ;			// Normalized Multiplicity for each Centrality class
  static Float_t  fMaxMult ;              		    	// Maximum expected multiplicity

 // Experimental Conditions
  static Double_t fMagneticField;	 		    	// magnetic field value (0.4 Tesla)
  static Double_t fCenterOfMassEnergy;   		    	// center of mass energy (5.5 TeV)
  static Short_t  fBeamMassNumberEast;   		    	// beam mass (Pb = 208)
  static Short_t  fBeamMassNumberWest;   		    	// beam mass (Pb = 208)

 // ALICE detector measures ...
  static Float_t fITSx ;	 		    		// inner ITS radial distance from the interaction point
  static Float_t fTPCx ;			    		// inner TPC radial distance from the interaction point 
  static Float_t fTRDx ;			    		// inner TRD radial distance from the interaction point 
  static Float_t fTOFx ;			    		// inner TOF radial distance from the interaction point 

  static Float_t  fMaxInt ;					// big number (to avoid infinity)

  ClassDef(Flow,1)               				// macro for rootcint
};

typedef Flow AliFlowConstants; //PK: Make rootcint happy, class name=file name

#endif

//////////////////////////////////////////////////////////////////////
