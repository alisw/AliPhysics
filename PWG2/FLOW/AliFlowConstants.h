//////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//
// Description: constants for the flow makers
//  bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla  
//  bla bla bla bla bla bla bla bla bla bla bla bla bla bla ...
//
// Original Authors:                 Art Poskanzer & Raimond Snellings
//
//////////////////////////////////////////////////////////////////////

#ifndef ALIFLOWCONSTANTS_H
#define ALIFLOWCONSTANTS_H

#include <TROOT.h>

class AliFlowConstants {

 public:

  enum {
    kHars        =   2, 					// number of harmonics (>= 2)
    kSels        =   2, 					// number of selections (for each harmonics) 
    kSubs        =   2, 					// number of sub-events
    kPhiBins     = 360, 					// number of phi bins in weighting histograms
    kEtaBins     = 100,						// number of eta bins in FlowAnalysis histograms
    kPtBins      = 100,						// number of pT bins in FlowAnalysis histograms
    kPtBinsPart  = 100,						// number of pT bins in flow histograms
    kCumulIntegOrders = 3, 					// order of the cumulant analysis
    kCumulIntegQmax   = 8,					// ..
    kCumulDiffOrders  = 2,					// ..
    kCumulDiffQmax    = 8,					// ..
    kCents       = 9,						// total number of centrality classes
    kPid         = 6						// total number of p.id. hypotesis
  };

  typedef Double_t PhiWgt_t[kSels][kHars][kPhiBins]; 		// intermediate type to import weights from histograms' file 

  static Int_t   fgMClabel ;	  				// checking the simulation: pTrack->Label()<fMClabel = primary track 
  static Bool_t  fgDebug ;  	  				// for more cout statements (debugging purpose)

 // Histograms limits
  static Float_t fgEtaMin ;					// eta lower limit for FlowAnalysis histograms
  static Float_t fgEtaMax ;					// eta upper limit for FlowAnalysis histograms
  static Float_t fgPtMin ;					// pT lower limit for FlowAnalysis histograms
  static Float_t fgPtMax ;					// pT upper limit for yield histograms
  static Float_t fgPtMaxPart ;					// pT upper limit for flow histograms
  static Float_t fgPtWgtSaturation ;				// flow(pT) saturation value
  static Float_t fgEtaMinTpcOnly ;				// eta lower limit of the full TPC acceptance (unused)
  static Float_t fgEtaMaxTpcOnly ;				// eta upper limit of the full TPC acceptance (unused)

 // Centrality Measurement
  static Float_t  fgEtaMid ;              		    	// Mid-Rapidity interval, used for Centrality measurement
  static Float_t  fgEtaGood ;             		    	// Good Rapidity interval (TPC acceptance)
  static Int_t    fgCent0[AliFlowConstants::kCents] ;		// Expected Multiplicity for each Centrality class
  static Double_t fgBayesian[AliFlowConstants::kPid] ;		// Expected particles' abundance
  static Float_t  fgCentNorm[AliFlowConstants::kCents] ;	// Normalized Multiplicity for each Centrality class
  static Float_t  fgMaxMult ;              		    	// Maximum expected multiplicity

 // Experimental Conditions
  static Double_t fgMagneticField;	 		    	// magnetic field value (0.4 Tesla)
  static Double_t fgCenterOfMassEnergy;   		    	// center of mass energy (5.5 TeV)
  static Short_t  fgBeamMassNumberEast;   		    	// beam mass (Pb = 208)
  static Short_t  fgBeamMassNumberWest;   		    	// beam mass (Pb = 208)

 // ALICE detector measures ...
  static Float_t fgITSx ;	 		    		// inner ITS radial distance from the interaction point
  static Float_t fgTPCx ;			    		// inner TPC radial distance from the interaction point 
  static Float_t fgTRDx ;			    		// inner TRD radial distance from the interaction point 
  static Float_t fgTOFx ;			    		// inner TOF radial distance from the interaction point 

  static Float_t fgMaxInt ;					// big number (to avoid infinity)

  ClassDef(AliFlowConstants,1)          			// macro for rootcint
};

#endif

//////////////////////////////////////////////////////////////////////
