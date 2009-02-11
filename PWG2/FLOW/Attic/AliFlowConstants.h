//////////////////////////////////////////////////////////////////////
//
// $Id: AliFlowConstants.h 18618 2007-05-16 15:38:22Z snelling $
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

namespace AliFlowConstants {

 // Enumerators
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

 // new type
  typedef Double_t PhiWgt_t[kSels][kHars][kPhiBins]; 		// 3d_array type,used to import phi weights from histograms 

 // Histograms limits
  extern Float_t  fgEtaMin ;	     			        // eta lower limit for FlowAnalysis histograms
  extern Float_t  fgEtaMax ;	     			        // eta upper limit for FlowAnalysis histograms
  extern Float_t  fgPtMin ;	     			        // pT lower limit for FlowAnalysis histograms
  extern Float_t  fgPtMax ;	     			        // pT upper limit for yield histograms
  extern Float_t  fgPtMaxPart ;      			        // pT upper limit for flow histograms
  extern Float_t  fgPtWgtSaturation ;			        // flow(pT) saturation value
  extern Float_t  fgEtaMinTpcOnly ;  			        // eta lower limit of the full TPC acceptance (-0.9)
  extern Float_t  fgEtaMaxTpcOnly ;  			        // eta upper limit of the full TPC acceptance (+0.9)

 // Centrality Measurement
  extern Float_t  fgEtaMid ;              		    	// Mid-Rapidity interval, used for Centrality measurement
  extern Float_t  fgEtaGood ;             		    	// Good Rapidity interval (TPC acceptance)
  extern Int_t    fgCent0[AliFlowConstants::kCents] ;		// Expected Multiplicity for each Centrality class
  extern Double_t fgBayesian[AliFlowConstants::kPid] ;		// Expected particles' abundance
  extern Float_t  fgCentNorm[AliFlowConstants::kCents] ;	// Normalized Multiplicity for each Centrality class
  extern Float_t  fgMaxMult ;              		    	// Maximum expected multiplicity

 // Experimental Conditions
  extern Double_t fgMagneticField;	 		    	// magnetic field value (0.4 Tesla)
  extern Double_t fgCenterOfMassEnergy;   		    	// center of mass energy (5.5 TeV)
  extern Short_t  fgBeamMassNumberEast;   		    	// beam mass (Pb = 208)
  extern Short_t  fgBeamMassNumberWest;   		    	// beam mass (Pb = 208)

 // ALICE detector measures
  extern Float_t  fgITSx ;      			        // inner ITS radial distance from the interaction point
  extern Float_t  fgTPCx ;      			        // inner TPC radial distance from the interaction point 
  extern Float_t  fgTRDx ;      			        // inner TRD radial distance from the interaction point 
  extern Float_t  fgTOFx ;      			        // inner TOF radial distance from the interaction point 

 // ...
  extern Int_t    fgMClabel ;				        // checking the simulation: pTrack->Label()<fMClabel = primary track 
  extern Bool_t   fgDebug ;  				        // for more cout statements (debugging purpose)
  extern Float_t  fgMaxInt ; 				        // big number (to avoid overflows)

}

#endif

//////////////////////////////////////////////////////////////////////
