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

//  // Get methods  (not working ...)
//   Int_t   MClabel() const	   { return fgMClabel ; }
//   Bool_t  Debug() const	   { return fgDebug ; } 
//   Float_t MaxInt() const	   { return fgMaxInt ; }	  
// 
//   Float_t EtaMin() const	   { return fgEtaMin ; }	         
//   Float_t EtaMax() const	   { return fgEtaMax ; }	         
//   Float_t PtMin() const	   { return fgPtMin ; } 	       
//   Float_t PtMax() const	   { return fgPtMax ; } 	       
//   Float_t PtMaxPart() const	   { return fgPtMaxPart ; }	          
//   Float_t PtWgtSaturation() const  { return fgPtWgtSaturation ; }         
// 
//   Float_t EtaMid() const	   { return fgEtaMid ; }	 
//   Float_t EtaGood() const	   { return fgEtaGood ; }	        	
// 
//   Float_t  MaxMult() const	   { return fgMaxMult ; }			
//   Int_t    Cent(Int_t i) const     { return (Int_t)(fgMaxMult * fgCentNorm[i]) ; }
//   Int_t    Cent0(Int_t i) const	   { return fgCent0[i] ; }
//   Double_t Bayesian(Int_t i) const { return fgBayesian[i] ; }
// 
//   Float_t  ITSx() const     	   { return fgITSx ; }      	   
//   Float_t  TPCx() const     	   { return fgTPCx ; }      	   
//   Float_t  TRDx() const     	   { return fgTRDx ; }      	   
//   Float_t  TOFx() const     	   { return fgTOFx ; }      	   
// 
//   Double_t MagneticField() const      { return fgMagneticField ; } 	      
//   Double_t CenterOfMassEnergy() const { return fgCenterOfMassEnergy ; }  
//   Short_t  BeamMassNumberEast() const { return fgBeamMassNumberEast ; }  
//   Short_t  BeamMassNumberWest() const { return fgBeamMassNumberWest ; }  
// 
//  // Set methods
//   void SetEtaGood(Float_t eta)	           { fgEtaGood = eta; }			
//   void SetMClabel(Int_t mc)		   { fgMClabel = mc ; }
//   void SetDebug(Int_t db)		   { fgDebug = db ; } 
//   void SetMaxMult(Float_t mult)   	   { fgMaxMult = mult ; }		
//   void SetCent0(Int_t i, Int_t mult)       { fgCent0[i] = mult ; }
//   void SetBayesian(Int_t i, Double_t mult) { fgBayesian[i] = mult ; }
// 
// 
//  private:

 // Histograms limits
  static Float_t  fgEtaMin ;	     			        // eta lower limit for FlowAnalysis histograms
  static Float_t  fgEtaMax ;	     			        // eta upper limit for FlowAnalysis histograms
  static Float_t  fgPtMin ;	     			        // pT lower limit for FlowAnalysis histograms
  static Float_t  fgPtMax ;	     			        // pT upper limit for yield histograms
  static Float_t  fgPtMaxPart ;      			        // pT upper limit for flow histograms
  static Float_t  fgPtWgtSaturation ;			        // flow(pT) saturation value
  static Float_t  fgEtaMinTpcOnly ;  			        // eta lower limit of the full TPC acceptance (-0.9)
  static Float_t  fgEtaMaxTpcOnly ;  			        // eta upper limit of the full TPC acceptance (+0.9)

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

 // ALICE detector measures
  static Float_t  fgITSx ;      			        // inner ITS radial distance from the interaction point
  static Float_t  fgTPCx ;      			        // inner TPC radial distance from the interaction point 
  static Float_t  fgTRDx ;      			        // inner TRD radial distance from the interaction point 
  static Float_t  fgTOFx ;      			        // inner TOF radial distance from the interaction point 

 // ...
  static Int_t    fgMClabel ;				        // checking the simulation: pTrack->Label()<fMClabel = primary track 
  static Bool_t   fgDebug ;  				        // for more cout statements (debugging purpose)
  static Float_t  fgMaxInt ; 				        // big number (to avoid overflows)

  ClassDef(AliFlowConstants,2)          			// macro for rootcint
};

#endif

//////////////////////////////////////////////////////////////////////
