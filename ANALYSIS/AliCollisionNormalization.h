
#ifndef ALICOLLISIONNORMALIZATION_H
#define ALICOLLISIONNORMALIZATION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                      Implementation of   Class AliCollisionNormalization
//
//  This class is used to store the vertex ditributions in the data
//  and in Monte Carlo, needed to compute the real number of
//  collisions a given sample is corresponding to.
//  The strategy matches what described in CERN-THESIS-2009-033 p 119.
//
//    Author:     Michele Floris, CERN
//-------------------------------------------------------------------------

#include "TH2F.h"

class TH1F;
class TH1I;
class AliMCEvent;

class AliCollisionNormalization : public TObject

{


public:
  enum { kNevBin0, kNevCollisions, kNevNbin };
  typedef enum { kProcSD, kProcDD, kProcND, kProcUnknown, kNProcs } ProcType_t; 
  AliCollisionNormalization();
  AliCollisionNormalization(Int_t nbinz, Float_t minz, Float_t maxz);
  AliCollisionNormalization(const char * dataFile, const char * dataListName, 
			    const char * mcFile,   const char * mcListName,
			    const char * eventStatFile);

  ~AliCollisionNormalization();
  
  void SetMC(Bool_t flag = kTRUE) { fIsMC = flag;}

  void BookAllHistos();
  TH1 * BookVzHisto(const char * name , const char * title, Bool_t vzOnly=kFALSE);

  void FillVzMCGen(Float_t vz, Int_t ntrk, AliMCEvent * mcEvt);      
  void FillVzMCRec(Float_t vz, Int_t ntrk, AliMCEvent * mcEvt);      
  void FillVzMCTrg(Float_t vz, Int_t ntrk, AliMCEvent * mcEvt);      
  void FillVzData (Float_t vz, Int_t ntrk)      {fHistVzData       ->Fill(vz,ntrk);}

  TH2F *   GetVzMCGen       (Int_t procType) ;
  TH2F *   GetVzMCRec       (Int_t procType) ;
  TH2F *   GetVzMCTrg       (Int_t procType) ;
  TH2F *   GetVzData        () { return fHistVzData       ; }
  TH1F *   GetStatBin0      () { return fHistStatBin0     ; }
  TH1F *   GetStat          () { return fHistStat         ; }
  TH1F *   GetHistProcTypes () { return fHistProcTypes    ; }
   

  Int_t GetProcessType(AliMCEvent * mcEvt) ;
  Double_t GetProcessWeight(Int_t proctype);

  void SetReferencsXS(Int_t ref) { fReferenceXS = ref;}
  
  Double_t ComputeNint();
  void SetZRange(Float_t zrange) { fZRange = zrange ;}

  void SetReferenceXS(Int_t ref) { fReferenceXS = ref ;}
  void GetRelativeFractions(Int_t origin, Float_t& ref_SD, Float_t& ref_DD, Float_t& ref_ND, Float_t& error_SD, Float_t& error_DD, Float_t& error_ND);

  void SetVerbose(Int_t lev) { fVerbose = lev ;}

  void SetEnergy(Float_t en) { fEnergy = en; }

  Long64_t Merge(TCollection* list);

  Double_t GetInputEvents(){return fInputEvents;}   // number of Input Events 
  Double_t GetPhysSelEvents(){return fPhysSelEvents;} // number of  Events after Physics Selection 
  Double_t GetBgEvents(){return fBgEvents;}      // number of background events 

  Double_t GetAllEvents(){return fAllEvents;}             // number of corrected events 
  Double_t GetAllEventsZRange(){return fAllEventsZRange;}       // number of corrected events in z range
  Double_t GetAllEventsZRangeMult1(){return fAllEventsZRangeMult1;}  // number of corrected events with multiplicity larger 1 in z range
  Double_t GetAllEventsInBin0ZRange(){return fAllEventsInBin0ZRange;} // number of corrected events in bin0 in z range
  Double_t GetTrigEffBin0(){return fTrigEffBin0;}                // trigger efficiency


protected:

  Int_t   fNbinsVz; // number of z bins in the vz histo
  Float_t fMinVz  ; // lowest Z
  Float_t fMaxVz  ; // highest Z

  Float_t fZRange; // max |Z| vertex to be considered

  Bool_t fIsMC; // True if processing MC
  
  Int_t fReferenceXS;                // index of reference cross section to be used to rescale process types in the calculation of the efficiency

  Int_t fVerbose;                    // Determines the ammount of printout

  Float_t fEnergy;                     // Beam energy in GeV. Defaults to 900.

  TH2F * fHistVzMCGen[kNProcs]    ;    // Vz distribution of generated events vs rec multiplicity
  TH2F * fHistVzMCRec[kNProcs]    ; 	// Vz distribution of reconstructed events vs rec multiplicity
  TH2F * fHistVzMCTrg[kNProcs]    ; 	// Vz distribution of triggered events vs rec multiplicity
  TH2F * fHistVzData              ; 	// Vz distribution of triggered events vs rec multiplicity    
  TH1F * fHistProcTypes           ;    // Number of evts for different Process types 

  TH1F * fHistStatBin0     ; // event stat histogram, created by physiscs selection; used in ComputeNint;
  TH1F * fHistStat         ; // event stat histogram, created by physiscs selection; used in ComputeNint;

  Double_t fInputEvents;   // number of Input Events 
  Double_t fPhysSelEvents; // number of  Events after Physics Selection 
  Double_t fBgEvents;            // number of background events 

  Double_t fAllEvents;     // number of corrected events 
  Double_t fAllEventsZRange;  // number of corrected events in z range
  Double_t fAllEventsZRangeMult1; // number of corrected events with multiplicity larger 1 in z range
  Double_t fAllEventsInBin0ZRange; // number of corrected events in bin0 in z range
  Double_t fTrigEffBin0;  // trigger efficiency in Bin0

  static const char * fProcLabel[] ; // labels of the different process types
  
  ClassDef(AliCollisionNormalization, 4);
    
private:
  AliCollisionNormalization(const AliCollisionNormalization&);
  AliCollisionNormalization& operator=(const AliCollisionNormalization&);
};

#endif
