#ifndef ALIPHOSTRIGGERUTILS_H
#define ALIPHOSTRIGGERUTILS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// class to test trigger information of PHOS trigger
// and to model it in MC simulations
// 
//*-- Author: Dmitri Peressounko (RRC "KI")

// --- ROOT system ---
#include "TNamed.h"
class TF1 ;
class TH2I;

// --- AliRoot header files ---
class AliVCluster ;
class AliVEvent ;
class AliPHOSGeoUtils;

class AliPHOSTriggerUtils : public TNamed {

public: 

  AliPHOSTriggerUtils() ;
  AliPHOSTriggerUtils(const Text_t* name, const Text_t* title="") ;
  AliPHOSTriggerUtils(const AliPHOSTriggerUtils & utils) ;
  
  virtual ~AliPHOSTriggerUtils(void){} ; 
  AliPHOSTriggerUtils & operator = (const AliPHOSTriggerUtils  & rvalue) ;

  void ForseUsingRun(Int_t run){fRun=run; fFixedRun=kTRUE; }  //Do not get run number from header, instead use this one
  
  void SetEvent(AliVEvent * event); //sets ref. to current event; inits class for new run if necessary
  
  Int_t IsFiredTrigger(AliVCluster * clu) ; //Returns bits if this cluster fired PHOS trigger in event: L0, L1low, L1med, L1high
  
  Int_t IsFiredTriggerMC(AliVCluster * clu) ;  //For MC simulations without detailed PHOS trigger 
                                                    //Returns (L0, L1low, L1med, L1high) if this cluster should fire PHOS trigger
                                                    //according to parameterization of turn-on curve and trigger bad map
                                                    
  void ReadTriggerParams(const char * filename); //Read trigger info (bad map, parameterizations) from specific file, ignore OADB        
  
  void SetTileOffset(Int_t left=-1, Int_t right=3){fbdrL=left;fbdrR=right ;}
  
  Bool_t TestBadMap(Int_t mod, Int_t ix,Int_t iz) ; //Check if trigger is in good/active region
  Int_t  WhichDDL(Int_t module, Int_t cellx) ;      //Calculates DDL number from (offline) module number and x cell coordinate
  
protected:
  
  void InitForRun(Int_t run) ; //read trigger bad map for this run from OADB. Should be called once per run
  Int_t FindBranch(Int_t nX, Int_t nZ) ; //Calculate number of PHOS branch
  Double_t TriggerProbabilityLHC13bcdef(Double_t eClu, Int_t module) ; //Parameterization of turn-on curve in LHC13bcdef
  Double_t TriggerProbabilityLHC16qrst(Double_t x, Int_t ddl);
  Double_t TriggerL1ProbabilityLHC16qrst(Double_t x, Int_t ddl);
  Double_t TriggerProbabilityLHC17pq(Double_t x, Int_t ddl);
  Double_t TriggerProbability(Double_t eClu, Int_t module, Int_t triggerBit) ; //Parameterization of turn-on curves
  
  
private:
  Int_t fbdrL ;  //Offset between cell and trigger 4x4 tile
  Int_t fbdrR ;  //in left/right (top/bottom) directions              
  
  Int_t fRun ;         //current run number (-1 not set yet, -2 use input file set with ReadTriggerParams)  
  Bool_t fFixedRun;    //do not read runnumber from header
  AliVEvent * fEvent ; //! Ref to current ESD/AOD event
  
  //Trigger bad map for 5 modules
  TH2I * fPHOSBadMap[5] ;
  AliPHOSGeoUtils* fGeom ;  //! PHOS geometry
  
  //Parameterization of turn-on curves
  TF1 * fTOCL0[5][8];       //parameterization of turn-on curves for 5 modules and 8 branches/module
  TF1 * fTOCL1low[5][8];    //parameterization of turn-on curves for 5 modules and 8 branches/module
  TF1 * fTOCL1med[5][8];    //parameterization of turn-on curves for 5 modules and 8 branches/module
  TF1 * fTOCL1high[5][8];   //parameterization of turn-on curves for 5 modules and 8 branches/module
  
  Int_t fNtrg4x4 ;       //! Number of triggers in current event
  Int_t fTrMod4x4[1000]; //! trigger modules
  Int_t fTrX4x4[1000] ;  //! trigger X coordinates
  Int_t fTrZ4x4[1000] ;  //! trigger Z coordinates
  Int_t fTrK4x4[1000] ;  //! trigger Kind,bits for L0, L1low, L1med, L1high
  
  
  ClassDef(AliPHOSTriggerUtils,2)       // PHOS trigger analysis class 

} ;

#endif // AliPHOSTRIGERUTILS_H

