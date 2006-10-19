#ifndef ALIPHOSTrigger_H
#define ALIPHOSTrigger_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id $ */
/* $Log $ */

//____________________________________________________________
//  Class for trigger analysis.
//  Digits are grouped in TRU's (Trigger Units). A TRU consist of 16x28 
//  crystals ordered fNTRUPhi x fNTRUZ. The algorithm searches all possible 
//  2x2 and nxn  (n multiple of 4) crystal combinations per each TRU, adding the 
//  digits amplitude and finding the maximum. Maxima are transformed in ADC 
//  time samples.  Each time bin is compared to the trigger threshold until it is larger 
//  and then, triggers are set. Thresholds need to be fixed. 
//  Usage:
//
//  //Inside the event loop
//  AliPHOSTrigger *tr = new AliPHOSTrigger();//Init Trigger
//  tr->SetL0Threshold(100);
//  tr->SetL1JetLowPtThreshold(1000);
//  tr->SetL1JetHighPtThreshold(20000);
//  tr->Trigger(); //Execute Trigger
//  tr->Print("");  //Print results
//
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, SUBATECH, CERN)
     
// --- ROOT system ---

class TClonesArray ;
#include "TMatrixD.h"

// --- AliRoot header files ---
#include "AliTriggerDetector.h"

class AliPHOSGeometry ;

class AliPHOSTrigger : public AliTriggerDetector {
  
 public: 
  
  AliPHOSTrigger() ; //  ctor
  AliPHOSTrigger(const AliPHOSTrigger & trig) ; // cpy ctor
  virtual ~AliPHOSTrigger() {}; //virtual dtor


  virtual void    CreateInputs(); //Define trigger inputs for Central Trigger Processor
  void            Print(const Option_t * opt ="") const ;  
  virtual void    Trigger();  //Make PHOS trigger

  //Getters
  Float_t  Get2x2MaxAmplitude()  const {return f2x2MaxAmp ; }
  Float_t  GetnxnMaxAmplitude()  const {return fnxnMaxAmp ; }
  Int_t    Get2x2CrystalPhi()    const {return f2x2CrystalPhi ; }
  Int_t    GetnxnCrystalPhi()    const {return fnxnCrystalPhi ; }
  Int_t    Get2x2CrystalEta()    const {return f2x2CrystalEta ; }
  Int_t    GetnxnCrystalEta()    const {return fnxnCrystalEta ; }
  Int_t    Get2x2SuperModule()   const {return f2x2SM ; }
  Int_t    GetnxnSuperModule()   const {return fnxnSM ; }

  Int_t *  GetADCValuesLowGainMax2x2Sum()  {return fADCValuesLow2x2; }
  Int_t *  GetADCValuesHighGainMax2x2Sum() {return fADCValuesHigh2x2; }
  Int_t *  GetADCValuesLowGainMaxnxnSum()  {return fADCValuesLownxn; }
  Int_t *  GetADCValuesHighGainMaxnxnSum() {return fADCValuesHighnxn; }

  void GetCrystalPhiEtaIndexInModuleFromTRUIndex(Int_t itru, Int_t iphitru, Int_t ietatru,Int_t &ietaMod,Int_t &iphiMod, const AliPHOSGeometry *geom) const ;

  Float_t  GetL0Threshold()          const {return fL0Threshold ; } 
  Float_t  GetL1JetLowPtThreshold()  const {return fL1JetLowPtThreshold ; }
  Float_t  GetL1JetHighPtThreshold() const {return fL1JetHighPtThreshold ; }

  Int_t    GetNTRU()    const  {return fNTRU ; }
  Int_t    GetNTRUZ()   const  {return fNTRUZ ; }
  Int_t    GetNTRUPhi() const  {return fNTRUPhi ; }
  
  Float_t  GetPatchSize() const  {return fPatchSize ; }
  Bool_t   IsSimulation() const {return fSimulation ; }

  //Setters

  void     SetDigitsList(TClonesArray * digits)          
   {fDigitsList  = digits ; }


  void     SetNTRU(Int_t ntru)             {fNTRU     = ntru ; }
  void     SetNTRUZ(Int_t ntru)            {fNTRUZ    = ntru ; }
  void     SetNTRUPhi(Int_t ntru)          {fNTRUPhi  = ntru ; } 

  void     SetL0Threshold(Int_t amp)         
    {fL0Threshold          = amp ; }
  void     SetL1JetLowPtThreshold(Int_t amp) 
    {fL1JetLowPtThreshold  = amp ; } 
  void     SetL1JetHighPtThreshold(Int_t amp)
    {fL1JetHighPtThreshold = amp ; }

  void SetPatchSize(Int_t ps)                {fPatchSize = ps ; }
  void SetSimulation(Bool_t sim )          {fSimulation = sim ; }

 private:

  AliPHOSTrigger & operator = (const AliPHOSTrigger & trig) ;//cpy assignment

  void FillTRU(const TClonesArray * digits, const AliPHOSGeometry * geom, TClonesArray * amptru, TClonesArray * timeRtru) const ;

  void MakeSlidingCell(const TClonesArray * amptrus, const TClonesArray * timeRtrus, Int_t mod, TMatrixD *ampmax2, TMatrixD *ampmaxn, const AliPHOSGeometry *geom) ;

  void SetTriggers(Int_t iMod, const TMatrixD *ampmax2,const TMatrixD *ampmaxn, const AliPHOSGeometry *geom) ;

 private: 

  Float_t f2x2MaxAmp ;     //! Maximum 2x2 added amplitude (not overlapped) 
  Int_t   f2x2CrystalPhi ; //! upper right cell, row(phi)   
  Int_t   f2x2CrystalEta ; //! and column(eta) 
  Int_t   f2x2SM ;         //! Module where maximum is found
  Float_t fnxnMaxAmp ;     //! Maximum nxn added amplitude (overlapped)
  Int_t   fnxnCrystalPhi ; //! upper right cell, row(phi)   
  Int_t   fnxnCrystalEta ; //! and column(eta)
  Int_t   fnxnSM ;         //! Module where maximum is found

  Int_t*   fADCValuesHighnxn ; //! Sampled ADC high gain values for the nxn crystals amplitude sum
  Int_t*   fADCValuesLownxn  ; //! " low gain  " 
  Int_t*   fADCValuesHigh2x2 ; //! " high gain " 2x2 "
  Int_t*   fADCValuesLow2x2  ; //! " low gaing " "

  TClonesArray* fDigitsList ;  //Array of digits 
 
  Float_t fL0Threshold ;          //! L0 trigger energy threshold
  Float_t fL1JetLowPtThreshold ;  //! L1 Low  pT trigger threshold
  Float_t fL1JetHighPtThreshold ; //! L1 High pT trigger threshold

  Int_t   fNTRU ;                 //! Number of TRUs per module
  Int_t   fNTRUZ ;                //! Number of crystal rows per Z in one TRU
  Int_t   fNTRUPhi ;              //! Number of crystal rows per Phi in one TRU
  Int_t fPatchSize;               //! Trigger patch factor, to be multiplied to 2x2 cells
                                          // 0 means 2x2, 1 means nxn, 2 means 8x8 ...
  Bool_t  fSimulation ;           //! Flag to do the trigger during simulation or reconstruction
  ClassDef(AliPHOSTrigger,4)
} ;


#endif //ALIPHOSTrigger_H
