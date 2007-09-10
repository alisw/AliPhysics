#ifndef ALIPHOSTrigger_H
#define ALIPHOSTrigger_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//____________________________________________________________
//  Class for trigger analysis.
//
//  -- Author: Gustavo Conesa & Yves Schutz (IFIC, SUBATECH, CERN)
//  Digits are grouped in TRU's (Trigger Units). A TRU consist of 16x28 
//  crystals ordered fNTRUPhi x fNTRUZ matrix. The algorithm searches all possible 
//  2x2 and nxn (n multiple of 4) crystal combinations per each TRU, adding the 
//  digits amplitude and  finding the maximum. Iti is found is maximum is isolated.
//  Maxima are transformed in ADC time samples. Each time bin is compared to the trigger 
//  threshold until it is larger and then, triggers are set. Thresholds need to be fixed. 
//  Usage:
//
//  //Inside the event loop
//  AliEMCALTrigger *tr = new AliEMCALTrigger();//Init Trigger
//  tr->SetL0Threshold(100);
//  tr->SetL1JetLowPtThreshold(1000);
//  tr->SetL1JetMediumPtThreshold(10000);
//  tr->SetL1JetHighPtThreshold(20000);
//  ....
//  tr->Trigger(); //Execute Trigger
//  tr->Print(""); //Print data members after calculation.
//     
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
  virtual ~AliPHOSTrigger();

  virtual void    CreateInputs(); //Define trigger inputs for Central Trigger Processor
  void            Print(const Option_t * opt ="") const ;  
  virtual void    Trigger();  //Make PHOS trigger
  void    Trigger(const char * fileName);  //Make PHOS trigger

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

  void     GetCrystalPhiEtaIndexInModuleFromTRUIndex(Int_t itru, Int_t iphitru, Int_t ietatru,Int_t &ietaMod,Int_t &iphiMod) const ;

  Float_t  GetL0Threshold()            const {return fL0Threshold ; } 
  Float_t  GetL1JetLowPtThreshold()    const {return fL1JetLowPtThreshold ; }
  Float_t  GetL1JetMediumPtThreshold() const {return fL1JetMediumPtThreshold ; }
  Float_t  GetL1JetHighPtThreshold()   const {return fL1JetHighPtThreshold ; }

  Int_t    GetNTRU()                   const {return fNTRU ; }
  Int_t    GetNTRUZ()                  const {return fNTRUZ ; }
  Int_t    GetNTRUPhi()                const {return fNTRUPhi ; }
  
  Int_t    GetPatchSize()              const {return fPatchSize ; }
  Int_t    GetIsolPatchSize()          const {return fIsolPatchSize ; }

  Float_t  Get2x2AmpOutOfPatch()       const {return  f2x2AmpOutOfPatch; }
  Float_t  GetnxnAmpOutOfPatch()       const {return  fnxnAmpOutOfPatch; }
  Float_t  Get2x2AmpOutOfPatchThres()  const {return  f2x2AmpOutOfPatchThres; }
  Float_t  GetnxnAmpOutOfPatchThres()  const {return  fnxnAmpOutOfPatchThres; }

  Bool_t   Is2x2Isol()                 const {return  fIs2x2Isol; }
  Bool_t   IsnxnIsol()                 const {return  fIsnxnIsol; }

  Bool_t   IsSimulation()              const {return fSimulation ; }
  Bool_t   IsIsolatedInModule()        const {return fIsolateInModule ; }

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
  void     SetL1JetMediumPtThreshold(Int_t amp) 
    {fL1JetMediumPtThreshold = amp; } 
  void     SetL1JetHighPtThreshold(Int_t amp)
    {fL1JetHighPtThreshold = amp ; }

  void SetPatchSize(Int_t ps)               { fPatchSize = ps ; }
  void SetIsolPatchSize(Int_t ps)           { fIsolPatchSize = ps ; }
  void Set2x2AmpOutOfPatchThres(Float_t th) { f2x2AmpOutOfPatchThres = th; }
  void SetnxnAmpOutOfPatchThres(Float_t th) { fnxnAmpOutOfPatchThres = th; }
  void SetSimulation(Bool_t sim )           { fSimulation = sim ; }
  void SetIsolateInModule(Bool_t isol )     { fIsolateInModule = isol ; }

 private:

  AliPHOSTrigger & operator = (const AliPHOSTrigger & trig) ;//cpy assignment

  void FillTRU(const TClonesArray * digits, const AliPHOSGeometry * geom, TClonesArray * amptru, TClonesArray * ampmod, TClonesArray * timeRtru) const ;

  Bool_t IsPatchIsolated(Int_t iPatchType, const TClonesArray * ampmods, const Int_t imod, const Int_t mtru, const Float_t maxamp, const Int_t maxphi, const Int_t maxeta) ;

  void MakeSlidingCell(const TClonesArray * amptrus, const TClonesArray * timeRtrus, Int_t mod, TMatrixD &ampmax2, TMatrixD &ampmaxn) ;

  void SetTriggers(const TClonesArray * amptrus, Int_t iMod, const TMatrixD &ampmax2,const TMatrixD &ampmaxn) ;

  void DoIt(const char * fileName) ; 
 
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
 
  Float_t fL0Threshold ;             //! L0 trigger energy threshold
  Float_t fL1JetLowPtThreshold ;     //! L1 Low  pT trigger threshold
  Float_t fL1JetMediumPtThreshold ;  //! L1 Medium  pT trigger threshold
  Float_t fL1JetHighPtThreshold ;    //! L1 High pT trigger threshold

  Int_t   fNTRU ;                //! Number of TRUs per module
  Int_t   fNTRUZ ;               //! Number of crystal rows per Z in one TRU
  Int_t   fNTRUPhi ;             //! Number of crystal rows per Phi in one TRU
  Int_t   fNCrystalsPhi;         //! Number of rows in a TRU
  Int_t   fNCrystalsZ;           //! Number of columns in a TRU
  
  Int_t fPatchSize;              //! Trigger patch factor, to be multiplied to 2x2 cells
                                 //  0 means 2x2, 1 means 4x4, 2 means 6x6 ...
  Int_t fIsolPatchSize ;         //  Isolation patch size, number of rows or columns to add to 
                                 //  the 2x2 or nxn maximum amplitude patch. 
                                 //  1 means a patch around max amplitude of 2x2 of 4x4 and around         
                                 //  max ampl patch of 4x4 of 8x8 
    
  Float_t f2x2AmpOutOfPatch;      // Amplitude in isolation cone minus maximum amplitude of the reference patch
  Float_t fnxnAmpOutOfPatch; 
  Float_t f2x2AmpOutOfPatchThres; // Threshold to select a trigger as isolated on f2x2AmpOutOfPatch value
  Float_t fnxnAmpOutOfPatchThres; 
  Float_t fIs2x2Isol;             //Patch is isolated if f2x2AmpOutOfPatchThres threshold is passed
  Float_t fIsnxnIsol ; 
  
  Bool_t  fSimulation ;           //! Flag to do the trigger during simulation or reconstruction
  Bool_t  fIsolateInModule;       //! Flag to isolate trigger patch in Module or in TRU acceptance

  ClassDef(AliPHOSTrigger,4)
} ;


#endif //ALIPHOSTrigger_H
