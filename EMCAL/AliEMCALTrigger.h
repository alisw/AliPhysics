#ifndef ALIEMCALTrigger_H
#define ALIEMCALTrigger_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//___________________________________________________________
//  Class for trigger analysis.
//
//  -- Author: Gustavo Conesa & Yves Schutz (IFIC, SUBATECH, CERN)
//  Digits are grouped in TRU's (Trigger Units). A TRU consist of 384 cells 
//  ordered fNTRUPhi x fNTRUEta matrix. The algorithm searches all possible 
//  2x2 and nxn (n multiple of 4) crystal combinations per each TRU, adding the 
//  digits amplitude and finding the maximum.  It is found is maximum is isolated. 
//  Maxima are transformed in adc time samples. Each time bin is compared to the 
//  trigger threshold until it is larger and then, triggers are set. 
//  Thresholds need to be fixed. 
//  Last 2 modules are half size in Phi, I considered that the number 
//  of TRU is maintained for the last modules but final decision has not 
//  been taken. If different, then this must to be changed. 
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
//  tr->Print("");  //Print data members after calculation.
//
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, SUBATECH, CERN)
     
// --- ROOT system ---

class TClonesArray ;
#include "TMatrixD.h"

// --- AliRoot header files ---
#include "AliTriggerDetector.h"

class AliEMCALGeometry ;

class AliEMCALTrigger : public AliTriggerDetector {
  
 public:   

  AliEMCALTrigger() ; //  ctor
  AliEMCALTrigger(const AliEMCALTrigger & trig) ; // cpy ctor
  virtual ~AliEMCALTrigger() {}; //virtual dtor


  virtual void    CreateInputs(); //Define trigger inputs for Central Trigger Processor
  void            Print(const Option_t * opt ="") const ;  
  virtual void    Trigger();  //Make EMCAL trigger

  //assignment operator for coding convention
  const AliEMCALTrigger & operator = (const AliEMCALTrigger & ) {return *this;}

  //Getters
  Float_t  Get2x2MaxAmplitude()  const { return f2x2MaxAmp ; }
  Float_t  GetnxnMaxAmplitude()  const { return fnxnMaxAmp ; }
  Int_t    Get2x2CellPhi()       const { return f2x2CellPhi ; }
  Int_t    GetnxnCellPhi()       const { return fnxnCellPhi ; }
  Int_t    Get2x2CellEta()       const { return f2x2CellEta ; }
  Int_t    GetnxnCellEta()       const { return fnxnCellEta ; }
  Int_t    Get2x2SuperModule()   const { return f2x2SM ; }
  Int_t    GetnxnSuperModule()   const { return fnxnSM ; }

  Int_t *  GetADCValuesLowGainMax2x2Sum()  { return fADCValuesLow2x2; }
  Int_t *  GetADCValuesHighGainMax2x2Sum() { return fADCValuesHigh2x2; }
  Int_t *  GetADCValuesLowGainMaxnxnSum()  { return fADCValuesLownxn; }
  Int_t *  GetADCValuesHighGainMaxnxnSum() { return fADCValuesHighnxn; }

  Float_t  GetL0Threshold() const            { return fL0Threshold ; } 
  Float_t  GetL1JetLowPtThreshold()    const { return fL1JetLowPtThreshold ; }
  Float_t  GetL1JetMediumPtThreshold() const { return fL1JetMediumPtThreshold ; }
  Float_t  GetL1JetHighPtThreshold()   const { return fL1JetHighPtThreshold ; }

  Int_t    GetPatchSize()              const { return fPatchSize ; }
  Int_t    GetIsolPatchSize()          const { return fIsolPatchSize ; }

  Float_t  Get2x2AmpOutOfPatch()       const { return  f2x2AmpOutOfPatch ; }
  Float_t  GetnxnAmpOutOfPatch()       const { return  fnxnAmpOutOfPatch ; }
  Float_t  Get2x2AmpOutOfPatchThres()  const { return  f2x2AmpOutOfPatchThres ; }
  Float_t  GetnxnAmpOutOfPatchThres()  const { return  fnxnAmpOutOfPatchThres ; } 

  Bool_t   Is2x2Isol()                 const { return  fIs2x2Isol ; }
  Bool_t   IsnxnIsol()                 const { return  fIsnxnIsol ; }

  Bool_t   IsSimulation()              const { return fSimulation ; }
  Bool_t   IsIsolatedInSuperModule()   const { return fIsolateInSuperModule ; }

  //Setters
  void     SetDigitsList(TClonesArray * digits)          
   {fDigitsList  = digits ; }

  void     SetL0Threshold(Int_t amp)     
    {fL0Threshold            = amp; }
  void     SetL1JetLowPtThreshold(Int_t amp) 
    {fL1JetLowPtThreshold    = amp; } 
  void     SetL1JetMediumPtThreshold(Int_t amp) 
    {fL1JetMediumPtThreshold = amp; } 
  void     SetL1JetHighPtThreshold(Int_t amp)
    {fL1JetHighPtThreshold   = amp; }

  void SetPatchSize(Int_t ps)                {fPatchSize = ps ; }
  void SetIsolPatchSize(Int_t ps)          {fIsolPatchSize = ps ; }
  void Set2x2AmpOutOfPatchThres(Float_t th) { f2x2AmpOutOfPatchThres = th; }
  void SetnxnAmpOutOfPatchThres(Float_t th) { fnxnAmpOutOfPatchThres = th; }
  void SetSimulation(Bool_t sim )          {fSimulation = sim ; }
  void SetIsolateInSuperModule(Bool_t isol )          {fIsolateInSuperModule = isol ; }

 private:

  Bool_t IsPatchIsolated(Int_t iPatchType, const TClonesArray * ampmods, const Int_t imod, const Int_t mtru, const Float_t maxamp, const Int_t maxphi, const Int_t maxeta) ;
  
  void MakeSlidingCell(const TClonesArray * amptrus, const TClonesArray * timeRtrus,const Int_t supermod, TMatrixD *ampmax2, TMatrixD *ampmaxn) ; 
  
  void SetTriggers(const TClonesArray * amptrus,const Int_t iSM, const TMatrixD *ampmax2, const TMatrixD *ampmaxn, const AliEMCALGeometry * geom) ;
  

 private: 

  Float_t f2x2MaxAmp ;         //! Maximum 2x2 added amplitude (not overlapped) 
  Int_t   f2x2CellPhi ;        //! upper right cell, row(phi)   
  Int_t   f2x2CellEta ;        //! and column(eta)  
  Int_t   f2x2SM ;             //! Super Module where maximum is found
  Float_t fnxnMaxAmp ;         //! Maximum nxn added amplitude (overlapped)
  Int_t   fnxnCellPhi ;        //! upper right cell, row(phi)   
  Int_t   fnxnCellEta ;        //! and column(eta)
  Int_t   fnxnSM ;             //! Super Module where maximum is found

  Int_t*   fADCValuesHighnxn ; //! Sampled ADC high gain values for the nxn crystals amplitude sum
  Int_t*   fADCValuesLownxn  ; //! " low gain  " 
  Int_t*   fADCValuesHigh2x2 ; //! " high gain " 2x2 "
  Int_t*   fADCValuesLow2x2  ; //! " low gaing " "

  TClonesArray* fDigitsList ;  //Array of digits 

  Float_t fL0Threshold ;            //! L0 trigger energy threshold
  Float_t fL1JetLowPtThreshold ;    //! L1 Low pT trigger energy threshold
  Float_t fL1JetMediumPtThreshold ; //! L1 Medium pT trigger energy threshold
  Float_t fL1JetHighPtThreshold ;   //! L1 High pT trigger energy threshold

  Int_t   fNTRU;             //! Number of TRU per SuperModule (3)
  Int_t   fNTRUEta ;         //! Number of crystal rows per Eta in one TRU (3)
  Int_t   fNTRUPhi ;         //! Number of crystal rows per Phi in one TRU (1)
  Int_t   fNCellsPhi;        //! Number of rows in a TRU (24)
  Int_t   fNCellsEta;        //! Number of columns in a TRU (16)

  Int_t fPatchSize;          //! Trigger patch factor, to be multiplied to 2x2 cells
                             //  0 means 2x2, 1 means 4x4, 2 means 6x6 ...
  Int_t fIsolPatchSize ;     //  Isolation patch size, number of rows or columns to add to 
                             //  the 2x2 or nxn maximum amplitude patch. 
                             //  1 means a patch around max amplitude of 2x2 of 4x4 and around         
                             //  max ampl patch of 4x4 of 8x8 
    
  Float_t f2x2AmpOutOfPatch;      //  Amplitude in isolation cone minus maximum amplitude of the reference patch
  Float_t fnxnAmpOutOfPatch; 
  Float_t f2x2AmpOutOfPatchThres; //  Threshold to select a trigger as isolated on f2x2AmpOutOfPatch value
  Float_t fnxnAmpOutOfPatchThres; 
  Float_t fIs2x2Isol;             //  Patch is isolated if f2x2AmpOutOfPatchThres threshold is passed
  Float_t fIsnxnIsol ; 

  Bool_t  fSimulation ;           //! Flag to do the trigger during simulation or reconstruction
  Bool_t  fIsolateInSuperModule;  //! Flag to isolate trigger patch in SuperModule or in TRU acceptance

  ClassDef(AliEMCALTrigger,1)
} ;
    
    
#endif //ALIEMCALTrigger_H
    
