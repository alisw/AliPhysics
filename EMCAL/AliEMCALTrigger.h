#ifndef ALIEMCALTRIGGER_H
#define ALIEMCALTRIGGER_H
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
//  tr->SetL1GammaLowPtThreshold(1000);
//  tr->SetL1GammaMediumPtThreshold(10000);
//  tr->SetL1GammaHighPtThreshold(20000);
//  ....
//  tr->Trigger();  //Execute Trigger
//  tr->Print("");  //Print data members after calculation.
//
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, SUBATECH, CERN)
//* -- Author: Aleksei Pavlinov, WSU, Detroit, USA
// Nov 2, 2007
// One TRU card receives 96 analogue sums from 12 FEE cards.
// One sum is correcponding output from on module.      
// This patch has size 12x8 modules (24x16 modules).
// Each SM has 3 TRU cards.

// --- ROOT system ---

class TClonesArray ;
class TTree;

#include <TMatrixD.h>
#include <TArrayF.h>

// --- AliRoot header files ---
#include "AliTriggerDetector.h"

class TBrowser;
class AliEMCALGeometry ;
class TH2F;

class AliEMCALTrigger : public AliTriggerDetector {
  
 public:   

  AliEMCALTrigger() ; //  ctor
  virtual ~AliEMCALTrigger(); //virtual dtor


  virtual void    CreateInputs(); //Define trigger inputs for Central Trigger Processor
  void            Print(const Option_t * opt ="") const ;  
  virtual void    Trigger();  //Make EMCAL trigger

  //Getters
  Float_t  Get2x2MaxAmplitude()  const { return f2x2MaxAmp ; }
  Float_t  GetnxnMaxAmplitude()  const { return fnxnMaxAmp ; }
  Int_t    Get2x2ModulePhi()     const { return f2x2ModulePhi ; }
  Int_t    GetnxnModulePhi()     const { return fnxnModulePhi ; }
  Int_t    Get2x2ModuleEta()     const { return f2x2ModuleEta ; }
  Int_t    GetnxnModuleEta()     const { return fnxnModuleEta ; }
  Int_t    Get2x2SuperModule()   const { return f2x2SM ; }
  Int_t    GetnxnSuperModule()   const { return fnxnSM ; }

  Int_t *  GetADCValuesLowGainMax2x2Sum()  const { return fADCValuesLow2x2; }
  Int_t *  GetADCValuesHighGainMax2x2Sum() const { return fADCValuesHigh2x2; }
  Int_t *  GetADCValuesLowGainMaxnxnSum()  const { return fADCValuesLownxn; }
  Int_t *  GetADCValuesHighGainMaxnxnSum() const { return fADCValuesHighnxn; }

  Float_t  GetL0Threshold()              const { return fL0Threshold ; } 
  Float_t  GetL1GammaLowPtThreshold()    const { return fL1GammaLowPtThreshold ; }
  Float_t  GetL1GammaMediumPtThreshold() const { return fL1GammaMediumPtThreshold ; }
  Float_t  GetL1GammaHighPtThreshold()   const { return fL1GammaHighPtThreshold ; }

  Int_t    GetPatchSize()              const { return fPatchSize ; }
  Int_t    GetIsolPatchSize()          const { return fIsolPatchSize ; }

  Float_t  Get2x2AmpOutOfPatch()       const { return  f2x2AmpOutOfPatch ; }
  Float_t  GetnxnAmpOutOfPatch()       const { return  fnxnAmpOutOfPatch ; }
  Float_t  Get2x2AmpOutOfPatchThres()  const { return  f2x2AmpOutOfPatchThres ; }
  Float_t  GetnxnAmpOutOfPatchThres()  const { return  fnxnAmpOutOfPatchThres ; } 

  Bool_t   Is2x2Isol()                 const { return  fIs2x2Isol ; }
  Bool_t   IsnxnIsol()                 const { return  fIsnxnIsol ; }

  Bool_t    IsSimulation()              const { return fSimulation ; }
  Bool_t    IsIsolatedInSuperModule()   const { return fIsolateInSuperModule ; }
  Bool_t    GetTimeKey()                const { return fTimeKey;}
  TH2F*     GetJetMatrixE()             const { return fJetMatrixE;}
  Double_t  GetEmcalSumAmp()            const;
  
  Int_t     GetNJetThreshold()   const {return fNJetThreshold;}
  Double_t* GetL1JetThresholds() const {return fL1JetThreshold;}
  TMatrixD  GetAmpJetMax()       const {return fAmpJetMax;}

  void PrintJetMatrix() const;                  // *MENU*
  void PrintAmpTruMatrix(Int_t ind) const;      // *MENU*
  void PrintAmpSmMatrix(Int_t ind) const;       // *MENU*
  void PrintMatrix(const TMatrixD &mat) const;  // *MENU*
  Bool_t CheckConsistentOfMatrixes(const Int_t pri=0); // *MENU*


  //Setters
  void     SetDigitsList(TClonesArray * digits)          
   {fDigitsList  = digits ; }

  void     SetL0Threshold(Int_t amp)     
    {fL0Threshold            = amp; }
  void     SetL1GammaLowPtThreshold(Int_t amp) 
    {fL1GammaLowPtThreshold    = amp; } 
  void     SetL1GammaMediumPtThreshold(Int_t amp) 
    {fL1GammaMediumPtThreshold = amp; } 
  void     SetL1GammaHighPtThreshold(Int_t amp)
    {fL1GammaHighPtThreshold   = amp; }

  void SetPatchSize(Int_t ps)                {fPatchSize = ps ; }
  void SetIsolPatchSize(Int_t ps)          {fIsolPatchSize = ps ; }
  void Set2x2AmpOutOfPatchThres(Float_t th) { f2x2AmpOutOfPatchThres = th; }
  void SetnxnAmpOutOfPatchThres(Float_t th) { fnxnAmpOutOfPatchThres = th; }
  void SetSimulation(Bool_t sim )          {fSimulation = sim ; }
  void SetIsolateInSuperModule(Bool_t isol )          {fIsolateInSuperModule = isol ; }
  void SetTimeKey(Bool_t timeKey) {fTimeKey = timeKey;}
  void SetJetPatchSize(const Int_t patchSize) {fNJetPatchPhi = fNJetPatchEta = patchSize;}
  void SetJetParameters(const Int_t patchSize, Double_t* jetThreshold)
  { // unused now
    fNJetPatchPhi = fNJetPatchEta = patchSize; 
    fL1JetThreshold = jetThreshold;
  }
  void SetVZER0Multiplicity(Double_t mult) {fVZER0Mult = mult;}

  //
  virtual void Browse(TBrowser* b);
  virtual Bool_t  IsFolder() const {return kTRUE;}

  // Name of Jet trigger(s)
  Char_t* GetNameOfJetTrigger(const Int_t i) {return Form("%s_Th_%2.2i",fgNameOfJetTriggers.Data(),i);}
  static TString GetNameOfJetTriggers() {return fgNameOfJetTriggers;}
  static TString fgNameOfJetTriggers; //Name of jet triggers
  // Estimation on EMCal energy from VZERO multiplicity
  // 0.0153 is coefficient from adc to energy
  // Dec 4, 2007
  // 1  p0           2.52248e-02   3.24364e-05   9.29319e-01  -2.34036e-06
  static Double_t GetMeanEmcalEnergy(const Int_t mult) {return 2.52248e-02*Double_t(mult);}
  static Double_t GetMeanEmcalPatchEnergy(const Int_t mult, Int_t patchSize) 
  {return GetMeanEmcalEnergy(mult)*Double_t(patchSize)*Double_t(patchSize)/208.;}
 private:

  void FillTRU(const TClonesArray * digits, TClonesArray * ampmatrix, TClonesArray * ampmatrixsmod, TClonesArray * timeRmatrix); 

  Bool_t IsPatchIsolated(Int_t iPatchType, const TClonesArray * ampmods, const Int_t imod, const Int_t mtru, const Float_t maxamp, const Int_t maxphi, const Int_t maxeta) ;
  
  void MakeSlidingTowers(const TClonesArray * amptrus, const TClonesArray * timeRtrus,
  const Int_t supermod, TMatrixD &ampmax2, TMatrixD &ampmaxn) ; 
  
  void SetTriggers(const TClonesArray * amptrus,const Int_t iSM, const TMatrixD &ampmax2, const TMatrixD &ampmaxn) ;
  void GetTriggerInfo(TArrayF &triggerPosition, TArrayF &triggerAmplitudes) const; 
  // Jet staff
  void FillJetMatrixFromSMs(TClonesArray *ampmatrixsmod, TMatrixD * const jetMat, AliEMCALGeometry * const g); 
  // no timing information here
  void MakeSlidingPatch(const TMatrixD &jm, const Int_t nPatchSize, TMatrixD &ampJetMax);

 private: 
  AliEMCALGeometry *fGeom;    //!

  Float_t f2x2MaxAmp ;         //! Maximum 2x2 added amplitude (not overlapped) 
  Int_t   f2x2ModulePhi ;      //! upper right cell, row(phi)   
  Int_t   f2x2ModuleEta ;      //! and column(eta)  
  Int_t   f2x2SM ;             //! Super Module where maximum is found
  Float_t fnxnMaxAmp ;         //! Maximum nxn added amplitude (overlapped)
  Int_t   fnxnModulePhi ;      //! upper right cell, row(phi)   
  Int_t   fnxnModuleEta ;      //! and column(eta)
  Int_t   fnxnSM ;             //! Super Module where maximum is found

  Int_t*   fADCValuesHighnxn ; //! Sampled ADC high gain values for the nxn crystals amplitude sum
  Int_t*   fADCValuesLownxn  ; //! " low gain  " 
  Int_t*   fADCValuesHigh2x2 ; //! " high gain " 2x2 "
  Int_t*   fADCValuesLow2x2  ; //! " low gaing " "

  TClonesArray* fDigitsList;   //! Array of digits 

  Float_t fL0Threshold ;              // L0 trigger energy threshold
  Float_t fL1GammaLowPtThreshold ;    // L1 gamma Low pT trigger energy threshold
  Float_t fL1GammaMediumPtThreshold ; // L1 gamma Medium pT trigger energy threshold
  Float_t fL1GammaHighPtThreshold ;   // L1 gamma High pT trigger energy threshold

  Int_t fPatchSize;          // Trigger patch factor, to be multiplied to 2x2 cells
                             //  0 means 2x2, 1 means 4x4 (max size 4x4 now)
  Int_t fIsolPatchSize ;     //  Isolation patch size, number of rows or columns to add to 
                             //  the 2x2 or nxn maximum amplitude patch. 
                             //  1 means a patch around max amplitude of 2x2 of 4x4 and around         
                             //  max ampl patch of 4x4 of 8x8 
    
  Float_t f2x2AmpOutOfPatch;      //  Amplitude in isolation cone minus maximum amplitude of the reference 2x2 patch
  Float_t fnxnAmpOutOfPatch;      //  Amplitude in isolation cone minus maximum amplitude of the reference nxn patch
  Float_t f2x2AmpOutOfPatchThres; //  Threshold to select a trigger as isolated on f2x2AmpOutOfPatch value
  Float_t fnxnAmpOutOfPatchThres; //  Threshold to select a trigger as isolated on fnxnAmpOutOfPatch value
  Float_t fIs2x2Isol;             //  2x2 Patch is isolated if f2x2AmpOutOfPatchThres threshold is passed
  Float_t fIsnxnIsol ;            //  nxn Patch is isolated if fnxnAmpOutOfPatchThres threshold is passed


  Bool_t  fSimulation ;           // Flag to do the trigger during simulation or reconstruction
  Bool_t  fIsolateInSuperModule;  // Flag to isolate trigger patch in SuperModule or in TRU acceptance
  Bool_t  fTimeKey;               // Flag to take into account the digits time information  
  // 
  TClonesArray *fAmpTrus;         //! Array of amplides of TRU matrixes
  TClonesArray *fTimeRtrus;       //! Array of recent times  (unused now)
  TClonesArray *fAmpSMods;        //! Array of amplides of SM  matrixes
  // Information for EMCAL ESD
  TArrayF fTriggerPosition;       // Triggered patch position
  TArrayF fTriggerAmplitudes;     // Triggered patch amplitude
  // Jet staf
  Int_t     fNJetPatchPhi;       // size of jet pathch in phi(row)   direction  (nJetPatchPhi*4 module) 
  Int_t     fNJetPatchEta;       // size of jet pathch in eta(column) direction (nJetPatchEta*4 module)
  Int_t     fNJetThreshold;      // number of jet threshold
  Double_t  *fL1JetThreshold;    //[fNJetThreshold] array of L1 jet energy threshold (this is not Et)
  Double_t  fJetMaxAmp;          // Max amp from patch (fNJetPatchPhi*fNJetPatchEta)
  TMatrixD* fAmpJetMatrix;       //-> Jet trigger matrix : (nphi(17), neta(12))
  TH2F*     fJetMatrixE;         //-> temporary solution for getting coordinate informatin
  TMatrixD  fAmpJetMax;          // 6 elements
  // VZER0 
  Double_t  fVZER0Mult;              // multiplicity (V0A+V0c)
  
  const AliEMCALTrigger & operator = (const AliEMCALTrigger & ) ;
  AliEMCALTrigger(const AliEMCALTrigger & trig) ; // cpy ctor
  
  ClassDef(AliEMCALTrigger, 2)
} ;
    
    
#endif //ALIEMCALTRIGGER_H
    
