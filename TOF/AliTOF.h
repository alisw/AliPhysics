////////////////////////////////////////////////
//					      //
//  Manager classe for TOF                    //
//  Interface :                               //
//  AliTOF                                    //
//  Associations between TOF related objects  //
//  are defined here                          //
// -- Authors: Pierella, Seganti, Vicinanza   //
//    (Bologna and Salerno University)        //
//                                            //
////////////////////////////////////////////////

#ifndef ALITOF_H
#define ALITOF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

class TFile;
class TDirectory;
class TString ;  
class TTask ;
class TFolder ;

#include "TObject.h"
#include "TTree.h" 

#include "AliDetector.h"

#include "AliTOFSDigitizer.h"
#include "AliTOFGeometry.h"


class AliTOF : public AliDetector {
public:
  AliTOF(); 
  AliTOF(const char *name, const char *title, Option_t *option="noTimeZero");
  virtual ~AliTOF() ;
// getters for AliTOF object status
  //Float_t GetTimeRes() const {return fTimeRes;};
  //Float_t GetChrgRes() const {return fChrgRes;};

  virtual void    SetTreeAddress();
  virtual void    AddHit(Int_t track, Int_t* vol, Float_t* hits);
  virtual void    AddT0Hit(Int_t track, Int_t* vol, Float_t* hits);
  virtual void    AddDigit(Int_t* tracks, Int_t* vol, Float_t* digits);
  virtual void    AddSDigit(Int_t tracknum, Int_t* vol, Float_t* digits);
  virtual void    CreateGeometry();
  virtual void    CreateMaterials();
  virtual void    Init();
  //virtual void    MakeBranch(Option_t* option, const char *file=0);
  virtual void    MakeBranch(Option_t *opt=" ");
  virtual void    Makehits(Bool_t hits=1);
  virtual void    FinishEvent();
  virtual Int_t   IsVersion() const =0;
  Int_t           DistancetoPrimitive(Int_t px, Int_t py) const;
  virtual void    StepManager()=0;
  virtual void    TOFpc(Float_t /*xtof*/, Float_t /*ytof*/, Float_t /*zlenC*/,
                        Float_t /*zlenB*/, Float_t /*zlenA*/, Float_t /*ztof0*/){}
  virtual void    DrawModule() const;
  virtual void    DrawDetectorModules()=0;
  virtual void    DrawDetectorStrips()=0;
  //virtual void   DrawDetectorModulesinFrame()=0;
  //virtual void   DrawDetectorStripsinFrame()=0;
          void    CreateTOFFolders();
  Bool_t    CheckOverlap(Int_t* vol, Float_t* digit, Int_t Track);
  //virtual void    Hits2Digits();   
  virtual void    Hits2SDigits();
  virtual void    Hits2SDigits(Int_t evNumber1, Int_t evNumber2);
  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const; 
  virtual void    Digits2Reco() {};
          void    Digits2Raw ();
          void    Raw2Digits (){};
  virtual void    ResetHits();
  virtual void    ResetDigits();
  virtual void    ResetSDigits();
  TClonesArray *SDigits() const {return fSDigits;}
  TClonesArray *ReconParticles() const {return fReconParticles;}
  void RecreateSDigitsArray();
  void CreateSDigitsArray();
  AliTOFGeometry *GetGeometry() const { return fTOFGeometry; }; 

protected:
  TFolder* fFGeom ;       //  Folder that holds the Geometry definition
  TTask*   fDTask ;       //  TOF Digitizer container
  TTask*   fReTask;       //  TOF Reconstructioner container
  TClonesArray* fSDigits; //! List of summable digits
  Int_t   fNSDigits;      //! Number of sdigits
  TClonesArray* fReconParticles; // List of reconstructed particles

  //Float_t fGapA;     //  Gap beetween tilted strip in A-type plate
  //Float_t fGapB;     //  Gap beetween tilted strip in B-type plate

  //Float_t fTimeRes;  //  Time resolution of the TOF
  //Float_t fChrgRes;  //  Charge resolution of ADC

  Int_t   fIdSens;     // The unique identifier for sensitive volume FPAD 

  Bool_t  fTZero;      // Flag indicating if T0 is used
  AliTOFGeometry *fTOFGeometry; //The TOF Geometry parameters
 
private:

  ClassDef(AliTOF,7)  // Time Of Flight base class
};
 
#endif /* ALITOF_H */
