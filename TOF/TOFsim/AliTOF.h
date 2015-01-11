#ifndef ALITOF_H
#define ALITOF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//					      //
//  Manager class for TOF                     //
//  Interface :                               //
//  AliTOF                                    //
//  Associations between TOF related objects  //
//  are defined here                          //
// -- Authors: Pierella, Seganti, Vicinanza   //
//    (Bologna and Salerno University)        //
//                                            //
////////////////////////////////////////////////

#include "AliDetector.h"

#include "AliTOFTrigger.h"
#include "AliTOFDDLRawData.h"

class TDirectory;
class TFile;
class TFolder ;
class TString ;  

class AliTOFGeometry;

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
  virtual void    AddDigit(Int_t* /*tracks*/, Int_t* /*vol*/) {};
  virtual void    AddDigit(Int_t* tracks, Int_t* vol, Int_t* digits);
  virtual void    AddSDigit(Int_t tracknum, Int_t* vol, Int_t* digits);
  virtual void    CreateGeometry();
  virtual void    CreateMaterials(){};
  virtual void    Init();
  //virtual void    MakeBranch(Option_t* option, const char *file=0);
  virtual void    MakeBranch(Option_t *opt=" ");
  virtual void    Makehits(Bool_t hits=1);
  virtual void    FinishEvent();
  virtual Int_t   IsVersion() const =0;
  virtual void    StepManager()=0;
  virtual void    TOFpc(Float_t /*xtof*/, Float_t /*ytof*/, Float_t /*zlenC*/,
                        Float_t /*zlenB*/, Float_t /*zlenA*/, Float_t /*ztof0*/){};
  virtual void    TOFpc(Float_t /*xtof*/,  Float_t /*ytof*/, Float_t /*zlenA*/,
			Float_t /*zlenB*/){};
  virtual void    TOFpc(Float_t /*xtof*/,  Float_t /*ytof*/, Float_t /*zlenA*/){};
          void    CreateTOFFolders();
  Bool_t    CheckOverlap(const Int_t * const vol, Int_t* digit, Int_t Track);
  //virtual void    Hits2Digits();   
  virtual void    Hits2SDigits();
  virtual void    Hits2SDigits(Int_t evNumber1, Int_t evNumber2);
  virtual AliDigitizer* CreateDigitizer(AliDigitizationInput* digInput) const; 
  virtual void    Digits2Reco () {};
          void    Digits2Raw  ();
	  void    Raw2Digits  () {};
	  void    Raw2Digits  (AliRawReader* rawReader);
	  Bool_t  Raw2SDigits (AliRawReader* rawReader);
  virtual void    ResetHits   ();
  virtual void    ResetDigits ();
  virtual void    ResetSDigits();
  TClonesArray *SDigits() const {return fSDigits;}
  TClonesArray *ReconParticles() const {return fReconParticles;}
  void RecreateSDigitsArray();
  void CreateSDigitsArray();
  virtual void   SetTOFSectors(Int_t * const sectors);
  virtual void   GetTOFSectors(Int_t *sectors) const;
  virtual void   SetTOFHoles(Bool_t holes) { fTOFHoles = holes; };
  virtual Bool_t GetTOFHoles() const { return fTOFHoles; };
  AliTOFGeometry *GetGeometry() const { return fTOFGeometry; }; 

  // Trigger
  virtual AliTriggerDetector* CreateTriggerDetector() const
  	{return new AliTOFTrigger();}

protected:
  TFolder* fFGeom ;       //  Folder that holds the Geometry definition
  TClonesArray* fSDigits; //! List of summable digits
  Int_t   fNSDigits;      //! Number of sdigits
  TClonesArray* fReconParticles; // List of reconstructed particles

  //Float_t fGapA;     //  Gap beetween tilted strip in A-type plate
  //Float_t fGapB;     //  Gap beetween tilted strip in B-type plate

  //Float_t fTimeRes;  //  Time resolution of the TOF
  //Float_t fChrgRes;  //  Charge resolution of ADC

  Int_t   fIdSens;     // The unique identifier for sensitive volume FPAD 

  Bool_t  fTZero;      // Flag indicating if T0 is used
  Int_t fTOFSectors[18]; // Selecting TOF Sectors to be simulated
  Bool_t fTOFHoles; // Selecting geometry with and w/o holes
  AliTOFGeometry *fTOFGeometry; //The TOF Geometry parameters

  AliTOFDDLRawData fTOFRawWriter; // AliTOFDDLRawData variable
 
private:
  AliTOF(const AliTOF &source); // copy constructor
  AliTOF& operator=(const AliTOF &source); // ass. op.

  ClassDef(AliTOF,12)  // Time Of Flight base class
};
 
#endif /* ALITOF_H */
