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
#include <iostream.h>

class AliTOF : public AliDetector {
public:
  AliTOF(); 
  AliTOF(const char *name, const char *title);
//  virtual        ~AliTOF() {} 
  virtual ~AliTOF() ;
// getters for AliTOF object status
  Int_t GetNStripA() const {return fNStripA;}
  Int_t GetNStripB() const {return fNStripB;}
  Int_t GetNStripC() const {return fNStripC;}
  Int_t GetNpadX()   const {return fNpadX;}
  Int_t GetNpadZ()   const {return fNpadZ;}
  Int_t GetPadXStr() const {return fPadXStr;}

  virtual void    AddHit(Int_t track, Int_t* vol, Float_t* hits);
  virtual void    AddDigit(Int_t* tracks, Int_t* vol, Float_t* digits);
  virtual void    CreateGeometry();
  virtual void    CreateMaterials();
  virtual void    Init();
  virtual void    MakeBranch(Option_t* option, const char *file=0);
  virtual void    Makehits(Bool_t hits=1);
  virtual void    FinishEvent();
  virtual Int_t   IsVersion() const =0;
  Int_t           DistancetoPrimitive(Int_t px, Int_t py) const;
  virtual void    StepManager()=0;
  virtual void    TOFpc(Float_t xtof, Float_t ytof, Float_t zlenC,
                        Float_t zlenB, Float_t zlenA, Float_t ztof0){}
  virtual void    DrawModule() const;
          void    CreateTOFFolders();
  virtual void    SDigits2Digits();
  virtual void    Hits2Digits();   
  virtual void    Hits2SDigits(){cout << "AliTOF::Hits2SDigits() dummy function called" << endl;}
  virtual void    Digits2Reco() {cout << "AliTOF::Digits2Reco()  dummy function called" << endl;}
          void    Digits2Raw (Int_t evNumber=0);
          void    Raw2Digits (Int_t evNumber=0);

protected:
  TFolder* fFGeom ;       //  Folder that holds the Geometry definition
  TTask*   fDTask ;       //  TOF Digitizer container
  TTask*   fReTask;       //  TOF Reconstructioner container
  TClonesArray* fSDigits; // List of summable digits
  Int_t   fNTof;  // number of TOF sectors
  Float_t fRmax;  // upper bound for radial extension of TOF detector
  Float_t fRmin;  // lower bound for radial extension of TOF detector
  Float_t fZlenA; // length along z-axis of type A module 
  Float_t fZlenB; // length along z-axis of type B module
  Float_t fZlenC; // length along z-axis of type C module
  Float_t fZtof;  // total semi-length of TOF detector
  Float_t fStripLn;  //  Strip Length
  Float_t fSpace;    //  Space Beetween the strip and the bottom of the plate
  Float_t fDeadBndZ; //  Dead Boundaries of a Strip along Z direction (width)
  Float_t fDeadBndX; //  Dead Boundaries of a Strip along X direction (length)
  Float_t fXpad;  //  X size of a pad
  Float_t fZpad;  //  Z size of a pad
  Float_t fGapA;  //  Gap beetween tilted strip in A-type plate
  Float_t fGapB;  //  Gap beetween tilted strip in B-type plate
  Float_t fOverSpc; // Space available for sensitive layers in radial direction
  Int_t   fNpadX; // Number of pads in a strip along the X direction
  Int_t   fNpadZ; // Number of pads in a strip along the Z direction
  Int_t   fPadXStr; // Number of pads per strip
  Int_t   fNStripA; // number of strips in A type module
  Int_t   fNStripB; // number of strips in B type module
  Int_t   fNStripC; // number of strips in C type module
  Float_t fTimeRes; // time resolution of the TOF
  Float_t fChrgRes; // charge resolution of ADC
  Int_t   fPadXSector; // number of pads per sector
  Int_t   fNRoc;       // number of ROC
  Int_t   fNFec;       // number of FEC
  Int_t   fNTdc;       // number of TDC
  Int_t   fNPadXRoc;   // number of pads for each ROC
  Int_t   fIdSens;     // the unique numeric identifier for sensitive volume FPAD 

private:
  Bool_t    CheckOverlap(Int_t* vol, Float_t* digit, Int_t Track);

  ClassDef(AliTOF,2)  // Time Of Flight base class
};
 
#endif /* ALITOF_H */
