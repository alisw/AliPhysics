////////////////////////////////////////////////
//					      //
//  Manager, hit and digit classes for TOF    //
//  Interfaces:
//  AliTOF                                    //
//  AliTOFhit                                 //
//  AliTOFdigit                               //
////////////////////////////////////////////////

#ifndef ALITOF_H
#define ALITOF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h"

class TFile;
 
#include "AliDetector.h"
#include "AliHit.h"
#include "AliDigit.h" 
//#include "AliTOFD.h"

class AliTOF : public AliDetector {

public:
  AliTOF();
  AliTOF(const char *name, const char *title);
  virtual        ~AliTOF() {}
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
  virtual void    MakeBranch(Option_t* option, char *file=0);
  virtual void    FinishEvent();
  virtual Int_t   IsVersion() const =0;
  Int_t           DistancetoPrimitive(Int_t px, Int_t py);
  virtual void    StepManager()=0;
  virtual void    TOFpc(Float_t xtof, Float_t ytof, Float_t zlenC,
                        Float_t zlenB, Float_t zlenA, Float_t ztof0){}
  virtual void    DrawModule();
  virtual void    SDigits2Digits();
          void    Hits2Digits(Int_t evNumber=0);
          void    Digits2Raw (Int_t evNumber=0);
          void    Raw2Digits (Int_t evNumber=0);

protected:
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
  Int_t   fIdSens;  // the unique numeric identifier for sensitive volume FPAD 
  
private:
	Bool_t    CheckOverlap(Int_t* vol, Float_t* digit, Int_t Track);

  ClassDef(AliTOF,1)  // Time Of Flight base class
};
 
//___________________________________________
 
class AliTOFhit : public AliHit {
  
public:
  AliTOFhit() {}
  AliTOFhit(Int_t shunt, Int_t track, Int_t* vol, 
            Float_t *hits);
  AliTOFhit(const AliTOFhit & hit) ;
  virtual ~AliTOFhit() {}
  // getters for AliTOFhit object
  Int_t   GetSector() const {return fSector;}
  Int_t   GetPlate() const {return fPlate;}
  Int_t   GetPadx() const {return fPadx;}
  Int_t   GetPadz() const {return fPadz;}
  Int_t   GetStrip() const {return fStrip;}
  Float_t GetTof() const {return fTof;}
  Float_t GetMom() const {return fPmom;}
  Float_t GetDx() const  {return fDx;}
  Float_t GetDz() const  {return fDz;}
  Float_t GetIncA() const {return fIncA;}
  Float_t GetEdep() const {return fEdep;}

protected:
  Int_t      fSector;  // number of sector 
  Int_t      fPlate;   // number of plate
  Int_t      fStrip;   // number of strip
  Int_t      fPadx;    // number of pad along x
  Int_t      fPadz;    // number of pad along z
// X, Y and Z coordinates of the hit are defined on mother class
// AliHit
  Float_t    fPx;      // px in TOF
  Float_t    fPy;      // py in TOF
  Float_t    fPz;      // pz in TOF
  Float_t    fPmom;    // P in TOF
  Float_t    fTof;     // Time of Flight
  Float_t    fDx;      // x of impact point in pad r.s.
  Float_t    fDy;      // y of impact point in pad r.s.
  Float_t    fDz;      // z of impact point in pad r.s.
  Float_t    fIncA;    // Incidence angle
  Float_t    fEdep;    // Energy lost in TOF sensitive layer

  ClassDef(AliTOFhit,1)  // Hits for Time Of Flight
};

//_______________________________________________________

class AliTOFdigit : public AliDigit {

 public:
  AliTOFdigit(){}
  AliTOFdigit(Int_t* tracks, Int_t* vol, Float_t* digit);
  AliTOFdigit(const AliTOFdigit & digit) ;
  virtual ~AliTOFdigit(){}
  void            GetLocation(Int_t* Loc);
  Int_t           GetTotPad();
  void            AddTrack(Int_t track);
  // getters for AliTOFdigit object 
  Float_t GetTdc()    const     {return fTdc;}
  Float_t GetAdc()    const     {return fAdc;}
  Int_t   GetSector() const     {return fSector;}
  // setters for AliTOFdigit object
  void    SetTdc(Float_t TDC){fTdc = TDC;}
  void    SetAdc(Float_t ADC){fAdc = ADC;}

protected:
  Int_t   fSector;  // number of sector
  Int_t   fPlate;   // number of plate
  Int_t   fStrip;   // number of strip
  Int_t   fPadx;    // number of pad along x
  Int_t   fPadz;    // number of pad along z
  Float_t fTdc;     // tdc values for digit
  Float_t fAdc;     // adc values for digit

  ClassDef(AliTOFdigit,2)  // Digits for Time Of Flight
};

#endif /* ALITOF_H */
