#ifndef TOF_H
#define TOF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set: TOF     //
////////////////////////////////////////////////
#include "TObject.h"

class TFile;
 
#include "AliDetector.h"
#include "AliHit.h"
#include "AliDigit.h" 
#include "AliTOFD.h"

class AliTOF : public AliDetector {

protected:
  Int_t   fIdSens;

public:
  Int_t   fNTof;
  Float_t fRmax;
  Float_t fRmin;
  Float_t fZlenA;
  Float_t fZlenB;
  Float_t fZlenC;
  Float_t fZtof;
   
  Float_t fStripLn;
  Float_t fSpace;
  Float_t fDeadBndZ;
  Float_t fDeadBndX;
  Float_t fXpad;
  Float_t fZpad;
  Float_t fGapA;
  Float_t fGapB;
  Float_t fOverSpc;
  Int_t   fNpadX;
  Int_t   fNpadZ;
  Int_t   fPadXStr;

  Int_t   fNStripA;
  Int_t   fNStripB;
  Int_t   fNStripC;

  Float_t fTimeRes;
  Float_t fChrgRes;
  
  Int_t   fPadXSector;
  Int_t   fNRoc;
  Int_t   fNFec;
  Int_t   fNTdc;
  Int_t   fNPadXRoc;


public:
  AliTOF();
  AliTOF(const char *name, const char *title);
  virtual        ~AliTOF() {}
  virtual void    AddHit(Int_t track, Int_t* vol, Float_t* hits);
  virtual void    AddDigit(Int_t*, Int_t*, Float_t*);
  virtual void    CreateGeometry();
  virtual void    CreateMaterials();
  virtual void    Init();
  virtual void    MakeBranch(Option_t*, char *file=0);
  virtual void    FinishEvent();
  virtual Int_t   IsVersion() const =0;
  Int_t           DistancetoPrimitive(Int_t px, Int_t py);
  virtual void    StepManager()=0;
  virtual void    TOFpc(Float_t, Float_t, Float_t,
                        Float_t, Float_t,Float_t) {}
  virtual void    DrawModule();
  virtual void    SDigits2Digits();
          void    Hits2Digits(Int_t evNumber=0);
          void    Digits2Raw (Int_t evNumber=0);
          void    Raw2Digits (Int_t evNumber=0);
  
private:
	Bool_t    CheckOverlap(Int_t*, Float_t*, Int_t);

  ClassDef(AliTOF,1)  // Time Of Flight base class
};
 
//___________________________________________
 
class AliTOFhit : public AliHit {
public:
  Int_t      fSector;  // number of sector 
  Int_t      fPlate;   // number of plate
  Int_t      fStrip;   // number of strip
  Int_t      fPad_x;   // number of pad along x
  Int_t      fPad_z;   // number of pad along z
  Float_t    fPx;      // px in TOF
  Float_t    fPy;      // py in TOF
  Float_t    fPz;      // pz in TOF
  Float_t    fPmom;    // P in TOF
  Float_t    fTof;     // Time of Flight
  Float_t    fDx;      // x of impact point in pad r.s.
  Float_t    fDy;      // y of impact point in pad r.s.
  Float_t    fDz;      // z of impact point in pad r.s.
  Float_t    fIncA;    // Incidence angle
  Float_t    fEdep;    // Energy lost in tof layer
 
public:
  AliTOFhit() {}
  AliTOFhit(Int_t shunt, Int_t track, Int_t* vol, 
            Float_t *hits);
  virtual ~AliTOFhit() {}

  inline  Int_t   GetSector() {return fSector;}
  inline  Int_t   GetPlate()  {return fPlate;}
  inline  Int_t   GetPad_x()  {return fPad_x;}
  inline  Int_t   GetPad_z()  {return fPad_z;}
  inline  Int_t   GetStrip()  {return (Int_t)(fPad_z*0.5);}
  inline  Float_t GetTof()    {return fTof;}
  inline  Float_t GetMom()    {return fPmom;}
  inline  Float_t GetDx()     {return fDx;}
  inline  Float_t GetDz()     {return fDz;}
  inline  Float_t GetIncA()   {return fIncA;}
  inline  Float_t GetEdep()   {return fEdep;}

  ClassDef(AliTOFhit,1)  // Hits for Time Of Flight
};

//_______________________________________________________

class AliTOFdigit : public AliDigit {

 public:
  Int_t   fSector;
  Int_t   fPlate;
  Int_t   fStrip;
  Int_t   fPad_x;
  Int_t   fPad_z;
  Float_t fTdc;
  Float_t fAdc;

 public:
  AliTOFdigit(){}
  AliTOFdigit(Int_t*, Int_t*, Float_t*);
  virtual ~AliTOFdigit(){}
  void            GetLocation(Int_t*);
  Int_t           GetTotPad();
  void            AddTrack(Int_t);

  inline  Float_t GetTdc()           {return fTdc;}
  inline  Float_t GetAdc()           {return fAdc;}
  inline  Int_t   GetSector()        {return fSector;}
  inline  void    SetTdc(Float_t TDC){fTdc = TDC;}
  inline  void    SetAdc(Float_t ADC){fAdc = ADC;}

  ClassDef(AliTOFdigit,2)  // Digits for Time Of Flight
};

#endif
