#ifndef ALIDETECTOR_H
#define ALIDETECTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <AliModule.h>
#include <AliHit.h>

class AliDetector : public AliModule {

public:

  // Creators - distructors
  AliDetector(const char* name, const char *title);
  AliDetector();
  AliDetector(const AliDetector &det) {det.Copy(*this);}
  virtual ~AliDetector();

  // Inline functions
          void  Copy(AliDetector &det) const;
  virtual int   GetNdigits() {return fNdigits;}
  virtual int   GetNhits()   {return fNhits;}
  TClonesArray *Digits() {return fDigits;}
  TClonesArray *Hits()   {return fHits;}
  TObjArray    *Points() {return fPoints;}
  Int_t         GetIshunt() {return fIshunt;}
  void          SetIshunt(Int_t ishunt) {fIshunt=ishunt;}
  AliDetector &operator=(const AliDetector &det) 
  {det.Copy(*this); return (*this);}
  
  // Other methods
  virtual void        Browse(TBrowser *b);
  virtual void        FinishRun();
  virtual void        LoadPoints(Int_t track);
  virtual void        MakeBranch(Option_t *opt=" ");
  virtual void        ResetDigits();
  virtual void        ResetHits();
  virtual void        ResetPoints();
  virtual void        SetTreeAddress();
  virtual void        SetTimeGate(Float_t gate) {fTimeGate=gate;}
  virtual Float_t     GetTimeGate() {return fTimeGate;}
  virtual void        StepManager() {}
  virtual void        DrawModule() {}
  virtual AliHit*     FirstHit(Int_t track);
  virtual AliHit*     NextHit();
  virtual void        SetBufferSize(Int_t bufsize=8000) {fBufferSize = bufsize;}  
 
  // Data members
protected:      
  
  Float_t       fTimeGate;    //Time gate in seconds

  Int_t         fIshunt;      //1 if the hit is attached to the primary
  Int_t         fNhits;       //Number of hits
  Int_t         fNdigits;     //Number of digits
  Int_t         fBufferSize;  //buffer size for Tree detector branches
  TClonesArray *fHits;        //List of hits for one track only
  TClonesArray *fDigits;      //List of digits for this detector
  TObjArray    *fPoints;      //Array of points for each track (all tracks in memory)

  ClassDef(AliDetector,1)  //Base class for ALICE detectors
};
#endif
