#ifndef ALIKINETRACKCUTS_H
#define ALIKINETRACKCUTS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisCuts.h"

class  TObject;
class  TList;

class AliKineTrackCuts : public AliAnalysisCuts
{

public:
                      AliKineTrackCuts(const Char_t* name = "AliKineTrackCuts", const Char_t* title = "");
  virtual            ~AliKineTrackCuts(){;}
  
              Bool_t  IsSelected(TObject* obj);
	      Bool_t  IsSelected(TList* /*list*/) {return kTRUE;}
            
              void    SetFinalParticles( Bool_t val=kTRUE )          { fOnlyFinalParticles = val; }
              void    SetPrimaryParticles( Bool_t val=kTRUE )        { fOnlyPrimary = val; }
  // track kinematic cut setters
	      void    SetPRange(Float_t r1=0, Float_t r2=1e10)       { fPMin=r1;   fPMax=r2;}
              void    SetPtRange(Float_t r1=0, Float_t r2=1e10)      { fPtMin=r1;  fPtMax=r2;}
              void    SetPxRange(Float_t r1=-1e10, Float_t r2=1e10)  { fPxMin=r1;  fPxMax=r2;}
              void    SetPyRange(Float_t r1=-1e10, Float_t r2=1e10)  { fPyMin=r1;  fPyMax=r2;}
              void    SetPzRange(Float_t r1=-1e10, Float_t r2=1e10)  { fPzMin=r1;  fPzMax=r2;}
              void    SetEtaRange(Float_t r1=-1e10, Float_t r2=1e10) { fEtaMin=r1; fEtaMax=r2;}
              void    SetRapRange(Float_t r1=-1e10, Float_t r2=1e10) { fRapMin=r1; fRapMax=r2;}
                  
protected:        
                  
           Bool_t    fOnlyFinalParticles;   // true => skip part with GetStatusCode()!=1
           Bool_t    fOnlyPrimary;          // Only Primary Particles
  // kinematics cuts
           Float_t   fPMin,   fPMax;        // definition of the range of the P
           Float_t   fPtMin,  fPtMax;       // definition of the range of the Pt
           Float_t   fPxMin,  fPxMax;       // definition of the range of the Px
           Float_t   fPyMin,  fPyMax;       // definition of the range of the Py
           Float_t   fPzMin,  fPzMax;       // definition of the range of the Pz
           Float_t   fEtaMin, fEtaMax;      // definition of the range of the eta
           Float_t   fRapMin, fRapMax;      // definition of the range of the y
  
  
  ClassDef(AliKineTrackCuts, 1)
};


#endif
