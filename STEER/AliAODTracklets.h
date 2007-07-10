/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD class to store tracklets
//     Author: Jan Fiete Grosse-Oetringhaus, CERN
//     Class created from AliMultiplicity
//-------------------------------------------------------------------------

#ifndef ALIAODTRACKLETS_H
#define ALIAODTRACKLETS_H

#include <TNamed.h>

class AliAODTracklets : public TNamed 
{
 public:
  AliAODTracklets();
  AliAODTracklets(const char* name, const char* title);

  virtual ~AliAODTracklets();

  void CreateContainer(Int_t nTracks);
  void DeleteContainer();

  Bool_t SetTracklet(Int_t pos, Float_t theta, Float_t phi, Float_t deltaPhi, Int_t label);

  Int_t GetNumberOfTracklets() const { return fNTracks; }
  inline Float_t GetTheta(Int_t i) const;
  inline Float_t GetPhi(Int_t i) const;
  inline Float_t GetDeltaPhi(Int_t i) const;
  inline Int_t   GetLabel(Int_t i) const;

 protected:
  Int_t    fNTracks;      // Number of tracklets
  Float_t *fTheta;        //[fNTracks] array with theta values
  Float_t *fPhi;          //[fNTracks] array with phi values
  Float_t *fDeltaPhi;     //[fNTracks] array with delta phi values
  Int_t   *fLabels;       //[fNTracks] array with labels of tracklets

 private:
  AliAODTracklets(const AliAODTracklets& evt); 
  AliAODTracklets& operator=(const AliAODTracklets& evt);

  ClassDef(AliAODTracklets, 1);
};

Float_t AliAODTracklets::GetTheta(Int_t i) const 
{ 
  if (i>=0 && i<fNTracks) 
  {
    return fTheta[i];
  }
  else 
    Error("GetTheta","Invalid track number %d",i); return -9999.;
}

Float_t AliAODTracklets::GetPhi(Int_t i) const 
{ 
  if (i>=0 && i<fNTracks) 
  {
    return fPhi[i];
  }
  else 
    Error("GetPhi","Invalid track number %d",i); return -9999.;
}

Float_t AliAODTracklets::GetDeltaPhi(Int_t i) const 
{
  if (i>=0 && i<fNTracks) 
  {
    return fDeltaPhi[i];
  }
  else 
    Error("GetDeltaPhi","Invalid track number %d",i); return -9999.;
}

Int_t AliAODTracklets::GetLabel(Int_t i) const 
{
  if (i>=0 && i<fNTracks) 
  {
    return fLabels[i];
  }
  else 
    Error("GetLabel","Invalid track number %d",i); return -9999;
}

#endif
