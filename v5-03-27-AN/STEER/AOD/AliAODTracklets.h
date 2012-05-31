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
  AliAODTracklets(const AliAODTracklets& evt); 
  AliAODTracklets& operator=(const AliAODTracklets& evt);

  virtual ~AliAODTracklets();

  void CreateContainer(Int_t nTracks);
  void DeleteContainer();

  Bool_t SetTracklet(Int_t pos, Double32_t theta, Double32_t phi, Double32_t deltaPhi, Int_t labelL1, Int_t labelL2);

  Int_t GetNumberOfTracklets() const { return fNTracks; }
  inline Double32_t GetTheta(Int_t i) const;
  inline Double32_t GetPhi(Int_t i) const;
  inline Double32_t GetDeltaPhi(Int_t i) const;
  inline Int_t   GetLabel(Int_t i, Int_t layer) const;
  inline void    SetLabel(Int_t i, Int_t layer,Int_t label);

 protected:
  Int_t      fNTracks;       // Number of tracklets
  Double32_t *fTheta;        //[fNTracks] array with theta values
  Double32_t *fPhi;          //[fNTracks] array with phi values
  Double32_t *fDeltaPhi;     //[fNTracks] array with delta phi values
  Int_t      *fLabels;       //[fNTracks] array with labels of cluster in L1 used for the tracklet
  Int_t      *fLabelsL2;     //[fNTracks] array with labels of cluster in L2 used for the tracklet


  ClassDef(AliAODTracklets, 3);
};

Double32_t AliAODTracklets::GetTheta(Int_t i) const 
{ 
  if (i>=0 && i<fNTracks) 
  {
    return fTheta[i];
  }
  else 
    Error("GetTheta","Invalid track number %d",i); return -9999.;
}

Double32_t AliAODTracklets::GetPhi(Int_t i) const 
{ 
  if (i>=0 && i<fNTracks) 
  {
    return fPhi[i];
  }
  else 
    Error("GetPhi","Invalid track number %d",i); return -9999.;
}

Double32_t AliAODTracklets::GetDeltaPhi(Int_t i) const 
{
  if (i>=0 && i<fNTracks) 
  {
    return fDeltaPhi[i];
  }
  else 
    Error("GetDeltaPhi","Invalid track number %d",i); return -9999.;
}

Int_t AliAODTracklets::GetLabel(Int_t i, Int_t layer) const 
{
  if (i>=0 && i<fNTracks) 
  {
    return (layer == 0) ? fLabels[i] : fLabelsL2[i];
  }
  else 
    Error("GetLabel","Invalid track number %d",i); return -9999;
}


void AliAODTracklets::SetLabel(Int_t i, Int_t layer,Int_t label)  
{
  if (i>=0 && i<fNTracks) 
  {
    if(layer == 0) fLabels[i] = label;
    else fLabelsL2[i] = label;
  }
}


#endif
