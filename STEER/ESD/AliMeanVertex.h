#ifndef ALIMEANVERTEX_H
#define ALIMEANVERTEX_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/*****************************************************************************
 *                                                                           *
 * This class contains the coordinates of the mean primary vertex position   *
 * computed by AliITSMeanVertex                                              *
 *                                                                           *
*****************************************************************************/
#include "AliESDVertex.h"

class AliMeanVertex : public AliESDVertex {
 public:
  AliMeanVertex();
  AliMeanVertex(Double_t pos[3],Double_t err[3],Double_t cov[6],Int_t nevents, Float_t notracklets, Float_t avertracklets, Float_t signotrackl);
  virtual ~AliMeanVertex() {}

  Int_t GetNumberOfContributingEvents() const { return GetNContributors(); }
  void GetErrorsOnPosition(Double_t err[3]) const;
  Float_t GetTotalNumbOfTracklets() const { return fTotTracklets; }
  Float_t GetAverageNumbOfTracklets() const { return fAverTracklets; }
  Float_t GetSigmaOnAvNumbOfTracks() const { return fSigmaOnAverTrack; }

 protected:
  Double32_t fErrW[3];       // errors on vertex coordinates (weighted average)
  Float_t      fTotTracklets;   // total number of tracklets used for M.V.
  Float_t    fAverTracklets;  // average number of tracklets per event
  Float_t fSigmaOnAverTrack;  // sigma on fAverTracklets

  ClassDef(AliMeanVertex,1)  // Class for mean Vertex   
}; 

#endif
