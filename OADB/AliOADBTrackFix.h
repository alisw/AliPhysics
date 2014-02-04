#ifndef ALIOADBTRACKFIX_H
#define ALIOADBTRACKFIX_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------------------------
//     OADB class for run dependent track fixing parameters
//     Convention for phi-dependent data: 0 : 2pi
//     Author: ruben.shahoyan@cern.ch
//-------------------------------------------------------------------------

#include <TNamed.h>
class TGraph;

class AliOADBTrackFix : public TNamed 
{
 public :
  enum CorMode_t {kCorModeGlob, kCorModeTPCInner, kNCorModes};
  //
  AliOADBTrackFix();
  AliOADBTrackFix(const char* name);
  virtual ~AliOADBTrackFix();
  //
  Double_t GetPtInvCorr(int mode, double sideAfrac, double phi=0) const;
  //
  TGraph*  GetPtInvCorrGraph(int mode,int side)             const {return (TGraph*)fPtInvCor[mode][side];}
  Double_t GetXIniPtInvCorr(int mode)                       const {return fXIniPtInvCorr[mode];}
  //
  void     SetPtInvCorr(int mode,int side, const TGraph* gr);
  void     SetXIniPtInvCorr(int mode, double x=0)                 {fXIniPtInvCorr[mode] = x;}
  //
 private:
  AliOADBTrackFix(const AliOADBTrackFix& cont); 
  AliOADBTrackFix& operator=(const AliOADBTrackFix& cont);

 protected:
  const TGraph   *fPtInvCor[kNCorModes][2];    // graphs with 1/pt correction vs phi for A,C sides
  Double_t        fXIniPtInvCorr[kNCorModes];  // if >0 use as the reper X for slope,position correction of corresponding mode
  //
  ClassDef(AliOADBTrackFix, 1);
};

#endif
