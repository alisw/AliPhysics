//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTD0CANDIDATE_H
#define ALIHLTD0CANDIDATE_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTD0Candidate.h
/// @author Gaute Ovrebekk
/// @date   2010-11-19
/// @brief  Class for storing the D0 candidates

#include "TObject.h"

class AliHLTD0Candidate : public TObject
{
 public:
  AliHLTD0Candidate();
  AliHLTD0Candidate(Double_t InvD0,Double_t InvD0bar,Double_t pt,Int_t lPos, Int_t lNeg, Double_t PtPos, Double_t PtNeg);
  ~AliHLTD0Candidate();
  
  // Get function for D0 ivariant mass
  Double_t GetInvMassD0()const {return fInvD0;} 
  // Get function for D0 bar ivariant mass
  Double_t GetInvMassD0bar()const {return fInvD0bar;} 
  // Get function for pt of D0
  Double_t GetD0Pt()const {return fpt;} 
  // Get function for the label of the positive track
  Int_t GetPosLabel()const {return flabelPos;}
  // Get function for the negative track
  Int_t GetNegLabel()const {return flabelNeg;}
  // Get function for the Pt of the positive track 
  Double_t GetPtPos()const {return fPtPos;} 
  // Get function for the Pt of the negative track
  Double_t GetPtNeg()const {return fPtNeg;}

 private:
  //Invariant mass D0
  Double_t fInvD0;                                     // Invariant mass D0
  //Invariant mass D0 bar
  Double_t fInvD0bar;                                  // Invariant mass D0bar
  // Pt for the D0 Candidate
  Double_t fpt;                                        // Pt for the D0 Candidate
  // Storng the lable of the positvie track
  Int_t flabelPos;                                     // Storng the lable of the positvie track
  // Storng the lable of the negative track
  Int_t flabelNeg;                                     // Storng the lable of the negative track
  // Storng the Pt of the positive track
  Double_t fPtPos;                                     // Storng the Pt of the positive track
  // Storng the Pt of the negative track
  Double_t fPtNeg;                                     // Storng the Pt of the negative track


  ClassDef(AliHLTD0Candidate, 1)

};
#endif
