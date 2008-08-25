/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


////////////////////////////////////////////////
// Class to define PID cuts on a pair of tracks
// The track pair object to use is AliCFPair
//
// author : renaud.vernet@cern.ch
////////////////////////////////////////////////

#ifndef ALICFPAIRPIDCUT_H
#define ALICFPAIRPIDCUT_H

#include "AliCFCutBase.h"
#include "AliCFTrackCutPid.h"

//__________________________________________________________________________________
// CUT ON TRACK PID FOR V0 DAUGHTERS
//__________________________________________________________________________________

class AliCFPairPidCut : public AliCFCutBase
{
  public :
  AliCFPairPidCut() ;
  AliCFPairPidCut(const Char_t* name, const Char_t* title) ;
  AliCFPairPidCut(const AliCFPairPidCut& c) ;
  AliCFPairPidCut& operator=(const AliCFPairPidCut& c) ;
  virtual ~AliCFPairPidCut() {delete fCutNeg; delete fCutPos; };

  void SetDetectors(TString detsNeg, TString detsPos) 
  {fCutNeg->SetDetectors(detsNeg); fCutPos->SetDetectors(detsPos);}     //sets the chosen detectors
  void SetPriors(Double_t r[AliPID::kSPECIES]) 
  {fCutNeg->SetPriors(r); fCutPos->SetPriors(r);}                       //sets the a priori concentrations
  void SetProbabilityCut(Double32_t cut1, Double32_t cut2)
  {fCutNeg->SetProbabilityCut(cut1); fCutPos->SetProbabilityCut(cut2);} //sets the prob cut
  void SetParticleType(Int_t iType1, Bool_t tocombine1, Int_t iType2, Bool_t tocombine2)        //sets the particle to be identified and the mode
  {fCutNeg->SetParticleType(iType1,tocombine1); fCutPos->SetParticleType(iType2,tocombine2);}   // (single detector kFALSE/ combined kTRUE)
  void SetMinDiffResp(Bool_t check1, Double_t mindiff1, Bool_t check2, Double_t mindiff2) 
  {fCutNeg->SetMinDiffResp(check1,mindiff1); fCutPos->SetMinDiffResp(check2,mindiff2);}     //set checking at det. response level
  void SetMinDiffProb(Bool_t check1, Double_t mindiff1, Bool_t check2, Double_t mindiff2)
  {fCutNeg->SetMinDiffProb(check1,mindiff1); fCutPos->SetMinDiffProb(check2,mindiff2);}  //set checking at probability level
  void SetAODmode(Bool_t mode) {fCutNeg->SetAODmode(mode); fCutPos->SetAODmode(mode);}

  Bool_t IsSelected(TObject *obj); //boolean for detectors
  Bool_t IsSelected(TList* /*list*/) {return kTRUE;}
 private:
  AliCFTrackCutPid* fCutNeg; // PID cut on negative daughter
  AliCFTrackCutPid* fCutPos; // PID cut on positive daughter

  ClassDef(AliCFPairPidCut,1);
};
#endif
    
