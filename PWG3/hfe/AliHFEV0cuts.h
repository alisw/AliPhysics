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
//
// Class for the V0 cuts - tuned to obtain clean eletron, pion and proton samples.
// NOT suitable for V0 analysis
//
#ifndef ALIHFEV0CUTS_H
#define ALIHFEV0CUTS_H

#include "AliHFEcollection.h"

class TList;

class AliMCParticle;
class AliVEvent;
class AliMCEvent;
class AliESDtrack;
class AliESDv0;
class AliKFVertex;
class AliKFParticle;
class AliVTrack;

class AliHFEV0cuts : public TObject {
 public:
  enum{ // Reconstructed V0
    kUndef = 0,
      kRecoGamma = 1,
      kRecoK0 = 2,
      kRecoPhi = 3,
      kRecoLambda = 4,
      kRecoALambda = -4
      
      };
  enum{ // Identified Daughter particles
    kRecoElectron = 0,
      kRecoPionK0 = 1,
      kRecoPionL = 2,
      kRecoKaon = 3,
      kRecoProton = 4
      };
  AliHFEV0cuts();
  ~AliHFEV0cuts();
  AliHFEV0cuts(const AliHFEV0cuts &ref);
  AliHFEV0cuts &operator=(const AliHFEV0cuts &ref);

  void Init(const char* name);
  
  void RunQA();
  void SetMCEvent(AliMCEvent* const mce)      { fMCEvent = mce; };
  void SetInputEvent(AliVEvent* const e)      { fInputEvent = e; };
  void SetPrimaryVertex(AliKFVertex* const v) { fPrimaryVertex = v; };
  
  TList* GetList()    { return fQA->GetList(); };
  TList* GetListMC()  { return fQAmc->GetList(); };
  
  Bool_t   TrackCutsCommon(AliESDtrack* track);
  Bool_t   V0CutsCommon(AliESDv0 *v0);
  Bool_t   GammaCuts(AliESDv0 *v0);
  Bool_t   K0Cuts(AliESDv0 *v0);
  Bool_t   LambdaCuts(AliESDv0 *v0, Bool_t &isLambda);
 
  Bool_t LooseRejectK0(AliESDv0 * const v0) const;
  Bool_t LooseRejectLambda(AliESDv0 * const v0) const;
  Bool_t LooseRejectGamma(AliESDv0 * const v0) const;

  void     Armenteros(AliESDv0 *v0, Float_t val[2]);

  Double_t OpenAngle(AliESDv0 *v0) const;//opening angle between V0 daughters; close to zero for conversions
  Double_t PsiPair(AliESDv0 *v0);
  
  Bool_t   CheckSigns(AliESDv0* const v0);
  Bool_t   GetConvPosXY(AliESDtrack * const ptrack, AliESDtrack * const ntrack, Double_t convpos[2]);
  Bool_t   GetHelixCenter(AliESDtrack * const track, Double_t b,Int_t charge, Double_t center[2]);


  // MC stuff
  void     SetCurrentV0id(Int_t id) { fCurrentV0id = id; };
  void     SetDaughtersID(Int_t d[2]) {fPdaughterPDG = d[0]; fNdaughterPDG = d[1]; };
  
  AliKFParticle *CreateMotherParticle(AliVTrack* const pdaughter, AliVTrack* const ndaughter, Int_t pspec, Int_t nspec);

 private:
  void Copy(TObject &ref) const;
      
 private:
  
  AliHFEcollection     *fQA;            // store QA cut histograms
  AliHFEcollection     *fQAmc;          // store 
  AliMCEvent           *fMCEvent;       // MC event
  AliVEvent            *fInputEvent;    // Input Event
  AliKFVertex          *fPrimaryVertex; // primary vertex

  Int_t                fCurrentV0id;   // MC flagged V0    
  Int_t                fPdaughterPDG;   // MC id of the positive daugeter
  Int_t                fNdaughterPDG;   // MC id of the negative daugeter


  ClassDef(AliHFEV0cuts, 1)
};
    

#endif
