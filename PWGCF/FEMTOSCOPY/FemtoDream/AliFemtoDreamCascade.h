/*
 * AliFemtoDreamCascade.h
 *
 *  Created on: Jan 11, 2018
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMCASCADE_H_
#define ALIFEMTODREAMCASCADE_H_
#include "Rtypes.h"
#include "TVector3.h"

#include "AliFemtoDreamBasePart.h"
#include "AliFemtoDreamTrack.h"
#include "AliAODEvent.h"
#include "AliAODcascade.h"
#include "AliESDEvent.h"
#include "AliESDcascade.h"

class AliFemtoDreamCascade : public AliFemtoDreamBasePart {
 public:
  AliFemtoDreamCascade();
  virtual ~AliFemtoDreamCascade();
  void SetCascade(AliAODEvent *evt,AliAODcascade *casc);
  void SetCascade(AliESDEvent *evt,AliESDcascade *casc);
  TString ClassName(){return "Cascade";};
  float GetXiMass() {return fXiMass;};
  float GetOmegaMass() {return fOmegaMass;};
  float GetDCAXiPrimVtx() {return fDCAXiPrimVtx;};
  float GetXiDCADaug() {return fDCAXiDaug;};
  float GetXiTransverseRadius() {return fTransRadius;};
  float GetXiRapidity() {return fRapXi;};
  float GetOmegaRapidity() {return fRapOmega;};
  float GetXiAlpha() {return fAlphaXi;};
  float GetPtArmXi() {return fPtArmXi;};
  float GetXiDecayLength() {return fXiLength;};
  float GetOmegaDecayLength() {return fOmegaLength;};
  float BachDCAPrimVtx() {return fDCABachPrimVtx;};
  float Getv0Mass() {return fMassv0;};
  float Getv0DCADaug() {return fv0DCADaug;};
  float Getv0DCAPrimVtx() {return fv0DCAPrimVtx;};
  float Getv0TransverseRadius() {return fv0TransRadius;};
  float Getv0CPA() {return fv0CPA;};
  float Getv0PosToPrimVtx() {return fv0PosToPrimVtx;};
  float Getv0NegToPrimVtx() {return fv0NegToPrimVtx;};
  float Getv0XiPointingAngle() {return fv0ToXiPointAngle;};
  float Getv0DecayLength() {return fv0Length;};
  float GetDCAv0Xi() {return fDCAv0Xi;};
  void SetPDGDaugPos(int PDGCode){fPosDaug->SetPDGCode(PDGCode);};
  void SetPDGDaugNeg(int PDGCode){fNegDaug->SetPDGCode(PDGCode);};
  void SetPDGDaugBach(int PDGCode){fBach->SetPDGCode(PDGCode);};
  TVector3 Getv0P() {return fv0Momentum;};
  float Getv0Pt() {return fv0Pt;};
  void Setv0PDGCode(int PDG) {fv0PDG=PDG;};
  AliFemtoDreamTrack *GetPosDaug()const{return fPosDaug;};
  AliFemtoDreamTrack *GetNegDaug()const{return fNegDaug;};
  AliFemtoDreamTrack *GetBach()const{return fBach;};
 private:
  void Reset();
  void SetMCMotherInfo(AliAODEvent *evt, AliAODcascade *casc);
  AliFemtoDreamTrack *fPosDaug;
  AliFemtoDreamTrack *fNegDaug;
  AliFemtoDreamTrack *fBach;
  float fXiMass;
  float fOmegaMass;
  float fDCAXiDaug;
  float fTransRadius;
  float fRapXi;
  float fRapOmega;
  float fAlphaXi;
  float fPtArmXi;
  float fDCAXiPrimVtx;
  float fXiLength;
  float fOmegaLength;
  float fMassv0;
  float fv0DCADaug;
  float fv0DCAPrimVtx;
  TVector3 fv0Momentum;
  float fv0Pt;
  float fDCABachPrimVtx;
  float fv0TransRadius;
  float fv0CPA;
  float fv0PosToPrimVtx;
  float fv0NegToPrimVtx;
  float fv0ToXiPointAngle;
  float fv0Length;
  float fDCAv0Xi;
  int fv0PDG;
  ClassDef(AliFemtoDreamCascade,2)
};

#endif /* ALIFEMTODREAMCASCADE_H_ */
