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
class AliFemtoDreamCascade : public AliFemtoDreamBasePart {
 public:
  AliFemtoDreamCascade();
  virtual ~AliFemtoDreamCascade();
  void SetCascade(AliAODEvent *evt,AliAODcascade *casc);
  TString ClassName(){return "Cascade";};
  double GetXiMass() {return fXiMass;};
  double GetOmegaMass() {return fOmegaMass;};
  double GetDCAXiPrimVtx() {return fDCAXiPrimVtx;};
  double GetXiDCADaug() {return fDCAXiDaug;};
  double GetXiTransverseRadius() {return fTransRadius;};
  double GetXiRapidity() {return fRapXi;};
  double GetOmegaRapidity() {return fRapOmega;};
  double GetXiAlpha() {return fAlphaXi;};
  double GetPtArmXi() {return fPtArmXi;};
  double GetXiDecayLength() {return fXiLength;};
  double GetOmegaDecayLength() {return fOmegaLength;};
  double BachDCAPrimVtx() {return fDCABachPrimVtx;};
  double Getv0Mass() {return fMassv0;};
  double Getv0DCADaug() {return fv0DCADaug;};
  double Getv0DCAPrimVtx() {return fv0DCAPrimVtx;};
  double Getv0TransverseRadius() {return fv0TransRadius;};
  double Getv0CPA() {return fv0CPA;};
  double Getv0PosToPrimVtx() {return fv0PosToPrimVtx;};
  double Getv0NegToPrimVtx() {return fv0NegToPrimVtx;};
  double Getv0XiPointingAngle() {return fv0ToXiPointAngle;};
  double Getv0DecayLength() {return fv0Length;};
  double GetDCAv0Xi() {return fDCAv0Xi;};
  void SetPDGDaugPos(int PDGCode){fPosDaug->SetPDGCode(PDGCode);};
  void SetPDGDaugNeg(int PDGCode){fNegDaug->SetPDGCode(PDGCode);};
  void SetPDGDaugBach(int PDGCode){fBach->SetPDGCode(PDGCode);};
  TVector3 Getv0P() {return fv0Momentum;};
  AliFemtoDreamTrack *GetPosDaug()const{return fPosDaug;};
  AliFemtoDreamTrack *GetNegDaug()const{return fNegDaug;};
  AliFemtoDreamTrack *GetBach()const{return fBach;};
 private:
  void Reset();
  void SetMCMotherInfo(AliAODEvent *evt, AliAODcascade *casc);
  AliFemtoDreamTrack *fPosDaug;
  AliFemtoDreamTrack *fNegDaug;
  AliFemtoDreamTrack *fBach;
  double fXiMass;
  double fOmegaMass;
  double fDCAXiDaug;
  double fTransRadius;
  double fRapXi;
  double fRapOmega;
  double fAlphaXi;
  double fPtArmXi;
  double fDCAXiPrimVtx;
  double fXiLength;
  double fOmegaLength;
  double fMassv0;
  double fv0DCADaug;
  double fv0DCAPrimVtx;
  TVector3 fv0Momentum;
  double fDCABachPrimVtx;
  double fv0TransRadius;
  double fv0CPA;
  double fv0PosToPrimVtx;
  double fv0NegToPrimVtx;
  double fv0ToXiPointAngle;
  double fv0Length;
  double fDCAv0Xi;
  ClassDef(AliFemtoDreamCascade,1)
};

#endif /* ALIFEMTODREAMCASCADE_H_ */
