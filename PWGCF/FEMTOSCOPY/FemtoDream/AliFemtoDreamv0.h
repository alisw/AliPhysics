/*
 * AliFemtoDreamv0.h
 *
 *  Created on: Dec 12, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMV0_H_
#define ALIFEMTODREAMV0_H_
#include "Rtypes.h"
#include "AliFemtoDreamBasePart.h"
#include "AliFemtoDreamTrack.h"
#include "AliAODEvent.h"
#include "AliAODv0.h"

class AliFemtoDreamv0 : public AliFemtoDreamBasePart {
 public:
  AliFemtoDreamv0();
  virtual ~AliFemtoDreamv0();
  void Setv0(AliAODEvent *evt, AliAODv0 *v0, const int multiplicity = -1);
  bool GetOnlinev0() const {
    return fOnlinev0;
  }
  ;
  bool GetHasDaughters() const {
    return fHasDaughter;
  }
  ;
  AliFemtoDreamTrack* GetPosDaughter() const {
    return fpDaug;
  }
  ;
  void SetPDGDaughterPos(int PDGCode) {
    fpDaug->SetPDGCode(PDGCode);
  }
  ;
  AliFemtoDreamTrack* GetNegDaughter() const {
    return fnDaug;
  }
  ;
  void SetPDGDaughterNeg(int PDGCode) {
    fnDaug->SetPDGCode(PDGCode);
  }
  ;
  float Getv0Mass() const {
    return fv0Mass;
  }
  ;
  void Setv0Mass(float mass) {
    fv0Mass = mass;
  }
  ;
  float GetDCAv0Vtx(int i) const {
    if (i > -1 && i < 3) {
      return fv0Vtx[i];
    } else {
      return 0.;
    };
  }
  ;
  float GetDaugDCA() const {
    return fdcav0Daug;
  }
  ;
  float GetDCAPrimVtx() const {
    return fdcaPrim;
  }
  ;
  float GetDCADaugPosVtx() const {
    return fdcaPrimPos;
  }
  ;
  float GetDCADaugNegVtx() const {
    return fdcaPrimNeg;
  }
  ;
  float GetDecayLength() const {
    return flenDecay;
  }
  ;
  float GetTransverseRadius() const {
    return fTransRadius;
  }
  ;
  TString ClassName() {
    return "v0Class";
  }
  ;
 private:
  AliFemtoDreamv0 &operator=(const AliFemtoDreamv0 &obj);
  AliFemtoDreamv0(const AliFemtoDreamv0&);
  void Reset();
  void SetDaughter(AliAODv0 *v0);
  void SetDaughterInfo(AliAODv0 *v0);
  void SetMotherInfo(AliAODEvent *evt, AliAODv0 *v0);
  void SetMCMotherInfo(AliAODEvent *evt, AliAODv0 *v0);
  bool fOnlinev0;
  bool fHasDaughter;
  AliFemtoDreamTrack *fpDaug;
  AliFemtoDreamTrack *fnDaug;
  float fv0Mass;
  double fv0Vtx[3];  // Decay Vertex in xyz
  float fdcav0Daug;  // Daugther to Daughter DCA
  float fdcaPrim;
  float fdcaPrimPos;  // rphi impact params w.r.t. Primary Vtx [cm]
  float fdcaPrimNeg;  // rphi impact params w.r.t. Primary Vtx [cm]
  float flenDecay;   // Decay Length
  float fTransRadius;   // Decay Length in xy
ClassDef(AliFemtoDreamv0,2)
};

#endif /* ALIFEMTODREAMV0_H_ */
