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
#include "AliESDEvent.h"
#include "AliAODv0.h"
#include "AliESDv0.h"

class AliFemtoDreamv0 : public AliFemtoDreamBasePart {
 public:
  AliFemtoDreamv0(const int nDaugh = 3);
  virtual ~AliFemtoDreamv0();
  void Setv0(AliAODEvent *evt, AliAODv0 *v0, const int multiplicity = -1);
  void Setv0(AliVEvent *evt, AliAODv0 *v0, const int multiplicity = -1);
  void Setv0(AliESDEvent *evt, AliMCEvent *mcEvent, AliESDv0 *v0,
             const int multiplicity = -1);
  // the last two are switches to ignore the first entry in the phi, eta, ...
  // vector in case one cases a V0 in this object, and only wants to keep the
  // daughter information
  void Setv0(const AliFemtoDreamBasePart &posDaughter,
             const AliFemtoDreamBasePart &negDaughter,
             const bool ignoreFirstPos = false, const bool ignoreFirstNeg =
                 false);
  void Setv0(const AliFemtoDreamBasePart &posDaughter,
             const AliFemtoDreamBasePart &negDaughter, AliAODEvent *evt,
             const bool ignoreFirstPos = false, const bool ignoreFirstNeg =
                 false,
             const bool setDaughter = true);
  void Setv0(const AliFemtoDreamBasePart &posDaughter,
             const AliFemtoDreamBasePart &negDaughter, AliVEvent *evt,
             const bool ignoreFirstPos = false, const bool ignoreFirstNeg =
                 false,
             const bool setDaughter = true);
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
    return GetInvMass();
  }
  ;
  void Setv0Mass(float mass) {
    SetInvMass(mass);
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
  double DecayLengthV0(const double *DecayVtx, const double *point) const;
  double CosPointingAngle(const double *DecayVtx, const double *point) const;
  double DecayLengthXY(double const *DecayVtx, double const *point) const;
  float GetArmenterosAlpha() const;
  float GetArmenterosQt() const;
 private:
  AliFemtoDreamv0 &operator=(const AliFemtoDreamv0 &obj);
  AliFemtoDreamv0(const AliFemtoDreamv0&);
  void Reset();
  void SetDaughter(AliAODv0 *v0);
  void SetDaughter(AliAODv0 *v0, AliVEvent *evt);
  void SetDaughter(const AliFemtoDreamBasePart &posDaughter, const AliFemtoDreamBasePart &negDaughter);
  void SetDaughter(const AliFemtoDreamBasePart &posDaughter,
                   const AliFemtoDreamBasePart &negDaughter,
                   AliVEvent *evt);
  void SetDaughter(AliESDEvent *evt, AliMCEvent *mcEvent, AliESDv0 *v0);
  void SetDaughterInfo(AliAODv0 *v0);
  void SetDaughterInfo(AliESDv0 *v0);
  void SetMotherInfo(AliAODEvent *evt, AliAODv0 *v0);
  void SetMotherInfo(AliESDEvent *evt, AliESDv0 *v0);
  void SetMCMotherInfo(AliAODEvent *evt, AliAODv0 *v0);
  void SetMCMotherInfo(AliVEvent *evt, AliAODv0 *v0);
  void SetMCMotherInfo(TClonesArray *mcarray, AliAODv0 *v0);
  bool fOnlinev0;
  bool fHasDaughter;
  AliFemtoDreamTrack *fpDaug;
  AliFemtoDreamTrack *fnDaug;
  double fv0Vtx[3];  // Decay Vertex in xyz
  float fdcav0Daug;  // Daugther to Daughter DCA
  float fdcaPrim;
  float fdcaPrimPos;  // rphi impact params w.r.t. Primary Vtx [cm]
  float fdcaPrimNeg;  // rphi impact params w.r.t. Primary Vtx [cm]
  float flenDecay;   // Decay Length
  float fTransRadius;   // Decay Length in xy
ClassDef(AliFemtoDreamv0, 4)
};

inline double AliFemtoDreamv0::DecayLengthV0(const double *DecayVtx,
                                             const double *point) const {
  return ::sqrt(
      ::pow(DecayVtx[0] - point[0], 2) + ::pow(DecayVtx[1] - point[1], 2)
          + ::pow(DecayVtx[2] - point[2], 2));
}

inline float AliFemtoDreamv0::GetArmenterosAlpha() const {
  TVector3 v0P = GetMomentum();
  TVector3 posP = fpDaug->GetMomentum();
  TVector3 negP = fnDaug->GetMomentum();
  return 1. - 2. / (1. + posP.Dot(v0P) / negP.Dot(v0P));
}

inline float AliFemtoDreamv0::GetArmenterosQt() const {
  TVector3 posP = fpDaug->GetMomentum();
  TVector3 v0P = GetMomentum();
  return posP.Perp(v0P);
}

#endif /* ALIFEMTODREAMV0_H_ */
