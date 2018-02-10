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
  void Setv0(AliAODEvent *evt,AliAODv0 *v0);
  bool GetOnlinev0() const {return fOnlinev0;};
  bool GetHasDaughters() const {return fHasDaughter;};
  AliFemtoDreamTrack* GetPosDaughter() const {return fpDaug;};
  void SetPDGDaughterPos(int PDGCode) {fpDaug->SetPDGCode(PDGCode);};
  AliFemtoDreamTrack* GetNegDaughter() const {return fnDaug;};
  void SetPDGDaughterNeg(int PDGCode) {fnDaug->SetPDGCode(PDGCode);};
  double Getv0Mass() const {return fv0Mass;};
  void Setv0Mass(double mass){fv0Mass=mass;};
  double GetDCAv0Vtx(int i)const {
    if(i>-1 && i<3){return fv0Vtx[i];}else{return 0.;};
  };
  double GetDaugDCA() const {return fdcav0Daug;};
  double GetDCAPrimVtx() const {return fdcaPrim;};
  double GetDCADaugPosVtx() const {return fdcaPrimPos;};
  double GetDCADaugNegVtx() const {return fdcaPrimNeg;};
  double GetDecayLength() const {return flenDecay;};
  double GetTransverseRadius() const {return fTransRadius;};
  TString ClassName(){return "v0Class";};
 private:
  void Reset();
  void SetDaughter(AliAODv0 *v0);
  void SetDaughterInfo(AliAODv0 *v0);
  void SetMotherInfo(AliAODEvent *evt,AliAODv0 *v0);
  void SetMCMotherInfo(AliAODEvent *evt,AliAODv0 *v0);
  bool fOnlinev0;
  bool fHasDaughter;
  AliFemtoDreamTrack *fpDaug;
  AliFemtoDreamTrack *fnDaug;
  double fv0Mass;
  double fv0Vtx[3];   // Decay Vertex in xyz
  double fdcav0Daug;  // Daugther to Daughter DCA
  double fdcaPrim;
  double fdcaPrimPos; // rphi impact params w.r.t. Primary Vtx [cm]
  double fdcaPrimNeg; // rphi impact params w.r.t. Primary Vtx [cm]
  double flenDecay;   // Decay Length
  double fTransRadius;// Decay Length in xy
  ClassDef(AliFemtoDreamv0,1)
};

#endif /* ALIFEMTODREAMV0_H_ */
