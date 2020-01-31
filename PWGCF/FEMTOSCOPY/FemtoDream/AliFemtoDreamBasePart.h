/*
 * AliFemtoPPbpbLamBaseParticle.h
 *
 *  Created on: Aug 29, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMBASEPART_H_
#define ALIFEMTODREAMBASEPART_H_
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"
#include "Rtypes.h"
#include "TVector3.h"
#include "AliAODConversionPhoton.h"

class AliFemtoDreamBasePart {
 public:
  AliFemtoDreamBasePart(const int part = 1);
  AliFemtoDreamBasePart(const AliFemtoDreamBasePart& part);
  AliFemtoDreamBasePart &operator=(const AliFemtoDreamBasePart &obj);
  AliFemtoDreamBasePart(const AliAODConversionPhoton *gamma,
                        const AliVTrack *pos, const AliVTrack *neg,
                        const AliVEvent *inputEvent);
  virtual ~AliFemtoDreamBasePart();
  enum PartOrigin {
    kPhysPrimary = 0,
    kWeak = 1,
    kMaterial = 2,
    kFake = 3,
    kContamination = 4,
    kUnknown = 5
  };
  void SetMCParticle(AliAODMCParticle *mcPart, AliMCEvent *evt);
  void SetMCParticleRePart(AliAODMCParticle *mcPart);
  void ResetMCInfo();
  void SetMomentum(unsigned int iEntr, TVector3 mom) {
    SetMomentum(iEntr, mom.X(), mom.Y(), mom.Z());
  }
  void SetMomentum(unsigned int iEntr, float px, float py, float pz) {
    if (iEntr < fP.size()) {
      fP.at(iEntr).SetXYZ(px, py, pz);
    } else {
      AliWarning("Trying to set momentum for a too small vector!");
    }
  }
  ;
  TVector3 GetMomentum(unsigned int iEntr = 0) const {
    if (iEntr > fP.size()) {
      std::cout << "Trying to get a momentum which is out of bounds. iEntr = "
                << iEntr << std::endl;
    } else {
      return fP.at(iEntr);
    }
    return TVector3(-999,-999,-999);
  }
  ;
  float GetP() const {
    return GetMomentum().Mag();
  }
  void SetMCMomentum(float px, float py, float pz) {
    fMCP.SetXYZ(px, py, pz);
  }
  ;
  TVector3 GetMCMomentum() const {
    return fMCP;
  }
  ;
  void SetPt(float pT) {
    fPt = pT;
  }
  ;
  float GetPt() const {
    return fPt;
  }
  ;
  void SetMCPt(float pT) {
    fMCPt = pT;
  }
  ;
  float GetMCPt() const {
    return fMCPt;
  }
  ;
  void SetMomTPC(float pTPC) {
    fP_TPC = pTPC;
  }
  ;
  float GetMomTPC() const {
    return fP_TPC;
  }
  ;
  void SetEta(float eta) {
    fEta.push_back(eta);
  }
  ;
  std::vector<float> GetEta() const {
    return fEta;
  }
  ;
  void SetTheta(float theta) {
    fTheta.push_back(theta);
  }
  ;
  std::vector<float> GetTheta() const {
    return fTheta;
  }
  ;
  void SetMCTheta(float theta) {
    fMCTheta.push_back(theta);
  }
  ;
  std::vector<float> GetMCTheta() const {
    return fMCTheta;
  }
  ;
  void SetPhi(float phi) {
    fPhi.push_back(phi);
  }
  ;
  std::vector<float> GetPhi() const {
    return fPhi;
  }
  ;
  void SetPhiAtRadius(std::vector<float> phiAtRad) {
    fPhiAtRadius.push_back(phiAtRad);
  }
  ;
  std::vector<std::vector<float>> GetPhiAtRaidius() const {
    return fPhiAtRadius;
  }
  ;
  void ResizePhiAtRadii(size_t i) {
    fPhiAtRadius.resize(i);
  }
  float GetAveragePhiAtRadius(size_t iPart) {
    if (iPart > fPhiAtRadius.size()) {
      std::cout << "ERROR - AliFemtoDreamBasePart::GetAveragePhiAtRadius\n";
      return 999.f;
    }
    float nCounts = fPhiAtRadius[iPart].size();
    float avPhi = 0.f;
    for (size_t i = 0; i < nCounts; ++i) {
      avPhi += fPhiAtRadius[iPart][i];
    }
    return avPhi / nCounts;
  }
  void SetXYZAtRadius(TVector3 XYZAtRad) {
    fXYZAtRadius.push_back(XYZAtRad);
  }
  ;
  std::vector<TVector3> GetXYZAtRadius() {
    return fXYZAtRadius;
  }
  ;
  void SetMCPhi(float phi) {
    fMCPhi.push_back(phi);
  }
  ;
  std::vector<float> GetMCPhi() const {
    return fMCPhi;
  }
  ;
  void SetIDTracks(int idTracks) {
    fIDTracks.push_back(idTracks);
  }
  ;
  std::vector<int> GetIDTracks() const {
    return fIDTracks;
  }
  ;
  void SetCharge(int charge) {
    fCharge.push_back(charge);
  }
  ;
  std::vector<int> GetCharge() const {
    return fCharge;
  }
  ;
  void SetCPA(float cpa) {
    fCPA = cpa;
  }
  ;
  float GetCPA() const {
    return fCPA;
  }
  ;
  void SetInvMass(float invMass) {
    fInvMass = invMass;
  }
  float GetInvMass() const {
    return fInvMass;
  }
  void SetParticleOrigin(PartOrigin org) {
    fOrigin = org;
  }
  ;
  PartOrigin GetParticleOrigin() const {
    return fOrigin;
  }
  ;
  void SetPDGCode(int pdgCode) {
    fPDGCode = pdgCode;
  }
  ;
  int GetPDGCode() const {
    return fPDGCode;
  }
  ;
  void SetMCPDGCode(int pdgCode) {
    fMCPDGCode = pdgCode;
  }
  ;
  int GetMCPDGCode() const {
    return fMCPDGCode;
  }
  ;
  void SetPDGMotherWeak(int PDGMother) {
    fPDGMotherWeak = PDGMother;
  }
  ;
  int GetMotherWeak() const {
    return fPDGMotherWeak;
  }
  ;
  void SetMotherID(int ID) {
    fMotherID = ID;
  }
  ;
  int GetMotherID() const {
    return fMotherID;
  }
  ;
  void SetID(int ID) {
    fID = ID;
  }
  int GetID() const {
    return fID;
  }
  void SetMotherPDG(int pdg) {
    fMotherPDG = pdg;
  }
  ;
  int GetMotherPDG() const {
    return fMotherPDG;
  }
  ;
  void SetEvtNumber(int evtNumb) {
    fEvtNumber = evtNumb;
  }
  ;
  int GetEvtNumber() const {
    return fEvtNumber;
  }
  ;
  void SetUseMCInfo(bool useMC) {
    fIsMC = useMC;
  }
  ;
  bool IsSet() const {
    return fIsSet;
  }
  ;
  void SetUse(bool use) {
    fUse = use;
  }
  ;
  bool UseParticle() const {
    return fUse;
  }
  ;
  void SetGlobalTrackInfo(AliAODTrack **GTI, Int_t size) {
    fGTI = GTI;
    fTrackBufferSize = size;
  }
  void SetGlobalTrackInfo(AliVTrack **VGTI, Int_t size) {
    fVGTI = VGTI;
    fTrackBufferSize = size;
  }
  void KillGlobalTrackArray() {
    fGTI = nullptr;
    fVGTI = nullptr;
  }
  int GetEventMultiplicity() const {
    return fEvtMultiplicity;
  }
  void SetEventMultiplicity(int evtMulti) {
    fEvtMultiplicity = evtMulti;
  }
  void DumpParticleInformation();
  TString ClassName() {
    return "AliFemtoDreamBasePart";
  }
  ;
 protected:
  bool fIsReset;
  AliAODTrack **fGTI;
  AliVTrack **fVGTI;
  int fTrackBufferSize;
  std::vector<TVector3> fP;
  TVector3 fMCP;
  float fPt;
  float fMCPt;
  float fP_TPC;
  std::vector<float> fEta;
  std::vector<float> fTheta;
  std::vector<float> fMCTheta;
  std::vector<float> fPhi;
  std::vector<std::vector<float>> fPhiAtRadius;
  std::vector<TVector3> fXYZAtRadius;
  std::vector<float> fMCPhi;
  std::vector<int> fIDTracks;
  std::vector<int> fCharge;
  float fCPA;
  float fInvMass;
  PartOrigin fOrigin;
// pdg code as set by the track cuts, used for invariant mass calculation/mc matching in v0s
  int fPDGCode;
// pdg code as obtained from the MC for this particle
  int fMCPDGCode;
  int fPDGMotherWeak;
  int fMotherID;
  int fID;
  int fMotherPDG;
  int fEvtNumber;
  bool fIsMC;
  bool fUse;    //passes cuts
  bool fIsSet;  //has all the attributes set properly
  int fEvtMultiplicity;
 private:
  void PhiAtRadii(const AliVTrack *track, const float bfield,
                  std::vector<float> &tmpVec);
//  AliFemtoDreamBasePart(const AliFemtoDreamBasePart&);
ClassDef(AliFemtoDreamBasePart, 6)
  ;
};

#endif /* ALIFEMTODREAMBASEPART_H_ */
