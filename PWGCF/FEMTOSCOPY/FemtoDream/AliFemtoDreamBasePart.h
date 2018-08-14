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
#include "Rtypes.h"
#include "TVector3.h"

class AliFemtoDreamBasePart {
 public:
  AliFemtoDreamBasePart();
  AliFemtoDreamBasePart(const AliFemtoDreamBasePart& part);
  AliFemtoDreamBasePart &operator=(const AliFemtoDreamBasePart &obj);
  virtual ~AliFemtoDreamBasePart();
  enum PartOrigin {
    kPhysPrimary=0,
    kWeak=1,
    kMaterial=2,
    kFake=3,
    kContamination=4,
    kUnknown=5
  };
  void SetMCParticle(AliAODMCParticle *mcPart,AliMCEvent *evt);
  void ResetMCInfo();
  void SetMomentum(float px,float py,float pz) {fP.SetXYZ(px,py,pz);};
  TVector3 GetMomentum() const {return fP;};
  void SetMCMomentum(float px,float py,float pz) {fMCP.SetXYZ(px,py,pz);};
  TVector3 GetMCMomentum() const {return fMCP;};
  void SetPt(float pT){fPt=pT;};
  float GetPt() const {return fPt;};
  void SetMCPt(float pT){fMCPt=pT;};
  float GetMCPt() const {return fMCPt;};
  void SetMomTPC(float pTPC){fP_TPC=pTPC;};
  float GetMomTPC() const {return fP_TPC;};
  void SetEta(float eta){fEta.push_back(eta);};
  std::vector<float>  GetEta()  const {return fEta;};
  void SetTheta(float theta){fTheta.push_back(theta);};
  std::vector<float>  GetTheta()  const {return fTheta;};
  void SetMCTheta(float theta){fMCTheta.push_back(theta);};
  std::vector<float>  GetMCTheta()  const {return fMCTheta;};
  void SetPhi(float phi){fPhi.push_back(phi);};
  std::vector<float>  GetPhi()  const {return fPhi;};
  void SetPhiAtRadius(std::vector<float> phiAtRad) {fPhiAtRadius.push_back(phiAtRad);};
  std::vector<std::vector<float>> GetPhiAtRaidius() {return fPhiAtRadius;};
  void SetMCPhi(float phi){fMCPhi.push_back(phi);};
  std::vector<float>  GetMCPhi()  const {return fMCPhi;};
  void SetIDTracks(int idTracks) {fIDTracks.push_back(idTracks);};
  std::vector<int> GetIDTracks() {return fIDTracks;};
  void SetCharge(int charge){fCharge.push_back(charge);};
  std::vector<int> GetCharge() const {return fCharge;};
  void SetCPA(float cpa) {fCPA=cpa;};
  float GetCPA() const {return fCPA;};
  void SetParticleOrigin(PartOrigin org){fOrigin=org;};
  PartOrigin GetParticleOrigin() const {return fOrigin;};
  void SetPDGCode(int pdgCode) {fPDGCode=pdgCode;};
  int GetPDGCode() const {return fPDGCode;};
  void SetMCPDGCode(int pdgCode) {fMCPDGCode=pdgCode;};
  int GetMCPDGCode() const {return fMCPDGCode;};
  void SetPDGMotherWeak(int PDGMother){fPDGMotherWeak=PDGMother;};
  int GetMotherWeak() const {return fPDGMotherWeak;};
  void SetMotherID(int ID) {fMotherID=ID;};
  int GetMotherID() const {return fMotherID;};
  void SetEvtNumber(int evtNumb){fEvtNumber=evtNumb;};
  int GetEvtNumber() const {return fEvtNumber;};
  void SetUseMCInfo(bool useMC) {fIsMC=useMC;};
  bool IsSet() const {return fIsSet;};
  void SetUse(bool use) {fUse=use;};
  bool UseParticle() const {return fUse;};
  void SetGlobalTrackInfo(AliAODTrack **GTI, Int_t size)
  {fGTI = GTI; fTrackBufferSize = size;};
 protected:
  bool fIsReset;
  AliAODTrack **fGTI;
  int fTrackBufferSize;
  TVector3 fP;
  TVector3 fMCP;
  float fPt;
  float fMCPt;
  float fP_TPC;
  std::vector<float> fEta;
  std::vector<float> fTheta;
  std::vector<float> fMCTheta;
  std::vector<float> fPhi;
  std::vector<std::vector<float>> fPhiAtRadius;
  std::vector<float> fMCPhi;
  std::vector<int> fIDTracks;
  std::vector<int> fCharge;
  float fCPA;
  PartOrigin fOrigin;
  // pdg code as set by the track cuts, used for invariant mass calculation/mc matching in v0s
  int fPDGCode;
  // pdg code as obtained from the MC for this particle
  int fMCPDGCode;
  int fPDGMotherWeak;
  int fMotherID;
  int fEvtNumber;
  bool fIsMC;
  bool fUse;    //passes cuts
  bool fIsSet;  //has all the attributes set properly
 private:
//  AliFemtoDreamBasePart(const AliFemtoDreamBasePart&);
  ClassDef(AliFemtoDreamBasePart,2);
};

#endif /* ALIFEMTODREAMBASEPART_H_ */
