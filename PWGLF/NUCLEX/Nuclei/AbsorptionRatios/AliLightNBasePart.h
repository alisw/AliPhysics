#ifndef ALILIGHTNBASEPART_H
#define ALILIGHTNBASEPART_H

/*
 * AliFemtoPPbpbLamBaseParticle.h
 *
 *  Created on: Aug 29, 2017
 *      Author: gu74req
 */


#include "AliAODTrack.h"
#include "Rtypes.h"
#include "TVector3.h"

class AliLightNBasePart {
public:
    AliLightNBasePart();
    virtual ~AliLightNBasePart();
    enum PartOrigin {
        kPhysPrimary=0,
        kWeak=1,
        kMaterial=2,
        kContamination=3,
        kUnknown=4
    };
    
    void SetMomentum(double px,double py,double pz) {fP.SetXYZ(px,py,pz);};
    TVector3 GetMomentum() const {return fP;};
    void SetMCMomentum(double px,double py,double pz) {fMCP.SetXYZ(px,py,pz);};
    TVector3 GetMCMomentum() const {return fMCP;};
    void SetTrackChi2perNDF(double Chi2perNDF){fTrackChi2perNDF=Chi2perNDF;};
    double GetTrackChi2perNDF() const {return fTrackChi2perNDF;};
    void SetP(double p){fMom=p;};
    double GetP() const {return fMom;};
    void SetPt(double pT){fPt=pT;};
    double GetPt() const {return fPt;};
    void SetMCPt(double pT){fMCPt=pT;};
    double GetMCPt() const {return fMCPt;};
    void SetMomTPC(double pTPC){fP_TPC=pTPC;};
    double GetMomTPC() const {return fP_TPC;};
    void SetEta(double eta){fEta.push_back(eta);};
    std::vector<double>  GetEta()  const {return fEta;};
    void SetTheta(double theta){fTheta.push_back(theta);};
    std::vector<double>  GetTheta()  const {return fTheta;};
    void SetMCTheta(double theta){fMCTheta.push_back(theta);};
    std::vector<double>  GetMCTheta()  const {return fMCTheta;};
    void SetPhi(double phi){fPhi.push_back(phi);};
    std::vector<double>  GetPhi()  const {return fPhi;};
    void SetMCPhi(double phi){fMCPhi.push_back(phi);};
    std::vector<double>  GetMCPhi()  const {return fMCPhi;};
    void SetIDTracks(int idTracks) {fIDTracks.push_back(idTracks);};
    std::vector<int> GetIDTracks() {return fIDTracks;};
    void SetCharge(int charge){fCharge.push_back(charge);};
    std::vector<int> GetCharge() const {return fCharge;};
    void SetCPA(double cpa) {fCPA=cpa;};
    double GetCPA() const {return fCPA;};
    void SetParticleOrigin(PartOrigin org){fOrigin=org;};
    PartOrigin GetParticleOrigin() const {return fOrigin;};
    void SetPDGCode(int pdgCode) {fPDGCode=pdgCode;};
    int GetPDGCode() const {return fPDGCode;};
    void SetMCPDGCode(int pdgCode) {fMCPDGCode=pdgCode;};
    int GetMCPDGCode() const {return fMCPDGCode;};
    void SetPDGMotherWeak(int PDGMother){fPDGMotherWeak=PDGMother;};
    int GetMotherWeak() const {return fPDGMotherWeak;};
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
    double fTrackChi2perNDF;
    double fMom;
    double fPt;
    double fMCPt;
    double fP_TPC;
    std::vector<double> fEta;
    std::vector<double> fTheta;
    std::vector<double> fMCTheta;
    std::vector<double> fPhi;
    std::vector<double> fMCPhi;
    std::vector<int> fIDTracks;
    std::vector<int> fCharge;
    double fCPA;
    PartOrigin fOrigin;
    // pdg code as set by the track cuts, used for invariant mass calculation/mc matching in v0s
    int fPDGCode;
    // pdg code as obtained from the MC for this particle
    int fMCPDGCode;
    int fPDGMotherWeak;
    bool fIsMC;
    bool fUse;
    bool fIsSet;
    ClassDef(AliLightNBasePart,1);
};

#endif /* ALILIGHTNBASEPART_H */
