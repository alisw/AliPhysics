#ifndef AliSigma0PhotonCuts_H
#define AliSigma0PhotonCuts_H

#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliSigma0ParticleBase.h"
#include "AliV0ReaderV1.h"
#include "AliVEvent.h"
#include "Riostream.h"
#include "TObject.h"

#include "AliStack.h"

#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TProfile.h"

class AliPIDResponse;

class AliSigma0PhotonCuts : public TObject {
 public:
  AliSigma0PhotonCuts();
  AliSigma0PhotonCuts(const AliSigma0PhotonCuts &);
  AliSigma0PhotonCuts &operator=(const AliSigma0PhotonCuts &);
  virtual ~AliSigma0PhotonCuts();

  static AliSigma0PhotonCuts *DefaultCuts();

  void SelectPhotons(AliVEvent *inputEvent, AliMCEvent *mcEvent,
                     std::vector<AliAODConversionPhoton> &photon);
  double getCosineOfPointingAngle(const AliConversionPhotonBase *photon,
                                  AliVEvent *event) const;
  bool IsPhotonCandidate(AliAODConversionPhoton *TruePhotonCandidate) const;

  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetV0ReaderName(TString name) { fV0ReaderName = name; }
  void SetDCArMax(float dcaR) { fDCAr = dcaR; }
  void SetDCAzMax(float dcaZ) { fDCAz = dcaZ; }
  void SetCosPA(float cosPA) { fCosPA = cosPA; }
  void SetPmax(float pmax) { fPmax = pmax; }
  void SetPtmax(float ptmax) { fPtmax = ptmax; }
  void SetPtElemax(float ptelemax) { fPtElemax = ptelemax; }

  void InitCutHistograms();
  TList *GetCutHistograms() const { return fHistograms; }

 protected:
  TList *fHistograms;
  TList *fHistogramsMC;

  AliV0ReaderV1
      *fV0Reader;         //! basic photon Selection Task !!! needs //! for grid
  TString fV0ReaderName;  //
  TClonesArray *
      fReaderGammas;  //! Array with conversion photons selected by V0Reader Cut

  bool fIsMC;

  AliVEvent *fInputEvent;  //!
  AliMCEvent *fMCEvent;    //!
  AliStack *fMCStack;      //!

  float fDCAr;
  float fDCAz;
  float fCosPA;
  float fPmax;
  float fPtmax;
  float fPtElemax;

  AliPIDResponse *fPIDResponse;  //! pid response

  // Histograms
  // =====================================================================
  TProfile *fHistCuts;  //

  TH1F *fHistPhotonCuts;        //
  TH1F *fHistNPhoton;           //
  TH1F *fHistPhotonPt;          //
  TH1F *fHistPhotonP;           //
  TH1F *fHistPhotonInvMass;     //
  TH2F *fHistPhotonInvMassPt;   //
  TH2F *fHistPhotonInvMassEta;  //
  TH2F *fHistPhotonCPAPt;       //
  TH2F *fHistPhotonR;           //
  TH2F *fHistPhotonArm;         //
  TH1F *fHistPhotonDCAz;        //
  TH1F *fHistPhotonDCAr;        //
  TH2F *fHistPhotonEtaPhi;      //
  TH1F *fHistPhotonConvPointX;  //
  TH1F *fHistPhotonConvPointY;  //
  TH1F *fHistPhotonConvPointZ;  //
  TH1F *fHistPhotonEleP;        //
  TH1F *fHistPhotonElePt;       //

  TH2F *fHistPhotonEleNsigmaTPC;      //
  TH2F *fHistPhotonEleNsigmaTPCPion;  //
  TH2F *fHistPhotonEleTPCsignal;      //

  TH1F *fHistTwoPhotonPt;          //
  TH1F *fHistTwoPhotonP;           //
  TH1F *fHistTwoPhotonInvMass;     //
  TH2F *fHistTwoPhotonInvMassPt;   //
  TH2F *fHistTwoPhotonInvMassEta;  //

  TH1F *fHistMCRecSigma0PhotonPt;          //
  TH1F *fHistMCRecSigma0PhotonP;           //
  TH1F *fHistMCRecSigma0PhotonInvMass;     //
  TH2F *fHistMCRecSigma0PhotonInvMassPt;   //
  TH2F *fHistMCRecSigma0PhotonInvMassEta;  //
  TH2F *fHistMCRecSigma0PhotonR;           //
  TH2F *fHistMCRecSigma0PhotonArm;         //
  TH1F *fHistMCRecSigma0PhotonDCAz;        //
  TH1F *fHistMCRecSigma0PhotonDCAr;        //
  TH1F *fHistMCRecSigma0PhotonConvPointX;  //
  TH1F *fHistMCRecSigma0PhotonConvPointY;  //
  TH1F *fHistMCRecSigma0PhotonConvPointZ;  //
  TH1F *fHistMCRecSigma0PhotonEleP;        //
  TH1F *fHistMCRecSigma0PhotonElePt;       //
  TH2F *fHistMCRecSigma0CPAPt;
  TH2F *fHistMCRecSigma0PhotonDCAzPt;
  TH2F *fHistMCRecSigma0PhotonDCArPt;
  TH2F *fHistMCRecSigma0PhotonPsiPairPt;
  TH2F *fHistMCRecSigma0PhotonChi2Pt;

  TH1F *fHistMCFakeSigma0PhotonPt;          //
  TH1F *fHistMCFakeSigma0PhotonP;           //
  TH1F *fHistMCFakeSigma0PhotonInvMass;     //
  TH2F *fHistMCFakeSigma0PhotonInvMassPt;   //
  TH2F *fHistMCFakeSigma0PhotonInvMassEta;  //
  TH2F *fHistMCFakeSigma0PhotonR;           //
  TH2F *fHistMCFakeSigma0PhotonArm;         //
  TH1F *fHistMCFakeSigma0PhotonDCAz;        //
  TH1F *fHistMCFakeSigma0PhotonDCAr;        //
  TH1F *fHistMCFakeSigma0PhotonConvPointX;  //
  TH1F *fHistMCFakeSigma0PhotonConvPointY;  //
  TH1F *fHistMCFakeSigma0PhotonConvPointZ;  //
  TH1F *fHistMCFakeSigma0PhotonEleP;        //
  TH1F *fHistMCFakeSigma0PhotonElePt;       //
  TH2F *fHistMCFakeSigma0CPAPt;
  TH2F *fHistMCFakeSigma0PhotonDCAzPt;
  TH2F *fHistMCFakeSigma0PhotonDCArPt;
  TH2F *fHistMCFakeSigma0PhotonPsiPairPt;
  TH2F *fHistMCFakeSigma0PhotonChi2Pt;

  TH1F *fHistMCFakePhotonPt;          //
  TH1F *fHistMCFakePhotonP;           //
  TH1F *fHistMCFakePhotonInvMass;     //
  TH2F *fHistMCFakePhotonInvMassPt;   //
  TH2F *fHistMCFakePhotonInvMassEta;  //
  TH2F *fHistMCFakePhotonR;           //
  TH2F *fHistMCFakePhotonArm;         //
  TH1F *fHistMCFakePhotonDCAz;        //
  TH1F *fHistMCFakePhotonDCAr;        //
  TH1F *fHistMCFakePhotonConvPointX;  //
  TH1F *fHistMCFakePhotonConvPointY;  //
  TH1F *fHistMCFakePhotonConvPointZ;  //
  TH1F *fHistMCFakePhotonEleP;        //
  TH1F *fHistMCFakePhotonElePt;       //
  TH2F *fHistMCFakePhotonCPAPt;
  TH2F *fHistMCFakePhotonDCAzPt;
  TH2F *fHistMCFakePhotonDCArPt;
  TH2F *fHistMCFakePhotonPsiPairPt;
  TH2F *fHistMCFakePhotonChi2Pt;

 private:
  ClassDef(AliSigma0PhotonCuts, 1)
};

#endif
