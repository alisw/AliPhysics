/*
 *****************************************************************************************/
#ifndef ALIANALYSISTASKSEPPBCORRELATIONSJETV2
#define ALIANALYSISTASKSEPPBCORRELATIONSJETV2


#include "AliAnalysisTask.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliVAD.h"
#include "AliVParticle.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF2.h"
#include "THnSparse.h"
#include "TString.h"
#include "AliEventCuts.h"
#include "AliTrigAssoPairST.h"


class TList;
class AliCFContainer;
// class AliTHn;
class AliAODEvent;
class AliEventPoolManager;
class AliVParticle;
class AliPIDResponse;
// class AliPID;
class AliAODv0;
// class AliAnalyseLeadingTrackUE;
class THnSparse;
class AliAODcascade;
class AliAODVertex;


#ifndef ALIANALYSISTASKSEH
#include "AliAnalysisTaskSE.h"
#endif

//---------------------------------------------------------------------------------------
class AliAnalysisTaskSEpPbCorrelationsJetV2 : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskSEpPbCorrelationsJetV2();
  AliAnalysisTaskSEpPbCorrelationsJetV2(const char *name);
  virtual ~AliAnalysisTaskSEpPbCorrelationsJetV2();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  // Analysis Correlation
  void MakeAna();

  // configration
  virtual void SetAnalysisMode(TString mode) { fAnaMode = mode; }
  virtual void SetAssociatedTrack(TString mode) { fasso = mode; }
  virtual void SetDatatype(Bool_t mode) { fDataType = mode; }
  virtual void SetRunType(Bool_t mode) { frun2 = mode; }
  virtual void SetFilterBit(Int_t mode) { ffilterbit = mode; }
  virtual void SetFMDcut(Bool_t mode) {fFMDcut=mode;}
  virtual void SetFMDcutpar(Int_t mode){fFMDcutmode=mode;}
  virtual void SetReduceDphi(Double_t mode){fReduceDphi=mode;}
  virtual void SetSymmetricFMD(Bool_t mode){fSymmetricFMD=mode;}
  virtual void SetLikeSign(Bool_t mode){fLikeSign=mode;}
  virtual void SetPtMax(Float_t mode)   {fPtMax=mode;}
  virtual void SetPtMin(Float_t mode)   {fPtMin=mode;}
  virtual void SetAssoCut(Float_t mode) {fAsscoptCut=mode;}
  virtual void SetAnalysisCent(TString mode) { fCentType = mode; }
  virtual void SetAnalysisCollisionType(TString mode) { fcollisiontype = mode; }
  virtual void SetmcprimFMD(Bool_t mode){fprimFMD=mode;}
  virtual void SetmcprimTPC(Bool_t mode){fprimTPC=mode;}
  virtual void SetCentCalib(Bool_t mode) { fcentcalib= mode; }


  void SetMaxNEventsInPool(Int_t events) { fPoolMaxNEvents = events; }
  void SetMinNTracksInPool(Int_t tracks) { fPoolMinNTracks = tracks; }
  void SetMinEventsToMix(Int_t events) { fMinEventsToMix = events; }
  void SetCentrality(Double_t cenMin, Double_t cenMax) {fCenMin = cenMin; fCenMax = cenMax;}
  void SetPoolPVzBinLimits(Int_t Nzvtxbins, const Double_t *ZvtxBins) {
    fNzVtxBins = Nzvtxbins;
    for (int ix = 0; ix < fNzVtxBins + 1; ix++) {
      fZvtxBins[ix] = ZvtxBins[ix];
    }
  }

  void SetPoolCentBinLimits(Int_t Ncentbins, const Double_t *CentBins) {
    fNCentBins = Ncentbins;
    for (int ix = 0; ix < fNCentBins + 1; ix++) {
      fCentBins[ix] = CentBins[ix];
    }
  }
 

  void DumpTObjTable(const char* note);

  
private:

  AliAnalysisTaskSEpPbCorrelationsJetV2(
      const AliAnalysisTaskSEpPbCorrelationsJetV2 &det);
  AliAnalysisTaskSEpPbCorrelationsJetV2 &
  operator=(const AliAnalysisTaskSEpPbCorrelationsJetV2 &det);

  void DefineGeneralOutput();
  void DefineVZEROOutput();
  void DefineCorrOutput();
  void DefinedQAHistos();

  TObjArray *GetAcceptedTracksLeading(AliAODEvent *faod,Bool_t leading,TObjArray*tracks);
  TObjArray *GetAcceptedTracksPID(AliAODEvent *faod);
  TObjArray *GetAcceptedV0Tracks(const AliAODEvent *faod);
  TObjArray *GetAcceptedCascadeTracks(AliAODEvent *faod);
  TObjArray *GetAcceptedTracksAssociated(AliAODEvent *faod);

  TObjArray* GetFMDhitsYS(Bool_t Aside);
  Bool_t IsAcceptedDaughterTrack(const AliAODTrack *itrack);
  Bool_t IsAcceptedPhiDaughterTrack(const AliAODTrack *itrack);
  Bool_t IsAcceptedTrack(const AliAODTrack *aodTrack);
  Bool_t IsAcceptedDecayLength(const AliAODv0 *aodv0,Double_t mass,Double_t maxctau);
  Bool_t IsAcceptedV0(const AliAODv0 *aodv0);
  Bool_t IsAcceptedCascade(const AliAODcascade *casc);
  Bool_t IsAcceptedCascadeOmega(const AliAODcascade *casc);

  TObjArray* CloneTrack(TObjArray* track);
  Double_t RangePhi(Double_t DPhi);
  Double_t RangePhi_FMD(Double_t DPhi);
  Double_t RangePhi2(Double_t DPhi);
  Bool_t NotSPDClusterVsTrackletBG() {return !fUtils.IsSPDClusterVsTrackletBG(this->InputEvent());};

  void FillCorrelationTracks(Double_t MultipOrCent, TObjArray *triggerArray, 
                             TObjArray *selectedTrackArray, TObjArray *selectedTrackArray_TPC, AliTHn *triggerHist, AliTHn *associateHist, Float_t bSign, Int_t step, TObjArray *select_TPCPairs);

  void FillCorrelationTracksMixing(Double_t MultipOrCentMix, Double_t pvxMix,
                                   TObjArray *triggerArray,
                                   TObjArray *selectedTrackArray,
                                   TObjArray *selectedTrackArray_TPC,
                                   AliTHn *triggerHist, AliTHn *associateHist,
                                   Float_t, Int_t, TObjArray *select_TPCPairs);



  inline Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1,
                             Float_t phi2, Float_t pt2, Float_t charge2,
                             Float_t radius, Float_t bSign);

  TString fcollisiontype;
  Bool_t fDataType;
  Bool_t frun2;
  Bool_t fQA;
  Bool_t fFMDcut;
  Int_t fFMDcutmode;
  Double_t fReduceDphi;
  TString fAnaMode;
  TString fasso;
  Bool_t fSymmetricFMD;
  Bool_t fLikeSign;

  TString fCentType;
  Bool_t fprimTPC;
  Bool_t fprimFMD;
  Bool_t fcentcalib;
  Int_t fNEntries;
  
  Double_t lCentrality;
  Double_t fCenMax;
  Double_t fCenMin;

  Float_t bSign;
  Double_t fZVertex;

  TList *fOutputList;  // Output list
  TList *fOutputList1; // Output list
  TList *fOutputList2; // Output list

  //AliPIDResponse *fPIDResponse; // PID Response

  Int_t ffilterbit;
  Float_t fnoClusters;
  Float_t fCutChargedDCAzMax;
  Float_t fCutChargedDCAxyMax;  
  Double_t fPtMin;
  Double_t fPtMax;
  Double_t fAsscoptCut;
  Double_t fEtaMax;
  Double_t fEtaMaxExtra;
  Double_t fEtaMinExtra;
  // V0 particles
  // V0 daughter cut
  
  AliEventCuts fEventCuts; 
  AliAnalysisUtils fUtils;
  AliAODEvent *fEvent; //  AOD Event
  AliMCEvent* mcEvent;
  AliAODVertex *lPrimaryBestVtx;
  Double_t fPrimaryZVtx;

  AliAODVZERO *fvzero;

  // Event Pool for mixing
  AliEventPoolManager *fPoolMgr;  //  event pool manager for Event Mixing
  Int_t fPoolMaxNEvents;   // set maximum number of events in the pool
  Int_t fPoolMinNTracks;   // set minimum number of tracks in the pool
  Int_t fMinEventsToMix;   // set the minimum number of events want to mix
  Int_t fNzVtxBins;        // number of z vrtx bins
  Double_t fZvtxBins[100]; // [fNzVtxBinsDim]
  Int_t fNCentBins;        // number of centrality bins
  Double_t fCentBins[100]; // [fNCentBinsDim]
  // Track cuts

  // Globaal Histograms
  TH1F *fHistzvertex;
  TH1F *fHistCentrality;
  TH1F *fHistCentrality_beforecut;
  TH1F *fHistCentrality_beforeFMDMulcut;
  TH2F* fHistCentzvertex;
  TH2F* fHistCentV0vsTracklets;
  TH2F* fHistCentvsTrackletsbefore;
  TH2F* mixedDist;
  TH2F* mixedDist2;
  
  
  //AliTHn *fHistLeadQA;

  AliTHn* fhistmcprim;
  TH2D*fhmcprimvzeta;

  TH1F*frefetac;
  TH1F*frefetaa;
  TH1F*frefvz;
  TH1D*fhcorr[10];

  TH1D*fhrefetaFMD[4];
  TH1D*fhrefphiFMD[4];

  TH2D*  fh2_pt_trig_asso;
  TH2D*  fh2_FMD_eta_phi_prim;
  TH2D*  fh2_FMD_acceptance;
  TH2F*  fh2_SPD_multcorr;
  TH2F*  fh2_SPDtrack_multcorr;
  TH1F*  fhtrackletsdphi;
  TH2D*  fh2_FMD_eta_phi;
  TH2D*  fh2_FMD_eta_phi_afterCut;
  TH1F* fHist_NeventRun;
  TH1F* fHist_V0AMultRun;
  TH1F* fHist_V0CMultRun;
  TH1F* fHist_FMDAMultRun;
  TH1F* fHist_FMDCMultRun;

  TH2D*  fhistfmdphiacc;
  TH2F* fhFMDmultchannel;
  TH2F* fhFMDmultchannel_actual;
  TH2D* fhFMDmult_runbyrun_cside[31];
  TH2D* fhFMDmult_runbyrun_aside[65];
  
  AliTHn* fhistfmd;
  THnSparseF* fhistits;
  AliTHn* fhSecFMD;
  //  const TH2D& d2Ndetadphi;
  TH2F*fFMDV0;
  TH2F*fFMDV0_post;
  TH2F*fFMDV0A;
  TH2F*fFMDV0A_post;
  TH2F*fFMDV0C;
  TH2F*fFMDV0C_post;

  TH2F*fFMDV0same;
  TH2F*fFMDV0same_post;
  TH2F*fFMDV0Asame;
  TH2F*fFMDV0Asame_post;
  TH2F*fFMDV0Csame;
  TH2F*fFMDV0Csame_post;


  TH1F *fHist_Stat;

  AliTHn *fHistTriggerTrack;
  AliTHn *fHistReconstTrack;
  AliTHn *fHistTriggerTrackMix;
  AliTHn *fHistReconstTrackMix;

  // only for test of dphi
//  TH2D *fTPCTPCdphi_deta_4_2;
/*
  TH1D *fTPCTPCdphi_deta_4_1;
  TH1D *fTPCTPCdphi_deta_4_0;
  TH1D *fTPCTPCdphi_deta_3_2;
  TH1D *fTPCTPCdphi_deta_3_1;
  TH1D *fTPCTPCdphi_deta_3_0;
  TH1D *fTPCTPCdphi_deta_2_2;
*/

  ClassDef(AliAnalysisTaskSEpPbCorrelationsJetV2, 2);
};
//---------------------------------------------------------------------------------------
Float_t AliAnalysisTaskSEpPbCorrelationsJetV2::GetDPhiStar(
    Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2,
    Float_t charge2, Float_t radius, Float_t bSign) {
  //
  // calculates dphistar (Copied from AliUEHistogram)
  //

  Float_t dphistar = phi1 - phi2 -
                     charge1 * bSign * TMath::ASin(0.075 * radius / pt1) +
                     charge2 * bSign * TMath::ASin(0.075 * radius / pt2);

  static const Double_t kPi = TMath::Pi();


  if (dphistar > kPi)
    dphistar = kPi * 2 - dphistar;
  if (dphistar < -kPi)
    dphistar = -kPi * 2 - dphistar;
  if (dphistar > kPi) // might look funny but is needed
    dphistar = kPi * 2 - dphistar;

  return dphistar;
}

#endif
