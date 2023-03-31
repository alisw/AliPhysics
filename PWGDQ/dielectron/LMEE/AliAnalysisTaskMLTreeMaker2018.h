#ifndef ALIANALYSISTASKMLTREEMAKER2018_H
#define ALIANALYSISTASKMLTREEMAKER2018_H
class TH1F;
class TList;
class TH2D;
class TH3D;
//class AliESDtrackCuts;
#endif

#include "AliAnalysisTaskSE.h"
#include "AliDielectronVarCuts.h"
#include "AliDielectronTrackCuts.h"
#include "AliDielectronCutGroup.h"
#include "AliDielectronPID.h"
#include "AliAnalysisFilter.h"
#include "AliDielectronEventCuts.h"
#include "AliDielectronSignalMC.h"
#include "AliDielectronPair.h"

#define PRECISION 1e-6

#ifndef ALIANALYSISTASKSE_H
#endif

// Authors: Jerome Jung (Uni Frankfurt) - jerome.jung@cern.ch




class AliAnalysisTaskMLTreeMaker2018 : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskMLTreeMaker2018(const char *name);
  AliAnalysisTaskMLTreeMaker2018();
  ~AliAnalysisTaskMLTreeMaker2018(); 

  AliDielectronEventCuts* eventCuts;
  AliAnalysisFilter* evfilter;

  AliDielectronVarCuts* trcuts;
  AliDielectronTrackCuts *trfilter;
  AliDielectronPID *pidcuts;
  AliDielectronCutGroup* cuts;
  AliAnalysisFilter* filter;

  // need this to use PID in dielectron framework
  AliDielectronVarManager* varManager;

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   FinishTaskOutput();
  virtual void   Terminate(Option_t *);
  

  void SetGeneratorName(TString name) {fGeneratorName = name;}
  void SetGeneratorMCSignalName(TString name) {fGeneratorMCSignalName = name;}
  void SetGeneratorULSSignalName(TString name) {fGeneratorULSSignalName = name;}

  void SetupTrackCuts(AliDielectronCutGroup* f);
  void SetupEventCuts(AliDielectronEventCuts* f);

  void runningOnAOD(Bool_t AOD){isAOD=AOD;}
  void isMC(Bool_t isMC){hasMC=isMC;}
  void doULSpairs(Bool_t ULS){doULS=ULS;}
  void doLSpairs(Bool_t LS){doLS=LS;}
  void doSingleLegMCSignals(Bool_t singleMC){doSingleLegMCSignal = singleMC;}

  void   AddSingleLegMCSignal(AliDielectronSignalMC signal1)         {fSingleLegMCSignal.push_back(signal1);}
  void   AddPairMCSignal(AliDielectronSignalMC signal1)              {fPairMCSignal.push_back(signal1);}
  void   AddMCSignalsWhereDielectronPairNotFromSameMother(std::vector<bool> vec) {fDielectronPairNotFromSameMother = vec;}


  void filterTrackQuality(Bool_t isQuality){filterQuality=isQuality;}
  void SetupPairing(Bool_t pairing) {doPairing=pairing;}
  void SetupTracks(Bool_t tracks) {doTracks=tracks;}

  void SetCentralityPercentileRange(Double_t min, Double_t max){
    fCentralityPercentileMin = min;
    fCentralityPercentileMax = max;

    Printf("Thresholds Set");
    Printf("cent = %f - %f",fCentralityPercentileMin,fCentralityPercentileMax);
  }

  void SetPtRange(Double_t min, Double_t max){
    fPtMin = min;
    fPtMax = max;

    Printf("Thresholds Set");
    Printf("pT = %f - %f",fPtMin,fPtMax);
  }

  void SetEtaRange(Double_t min, Double_t max){
    fEtaMin = min;
    fEtaMax = max;

    Printf("Thresholds Set");
    Printf("eta = %f - %f",fEtaMin,fEtaMax);
  }

  void SetFilterBit(Int_t filterBit){
    fFilterBit = filterBit;
  }

  class Particle{
    public:
      Particle() :
        fPt(-99), fEta(-99), fPhi(-99), fCharge(-99), fXv(-99), fYv(-99), fZv(-99), fDCAxy(-99), fDCAz(-99), fDCAxy_res(-99), fDCAz_res(-99), fTrackID(0), fTrackLabel(0), fMotherID(0), fMCSignalPair(false), fULSSignalPair(false), isMCSignal(), DielectronPairFromSameMother() {}
      Particle(double pt, double eta, double phi, short charge, float x, float y, float z) :
        fPt(pt), fEta(eta), fPhi(phi), fCharge(charge), fXv(x), fYv(y), fZv(z), fDCAxy(-99), fDCAz(-99), fDCAxy_res(-99), fDCAz_res(-99), fTrackID(0), fTrackLabel(0), fMotherID(0), fMCSignalPair(false), fULSSignalPair(false), isMCSignal(), DielectronPairFromSameMother() {}

      void SetTrackID(int id) {fTrackID = id;}
      void SetTrackLabel(int label) {fTrackLabel = label;}

      void SetMotherID(int id) {fMotherID = id;}
      void SetMCSignalPair (bool value) {fMCSignalPair = value;}
      void SetULSSignalPair(bool value) {fULSSignalPair = value;}
      void SetDielectronPairFromSameMother(std::vector<Bool_t> vec){DielectronPairFromSameMother = vec;}

      void SetDCA(double xy, double z) {fDCAxy=xy; fDCAz=z;}
      void SetDCAres(double xy, double z) {fDCAxy_res=xy; fDCAz_res=z;}

      int  GetTrackID() {return fTrackID;}
      int  GetMotherID() {return fMotherID;}
      int  GetTrackLabel() {return fTrackLabel;}
      bool GetMCSignalPair() {return fMCSignalPair;}
      bool GetULSSignalPair() {return fULSSignalPair;}


      float fPt;
      float fEta;
      float fPhi;
      short fCharge;
   
      float fXv;
      float fYv;
      float fZv;

      float fDCAxy;
      float fDCAz;
      float fDCAxy_res;
      float fDCAz_res;

      int     fTrackID;
      int     fTrackLabel;

      int     fMotherID;
      bool    fMCSignalPair;
      bool    fULSSignalPair;
      std::vector<Bool_t> isMCSignal;
      std::vector<Bool_t> DielectronPairFromSameMother;
  };



 private:

  AliPIDResponse *fPIDResponse;     //! PID response object

  //AliPIDCombined *fPIDCombined;

  std::vector<Float_t> eta;
  std::vector<Float_t> phi;
  std::vector<Float_t> pt;
  std::vector<Int_t> charge;

  std::vector<Int_t> NClustersITS;
  std::vector<Float_t> NCrossedRowsTPC;
  std::vector<Int_t> NClustersTPC;
  std::vector<Bool_t> HasSPDfirstHit;
  std::vector<Float_t> RatioCrossedRowsFindableClusters;
  std::vector<Int_t> NTPCSignal;
  //Bool_t loCuts;        //loose cuts?

//  std::vector<Int_t> IsBG;

  Int_t runNumber;
  Int_t nTracks;
  std::vector<Int_t> nPairs;

  Double_t cent;

  AliAnalysisManager *man;

//  Double_t IsEventAccepted(AliVEvent *event);
  Int_t GetAcceptedTracks(AliVEvent *event, Double_t gCentrality);
  std::tuple<int, int, int> GetAcceptedPairs(AliVEvent *event, Double_t gCentrality);
  Int_t FillSimplePairs(std::vector<Particle> parts1, std::vector<Particle> parts2, AliVEvent *event, Int_t index);
  Int_t FillSignalPairs(AliVEvent *event);
  Bool_t GetDCA(const AliVEvent* event, const AliAODTrack *track, Double_t* d0z0, Double_t* covd0z0);
  AliAnalysisTaskMLTreeMaker2018(const AliAnalysisTaskMLTreeMaker2018&); // not implemented


  void    CheckSingleLegMCsignals(std::vector<Bool_t>& vec, const int track);
  void    CheckPairMCsignals(std::vector<Bool_t>& vec, AliVParticle* part1, AliVParticle* part2);
  bool    CheckGenerator(int trackID, std::vector<unsigned int> vecHashes, Bool_t isGen);

  void CheckIfFromMotherWithDielectronAsDaughter(Particle& part);
  Bool_t CheckIfOneIsTrue(std::vector<Bool_t>& vec);
  Particle CreateParticle(AliVParticle* mcPart1);

  AliAnalysisTaskMLTreeMaker2018& operator=(const AliAnalysisTaskMLTreeMaker2018&); // not implemented


  TList *fList;//output list for QA histograms

  Float_t fCentralityPercentileMin;// minimum centrality threshold (default = 0)
  Float_t fCentralityPercentileMax;// maximum centrality threshold (default = 80)

  Float_t fPtMin;// minimum pT threshold (default = 0)
  Float_t fPtMax;// maximum pT threshold (default = 1000)
  Float_t fEtaMin;// minimum eta threshold (default = -10)
  Float_t fEtaMax;// maximum eta threshold (default = 10)

  Int_t fFilterBit;// track cut bit from track selection (default = kTPCqualSPDany)

//  AliESDtrackCuts* fESDTrackCuts;

  AliMCEvent* mcEvent;
  AliVEvent*   fEvent;
  Bool_t isAOD;
  Bool_t hasMC;
  Bool_t filterQuality;
  Bool_t doPairing;
  Bool_t doTracks;

  Bool_t doULS;
  Bool_t doLS;
  Bool_t doSingleLegMCSignal;

  std::vector<AliDielectronSignalMC> fSingleLegMCSignal;
  std::vector<AliDielectronSignalMC> fPairMCSignal;
  std::vector<bool> fDielectronPairNotFromSameMother; // this is used to get electrons from charmed mesons in a environment where GEANT is doing the decay of D mesons, like in LHC18b5a
  std::vector<Particle> fNegPart;
  std::vector<Particle> fPosPart;

  TString fGeneratorName;
  TString fGeneratorMCSignalName;
  TString fGeneratorULSSignalName;
  std::vector<unsigned int> fGeneratorHashs;
  std::vector<unsigned int> fGeneratorMCSignalHashs;
  std::vector<unsigned int> fGeneratorULSSignalHashs;

  std::vector<Float_t> MCpt;
  std::vector<Float_t> MCeta;
  std::vector<Float_t> MCphi;

  std::vector<Float_t> MCvertx;
  std::vector<Float_t> MCverty;
  std::vector<Float_t> MCvertz;


  std::vector<Int_t> glabel ;
  std::vector<Int_t> gLabelFirstMother ;
  std::vector<Int_t> gLabelMinFirstMother ;
  std::vector<Int_t> gLabelMaxFirstMother ;
  std::vector<Int_t> iGenIndex ;
  std::vector<Int_t> iPdgFirstMother ;

  std::vector<Float_t> dcaXY;    //DCA
  std::vector<Float_t> dcaZ;

  std::vector<Float_t> dcaXY_res;    //DCA
  std::vector<Float_t> dcaZ_res;

  Double_t vertx;
  Double_t verty;
  Double_t vertz;

  std::vector<Int_t> nTPC;
  std::vector<Int_t> nITS;
  std::vector<Double_t> nITSshared;

  std::vector<Int_t> ITS1S;
  std::vector<Int_t> ITS2S;
  std::vector<Int_t> ITS3S;
  std::vector<Int_t> ITS4S;
  std::vector<Int_t> ITS5S;
  std::vector<Int_t> ITS6S;

  std::vector<Float_t> chi2ITS;
  std::vector<Float_t> chi2TPC;
  std::vector<Float_t> chi2GlobalPerNDF;
  std::vector<Float_t> chi2GlobalvsTPC;

  std::vector<Int_t> pdg;
  std::vector<Int_t> pdgmother;
  std::vector<Int_t> hasmother;
  std::vector<Int_t> label;
  std::vector<Int_t> motherlabel;

  //TBits*            fUsedVars;                // used variables by AliDielectronVarManager

//  Double_t probs[AliPID::kSPECIESC];

  TTree* fTreeTracks;
  TTree* fTreeULS;
  TTree* fTreeLSmm;
  TTree* fTreeLSpp;
  TTree* fTreePairSignal;
  std::vector<TTree*> fTreePairs;
  std::vector<TTree*> fTreePairSignals;
  TH1F* fQAHistEvents;
  TH1F* fQAHistTracks;

  Int_t fFillPairSignalOffset;

  std::vector<std::vector<Float_t>> eta_tracks1;
  std::vector<std::vector<Float_t>> phi_tracks1;
  std::vector<std::vector<Float_t>> pt_tracks1;
  std::vector<std::vector<Int_t>> charge_tracks1;

  std::vector<std::vector<Float_t>> MCpt_tracks1;
  std::vector<std::vector<Float_t>> MCeta_tracks1;
  std::vector<std::vector<Float_t>> MCphi_tracks1;

  std::vector<std::vector<Float_t>> MCvertx_tracks1;
  std::vector<std::vector<Float_t>> MCverty_tracks1;
  std::vector<std::vector<Float_t>> MCvertz_tracks1;

  std::vector<std::vector<Float_t>> dcaXY_tracks1;    //DCA
  std::vector<std::vector<Float_t>> dcaZ_tracks1;

  std::vector<std::vector<Float_t>> dcaXY_res_tracks1;    //DCA
  std::vector<std::vector<Float_t>> dcaZ_res_tracks1;

  std::vector<std::vector<Float_t>> verticesX_tracks1;
  std::vector<std::vector<Float_t>> verticesY_tracks1;
  std::vector<std::vector<Float_t>> verticesZ_tracks1;

  std::vector<std::vector<Float_t>> eta_tracks2;
  std::vector<std::vector<Float_t>> phi_tracks2;
  std::vector<std::vector<Float_t>> pt_tracks2;
  std::vector<std::vector<Int_t>> charge_tracks2;

  std::vector<std::vector<Float_t>> MCpt_tracks2;
  std::vector<std::vector<Float_t>> MCeta_tracks2;
  std::vector<std::vector<Float_t>> MCphi_tracks2;

  std::vector<std::vector<Float_t>> MCvertx_tracks2;
  std::vector<std::vector<Float_t>> MCverty_tracks2;
  std::vector<std::vector<Float_t>> MCvertz_tracks2;

  std::vector<std::vector<Float_t>> dcaXY_tracks2;    //DCA
  std::vector<std::vector<Float_t>> dcaZ_tracks2;

  std::vector<std::vector<Float_t>> dcaXY_res_tracks2;    //DCA
  std::vector<std::vector<Float_t>> dcaZ_res_tracks2;

  std::vector<std::vector<Float_t>> verticesX_tracks2;
  std::vector<std::vector<Float_t>> verticesY_tracks2;
  std::vector<std::vector<Float_t>> verticesZ_tracks2;



  std::vector<std::vector<Float_t>> primVerticesX;
  std::vector<std::vector<Float_t>> primVerticesY;
  std::vector<std::vector<Float_t>> primVerticesZ;

  std::vector<std::vector<Float_t>> pairMasses;
  std::vector<std::vector<Float_t>> pairPtees;
  std::vector<std::vector<Float_t>> pairOpAngs;
  std::vector<std::vector<Float_t>> pairDCAee;
  
  std::vector<std::vector<Float_t>> pairCosPointAngs;
  std::vector<std::vector<Float_t>> pairDecayLengths;
  std::vector<std::vector<Float_t>> pairRs;

//  TH2D* fHistTrackStats;//QA histogram for track filter bit statistics vs. centrality

//  TH3D* fHistEtaPhiPt;//QA histogram for eta/phi/pt distribution


  ClassDef(AliAnalysisTaskMLTreeMaker2018, 3);

};
