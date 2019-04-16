#ifndef ALIANALYSISTASKSEHFSYSTNSIGMAPID_H
#define ALIANALYSISTASKSEHFSYSTNSIGMAPID_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. */

/////////////////////////////////////////////////////////////////////////////////////////
// \class AliAnalysisTaskSEHFSystPID                                                   //
// \brief analysis task for the study of PID systematic uncertainties of HF particles  //
// \author: A. M. Barbano, anastasia.maria.barbano@cern.ch                             //
// \author: F. Grosa, fabrizio.grosa@cern.ch                                           //
/////////////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>

#include "AliAnalysisTaskSE.h"
#include "AliAODv0KineCuts.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"

using namespace std;

class AliAnalysisTaskSEHFSystPID : public AliAnalysisTaskSE {
  
public:

  enum tagflags {
    kIsPionFromK0s       = BIT(0),
    kIsPionFromL         = BIT(1),
    kIsProtonFromL       = BIT(2),
    kIsElectronFromGamma = BIT(3),
    kIsKaonFromKinks     = BIT(4),
    kIsKaonFromTOF       = BIT(5),
    kIsKaonFromTPC       = BIT(6)
  };

  enum centest {
    kCentOff,
    kCentV0M,
    kCentV0A,
    kCentZNA,
    kCentCL0,
    kCentCL1
  };

  enum SystemForNsigmaDataCorr {
    kNone=-1,
    kPbPb010,
    kPbPb3050
  };

  AliAnalysisTaskSEHFSystPID();
  AliAnalysisTaskSEHFSystPID(const char *name, int system=0);
  virtual ~AliAnalysisTaskSEHFSystPID();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);

  void SetReadMC(bool flag = true)                                            {fIsMC = flag;}
  void SetCentralityLimits(int mincent, int maxcent)                          {fCentMin = mincent; fCentMax = maxcent;}
  void SetCentralityEstimator(int centest=kCentV0M)                           {fCentEstimator=centest;}
  void SetESDtrackCuts(AliESDtrackCuts * trackCuts)                           {fESDtrackCuts = trackCuts;}
  void SetTriggerInfo(TString trigClass, unsigned long long mask=0)           {fTriggerClass = trigClass; fTriggerMask = mask;}
  void SetNsigmaKaonForTagging(float nsigmamax = 0.02)                        {fNsigmaMaxForTag=nsigmamax;}
  void SetKinksSelections(float qtmin=0.15, float Rmin=120, float Rmax=210)   {fQtMinKinks=qtmin; fRMinKinks=Rmin; fRMaxKinks=Rmax;}
  void SetfFillTreeWithNsigmaPIDOnly(bool fillonlyNsigma=true)                {fFillTreeWithNsigmaPIDOnly=fillonlyNsigma;}
  void EnableDownSampling(double fractokeep=0.1, double ptmax=1.5)            {fEnabledDownSampling=true; fFracToKeepDownSampling=fractokeep; fPtMaxDownSampling=ptmax;}
  void SetAODMismatchProtection(int opt=1)                                    {fAODProtection=opt;}
  
  void EnableNsigmaDataDrivenCorrection(int syst) {
    fEnableNsigmaTPCDataCorr = true;
    fSystNsigmaTPCDataCorr = syst;
  }

private:

  bool IsVertexAccepted();
  bool IsCentralitySelected();
  void GetTaggedV0s(vector<short> &idPionFromK0s, vector<short> &idPionFromL, vector<short> &idProtonFromL, vector<short> &idElectronFromGamma);
  short GetPDGcodeFromMC(AliAODTrack* track, TClonesArray* arrayMC);
  AliAODTrack* IsKinkDaughter(AliAODTrack* track);
  void GetTaggedKaonsFromKinks(vector<short> &idKaonFromKinks);
  float MaxOpeningAngleKnu(float p);
  float GetTOFmomentum(AliAODTrack* track);
  short ConvertFloatToShort(float num);
  unsigned short ConvertFloatToUnsignedShort(float num);
  void GetNsigmaTPCMeanSigmaData(float &mean, float &sigma, AliPID::EParticleType species, float pTPC);
  void SetNsigmaTPCDataCorr(int run);

  enum hypos{kPion,kKaon,kProton};
  static const int kNHypo = 3;
  const TString hyponames[kNHypo] = {"Pion","Kaon","Proton"};

  const float kCSPEED = 2.99792457999999984e-02; // cm / ps

  TList *fOutputList;                              //!<! output list for histograms

  TH1F *fHistNEvents;                              //!<! histo with number of events
  TH2F *fHistArmenteroPlot[5];                     //!<! histo for armenteros-podolanski plot
  TH2F *fHistQtVsMassKinks;                        //!<! histo for mother-kink qt vs. mass distribution
  TH2F *fHistPDaughterVsMotherKink;                //!<! histo for pT daughter vs. pT mother kink
  TH2F *fHistdEdxVsPMotherKink;                    //!<! histo for mother kink TPC dEdx vs. p
  TH2F *fHistOpeningAngleVsPMotherKink;            //!<! histo for opening angle vs. pT mother kink
  TH2F *fHistNTPCclsVsRadius;                      //!<! histo for nTPC clusters vs. R mother kink
  TH2F *fHistNsigmaTPCvsPt[kNHypo];                //!<! array of histos for nsigmaTPC vs pt (MC truth)
  TH2F *fHistNsigmaTOFvsPt[kNHypo];                //!<! array of histos for nsigmaTPC vs pt (MC truth)
  TTree* fPIDtree;                                 //!<! tree with PID info

  short fPIDNsigma[6];                             /// Nsigma PID to fill the tree
  unsigned short fPTPC;                            /// TPC momentum to fill the tree
  unsigned short fPTOF;                            /// TOF momentum to fill the tree
  unsigned short fdEdxTPC;                         /// TPC dEdX to fill the tree
  unsigned short fToF;                             /// ToF signal to fill the tree
  unsigned short fPt;                              /// transverse momentum to fill the tree
  unsigned char fTPCNcls;                          /// number of clusters in TPC to fill the tree
  unsigned char fTPCNclsPID;                       /// number of PID clusters in TPC to fill the tree
  unsigned short fTrackLength;                     /// track length for TOF PID
  unsigned short fStartTimeRes;                    /// start time resolution for TOF PID
  short fPDGcode;                                  /// PDG code in case of MC to fill the tree
  unsigned char fTag;                              /// bit map for tag (see enum above)
  float fNsigmaMaxForTag;                          /// max nSigma value to tag kaons
  float fQtMinKinks;                               /// min qt for kinks
  float fRMinKinks;                                /// min radius in XY for kinks
  float fRMaxKinks;                                /// max radius in XY for kinks

  float fCentMin;                                  /// min centrality
  float fCentMax;                                  /// max centrality
  int fCentEstimator;                              /// centrality estimator
  TString fTriggerClass;                           /// trigger class
  unsigned long long fTriggerMask;                 /// trigger mask
  bool fIsMC;                                      /// flag to switch on the MC analysis for the efficiency estimation
  int fSystem;                                     /// system: 0->pp,pPb 1->PbPb

  AliESDtrackCuts * fESDtrackCuts;                 /// single-track cut set 
  AliAODEvent *fAOD;                               /// AOD object
  AliPIDResponse *fPIDresp;                        /// basic pid object
  AliAODv0KineCuts *fV0cuts;                       /// AOD V0 cuts

  bool fFillTreeWithNsigmaPIDOnly;                 /// flag to enable filling of the tree with only Nsigma variables for the PID
  bool fEnabledDownSampling;                       /// flag to enable/disable downsampling
  double fFracToKeepDownSampling;                  /// fraction to keep when downsampling activated
  double fPtMaxDownSampling;                       /// pT max of tracks to downsample
  
  int fAODProtection;                              /// flag to activate protection against AOD-dAOD mismatch

  int fRunNumberPrevEvent;                         /// run number of previous event
  bool fEnableNsigmaTPCDataCorr;                   /// flag to enable data-driven NsigmaTPC correction
  int fSystNsigmaTPCDataCorr;                      /// system for data-driven NsigmaTPC correction
  float fMeanNsigmaTPCPionData[100];               /// array of NsigmaTPC pion mean in data 
  float fMeanNsigmaTPCKaonData[100];               /// array of NsigmaTPC kaon mean in data 
  float fMeanNsigmaTPCProtonData[100];             /// array of NsigmaTPC proton mean in data 
  float fSigmaNsigmaTPCPionData[100];              /// array of NsigmaTPC pion mean in data 
  float fSigmaNsigmaTPCKaonData[100];              /// array of NsigmaTPC kaon mean in data 
  float fSigmaNsigmaTPCProtonData[100];            /// array of NsigmaTPC proton mean in data 
  float fPlimitsNsigmaTPCDataCorr[101];            /// array of p limits for data-driven NsigmaTPC correction
  int fNPbinsNsigmaTPCDataCorr;                    /// number of p bins for data-driven NsigmaTPC correction

  ClassDef(AliAnalysisTaskSEHFSystPID, 3);
};

#endif
