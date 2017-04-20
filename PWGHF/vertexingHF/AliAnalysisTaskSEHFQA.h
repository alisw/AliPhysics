#ifndef ALIANALYSISTASKSEHFQUALITYASSURANCE_H
#define ALIANALYSISTASKSEHFQUALITYASSURANCE_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
/// \class Class AliAnalysisTaskSEHFQA
/// \brief AliAnalysisTaskSE for HF quality assurance
/// \author Authors: C.Bianchin, chiara.bianchin@pd.infn.it
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1I.h"
#include "TProfile2D.h"

#include "AliAnalysisTaskSE.h"

class AliRDHFCuts;
class TH1F;
class AliAODEvent;
class AliFlowEvent;
class AliFlowTrackCuts;

class AliAnalysisTaskSEHFQA : public AliAnalysisTaskSE
{

 public:

  enum DecChannel {kDplustoKpipi,kD0toKpi,kDstartoKpipi,kDstoKKpi,kD0toKpipipi,kLambdactopKpi,kLambdactoV0};

  AliAnalysisTaskSEHFQA();
  AliAnalysisTaskSEHFQA(const char *name, DecChannel ch, AliRDHFCuts* cuts);
  virtual ~AliAnalysisTaskSEHFQA();

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  /// setters
  void SetReadMC(Bool_t mcflag){fReadMC=mcflag;}
  void SetSimpleMode(Bool_t flag){fSimpleMode=flag;}
  void SetTrackOn(Bool_t trackon=kTRUE){fOnOff[0]=trackon;}
  void SetPIDOn(Bool_t pidon=kTRUE){fOnOff[1]=pidon;}
  void SetCentralityOn(Bool_t centron=kTRUE){fOnOff[2]=centron;}
  void SetEvSelectionOn(Bool_t evselon=kTRUE){fOnOff[3]=evselon;}
  void SetFlowObsOn(Bool_t flowobson=kTRUE){fOnOff[4]=flowobson;}
  void SetUseSelectionBit(Bool_t selectionbiton=kTRUE){fUseSelectionBit=selectionbiton;}
  void SetSecondCentralityEstimator(AliRDHFCuts::ECentrality est){fEstimator = est;}
  void SetFillDistributionsForTrackEffChecks(Bool_t filldistrtrackeffcheckson=kFALSE){fFillDistrTrackEffChecks = filldistrtrackeffcheckson;}
  void SetAODMismatchProtection(Int_t opt=1) {fAODProtection=opt;}

  /// getters
  AliRDHFCuts* GetCutObject() const {return fCuts;}
  DecChannel GetDecayChannel()const {return fDecayChannel;}
  Bool_t GetTrackStatus() const {return fOnOff[0];}
  Bool_t GetPIDStatus() const {return fOnOff[1];}
  Bool_t GetCentralityStatus() const {return fOnOff[2];}
  Bool_t GetEvSelStatus() const {return fOnOff[3];}
  Bool_t GetFlowObsStatus() const {return fOnOff[4];}
  Bool_t GetUseSelectionBit() const {return fUseSelectionBit;}
  AliRDHFCuts::ECentrality GetSecondCentralityEstimator()const {return fEstimator;}
  Bool_t GetFillDistributionsForTrackEffChecks()const {return fFillDistrTrackEffChecks;}

 private:
  AliAnalysisTaskSEHFQA(const AliAnalysisTaskSEHFQA &source);
  AliAnalysisTaskSEHFQA operator=(const AliAnalysisTaskSEHFQA &source);
  void FillFlowObs(AliAODEvent *aod);

 TList* fOutputEntries;    //!<! list sent on output slot 1
 TList* fOutputPID;        //!<! list sent on output slot 2
 TList* fOutputTrack;      //!<! list sent on output slot 3
 TList* fOutputCounters;   //!<! list sent on output slot 5
 TList* fOutputCheckCentrality;   //!<! list sent on output slot 6
 TList* fOutputEvSelection; //!<! list sent on output slot 7
 TList* fOutputFlowObs;    //!<! list sent on output slot 8
 DecChannel fDecayChannel; //identify the decay channel
 AliRDHFCuts* fCuts;       // object containing cuts 
 AliFlowEvent *fFlowEvent; //!<! to handle the reusage of the flowEvent object
 AliFlowTrackCuts *fRFPcuts; //!<! reference flow particle cuts
 AliRDHFCuts::ECentrality fEstimator; //2nd estimator for centrality
 Bool_t fReadMC;           /// flag to read MC
 Bool_t fSimpleMode;       /// if true, don't do candidates (much faster in PbPb)
 Bool_t fUseSelectionBit;  /// flag to use or not the selection bit
 Bool_t fOnOff[5];         /// on-off the QA on tracks (0), PID (1), centrality (2), event selection -- default is {kTRUE,kTRUE,kTRUE,kTRUE}
 Bool_t fFillDistrTrackEffChecks;
 Int_t fAODProtection;     /// flag to activate protection against AOD-dAOD mismatch.
                           /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names

 //___ Declaration of histograms
 TH1F* fHisNentries;                         //!<!  Histo. of output slot #1 (fOutputEntries)
 TH2F* fHisNentriesSelBit;                   //!<!  Histo. of output slot #1 (fOutputEntries)
 TH1F* fHisTOFflags;                         //!<!  Histo. of output slot #2 (fOutputPID)
 TH1F* fHisTOFsig;                           //!<!  Histo. of output slot #2 (fOutputPID)
 TH1F* fHisTOFstartTimeMask;                 //!<!  Histo. of output slot #2 (fOutputPID)
 TH1F* fHisTOFstartTimeRes;                  //!<!  Histo. of output slot #2 (fOutputPID)
 TH1F* fHisTOFstartTimeDistrib;              //!<!  Histo. of output slot #2 (fOutputPID)
 TH1F* fHisTOFtime;                          //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTOFtimeKaonHyptime;               //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTOFtimeKaonHyptimeAC;             //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTOFsigmaKSigPid;                  //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTOFsigmaPionSigPid;               //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTOFsigmaProtonSigPid;             //!<!  Histo. of output slot #2 (fOutputPID)
 TH1F* fHisTOFsigPid3sigPion;                //!<!  Histo. of output slot #2 (fOutputPID)
 TH1F* fHisTOFsigPid3sigKaon;                //!<!  Histo. of output slot #2 (fOutputPID)
 TH1F* fHisTOFsigPid3sigProton;              //!<!  Histo. of output slot #2 (fOutputPID)
 TH1F* fHisTPCsig;                           //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTPCsigvsp;                        //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTPCsigvspAC;                      //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTPCsigmaK;                        //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTPCsigmaPion;                     //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTPCsigmaProton;                   //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTPCsigNvsPtAllTracks;             //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTPCsigNvsPhiAllTracks;            //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTPCsigNvsEtaAllTracks;            //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTPCsigNvsPtDaughters;             //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTPCsigNvsPhiDaughters;            //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTPCsigNvsEtaDaughters;            //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTOFsigmaMCKSigPid;                //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTOFsigmaMCPionSigPid;             //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTOFsigmaMCProtonSigPid;           //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTPCsigmaMCK;                      //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTPCsigmaMCPion;                   //!<!  Histo. of output slot #2 (fOutputPID)
 TH2F* fHisTPCsigmaMCProton;                 //!<!  Histo. of output slot #2 (fOutputPID)
 TH1F* fHisnClsITS;                          //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisnClsITSselTr;                     //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisnClsITSSA;                        //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisnClsITSSAspdAny;                  //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisnClsITSSAspdIn;                   //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisnClsITSSAspdOut;                  //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisnLayerITS;                        //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisnLayerITSselTr;                   //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisnLayerITSsa;                      //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisnClsSPD;                          //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisptGoodTr;                         //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisptGoodTrFromDaugh;                //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisptGoodTrFromDaugh_filt;           //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisdistrGoodTr;                      //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisdistrSelTr;                       //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisd0dau;                            //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisd0dau_filt;                       //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisd0dauphi;                         //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisd0dauphi_filt;                    //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisd0zdau;                           //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisd0zdau_filt;                      //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisd0zdauphi;                        //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisd0zdauphi_filt;                   //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisd0TracksSPDin;                    //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisd0TracksSPDany;                   //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisd0TracksFilterBit4;               //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisd0TracksTPCITSSPDany;             //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisPtDaughters;                      //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisPhiDaughters;                     //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisEtaDaughters;                     //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisEtavsPhiDaughters;                //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNTPCclsvsPtDaughters;             //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNTPCclsvsPhiDaughters;            //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNTPCclsvsEtaDaughters;            //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNTPCCrossedRowsvsPtDaughters;     //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNTPCCrossedRowsvsPhiDaughters;    //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNTPCCrossedRowsvsEtaDaughters;    //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisRatioCRowsOverFclsvsPtDaughters;  //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisRatioCRowsOverFclsvsPhiDaughters; //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisRatioCRowsOverFclsvsEtaDaughters; //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNITSclsvsPtDaughters;             //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNITSclsvsPhiDaughters;            //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNITSclsvsEtaDaughters;            //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisSPDclsDaughters;                  //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisPtAllTracks;                      //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisPhiAllTracks;                     //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisEtaAllTracks;                     //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisEtavsPhiAllTracks;                //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNTPCclsvsPtAllTracks;             //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNTPCclsvsPhiAllTracks;            //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNTPCclsvsEtaAllTracks;            //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNTPCCrossedRowsvsPtAllTracks;     //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNTPCCrossedRowsvsPhiAllTracks;    //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNTPCCrossedRowsvsEtaAllTracks;    //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisRatioCRowsOverFclsvsPtAllTracks;  //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisRatioCRowsOverFclsvsPhiAllTracks; //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisRatioCRowsOverFclsvsEtaAllTracks; //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNITSclsvsPtAllTracks;             //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNITSclsvsPhiAllTracks;            //!<!  Histo. of output slot #3 (fOutputTrack)
 TH2F* fHisNITSclsvsEtaAllTracks;            //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisSPDclsAllTracks;                  //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisdistrFakeTr;                      //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisd0f;                              //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisd0f_filt;                         //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisptFakeTr;                         //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisptFakeTrFromDaugh;                //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisptFakeTrFromDaughFilt;            //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisNtracklets;                       //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisNtracklets01;                     //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisNtracklets01AllEv;                //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisMult;                             //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisMultFBit4;                        //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisMultComb05;                       //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisMultComb08;                       //!<!  Histo. of output slot #3 (fOutputTrack)
 TH1F* fHisNtrackletsIn;                     //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH1F* fHisMultIn;                           //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH1F* fHisNtrackletsOut;                    //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH1F* fHisMultOut;                          //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH2F* fHisMultvsPercentile;                 //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH2F* fHisntrklvsPercentile;                //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH2F* fHisntrklvsPercentile01;              //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH2F* fHisntrklvsPercentile01AllEv;         //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH2F* fHisnTPCTracksvsPercentile;           //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH2F* fHisnTPCITSTracksvsPercentile;        //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH2F* fHisnTPCITS1SPDTracksvsPercentile;    //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH2F* fHisStdEstimSignalPercentile;         //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH2F* fHisStdEstimSignalNtrackletsIn;       //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH1F* fHisStdEstimSignal;                   //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH2F* fHisStdPercentileSecondPercentile;    //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH2F* fHisStdSignalSecondSignal;            //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH2F* fHisStdPercentileOldFrwPercentile;    //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH1F* fHisStdPercentileOldFrwPercentileDev; //!<!  Histo. of output slot #6 (fOutputCheckCentrality)
 TH1F* fHisxvtx;                             //!<!  Histo. of output slot #7 (fOutputEvSelection)
 TH1F* fHisyvtx;                             //!<!  Histo. of output slot #7 (fOutputEvSelection)
 TH1F* fHiszvtx;                             //!<!  Histo. of output slot #7 (fOutputEvSelection)
 TH1F* fHisxvtxSelEv;                        //!<!  Histo. of output slot #7 (fOutputEvSelection)
 TH1F* fHisyvtxSelEv;                        //!<!  Histo. of output slot #7 (fOutputEvSelection)
 TH1F* fHiszvtxSelEv;                        //!<!  Histo. of output slot #7 (fOutputEvSelection)
 TH1F* fHisWhichVert;                        //!<!  Histo. of output slot #7 (fOutputEvSelection)
 TH1F* fHisWhichVertSelEv;                   //!<!  Histo. of output slot #7 (fOutputEvSelection)
 TH2F* fHisTrigCent;                         //!<!  Histo. of output slot #7 (fOutputEvSelection)
 TH2F* fHisTrigMul;                          //!<!  Histo. of output slot #7 (fOutputEvSelection)
 TH2F* fHisTrigCentSel;                      //!<!  Histo. of output slot #7 (fOutputEvSelection)
 TH2F* fHisTrigMulSel;                       //!<!  Histo. of output slot #7 (fOutputEvSelection)
 TH1F* fHisWhyEvRejected;                    //!<!  Histo. of output slot #7 (fOutputEvSelection)
 TH2F* fHisnClsITSvsNtrackletsSel;           //!<!  Histo. of output slot #7 (fOutputEvSelection)
 TH2F* fHiszvtxvsSPDzvtx;                    //!<!  Histo. of output slot #7 (fOutputEvSelection)
 TH2F* fHiszvtxvsSPDzvtxSel;                 //!<!  Histo. of output slot #7 (fOutputEvSelection)
 TH2F* fHisFEvents;                          //!<!  Histo. of output slot #8 (fOutputFlowObs)
 TH3F* fHisTPCVZE_AngleQ;                    //!<!  Histo. of output slot #8 (fOutputFlowObs)
 TH2F* fHisCentVsMultRPS;                    //!<!  Histo. of output slot #8 (fOutputFlowObs)
 TH2F* fHisAngleQ[3];                        //!<!  Histo. of output slot #8 (fOutputFlowObs)
 TH3F* fHisPhiEta[3];                        //!<!  Histo. of output slot #8 (fOutputFlowObs)
 TProfile2D *fHisQ[3];                       //!<!  Histo. of output slot #8 (fOutputFlowObs)

 /// \cond CLASSIMP
 ClassDef(AliAnalysisTaskSEHFQA,17); ///AnalysisTaskSE for the quality assurance of HF in hadrons
 /// \endcond
};

#endif

