/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIANALYSISTASKHAHFECORREL_H
#define ALIANALYSISTASKHAHFECORREL_H 

class THnSparse;
class TH2F;
class TProfile;
class TLorentzVector;

class AliEMCALTrack;
class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliAODMCParticle;
class AliAODMCHeader;
class AliGenEventHeader;
class AliVEvent;
class AliEMCALGeometry;
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliESDtrack;
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliCFManager;
class AliPIDResponse;
class AliMultSelection;
class AliEventPoolManager;
class AliAODv0KineCuts;
class AliESDv0KineCuts;
class AliVertexingHFUtils;
//class AliExternalTrackParam;
#include "AliAODv0KineCuts.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "TDatabasePDG.h"
#include <vector>


class AliAnalysisTaskHaHFECorrel : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskHaHFECorrel();
    AliAnalysisTaskHaHFECorrel(const char *name);
    ~AliAnalysisTaskHaHFECorrel();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    
    //******************** ANALYSIS
    AliVTrack* FindLPAndHFE(TObjArray* RedTracksHFE, const AliVVertex *pVtx, Int_t nMother, Int_t listMother[], Double_t mult, Bool_t &EvContTP, Bool_t &EvContNTP, Double_t EventWeight);
    void FindPhotonicPartner(Int_t iTracks, AliVTrack* track,  const AliVVertex *pVtx, Int_t nMother, Int_t listMother[], Int_t &LSPartner, Int_t &ULSPartner, Int_t *LSPartnerID, Int_t *ULSPartnerID,  Float_t *LSPartnerWeight, Float_t *ULSPartnerWeight, Bool_t &trueULSPartner, Bool_t &isPhotonic, Float_t &MCPartnerPt, Float_t &RecPartnerPt, Double_t EventWeight, Double_t mult);
    void CheckPhotonicPartner(AliVTrack* Vtrack, Bool_t Tagged, Float_t& MCPartnerPt, Float_t RecPartnerPt, Double_t EventWeight);
    void CorrelateElectron(TObjArray* RedTracksHFE);

    void CorrelateLP(AliVTrack* LPtrack,  const AliVVertex* pVtx, Int_t nMother, Int_t listMother[], TObjArray* RedTracksHFE, Double_t EventWeight);
    void CorrelateLPMixedEvent(AliVTrack* LPtrack, Float_t mult, Float_t zVtx, Float_t maxPt, Bool_t EvContTP, Bool_t EvContNTP, Double_t EventWeight);
    
    void CorrelateHadron(TObjArray* RedTracksHFE,  const AliVVertex* pVtx, Int_t nMother, Int_t listMother[], Float_t mult, Float_t maxPt, Double_t EventWeight);
    void CorrelateHadronMixedEvent(Float_t mult, const AliVVertex* zVtx, Float_t maxPt, Int_t nMother, Int_t listMother[], Bool_t EvContTP, Bool_t EvContNTP, Double_t EventWeight);

    void CorrelateWithHadrons(AliVTrack* TriggerTrack, const AliVVertex* pVtx, Int_t nMother, Int_t listMother[], Bool_t FillHadron, Bool_t FillLP,Bool_t** NonElecIsTrigger, Double_t *NonElecIsTriggerPt, Double_t *NonElecIsTriggerWeight, Int_t NumElectronsInEvent, Double_t EventWeight); 

    void MCTruthCorrelation(TObjArray* MCRedTracks, Bool_t AfterEventCuts, Int_t RecLPLabel, Float_t pVtxZ, Float_t mult, Int_t &LPinAcceptance, Int_t &LP, Double_t EventWeight);

    //********************MC
    void  MCEfficiencyCorrections(const AliVVertex * RecVertex, Double_t EventWeight);
    Int_t HFEisCharmOrBeauty(Int_t ElectronIndex);

    //*********************ANALYSIS Helper
    Bool_t ChargedHadronTrackCuts(const AliVVertex *pVtx,AliVTrack *Htrack, Int_t nMother, Int_t listMother[], Double_t EventWeight, Bool_t fillHists=kFALSE);
    Bool_t ChargedHadronPIDCuts(AliVTrack *Htrack, Double_t EventWeight);;

    Bool_t AssoHadronPIDCuts(AliVTrack *Htrack, Double_t EventWeight);

    Bool_t InclElecTrackCuts(const AliVVertex *pVtx,AliVTrack *ietrack, Int_t nMother, Int_t listMother[], Double_t EventWeight, Bool_t fillHists=kFALSE);
    Bool_t InclElecPIDCuts(AliVTrack *track,  Double_t EventWeight, Bool_t fillHists=kFALSE);

    Bool_t PhotElecPIDCuts(AliVTrack *track, Double_t EventWeight);
    Bool_t PhotElecTrackCuts(const AliVVertex *pVtx,AliVTrack *aetrack, Int_t nMother, Int_t listMother[], Double_t EventWeight);
    void   PhotULSLSElectronAcceptance(const AliVVertex *pVtx, Float_t mult,  Int_t nMother, Int_t listMother[], Double_t EventWeight);
    
    void EvaluateTaggingEfficiency(AliVTrack * track, Int_t LSPartner, Int_t ULSPartner, Bool_t trueULSPartner, Double_t EventWeight, Double_t mult, Double_t recEffE); 
    Bool_t CloneAndReduceTrackList(TObjArray* RedTracks, AliVTrack* track, Int_t LSPartner, Int_t ULSPartner, Int_t *LSPartnerID, Int_t *ULSPartnerID, Float_t *LSPartnerWeight, Float_t *ULSPartnerWeight, Bool_t trueULSPartner, Float_t MCPartnerPt, Float_t RecPartnerPt, Bool_t isPhotonic, Bool_t isHadron);

    void BinLogX(TAxis *axis);
    void SetPDGAxis(TAxis *axis, std::vector<TString> PDGLabel);
    void SetTriggerAxis(TAxis *axis);
    void CheckHadronIsTrigger(Double_t ptE, Bool_t *HadronIsTrigger);
    void CheckElectronIsTrigger(Double_t ptH, Bool_t *ElectronIsTrigger) ;
    Bool_t PassEventBias( const AliVVertex *pVtx, Int_t nMother, Int_t *listMother, Double_t EventWeight);    


    //**************  SETTINGS
    void SetMC (Bool_t IsMC) { fIsMC=IsMC;};
    void SetAODanalysis(Bool_t IsAOD) {fIsAOD = IsAOD;};
    void SetTender (Bool_t UseTender) {fUseTender = UseTender;};
    void SetPeriod (Double_t period) {fWhichPeriod = period;};
    void SetEpos (Bool_t UseEpos) {fUseEPOS = UseEpos;};
    
    void SetOnlyEfficiency() {
      fTRDQA = kFALSE;
      fCorrHadron = kFALSE;
      fCorrLParticle = kFALSE;
      fMixedEvent = kFALSE;
      fMCTrueCorrelation = kFALSE;
    }

    
    void SetEleVarOpt(Int_t VarOption);
    void SetHadVarOpt(Int_t VarOption);
    void SetPhotVarOpt(Int_t PhotVarOpt);
    void SetVtxVarOpt(Int_t VtxVarOpt) {fVarZVTXCut=VtxVarOpt;};
    
    void SetTRDQA(Bool_t TRDQA) {fTRDQA=TRDQA;};
    void SetPtMinEvent(Double_t PtMin) {fMinPtEvent=PtMin;};
    void SetPtMaxEvent(Double_t PtMax) {fMaxPtEvent=PtMax;};

    void SetMinNTr(Double_t MinNTr) {fMinNTr=MinNTr;};
    void SetMaxNTr(Double_t MaxNTr) {fMaxNTr=MaxNTr;};


    void SetEtaMax(Double_t EtaMax) {
      fMaxElectronEta = TMath::Abs(EtaMax);
      fMinElectronEta = -TMath::Abs(EtaMax);
    };

    void SetTPCnCut(Int_t TPCnCut) {fTPCnCut = TPCnCut;};
    void SetTPCnCutdEdx(Int_t TPCnCutdEdx) {fTPCndEdxCut = TPCnCutdEdx;};
    void SetITSnCut (Int_t ITSnCut) {fITSnCut = ITSnCut;};
    void SetITSSharedClusterCut (Float_t ITSSharedCluster) {fITSSharedClusterCut=ITSSharedCluster;};
   

    void SetPhotElecPtCut (Double_t AssPtCut) {fPhotElecPtCut = AssPtCut;};
    void SetPhotElecTPCnCut (Int_t AssTPCnCut) {fPhotElecTPCnCut = AssTPCnCut;};
    void SetPhotElecITSrefitCut(Bool_t AssITSrefitCut) {fPhotElecITSrefitCut = AssITSrefitCut;};  
    void SetPhotCorrCase(Int_t PhotCorrCase) {fPhotCorrCase = PhotCorrCase;};

    void SetHTPCnCut(Int_t HTPCnCut) {fHTPCnCut = HTPCnCut;}
    void SetHITSrefitCut(Bool_t HITSrefitCut) {fHITSrefitCut = HITSrefitCut;};   
    void SetHTPCrefitCut(Bool_t HTPCrefitCut) {fHTPCrefitCut = HTPCrefitCut;};

   
    void SetUseTRD(Bool_t UseTRD) {fUseTRD = UseTRD;}
    void SetUseITSsa(Bool_t UseITSsa) {fUseITSsa = UseITSsa;}
    void SetSigmaITScut(Double_t SigmaITScut) {fSigmaITScut = SigmaITScut;};
    void SetSigmaTOFcut(Double_t SigmaTOFcut) {fSigmaTOFcut = SigmaTOFcut;};
    void SetSigmaTPCcut(Double_t SigmaTPCcut) {fSigmaTPCcutLow = SigmaTPCcut;};

  
    void SetRecEff(Bool_t RecEff) { 
      fRecEff=RecEff;
    }
    void SetTagEff(Bool_t TagEff) {
      fTagEff=TagEff;
    }
    void SetOneTimeCheck(Bool_t OneTimeCheck) {
      fOneTimeCheck = OneTimeCheck;
    }
    void SetHadronCorrelation(Bool_t CorrHadron) {
      fCorrHadron = CorrHadron;
    };
    void SetLPCorrelation(Bool_t CorrLP) {
      fCorrLParticle = CorrLP;
    };
    void SetMCTruthCorrelation(Bool_t MCTruthCorr) {
      fMCTrueCorrelation = MCTruthCorr;
    };
    void SetUseEventWeights(Bool_t UseEventWeights) {
      fUseEventWeights = UseEventWeights;
    };

    void SetOpeningAngleCut(Bool_t OpeningAngleCut) {fOpeningAngleCut=OpeningAngleCut;};
    void SetInvmassCut(Double_t InvmassCut) {fInvmassCut=InvmassCut;};

    void SetPi0WeightToData(TH1F &  WPion) {fCorrectPiontoData = WPion; fCorrectPiontoData.SetName("fCorrectPiontoData");};
    void SetEtaWeightToData(TH1F &  WEta)  {fCorrectEtatoData  = WEta; fCorrectEtatoData.SetName("fCorrectEtatoData");};
    void SetBGWeight(TH2F & BGWeight) {fBgWeight = BGWeight; fBgWeight.SetName("fBgWeight");};
    void SetHadRecEff(TH3F & HadRecEff) {fHadRecEff = HadRecEff; fHadRecEff.SetName("fHadRecEff");};
    void SetEleRecEff(TH3F & EleRecEff) {fEleRecEff = EleRecEff; fEleRecEff.SetName("fEleRecEff");};
    //  void SetSPDnTrAvg(TProfile & SPDnTrAvg) {fSPDnTrAvg = SPDnTrAvg; fSPDnTrAvg.SetName("fSPDnTrAvg");}
    void SetSPDConfigHist(TH1I & SPDConfigHist) {fSPDConfigHist = SPDConfigHist; fSPDConfigHist.SetName("SPDConfigHist");
      /* for (Int_t i=1; i<300; i++) { */
      /* 	printf("%i, %s, %10.2f, %s, %10.2f", i, fSPDConfigHist.GetXaxis()->GetBinLabel(i), fSPDConfigHist.GetBinContent(i),  SPDConfigHist.GetXaxis()->GetBinLabel(i) ,  SPDConfigHist.GetBinContent(i)); */
      /* } */

    };
    void SetSPDConfigProfiles(TH3F & SPDConfigProfiles) {fSPDConfigProfiles = SPDConfigProfiles; fSPDConfigProfiles.SetName("fSPDConfigProfiles");}
    void SetNonTagCorr(TH1F & NonTagCorr) {fNonTagCorr = NonTagCorr; fNonTagCorr.SetName("fNonTagCorr");}
    void SetTriggerWeight(TH3F & TriggerWeight){fTriggerWeight = TriggerWeight; fTriggerWeight.SetName("fTriggerWeight");}
    void SetVtxWeight(TH2F & VtxWeight) {fVtxWeight = VtxWeight; fVtxWeight.SetName("fVtxWeight");};

    Bool_t   ESDkTrkGlobalNoDCA(AliVTrack* Vtrack);

  


 private:
    
    Bool_t                IsPhotonicElectron(Int_t Label1) const;
    Bool_t                HaveSameMother(Int_t Label1, Int_t Label2) const;
    Double_t              GetDeltaPhi(Double_t phiA,Double_t phiB) const;
    Double_t              GetDeltaEta(Double_t etaA,Double_t etaB) const;
    Double_t              Eta2y(Double_t pt, Double_t m, Double_t eta) const;
    Double_t              GetHadronRecEff(Int_t run, Double_t pt, Double_t phi, Double_t eta, Double_t zVtx);
    Double_t              GetElectronRecEff(Int_t run, Double_t pt, Double_t phi, Double_t eta, Double_t zVtx);
    Double_t              GetTriggerWeight(Int_t run, Double_t minV0, Double_t nTrAcc);
    Double_t              GetVtxWeight(Int_t run, Double_t nTrAcc);
    Double_t              GetNonTagCorr(Double_t ptTrack, Double_t ptAsso);

    Double_t              Sphericity(const TObjArray* tracks, Double_t MaxEta, Double_t MinPt);
    Bool_t                Thrust(const TObjArray* tracks, Double_t t[2], Double_t MaxEta, Double_t MinPt);

    Int_t                 CheckParticleOrigin(Int_t Label);

    Int_t                 fRunNumber;               //
    Bool_t                fUseTender;               // Use tender
    Int_t                 fWhichPeriod;             // period
    Bool_t                fUseEPOS;
    Bool_t                fUseKFforPhotonicPartner; //default ist DCA

    Float_t               fMaxPtEvent;              //
    Float_t               fMinPtEvent;              //
    Int_t                 fMaxNTr;                  //
    Int_t                 fMinNTr;                  //
    Int_t                 fVarZVTXCut;             //

    
    Double_t              fMaxElectronEta;          //
    Double_t              fMinElectronEta;          //
    Double_t              fMaxHadronEta;            //
    Double_t              fMinHadronEta;            //

    // HFECuts
    Int_t                 fVarEleOpt;               //
    Bool_t                fElectronkAny;            // True: kAny, False: kBoth
    Bool_t                fElectronkFirst;          // True: kFirst, False: kBoth
    Int_t                 fTPCnCut;                 // TPC number of clusters for tagged electron
    Int_t                 fTPCndEdxCut;             //
    Int_t                 fITSnCut;                 // ITs number of clusters for tagged electrons 
    Float_t               fITSSharedClusterCut;     //
    Double_t              fEleDCAr;                 //
    Double_t              fEleDCAz;                 //

    Bool_t                fUseTRD;                  //
    Bool_t                fUseITSsa;                // Use ITSsa tracks
    Double_t              fSigmaITScut;             // ITS nSigma cut
    Double_t              fSigmaTOFcut;             // TOF nSigma cut
    Double_t              fSigmaTPCcutLow;          // lower TPC nSigma cut
    Double_t              fSigmaTPCcutHigh;         //

    // Photonic  Electrons
    Int_t                 fVarPhotOpt;                   //
    Double_t              fPhotElecPtCut;                // pt cut for associated electron
    Double_t              fPhotElecSigmaTPCcut;          //
    Int_t                 fPhotElecTPCnCut;              // TPC number of clusters for associated electron
    Bool_t                fPhotElecITSrefitCut;          // ITS refit for associated electron
    Int_t                 fPhotCorrCase;           //

    // Associate Hadron (non Electron)
    Double_t              fAssNonEleTPCcut;         //  

    // Hadron Cut
    Int_t                 fVarHadOpt;               //
    Int_t                 fHTPCnCut;                // TPC number of clusters for trigger hadron
    Bool_t                fHITSrefitCut;            // ITS refit for trigger hadron
    Bool_t                fHTPCrefitCut;            // TPC refit for trigger hadron
    Double_t              fHadDCAr;                 //
    Double_t              fHadDCAz;                 //
    Bool_t                fHadkAny;                 //
    Bool_t                fHadTOFmatch;             // matching to TOF bunch crossing ID to suppress pileup
    
    Double_t              fOpeningAngleCut;         // openingAngle cut for non-HFE selection
    Double_t              fInvmassCut;              // invariant mass cut  for non-HFE selection
    Double_t              fChi2Cut;                 //! used?? Chi2 cut  for non-HFE selection
    Double_t              fDCAcut;                  //! used?? DCA cut  for non-HFE selection
  
    // ******* Switch for analysis modes
    Bool_t                fTRDQA;                   // TRDQA
    Bool_t                fMCTrueCorrelation;       //
    Bool_t                fUseEventWeights;         //
    Bool_t                fCorrHadron;              // Choose Hadron-HFE Correl
    Bool_t                fCorrLParticle;           // Choose LP-HFE Correl
    Bool_t                fMixedEvent;              // Fill Mixed Event for the cases chosen above
    Bool_t                fPionEtaProduction;       //
    Bool_t                fRecEff;                  //
    Bool_t                fTagEff;                  //
    Bool_t                fHadCont;                 //
    Bool_t                fOneTimeCheck;            //
    Bool_t                fLParticle;               // Is LP found?

    AliESDEvent           *fESD;                    //! ESD object
    AliESDtrackCuts       *fesdTrackCuts;           //!
    AliAODEvent           *fAOD;                    //! AOD object
    AliVEvent             *fVevent;                 //! VEvent
    AliPIDResponse        *fpidResponse;            //! PID response
    AliMultSelection      *fMultSelection;          //! MulSelection
    AliCentrality         *fCentrality;             //! Centrality

    AliEventPoolManager   *fPoolMgr;                //! event pool manager
    TH3F                  *fPoolIsFilled;           //! check if pool is filled
    
    AliMCEvent            *fMC;                     //! MC object
    AliStack              *fStack;                  //! stack
    AliAODMCParticle      *fMCparticle;             //! MC particle
    TClonesArray          *fMCarray;                //! MC array
    AliAODMCHeader        *fMCheader;               //! MC header
    TDatabasePDG          *PdgTable;                //!
    std::map<Int_t, Int_t>     PDGMap;                   //!

    TClonesArray          *fTracks_tender;          //Tender tracks
    TClonesArray          *fCaloClusters_tender;    //Tender clusters
      
    AliEventCuts          fEventCuts;               //! Test AliEventCuts
    AliHFEcuts            *fCuts;                   //! Cut Collection
    
     
    Bool_t                fIsMC;                    // flag for MC analysis
    Bool_t                fIsAOD;                   // flag for AOD analysis
    AliCFManager          *fCFM;                    //! Correction Framework Manager
    AliHFEpid             *fPID;                    //! PID
    AliHFEpidQAmanager    *fPIDqa;                  //! PID QA manager
        
    TList                 *fOutputList;             //! output list
    TList                 *fOutputListMain;         //!
    TList                 *fOutputListLP;           //!
    TList                 *fOutputListHadron;       //!
    TList                 *fOutputListQA;           //!
    TH1F                  *fNoEvents;               //! no of events for different cuts
    TH2F                  *fNoEventsNTr;               //! no of events for different cuts
    TH2F                  *fMCNoEvents;             //! no of events for different cuts
    TH2F                  *fHFENoEvents;            //! no of events for different cuts
    TH3F                  *fDiffractiveType;        //!
    TH2F                  *fV0ACTrueInel;           //!
    TH2F                  *fV0TrueMinInel;          //!
    TH3F                  *fV0TrueMinInelNTr;          //!
    TH2F                  *fV0ACTriggered;          //!
    TH2F                  *fV0MinTriggered;         //!
    TH3F                  *fV0MinTriggeredNTr;         //!
    TH3F                  fTriggerWeight;
    TH2F                  *fVtxEtaNTr;              //!
    TH2F                  *fVtxBeforeNTrAcc;        //!
    TH2F                  *fVtxAfterNTrAcc;         //!
    TH1F                  *fVtxRecBeforeNTr;        //!
    TH2F                  *fVtxRecAfterNTr;         //!
    TH2F                  fVtxWeight;
    TH2F                  *fTrkpt;                  //! track pt for different cuts
    TH2F                  *fEtaVtxZ;                //! Eta vs Vtx z (check for ITS acceptance problem)

    TH2F                  *fSPDVtxRes;              //!
    TH2F                  *fDiffSPDPrimVtx;         //!
    TH2F                  *fSPDnTrAcc;              //!
    TH2F                  *fSPDnTrCorrMax;          //!
    TH2F                  *fSPDnTrGen;              //!
    TH2F                  *fDiffSPDMCVtx;           //!
    THnSparseF            *fnTrAccMaxGen;           //!
    THnSparseF            *fnTrAccGen;           //!
    TH2F                  *fnTrAccGenTrueInel;   //!
    //TH2F                  *fnTrAccGenTrueInelTrig;   //!
    //TH2F                  *fnTrAccGenTrueInelVtxQA;   //!
    //TH2F                  *fnTrAccGenTrueInelVtxEx;   //!
    THnSparseF            *fnTrAccMinGen;           //!
    THnSparseF            *fnTrAccMeanGen;          //!
    THnSparseF            *fnTrAccMax;              //!
    THnSparseF            *fnTrAccMin;              //!
    THnSparseF            *fnTrAccMean;             //!
    TH3F                  *fMCThrustTagged; //!
    TH3F                  *fMCSpherTagged; //! 
    TH3F                  *fRecLPTagged; //!
    TH3F                  *fMultCorrTagged; //!
    TH3F                  *fNHadTagged; //!
    TH3F                  *fNHadTaggedA; //!
    TH3F                  *fNHadTaggedB; //!
    TH3F                  *fNHadTaggedC; //!
    TH3F                 *fMeanPtTagged; //!
    TH3F                 *fMeanPtTaggedA; //!
    TH3F                 *fMeanPtTaggedB; //!
    TH3F                 *fMeanPtTaggedC; //!
    TH3F                  *fMCThrustNTagged; //!
    TH3F                   *fMCSpherNTagged; //! 
    TH3F                  *fRecLPNTagged; //!
    TH3F                  *fMultCorrNTagged; //!
    TH3F                  *fNHadNTagged; //!
    TH3F                  *fNHadNTaggedA; //!
    TH3F                  *fNHadNTaggedB; //!
    TH3F                  *fNHadNTaggedC; //!
    TH3F                 *fMeanPtNTagged; //!
    TH3F                 *fMeanPtNTaggedA; //!
    TH3F                 *fMeanPtNTaggedB; //!
    TH3F                 *fMeanPtNTaggedC; //!

    TH3F *fPt2Tagged; //!
    TH3F *fPt2NTagged; //!
   TH2F                  *fMothMCThrustTagged; //!
    TH2F                  *fMothMCSpherTagged; //! 
    TH2F                  *fMothRecLPTagged; //!
    TH2F                  *fMothMultCorrTagged; //!
    TH2F                  *fMothNHadTagged; //!
    TH2F                 *fMothMeanPtTagged; //!
    TH2F                  *fMothMCThrustNTagged; //!
    TH2F                   *fMothMCSpherNTagged; //! 
    TH2F                  *fMothRecLPNTagged; //!
    TH2F                  *fMothMultCorrNTagged; //!
    TH2F                  *fMothNHadNTagged; //!
    TH2F                 *fMothMeanPtNTagged; //!

    TH2F                  *fMCThrustTaggedH; //!
    TH2F                  *fMCSpherTaggedH; //! 
    TH2F                  *fRecLPTaggedH; //!
    TH2F                  *fMultCorrTaggedH; //!
    TH2F                  *fNHadTaggedH; //!
    TH2F                 *fMeanPtTaggedH; //!
    TH2F                  *fMCThrustNTaggedH; //!
    TH2F                   *fMCSpherNTaggedH; //! 
    TH2F                  *fRecLPNTaggedH; //!
    TH2F                  *fMultCorrNTaggedH; //!
    TH2F                  *fNHadNTaggedH; //!
    TH2F                 *fMeanPtNTaggedH; //!

   TH2F                  *fMothMCThrustTaggedH; //!
    TH2F                  *fMothMCSpherTaggedH; //! 
    TH2F                  *fMothRecLPTaggedH; //!
    TH2F                  *fMothMultCorrTaggedH; //!
    TH2F                  *fMothNHadTaggedH; //!
    TH2F                 *fMothMeanPtTaggedH; //!
    TH2F                  *fMothMCThrustNTaggedH; //!
    TH2F                   *fMothMCSpherNTaggedH; //! 
    TH2F                  *fMothRecLPNTaggedH; //!
    TH2F                  *fMothMultCorrNTaggedH; //!
    TH2F                  *fMothNHadNTaggedH; //!
    TH2F                 *fMothMeanPtNTaggedH; //!


    
    THnSparse             *fMultiplicity;	    //! multiplicity distribution
    TH3F                  *fSPDMultiplicity;        //!
    Int_t                 *fRunList;                //!

    TH2F                  *fElectronTrackCuts;      //!
    TH2F                  *fElectronTrackTPCChi2;   //!
    TH2F                  *fElectronTrackTPCCrossedRows; //!
    TH2F                  *fElectronTrackTPCNcls;   //! 
    TH2F                  *fElectronTrackTPCNclsdEdx; //! 
    TH2F                  *fElectronTrackTPCFrac;   //! 
    TH2F                  *fElectronTrackITSNcls;   //!
    TH2F                  *fElectronTrackITSChi2;   //!
    TH3F                  *fElectronTrackITSLayer;  //!
    TH3F                  *fElectronTrackRefit;     //!
    TH3F                  *fElectronTrackDCA;       //! 
    THnSparseF*           fElectronTrackITSCuts;    //!
    THnSparseF*           fPhotTrackITSCuts;        //!
    TH2F                  *fHadronTrackCuts;        //!
    TH2F                  *fHadronTrackTPCNcls;     //! 
    TH3F                  *fHadronTrackRefit;       //!
    TH3F                  *fHadronTrackDCA;         //! 
    TH3F                  *fHadronTrackDCA_woITSAny;//!
    TH3F                  *fHadronTrackDCA_wITSAny; //!


   

    TH2F                  *fHistITSnSig;            //! ITS sigma vs p
    TH2F                  *fHistTOFnSig;            //! TOF sigma vs p
    TH2F                  *fHistTPCnSig;            //! TPC sigma vs p
    TH2F                  *fHistTPCnSigITScut;      //! TPC sigma vs p (ITS cut)
    TH2F                  *fHistTPCnSigTOFcut;      //! TPC sigma vs p (TOF cut)
    TH2F                  *fHistTPCnSigITSTOFcut;   //! TPC sigma vs p (ITS+TOF cuts)
    TH2F                  *fHistITSnSigTOFTPCcut;   //! ITS sigma vs p (TPC+TOF cuts)

    THnSparse             *fCheckNHadronScaling;    //!
    THnSparse             *fCheckNPhotHadScaling;   //!
    TH3F                  *fCheckTaggedEvent;       //!

    TH2F                  *fHadContPvsPt;           //!
    TH3F                  *fHadContEtaPhiPt;        //!
    TH3F                  *fHadContTPCEtaPhiPt;     //!
    THnSparse             *fHadContPPhiEtaTPC;      //!
    THnSparse             *fHadContamination;       //! HadronicContaminationTOF
    THnSparse             *fHadContaminationPt;     //! HadronicContaminationTOF
    THnSparse             *fHadContMC;              //!
    THnSparse             *fHadContMCPt;            //!
    

  
    TH3F                  *fInclElecPtEta;          //! inclusive electron p
    TH3F                  *fInclElecPtEtaWRecEff;        //! inclusive electron p
    TH1F                  *fInclElecP;              //! inclusive electron p
    TH2F                  *fULSElecPt;              //! ULS electron pt (after IM cut)
    TH2F                  *fULSElecPtWRecEff;            //! ULS electron pt (after IM cut)
    TH2F                  *fLSElecPt;               //! LS electron pt (after IM cut)
    TH2F                  *fLSElecPtWRecEff;             //! LS electron pt (after IM cut)
    TH2F                  *fInvmassLS;              //! Inv mass of LS (e,e)
    TH2F                  *fInvmassULS;             //! Inv mass of ULS (e,e)
    TH3F                  *fInvmassMCTrue;          //! Inv mass of ULS (e,e)
    TH3F                  *fRecMCInvMass;           //!
    TH1F                  *fPhotMixULS;             //!
    TH2F                  *fPhotMixLS;              //!
    TH2F                  *fPhotPt1PtMTag;          //!
    TH2F                  *fPhotPt1PtMNTag;         //!
    THnSparseF            *fPhotPt1Pt2;             //!
    TH2F                  *fPhotPt1Pt2Only;         //!
    THnSparseF            *fPhotPt1Pt2Corr;         //!
    THnSparseF            *fPhotPt1Pt2MC;           //!
    THnSparseF            *fPhotPt1RecPt2;         //!
    THnSparseF            *fPhotPt1RecPt2Corr;         //!
    THnSparseF            *fPhotPt1RecPt2Rec;         //!
    THnSparseF            *fPhotPt1RecPt2RecCorr;         //!
    THnSparseF            *fPhotPt1Pt2Rec;         //!
    THnSparseF            *fPhotPt1Pt2RecCorr;         //!
    TH2F                  *fPhotPt2MCRec;           //!
    THnSparseF            *fPhotPt1E;               //!
    THnSparseF            *fPhotPt1Pt2E;            //!
    THnSparseF            *fPhotPt1ECorr;           //!
    THnSparseF            *fPhotPt1Pt2ECorr;         //!
    THnSparseF            *fPhotPt1Mass;            //!
    THnSparseF            *fPhotPt1Mom;             //!
    TH2F                  *fOpeningAngleLS;         //! opening angle for LS pairs
    TH2F                  *fOpeningAngleULS;        //! opening angle for ULS pairs
    TH2F                  *fCheckLSULS;             //! check no of LS/ULS partner per electron
    TH3F                  *fTagEtaPt1Pt2;            //!
    TH3F                  *fTagEtaPhiPt;            //!
    TH3F                  *fTagEtaZvtxPt;           //!
    TH3F                  *fTagEtaPhiPtwW;          //!
    TH3F                  *fTagEtaZvtxPtwW;         //!
    TH3F                  *fNonTagEtaPt1Pt2;        //!
    TH3F                  *fNonTagEtaPhiPt;         //!
    TH3F                  *fNonTagEtaZvtxPt;        //!
    TH3F                  *fNonTagEtaPhiPtwW;       //!
    TH3F                  *fNonTagEtaZvtxPtwW;      //!



    THnSparse             *fTagMotherPt;            //!
    TH2F                  *fTagEffInclMult;         //!
    TH2F                  *fTagEffULSMult;          //!
    TH3F                  *fTagEffInclBGMult;       //!
    TH3F                  *fTagEffULSBGMult;        //!
    TH2F                  *fTagTruePairsMult;       //!
    TH2F                  *fTagEffInclMultWoW;         //!
    TH2F                  *fTagEffULSMultWoW;          //!
    TH3F                  *fTagEffInclBGMultWoW;       //!
    TH3F                  *fTagEffULSBGMultWoW;        //!
    TH2F                  *fTagTruePairsMultWoW;       //!
    TH2F                  *fTagEffInclMultWoWS;         //!
    TH2F                  *fTagEffULSMultWoWS;          //!
    TH3F                  *fTagEffInclBGMultWoWS;       //!
    TH3F                  *fTagEffULSBGMultWoWS;        //!
    TH2F                  *fTagTruePairsMultWoWS;       //!
    THnSparse             *fTagEffIncl;             //! 
    THnSparse             *fTagEffLS;               //!
    THnSparse             *fTagEffULS;              //!
    THnSparse             *fTagTruePairs;           //!
    THnSparse             *fTagEffInclWoWeight;     //! 
    THnSparse             *fTagEffLSWoWeight;       //!
    THnSparse             *fTagEffULSWoWeight;      //!
    THnSparse             *fTagTruePairsWoWeight;   //!

    TH1F                  fCorrectPiontoData;   
    Double_t              GetPionWeight(Double_t pt);
    TH1F                  fCorrectEtatoData;       
    Double_t              GetEtaWeight(Double_t pt);
    TH2F                  fBgWeight;       
    Double_t              GetBackgroundWeight(Int_t PDGMother, Double_t pt);
   

    TH3F                  fHadRecEff;
    TH3F                  fEleRecEff;
    Int_t                 fSPDConfig;
    TH1I                  fSPDConfigHist;
    TH3F                  fSPDConfigProfiles;
    TProfile*             fSPDnTrAvg;               //!
    TH1F                  fNonTagCorr;

    Int_t                 fAssPtHad_Nbins;
    TArrayF               fAssPtHad_Xmin;
    TArrayF               fAssPtHad_Xmax;

    Int_t                 fAssPtElec_Nbins;
    TArrayF               fAssPtElec_Xmin;
    TArrayF               fAssPtElec_Xmax;



    // HFE HFE
    TH1F                  *fElecTrigger;            //! trigger electron vs pt
    TH2F                  *fInclElecPhi;            //! electron (trigger): phi vs pt
    TH2F                  *fInclElecEta;            //! electron (trigger): phi vs pt
    TH2F                  *fInclElecPhiEta;            //! electron (trigger): phi vs pt
    TH2F                  *fULSElecPhi;             //! phi vs pt for electrons from ULS pairs
    TH2F                  *fLSElecPhi;              //! phi vs pt for electrons from LS pairs
    TH2F                  *fElecDphi;               //! inlcusive electron: dPhi vs pt of triggered electron
    TH2F                  *fULSElecDphi;            //! electron from ULS pairs: dPhi vs pt of triggered electron
    TH2F                  *fLSElecDphi;             //! electron from LS pairs: dPhi vs pt of triggered electron
    TH2F                  *fULSElecDphiDiffMethod;  //! electron from ULS pairs: dPhi vs pt of triggered electron
    TH2F                  *fLSElecDphiDiffMethod;   //! electron from LS pairs: dPhi vs pt of triggered electron

   
    TH1F                  *fNoPartnerNoT; //!
    THnSparse             *fNoPartnerNoTPt2; //!
    TH1F                  *fTPartnerNoT; //!
    THnSparse             *fTPartnerNoTPt2; //!
    TH3F                  *fElecHadTrigger;         //!
    TH2F                  *fElecHadTriggerLS;         //!
    TH2F                  *fElecHadTriggerULS;         //!
    TH2F                  *fElecHadTriggerLSNoP;         //!
    TH2F                  *fElecHadTriggerULSNoP;         //!
    TH2F                  *fElecHadTriggerLSNoPCorr;         //!
    TH2F                  *fElecHadTriggerULSNoPCorr;         //!
    TH2F                  *fElecHadTriggerULSNoPCorrTrue;         //!

    TH2F                  *fHadContTrigger;         //!
    TH2F                  *fHadElecTrigger;         //!
    TH2F                  *fNonElecHadTrigger;      //!
    TH2F                  *fHadNonElecTrigger;      //!
    THnSparse             *fInclElecHa;             //!
    THnSparse             *fLSElecHa;               //!
    THnSparse             *fULSElecHa;              //!
    THnSparse             *fULSElecHaTrue;          //!
    THnSparse             *fSignalElecHa;           //!
    THnSparse             *fBackgroundElecHa;       //!
    THnSparse             *fMCElecHaHadron;         //!
    THnSparse             *fElecHaHa;               //!
    THnSparse             *fElecHaLSNoPartner;      //!
    THnSparse             *fElecHaULSNoPartner;     //!

    THnSparse             *fElecHaLSNoPartnerCorrTrue; //!
    THnSparse             *fElecHaULSNoPartnerCorrTrue; //!
    THnSparse             *fElecHaLSNoPartnerCorr; //!
    THnSparse             *fElecHaULSNoPartnerCorr; //!


    THnSparse             *fMCElecHaTruePartner;    //!
    THnSparse             *fMCElecHaNoPartner;      //!
    THnSparse             *fMCElecHaRemovedPartner; //!
    TH2F                  *fMCElecHaTruePartnerTrigger;        //!
    TH2F                  *fMCElecHaTruePartnerTriggerWW;        //!
    TH2F                  *fMCElecHaNoPartnerTrigger;          //!
    TH2F                  *fMCElecHaNoPartnerTriggerWW;          //!
    TH2F                  *fMCElecHaRemovedPartnerTrigger; //!
    THnSparse             *fElecHaMixedEvent;       //!
    THnSparse             *fLSElecHaMixedEvent;     //!
    THnSparse             *fULSElecHaMixedEvent;    //!
    THnSparse             *fTagHaMixedEvent; //!
    THnSparse             *fNonTagHaMixedEvent; //!


    TH3F                  *fElecLPTrigger;          //!
    TH2F                  *fElecLPTriggerLS;         //!
    TH2F                  *fElecLPTriggerULS;        //!
    TH2F                  *fElecLPTriggerLSNoP;         //!
    TH2F                  *fElecLPTriggerULSNoP;        //!
    TH2F                  *fElecLPTriggerULSNoPCorr; //!
    TH2F                  *fHadContLPTrigger;       //!
    TH2F                  *fLPElecTrigger;          //!
    TH2F                  *fLPNonElecTrigger;       //!
    TH2F                  *fNonElecLPTrigger;       //!
    THnSparse             *fInclElecLP;             //!
    THnSparse             *fLSElecLP;               //! 
    THnSparse             *fULSElecLP;              //! 
    THnSparse             *fMCElecLPHadron;         //! 
    THnSparse             *fElecLPHa;               //!
    THnSparse             *fElecLPLSNoPartner;      //! 
    THnSparse             *fElecLPULSNoPartner;     //!
    THnSparse             *fElecLPLSNoPartnerCorrTrue; //!
    THnSparse             *fElecLPULSNoPartnerCorrTrue; //!
    THnSparse             *fElecLPLSNoPartnerCorr; //!
    THnSparse             *fElecLPULSNoPartnerCorr; //!

    THnSparse             *fMCElecLPTruePartner;    //! 
    THnSparse             *fMCElecLPNoPartner;      //! 
    THnSparse             *fMCElecLPRemovedPartner; //!
    TH2F                  *fMCElecLPTruePartnerTrigger;        //!
    TH2F                  *fMCElecLPNoPartnerTrigger;          //!
    TH2F                  *fMCElecLPRemovedPartnerTrigger; //!
    THnSparse             *fElecLPMixedEvent;       //!
    THnSparse             *fLSElecLPMixedEvent;     //!
    THnSparse             *fULSElecLPMixedEvent;    //!
    THnSparse             *fTagLPMixedEvent; //!
    THnSparse             *fNonTagLPMixedEvent; //!
   
    
    TH2F                  *fCheckMCVertex;           //!
    
    TH2F                  *fCheckMCPtvsRecPtHad;     //!
    TH2F                  *fCheckMCEtavsRecEtaHad;   //!
    TH2F                  *fCheckMCPhivsRecPhiHad;   //!
    THnSparse             *fMCHadPtEtaPhiVtx;        //!
    TH2F                  *fRecHadMCSecondaryCont;   //!
    THnSparse             *fRecHadMCPtEtaPhiVtx;     //!
    THnSparse             *fRecHadPtEtaPhiVtx;       //!
    THnSparse             *fRecHadPtEtaPhiVtxWRecEff; //!

    TH2F                  *fCheckMCPtvsRecPtEle;     //!
    TH1F                  *fRecHFE; //!
    THnSparse             *fMCElecPtEtaPhiVtx;       //!
    TH2F                  *fRecElecMCSecondaryCont;   //!
    THnSparse             *fRecElecPtEtaPhiVtx;      //!
    THnSparse             *fRecElecPtEtaPhiVtxWRecEff;  //!
    THnSparse             *fRecElecMCPtEtaPhiVtx;    //!
    TH1F                  *fMCElecPDG;               //!
    THnSparse             *fMCElecPtEtaPhiStrictVtx; //!

    THnSparse             *fMCPi0Prod;               //!
    THnSparse             *fMCEtaProd;               //!
    THnSparse             *fMCPiPlusProd;            //!
    THnSparse             *fMCPiPlusProdV2;          //!
    THnSparse             *fMCBGProd;                //!
    THnSparse             *fMCLeadingParticle;       //!
    TH3F                  *fCompareLPRecCheck;       //!

    AliEventPoolManager   *fMCTruePoolMgr;            //! event pool manager
    THnSparse             *fTrueMCHadronEventCuts;    //!
    THnSparse             *fTrueMCHadronEventCutsZvtx;//!
    THnSparse             *fTrueMCHadronEventCutsZvtxMEv; //!
    THnSparse             *fTrueMCHadron;             //!
    TH3F                  *fTrueMCElecHaTriggerEventCuts; //!
    TH3F                  *fTrueMCElecHaTrigger;      //!
    THnSparse             *fTrueMCLPEventCuts;        //!
    THnSparse             *fTrueMCLPEventCutsZvtx;    //!
    THnSparse             *fTrueMCLPEventCutsZvtxMEv; //!
    THnSparse             *fTrueMCLP;                 //!
    TH3F                  *fTrueMCElecLPTriggerEventCuts; //!
    TH3F                  *fTrueMCElecLPTrigger;      //!
    TH3F                  *fTrueElectronEta;          //!
    TH2F                  *fRecHFEEtaWRecEff;         //!
    TH2F                  *fTrueLPinAcceptanceEta;    //!
    TH2F                  *fTrueLPEta;                //!
    TH2F                  *fRecLPEta;                 //!
    TH2F                  *fTrueHadronEta;            //!
    TH2F                  *fRecHadronEtaWRecEff;      //!
    TH3F                  *fCompareLP;                //!


    AliESDv0KineCuts *fV0cutsESD;        //! ESD V0 cuts
    AliAODv0KineCuts *fV0cuts;           //! AOD V0 cuts
    TObjArray *fV0electrons;             //! array with pointer to identified particles from V0 decays (electrons)
    TObjArray *fV0pions;                 //! array with pointer to identified particles from V0 decays (pions)
    TObjArray *fV0protons;               //! array with pointer to identified particles from V0 decays (ptotons)
    TH2F      *fhArmenteros;             //!
    TH1F      *fEventsPerRun;            //!
    TH2F      *fTRDnTrackRun;            //!
    Int_t     *fV0tags;                  //!
    void      FindV0CandidatesAOD(AliAODEvent *Event);
    void      FindV0CandidatesESD(AliESDEvent *Event);
    void      ClearV0PIDList();
    void      TRDQA(Int_t RunNumber, const AliVVertex *pVtx, Int_t nMother, Int_t listMother[], Double_t EventWeight);
    void      FillV0Histograms(AliVTrack* track, Int_t Species, Int_t RunNumber);
    THnSparse *fTRDEtaPhi;               //!
    THnSparse *fTRDNTracklets;           //!
    THnSparse *fTRDV0NTracklets;         //!
    THnSparse *fTRDSpectra;              //!
    THnSparse *fTRDV0Spectra;            //!
    THnSparse *fTRDMCSpectra;            //!

    
    AliAnalysisTaskHaHFECorrel(const AliAnalysisTaskHaHFECorrel&);
    AliAnalysisTaskHaHFECorrel& operator=(const AliAnalysisTaskHaHFECorrel&);
    
    ClassDef(AliAnalysisTaskHaHFECorrel, 5);
};


// class storing reduced track information for mixed event pool

class AliBasicParticleHaHFE : public AliVParticle
{
 public:
 AliBasicParticleHaHFE() 
   : fID(0), fEta(0), fPhi(0), fpT(0), fCharge(0), fULSpartner(0), fLSpartner(0) , fIDLSPartner(0), fIDULSPartner(0), fWeightLSPartner(0), fWeightULSPartner(0), fTrueULSPartner(kFALSE), fTruePartnerMCPt(-999), fTruePartnerRecPt(-999), fIsPhotonic(kFALSE), fIsHadron(kFALSE), fLabel(0)
    {
      fExtTrackParam = AliExternalTrackParam();
    }
 AliBasicParticleHaHFE(Int_t id, Float_t eta, Float_t phi, Float_t pt, Short_t charge, Short_t LS, Short_t ULS, Int_t *LSPartner, Int_t *ULSPartner, Float_t *LSPartnerWeight, Float_t *ULSPartnerWeight, Bool_t trueULSPartner, Float_t truePartnerMCPt, Float_t truePartnerRecPt, Bool_t isPhotonic, Bool_t isHadron, Int_t label, AliExternalTrackParam & ExtTrackParam)
   : fID(id), fEta(eta), fPhi(phi), fpT(pt), fCharge(charge), fULSpartner(ULS), fLSpartner(LS), fIDLSPartner(0), fIDULSPartner(0),  fWeightLSPartner(0), fWeightULSPartner(0),  fTrueULSPartner(trueULSPartner), fTruePartnerMCPt(truePartnerMCPt), fTruePartnerRecPt(truePartnerRecPt), fIsPhotonic(isPhotonic), fIsHadron(isHadron), fLabel(label), fExtTrackParam(ExtTrackParam)
  {
 
    
    fIDLSPartner = new Int_t[LS];
    fIDULSPartner = new Int_t[ULS];
    fWeightLSPartner = new Float_t[LS];
    fWeightULSPartner = new Float_t[ULS];
    for (Int_t i=0; i<LS; i++) {
      fIDLSPartner[i]=LSPartner[i];
      fWeightLSPartner[i]=LSPartnerWeight[i];
    }
    for (Int_t i=0; i<ULS; i++) {
      fIDULSPartner[i]=ULSPartner[i];
      fWeightULSPartner[i]=ULSPartnerWeight[i];
    }
    fExtTrackParam = ExtTrackParam;
  }
  virtual ~AliBasicParticleHaHFE() {
    if (fIDLSPartner)  delete[] fIDLSPartner;
    if (fIDULSPartner) delete[] fIDULSPartner;
    if (fWeightLSPartner)  delete[] fWeightLSPartner;
    if (fWeightULSPartner) delete[] fWeightULSPartner;
  }
  AliBasicParticleHaHFE(const AliBasicParticleHaHFE &CopyClass) 
    : fID(CopyClass.fID), fEta(CopyClass.fEta), fPhi(CopyClass.fPhi), fpT(CopyClass.fpT), fCharge(CopyClass.fCharge), fULSpartner(CopyClass.fULSpartner), fLSpartner(CopyClass.fLSpartner), fIDLSPartner(0), fIDULSPartner(0), fWeightLSPartner(0), fWeightULSPartner(0), fTrueULSPartner(CopyClass.fTrueULSPartner), fTruePartnerMCPt(CopyClass.fTruePartnerMCPt),fTruePartnerRecPt(CopyClass.fTruePartnerRecPt),fIsPhotonic(CopyClass.fIsPhotonic), fIsHadron(CopyClass.fIsHadron), fLabel(CopyClass.fLabel), fExtTrackParam(CopyClass.fExtTrackParam)
    {
      fIDLSPartner = new Int_t[CopyClass.fLSpartner];
      fIDULSPartner = new Int_t[CopyClass.fULSpartner];
      for (Int_t i=0; i<fLSpartner; i++) {fIDLSPartner[i]=CopyClass.fIDLSPartner[i];}
      for (Int_t i=0; i<fULSpartner; i++) {fIDULSPartner[i]=CopyClass.fIDULSPartner[i];}
      fWeightLSPartner = new Float_t[CopyClass.fLSpartner];
      fWeightULSPartner = new Float_t[CopyClass.fULSpartner];
      for (Int_t i=0; i<fLSpartner; i++) {fWeightLSPartner[i]=CopyClass.fWeightLSPartner[i];}
      for (Int_t i=0; i<fULSpartner; i++) {fWeightULSPartner[i]=CopyClass.fWeightULSPartner[i];}
    }
  AliBasicParticleHaHFE& operator=(const AliBasicParticleHaHFE &CopyClass) 
    {
      if (this==&CopyClass) return *this;
      if (fIDLSPartner) delete[] fIDLSPartner;
      if (fIDULSPartner) delete[] fIDULSPartner;
      if (fWeightLSPartner) delete[] fWeightLSPartner;
      if (fWeightULSPartner) delete[] fWeightULSPartner;
      fID=CopyClass.fID;
      fEta=CopyClass.fEta;
      fPhi=CopyClass.fPhi;
      fpT=CopyClass.fpT;
      fCharge=CopyClass.fCharge;
      fULSpartner=CopyClass.fULSpartner;
      fLSpartner=CopyClass.fLSpartner;
      fTrueULSPartner=CopyClass.fTrueULSPartner;
      fTruePartnerMCPt=CopyClass.fTruePartnerMCPt;
      fTruePartnerRecPt=CopyClass.fTruePartnerRecPt;
      fIsPhotonic=CopyClass.fIsPhotonic;
      fIsHadron=CopyClass.fIsHadron;
      fIDLSPartner = new Int_t[CopyClass.fLSpartner];
      fIDULSPartner = new Int_t[CopyClass.fULSpartner];
      for (Int_t i=0; i<fLSpartner; i++) {fIDLSPartner[i]=CopyClass.fIDLSPartner[i];}
      for (Int_t i=0; i<fULSpartner; i++) {fIDULSPartner[i]=CopyClass.fIDULSPartner[i];}
      fWeightLSPartner = new Float_t[CopyClass.fLSpartner];
      fWeightULSPartner = new Float_t[CopyClass.fULSpartner];
      for (Int_t i=0; i<fLSpartner; i++) {fWeightLSPartner[i]=CopyClass.fWeightLSPartner[i];}
      for (Int_t i=0; i<fULSpartner; i++) {fWeightULSPartner[i]=CopyClass.fWeightULSPartner[i];}
      fLabel=CopyClass.fLabel;
      fExtTrackParam = CopyClass.fExtTrackParam;
      return *this;
    }

  // kinematics
  virtual Double_t Px() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Py() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Pz() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Pt() const { return fpT; }
  virtual Double_t P() const { AliFatal("Not implemented"); return 0; }
  virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

  virtual Double_t Xv() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Yv() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Zv() const { AliFatal("Not implemented"); return 0; }
  virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

  virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Phi()        const { return fPhi; }
  virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }


  virtual Double_t E()          const { AliFatal("Not implemented"); return 0; }
  virtual Double_t M()          const { AliFatal("Not implemented"); return 0; }
 
  virtual Double_t Eta()        const { return fEta; }
  virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }

  virtual Short_t Charge()      const { return fCharge; }
  virtual Int_t   GetLabel()    const { return fLabel; }
    // PID
  virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }
  virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }
  
  virtual Short_t LS()          const {return fLSpartner; }
  virtual Short_t ULS()         const {return fULSpartner; }
  virtual Short_t ID()          const {return fID;}
  virtual Int_t  LSPartner(Int_t i)   const {return fIDLSPartner[i];}
  virtual Int_t  ULSPartner(Int_t i)  const {return fIDULSPartner[i];}
  virtual Float_t  LSPartnerWeight(Int_t i)   const {return fWeightLSPartner[i];}
  virtual Float_t  ULSPartnerWeight(Int_t i)  const {return fWeightULSPartner[i];}
  virtual Bool_t TruePartner() const {return fTrueULSPartner;}
  virtual Float_t TruePartnerMCPt() const {return fTruePartnerMCPt;}
  virtual Float_t TruePartnerRecPt() const {return fTruePartnerRecPt;}
  virtual Bool_t IsPhotonic() const {return fIsPhotonic;}
  virtual Bool_t IsHadron() const {return fIsHadron;}
  AliExternalTrackParam GetExtTrackParam() {return fExtTrackParam;};
  //  virtual Int_t  PDG() const {return fPDG;}


 private:
  Int_t   fID;             // particle id
  Float_t fEta;            // eta
  Float_t fPhi;            // phi
  Float_t fpT;             // pT
  Short_t fCharge;         // charge
  Short_t fULSpartner;     // no of ULS partner
  Short_t fLSpartner;      // no of LS partner

  Int_t*  fIDLSPartner;    //! particle id of LS Partner
  Int_t*  fIDULSPartner;   //! partilce id of ULS partner
  Float_t*  fWeightLSPartner;    //! CorrWeight for NoPartner
  Float_t*  fWeightULSPartner;   //! CorrWeight for NoPartner
  Bool_t  fTrueULSPartner; // check if true partner was tagged
  Float_t fTruePartnerMCPt; // only for
  Float_t fTruePartnerRecPt; //
  Bool_t  fIsPhotonic;     //
  Bool_t  fIsHadron;       // only for MC
  Int_t   fLabel;         //
  AliExternalTrackParam fExtTrackParam; //

  ClassDef(AliBasicParticleHaHFE, 6); // class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};



#endif
