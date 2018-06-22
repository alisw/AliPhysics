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
    AliVTrack* FindLPAndHFE(TObjArray* RedTracksHFE, const AliVVertex *pVtx, Int_t nMother, Int_t listMother[], Double_t mult, Bool_t &EvContTP, Bool_t &EvContNTP);
    void FindPhotonicPartner(Int_t iTracks, AliVTrack* track,  const AliVVertex *pVtx, Int_t nMother, Int_t listMother[], Int_t &LSPartner, Int_t&ULSPartner, Int_t *LSPartnerID, Int_t *ULSPartnerID, Bool_t &trueULSPartner, Bool_t &iHsPhotonic);
    void CorrelateElectron(TObjArray* RedTracksHFE);

    void CorrelateLP(AliVTrack* LPtrack,  const AliVVertex* pVtx, Int_t nMother, Int_t listMother[], TObjArray* RedTracksHFE);
    void CorrelateLPMixedEvent(AliVTrack* LPtrack, Float_t mult, Float_t zVtx, Float_t maxPt, Bool_t EvContTP, Bool_t EvContNTP);
    
    void CorrelateHadron(TObjArray* RedTracksHFE,  const AliVVertex* pVtx, Int_t nMother, Int_t listMother[], Float_t mult, Float_t maxPt);
    void CorrelateHadronMixedEvent(Float_t mult, const AliVVertex* zVtx, Float_t maxPt, Int_t nMother, Int_t listMother[], Bool_t EvContTP, Bool_t EvContNTP);

    void CorrelateWithHadrons(AliVTrack* TriggerTrack, const AliVVertex* pVtx, Int_t nMother, Int_t listMother[], Bool_t FillHadron, Bool_t FillLP,Bool_t** NonElecIsTrigger, Double_t *NonElecIsTriggerPt, Double_t *NonElecIsTriggerWeight, Int_t NumElectronsInEvent); 

    void MCTruthCorrelation(Bool_t AfterEventCuts, Int_t RecLPLabel, Float_t pVtxZ, Float_t mult, Int_t &LPinAcceptance, Int_t &LP);

    //********************MC
    void  MCEfficiencyCorrections(const AliVVertex * RecVertex);
    Int_t HFEisCharmOrBeauty(Int_t ElectronIndex);

    //*********************ANALYSIS Helper
    Bool_t ChargedHadronTrackCuts(const AliVVertex *pVtx,AliVTrack *Htrack, Int_t nMother, Int_t listMother[]);
    Bool_t ChargedHadronPIDCuts(AliVTrack *Htrack);

    Bool_t AssoHadronPIDCuts(AliVTrack *Htrack);

    Bool_t InclElecTrackCuts(const AliVVertex *pVtx,AliVTrack *ietrack, Int_t nMother, Int_t listMother[]);
    Bool_t InclElecPIDCuts(AliVTrack *track, Bool_t IsPrimary);

    Bool_t PhotElecPIDCuts(AliVTrack *track);
    Bool_t PhotElecTrackCuts(const AliVVertex *pVtx,AliVTrack *aetrack, Int_t nMother, Int_t listMother[]);
    
    void EvaluateTaggingEfficiency(AliVTrack * track, Int_t LSPartner, Int_t ULSPartner, Bool_t trueULSPartner); 
    Bool_t CloneAndReduceTrackList(TObjArray* RedTracks, AliVTrack* track, Int_t LSPartner, Int_t ULSPartner, Int_t *LSPartnerID, Int_t *ULSPartnerID, Bool_t trueULSPartner, Bool_t isPhotonic, Bool_t isHadron);

    void BinLogX(TAxis *axis);
    void SetPDGAxis(TAxis *axis, std::vector<TString> PDGLabel);
    void CheckHadronIsTrigger(Double_t ptE, Bool_t *HadronIsTrigger);
    void CheckElectronIsTrigger(Double_t ptH, Bool_t *ElectronIsTrigger) ;
    Bool_t PassEventBias( const AliVVertex *pVtx, Int_t nMother, Int_t *listMother);    


    //**************  SETTINGS
    void SetMC (Bool_t IsMC) { fIsMC=IsMC;};
    void SetAODanalysis(Bool_t IsAOD) {fIsAOD = IsAOD;};
    void SetTender (Bool_t UseTender) {fUseTender = UseTender;};
    void SetPeriod (Double_t period) {fWhichPeriod = period;};

    void SetOnlyEfficiency() {
      fTRDQA = kFALSE;
      fCorrHadron = kFALSE;
      fCorrLParticle = kFALSE;
      fMixedEvent = kFALSE;
      fMCTrueCorrelation = kFALSE;
    }

    

    void SetTRDQA(Bool_t TRDQA) {fTRDQA=TRDQA;};
    void SetPtMinEvent(Double_t PtMin) {fMinPtEvent=PtMin;};
    void SetPtMaxEvent(Double_t PtMax) {fMaxPtEvent=PtMax;};


    void SetEtaMax(Double_t EtaMax) {
      fMaxElectronEta = TMath::Abs(EtaMax);
      fMinElectronEta = -TMath::Abs(EtaMax);
    };

    void SetTPCnCut(Int_t TPCnCut) {fTPCnCut = TPCnCut;};
    void SetTPCnCutdEdx(Int_t TPCnCutdEdx) {fTPCndEdxCut = TPCnCutdEdx;};
    void SetITSnCut (Int_t ITSnCut) {fITSnCut = ITSnCut;};
   

    void SetPhotElecPtCut (Double_t AssPtCut) {fPhotElecPtCut = AssPtCut;};
    void SetPhotElecTPCnCut (Int_t AssTPCnCut) {fPhotElecTPCnCut = AssTPCnCut;};
    void SetPhotElecITSrefitCut(Bool_t AssITSrefitCut) {fPhotElecITSrefitCut = AssITSrefitCut;};  

    void SetHTPCnCut(Int_t HTPCnCut) {fHTPCnCut = HTPCnCut;}
    void SetHITSrefitCut(Bool_t HITSrefitCut) {fHITSrefitCut = HITSrefitCut;};   
    void SetHTPCrefitCut(Bool_t HTPCrefitCut) {fHTPCrefitCut = HTPCrefitCut;};

   
    void SetUseTRD(Bool_t UseTRD) {fUseTRD = UseTRD;}
    void SetUseITS(Bool_t UseITS) {fUseITS = UseITS;}
    void SetSigmaITScut(Double_t SigmaITScut) {fSigmaITScut = SigmaITScut;};
    void SetSigmaTOFcut(Double_t SigmaTOFcut) {fSigmaTOFcut = SigmaTOFcut;};
    void SetSigmaTPCcut(Double_t SigmaTPCcut) {fSigmaTPCcut = SigmaTPCcut;};

  
    void SetRecEff(Bool_t RecEff) { 
      fRecEff=RecEff;
    }
    void SetTagEff(Bool_t TagEff) {
      fTagEff=TagEff;
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

    void SetOpeningAngleCut(Bool_t OpeningAngleCut) {fOpeningAngleCut=OpeningAngleCut;};
    void SetInvmassCut(Double_t InvmassCut) {fInvmassCut=InvmassCut;};

    void SetPi0WeightToData(TH1F &  WPion) {fCorrectPiontoData = WPion; fCorrectPiontoData.SetName("fCorrectPiontoData");}
    void SetEtaWeightToData(TH1F &  WEta)  {fCorrectEtatoData  = WEta; fCorrectEtatoData.SetName("fCorrectEtatoData");}
    void SetHadRecEff(TH3F & HadRecEff) {fHadRecEff = HadRecEff; fHadRecEff.SetName("fHadRecEff");}
    void SetEleRecEff(TH2F & EleRecEff) {fEleRecEff = EleRecEff; fEleRecEff.SetName("fEleRecEff");}
    void SetSPDnTrAvg(TProfile & SPDnTrAvg) {fSPDnTrAvg = SPDnTrAvg; fSPDnTrAvg.SetName("fSPDnTrAvg");}
    

    Bool_t   ESDkTrkGlobalNoDCA(AliVTrack* Vtrack);

  


 private:
    
    Bool_t                IsPhotonicElectron(Int_t Label1) const;
    Bool_t                HaveSameMother(Int_t Label1, Int_t Label2) const;
    Double_t              GetDeltaPhi(Double_t phiA,Double_t phiB) const;
    Double_t              GetDeltaEta(Double_t etaA,Double_t etaB) const;
    Double_t              Eta2y(Double_t pt, Double_t m, Double_t eta) const;
    Double_t              GetHadronRecEff(Double_t pt, Double_t phi, Double_t eta, Double_t zVtx);
    Double_t              GetElectronRecEff(Double_t pt, Double_t phi, Double_t eta, Double_t zVtx);

    
    Bool_t                fUseTender;               // Use tender
    Int_t                 fWhichPeriod;             // period
    Bool_t                fUseKFforPhotonicPartner; //default ist DCA

    Float_t               fMaxPtEvent;              //
    Float_t               fMinPtEvent;

    Double_t              fMaxElectronEta;          //
    Double_t              fMinElectronEta;          //
    Double_t              fMaxHadronEta;            //
    Double_t              fMinHadronEta;            //

    // HFECuts
    Int_t                 fTPCnCut;                 // TPC number of clusters for tagged electron
    Int_t                 fTPCndEdxCut;             //
    Int_t                 fITSnCut;                 // ITs number of clusters for tagged electrons 

    Bool_t                fUseTRD;                  //
    Bool_t                fUseITS;                  //
    Double_t              fSigmaITScut;             // ITS nSigma cut
    Double_t              fSigmaTOFcut;             // TOF nSigma cut
    Double_t              fSigmaTPCcut;             // lower TPC nSigma cut 

    // Photonic  Electrons
    Double_t              fPhotElecPtCut;                // pt cut for associated electron
    Double_t              fPhotElecSigmaTPCcut;          //
    Int_t                 fPhotElecTPCnCut;              // TPC number of clusters for associated electron
    Bool_t                fPhotElecITSrefitCut;          // ITS refit for associated electron
    

    // Associate Hadron (non Electron)
    Double_t              fAssNonEleTPCcut;         //  

    // Hadron Cut
    Int_t                 fHTPCnCut;                // TPC number of clusters for trigger hadron
    Bool_t                fHITSrefitCut;            // ITS refit for trigger hadron
    Bool_t                fHTPCrefitCut;            // TPC refit for trigger hadron 
    
    Double_t              fOpeningAngleCut;         // openingAngle cut for non-HFE selection
    Double_t              fInvmassCut;              // invariant mass cut  for non-HFE selection
    Double_t              fChi2Cut;                 //! used?? Chi2 cut  for non-HFE selection
    Double_t              fDCAcut;                  //! used?? DCA cut  for non-HFE selection
  
    // ******* Switch for analysis modes
    Bool_t                fTRDQA;                   // TRDQA
    Bool_t                fMCTrueCorrelation;       //
    Bool_t                fCorrHadron;              // Choose Hadron-HFE Correl
    Bool_t                fCorrLParticle;           // Choose LP-HFE Correl
    Bool_t                fMixedEvent;              // Fill Mixed Event for the cases chosen above
    Bool_t                fPionEtaProduction;       //
    Bool_t                fRecEff;                  //
    Bool_t                fTagEff;                  //
    Bool_t                fHadCont;                 //
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
    TH2F                  *fTrkpt;                  //! track pt for different cuts
    TH2F                  *fEtaVtxZ;                //! Eta vs Vtx z (check for ITS acceptance problem)

    TH2F                  *fSPDVtxRes;              //!
    TH2F                  *fDiffSPDPrimVtx;         //!
    TH2F                  *fSPDnTrAcc;              //!
    TH2F                  *fSPDnTrCorrMax;          //!
    TH2F                  *fSPDnTrGen;              //!
    TH2F                  *fDiffSPDMCVtx;           //!
    THnSparseF            *fnTrAccMaxGen;           //!
    THnSparseF            *fnTrAccMinGen;           //!
    THnSparseF            *fnTrAccMeanGen;          //!
    THnSparseF            *fnTrAccMax;              //!
    THnSparseF            *fnTrAccMin;              //!
    THnSparseF            *fnTrAccMean;             //!


    
    THnSparse             *fMultiplicity;	    //! multiplicity distribution
    TH3F                  *fSPDMultiplicity;        //!
    Int_t                 *fRunList;                //!

    TH2F                  *fElectronTrackCuts;      //! 
    TH2F                  *fElectronTrackTPCNcls;   //! 
    TH2F                  *fElectronTrackTPCNclsdEdx; //! 
    TH2F                  *fElectronTrackTPCFrac;    //! 
    TH2F                  *fElectronTrackITSNcls;   //!
    TH3F                  *fElectronTrackITSLayer;   //!
    TH3F                  *fElectronTrackRefit;     //!
    TH3F                  *fElectronTrackDCA;       //! 
    
    TH2F                  *fHadronTrackTPCNcls;   //! 
    TH3F                  *fHadronTrackRefit;     //!
    TH3F                  *fHadronTrackDCA;       //! 



   

    TH2F                  *fHistITSnSig;            //! ITS sigma vs p
    TH2F                  *fHistTOFnSig;            //! TOF sigma vs p
    TH2F                  *fHistTPCnSig;            //! TPC sigma vs p
    TH2F                  *fHistTPCnSigITScut;      //! TPC sigma vs p (ITS cut)
    TH2F                  *fHistTPCnSigTOFcut;      //! TPC sigma vs p (TOF cut)
    TH2F                  *fHistTPCnSigITSTOFcut;   //! TPC sigma vs p (ITS+TOF cuts)

    THnSparse             *fCheckNHadronScaling;    //!
    THnSparse             *fCheckNPhotHadScaling;  //!
    TH3F                  *fCheckTaggedEvent;       //!

    TH2F                  *fHadContPvsPt;           //!
    TH3F                  *fHadContEtaPhiPt;        //!
    TH3F                  *fHadContTPCEtaPhiPt;     //!
    THnSparse             *fHadContPPhiEtaTPC;      //!
    THnSparse             *fHadContamination;       //! HadronicContaminationTOF
    THnSparse             *fHadContaminationPt;       //! HadronicContaminationTOF
    THnSparse             *fHadContMC;              //!
    THnSparse             *fHadContMCPt;              //!
    

  
    TH2F                  *fInclElecPtEta;             //! inclusive electron p
    TH2F                  *fInclElecPtEtaWW;             //! inclusive electron p
    TH1F                  *fInclElecP;             //! inclusive electron p
    TH1F                  *fULSElecPt;              //! ULS electron pt (after IM cut)
    TH1F                  *fULSElecPtWW;              //! ULS electron pt (after IM cut)
    TH1F                  *fLSElecPt;               //! LS electron pt (after IM cut)
    TH1F                  *fLSElecPtWW;               //! LS electron pt (after IM cut)
    TH2F                  *fInvmassLS;              //! Inv mass of LS (e,e)
    TH2F                  *fInvmassULS;             //! Inv mass of ULS (e,e)
    TH3F                  *fInvmassMCTrue;           //! Inv mass of ULS (e,e)
    TH2F                  *fOpeningAngleLS;         //! opening angle for LS pairs
    TH2F                  *fOpeningAngleULS;        //! opening angle for ULS pairs
    TH2F                  *fCheckLSULS;             //! check no of LS/ULS partner per electron
    TH3F                  *fTagEtaPhiPt;            //!
    TH3F                  *fTagEtaZvtxPt;            //!
    TH3F                  *fTagEtaPhiPtwW;           //!
    TH3F                  *fTagEtaZvtxPtwW;           //!
    TH3F                  *fNonTagEtaPhiPt;         //!
    TH3F                  *fNonTagEtaZvtxPt;         //!
    TH3F                  *fNonTagEtaPhiPtwW;       //!
    TH3F                  *fNonTagEtaZvtxPtwW;        //!



    THnSparse             *fTagMotherPt;              //!
    THnSparse             *fTagEffIncl;             //! 
    THnSparse             *fTagEffLS;               //!
    THnSparse             *fTagEffULS;              //!
    THnSparse             *fTagTruePairs;           //!
    THnSparse             *fTagEffInclWoWeight;      //! 
    THnSparse             *fTagEffLSWoWeight;        //!
    THnSparse             *fTagEffULSWoWeight;        //!
    THnSparse             *fTagTruePairsWoWeight;     //!

    TH1F                  fCorrectPiontoData;   
    Double_t              GetPionWeight(Double_t pt);
    TH1F                  fCorrectEtatoData;       
    Double_t              GetEtaWeight(Double_t pt);
    TH3F                  fHadRecEff;
    TH2F                  fEleRecEff;
    TProfile              fSPDnTrAvg;

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
    TH1F                  *fTPartnerNoT; //!
    TH3F                  *fElecHadTrigger;         //!
    TH2F                  *fElecHadTriggerLS;         //!
    TH2F                  *fElecHadTriggerULS;         //!
    TH2F                  *fElecHadTriggerLSNoP;         //!
    TH2F                  *fElecHadTriggerULSNoP;         //!
    TH2F                  *fElecHadTriggerLSNoPCorr;         //!
    TH2F                  *fElecHadTriggerULSNoPCorr;         //!

    TH2F                  *fHadContTrigger;         //!
    TH2F                  *fHadElecTrigger;         //!
    TH2F                  *fNonElecHadTrigger;      //!
    TH2F                  *fHadNonElecTrigger;      //!
    THnSparse             *fInclElecHa;             //!
    THnSparse             *fLSElecHa;               //!
    THnSparse             *fULSElecHa;              //!
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
    TH2F                  *fMCElecHaNoPartnerTrigger;          //!
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
    THnSparse             *fRecHadMCPtEtaPhiVtx;     //!
    THnSparse             *fRecHadPtEtaPhiVtx;       //!
    THnSparse             *fRecHadPtEtaPhiVtxwW;       //!

    TH2F                  *fCheckMCPtvsRecPtEle;     //!
    THnSparse             *fMCElecPtEtaPhiVtx;       //!
    THnSparse             *fRecElecPtEtaPhiVtx;      //!
    THnSparse             *fRecElecPtEtaPhiVtxwW;    //!
    THnSparse             *fRecElecMCPtEtaPhiVtx;    //!
    TH1F                  *fMCElecPDG;               //!
    THnSparse             *fMCElecPtEtaPhiStrictVtx; //!

    THnSparse             *fMCPi0Prod;               //!
    THnSparse             *fMCEtaProd;               //!
    THnSparse             *fMCPiPlusProd;            //!
    THnSparse             *fMCPiPlusProdV2;          //!
    THnSparse             *fMCLeadingParticle;       //!

    AliEventPoolManager   *fMCTruePoolMgr;           //! event pool manager
    THnSparse             *fTrueMCHadronEventCuts;   //!
    THnSparse             *fTrueMCHadronEventCutsZvtx; //!
    THnSparse             *fTrueMCHadronEventCutsZvtxMEv; //!
    THnSparse             *fTrueMCHadron;            //!
    TH3F                  *fTrueMCElecHaTriggerEventCuts; //!
    TH3F                  *fTrueMCElecHaTrigger; //!
    THnSparse             *fTrueMCLPEventCuts;   //!
    THnSparse             *fTrueMCLPEventCutsZvtx; //!
    THnSparse             *fTrueMCLPEventCutsZvtxMEv; //!
    THnSparse             *fTrueMCLP;   //!
    TH3F                  *fTrueMCElecLPTriggerEventCuts; //!
    TH3F                  *fTrueMCElecLPTrigger; //!
    TH2F                  *fTrueElectronEta; //!
    TH2F                  *fRecElectronEta; //!
    TH2F                  *fTrueLPinAcceptanceEta; //!
    TH2F                  *fTrueLPEta; //!
    TH2F                  *fRecLPEta; //!
    TH2F                  *fTrueHadronEta; //!
    TH2F                  *fRecHadronEta; //!
    TH3F                  *fCompareLP; //!


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
    void      TRDQA(Int_t RunNumber, const AliVVertex *pVtx, Int_t nMother, Int_t listMother[]);
    void      FillV0Histograms(AliVTrack* track, Int_t Species, Int_t RunNumber);
    THnSparse *fTRDEtaPhi;               //!
    THnSparse *fTRDNTracklets;           //!
    THnSparse *fTRDV0NTracklets;         //!
    THnSparse *fTRDSpectra;              //!
    THnSparse *fTRDV0Spectra;            //!
    THnSparse *fTRDMCSpectra;            //!

    
    AliAnalysisTaskHaHFECorrel(const AliAnalysisTaskHaHFECorrel&);
    AliAnalysisTaskHaHFECorrel& operator=(const AliAnalysisTaskHaHFECorrel&);
    
    ClassDef(AliAnalysisTaskHaHFECorrel, 4);
};


// class storing reduced track information for mixed event pool

class AliBasicParticleHaHFE : public AliVParticle
{
 public:
 AliBasicParticleHaHFE() 
   : fID(0), fEta(0), fPhi(0), fpT(0), fCharge(0), fULSpartner(0), fLSpartner(0) , fIDLSPartner(0), fIDULSPartner(0), fTrueULSPartner(kFALSE), fIsPhotonic(kFALSE), fIsHadron(kFALSE), fLabel(0)
    {}
 AliBasicParticleHaHFE(Int_t id, Float_t eta, Float_t phi, Float_t pt, Short_t charge, Short_t LS, Short_t ULS, Int_t *LSPartner, Int_t *ULSPartner, Bool_t trueULSPartner, Bool_t isPhotonic, Bool_t isHadron, Int_t label)
   : fID(id), fEta(eta), fPhi(phi), fpT(pt), fCharge(charge), fULSpartner(ULS), fLSpartner(LS), fIDLSPartner(0), fIDULSPartner(0), fTrueULSPartner(trueULSPartner), fIsPhotonic(isPhotonic), fIsHadron(isHadron), fLabel(label)
  {
    fIDLSPartner = new Int_t[LS];
    fIDULSPartner = new Int_t[ULS];
    for (Int_t i=0; i<LS; i++)  fIDLSPartner[i]=LSPartner[i];
    for (Int_t i=0; i<ULS; i++) fIDULSPartner[i]=ULSPartner[i];
  }
  virtual ~AliBasicParticleHaHFE() {
    if (fIDLSPartner)  delete fIDLSPartner;
    if (fIDULSPartner) delete fIDULSPartner;
  }
  AliBasicParticleHaHFE(const AliBasicParticleHaHFE &CopyClass) 
    : fID(CopyClass.fID), fEta(CopyClass.fEta), fPhi(CopyClass.fPhi), fpT(CopyClass.fpT), fCharge(CopyClass.fCharge), fULSpartner(CopyClass.fULSpartner), fLSpartner(CopyClass.fLSpartner), fIDLSPartner(0), fIDULSPartner(0), fTrueULSPartner(CopyClass.fTrueULSPartner), fIsPhotonic(CopyClass.fIsPhotonic), fIsHadron(CopyClass.fIsHadron), fLabel(CopyClass.fLabel)
    {
      fIDLSPartner = new Int_t[CopyClass.fLSpartner];
      fIDULSPartner = new Int_t[CopyClass.fULSpartner];
      for (Int_t i=0; i<fLSpartner; i++) {fIDLSPartner[i]=CopyClass.fIDLSPartner[i];}
      for (Int_t i=0; i<fULSpartner; i++) {fIDULSPartner[i]=CopyClass.fIDULSPartner[i];}
    }
  AliBasicParticleHaHFE& operator=(const AliBasicParticleHaHFE &CopyClass) 
    {
      if (this==&CopyClass) return *this;
      if (fIDLSPartner) delete[] fIDLSPartner;
      if (fIDULSPartner) delete[] fIDULSPartner;
      fID=CopyClass.fID;
      fEta=CopyClass.fEta;
      fPhi=CopyClass.fPhi;
      fpT=CopyClass.fpT;
      fCharge=CopyClass.fCharge;
      fULSpartner=CopyClass.fULSpartner;
      fLSpartner=CopyClass.fLSpartner;
      fTrueULSPartner=CopyClass.fTrueULSPartner;
      fIsPhotonic=CopyClass.fIsPhotonic;
      fIsHadron=CopyClass.fIsHadron;
      fIDLSPartner = new Int_t[CopyClass.fLSpartner];
      fIDULSPartner = new Int_t[CopyClass.fULSpartner];
      for (Int_t i=0; i<fLSpartner; i++) {fIDLSPartner[i]=CopyClass.fIDLSPartner[i];}
      for (Int_t i=0; i<fULSpartner; i++) {fIDULSPartner[i]=CopyClass.fIDULSPartner[i];}
      fLabel=CopyClass.fLabel;
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
  virtual Bool_t TruePartner() const {return fTrueULSPartner;}
  virtual Bool_t IsPhotonic() const {return fIsPhotonic;}
  virtual Bool_t IsHadron() const {return fIsHadron;}
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
  Bool_t  fTrueULSPartner; // check if true partner was tagged
  Bool_t  fIsPhotonic;     //
  Bool_t  fIsHadron;       // only for MC
  Int_t   fLabel;

  ClassDef(AliBasicParticleHaHFE, 4); // class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};



#endif
