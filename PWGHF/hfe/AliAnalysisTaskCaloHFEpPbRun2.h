/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

/////////////////////////////////////////////////////////////////////////////
//   Task for Heavy-flavour electrons analysis in p-Pb collisions
//   at 8.16 TeV with TPC + EMCal + DCal
//
//   - [Preliminary] Nuclear modification factor RpPb
//   - Multiplicity dependence (QpPb)
//   - HFe-h correlation
//
//   Author : Daichi Kawana (daichi.kawana@cern.ch)
//
/////////////////////////////////////////////////////////////////////////////

#ifndef AliAnalysisTaskCaloHFEpPbRun2_H
#define AliAnalysisTaskCaloHFEpPbRun2_H

class AliMultSelection;
// classes for MC
class AliAODMCHeader;
class AliAODMCParticle;
// classed for correlation
class AliEventPoolManager;

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

class AliAnalysisTaskCaloHFEpPbRun2 : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskCaloHFEpPbRun2();
                                AliAnalysisTaskCaloHFEpPbRun2(const char *name);
        virtual                 ~AliAnalysisTaskCaloHFEpPbRun2();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        //--- MC Switch ---//
        void SetMC(Bool_t flagMC) {fFlagMC = flagMC;};
        void SetRunPeriod(TString period) {fPeriodName = period;};
        //---EMCal/DCal switch---//
        void SetClusterTypeEMC(Bool_t flagClsEMCal) {fFlagClsTypeEMCal = flagClsEMCal;};
        void SetClusterTypeDCAL(Bool_t flagClsDCal) {fFlagClsTypeDCal = flagClsDCal;};
        //---EMCal Correction switch---//
        void SetEMCalCorrection(Bool_t flagEMCalCorrection) {fFlagEMCalCorrection = flagEMCalCorrection;};
        //---Trigger Selection switch---//
        void SetEG1(Bool_t flagEG1) {fFlagEG1 = flagEG1;};
        void SetEG2(Bool_t flagEG2) {fFlagEG2 = flagEG2;};
        void SetDG1(Bool_t flagDG1) {fFlagDG1 = flagDG1;};
        void SetDG2(Bool_t flagDG2) {fFlagDG2 = flagDG2;};
        //---Centrality selection---//
        void SetCentrality(Double_t min, Double_t max) {fCentMin = min, fCentMax = max;};
        //---ch particle eff---//
        void SetEffHadron(TH1D *hist) {if(fEffHadron) delete fEffHadron;fEffHadron = (TH1D*)hist->Clone();};

        //#########################//
        //Systematic uncertainties //
        //#########################//
        void SetTrackEta(Double_t low, Double_t high) {TrackEtaLow = low, TrackEtaHigh = high;};
        void SetTrackClust(Int_t TPC, Int_t ITS, Int_t Crossed) {NTPCClust = TPC, NITSClust = ITS, NCrossedRow = Crossed;};
        void SetDCA(Double_t xy, Double_t z) {DCAxy = xy, DCAz = z;};
        void SetTrackMatch(Double_t phi, Double_t eta) {TrackMatchPhi = phi, TrackMatchEta = eta;};
        void SetNsigma(Double_t low, Double_t high) {NsigmaLow = low, NsigmaHigh = high;};
        void SetM20(Double_t low, Double_t high) {M20Low = low, M20High = high;};
        void SetEop(Double_t low, Double_t high) {EopLow = low, EopHigh = high;};
        void SetMassPara(Double_t mimpt, Int_t TPC, Double_t mass) {AssoMinpT = mimpt, AssoNTPCClust = TPC, MassCut = mass;};

        //--- centrality, Multiplicity ---//
        void CheckCentrality(AliAODEvent* fAOD, Bool_t &IsCentPass);

        //--- two particles correlation ---//
        void ElectronHadCorrel(Int_t itrack, AliAODTrack *Trigtrack, THnSparse *SparseEHCorrl);
        void ElectronHadCorrelNopartner(Int_t itrack, Int_t pairtrack, AliAODTrack *Trigtrack, THnSparse *SparseEHCorrl);
        void MixedEvent(AliAODTrack *Trigtrack, THnSparse *SparseMixEHCorrl);
        TObjArray* CloneAndReduceTrackList();
        Bool_t PassHadronCuts(AliAODTrack *HadTrack);
        Double_t GetHadronEfficiency(Double_t pT);

        void FindMother(AliAODMCParticle* part, Int_t &label, Int_t &pid);
        Bool_t ConversionCheck(Int_t &pidM);
        Bool_t DalitzCheck(Int_t &pidM);
        Bool_t CharmCheck(Int_t &pidM, Int_t &pidGM);
        Bool_t EnhanceCharmCheck(Int_t &pdgM, Int_t &pdgGM, Int_t &NpureMC, Int_t &labelM, Int_t &labelGM);
        Bool_t BeautyCheck(Int_t &pidM, Int_t &pidGM);
        Bool_t EnhanceBeautyCheck(Int_t &pdgM, Int_t &pdgGM, Int_t &NpureMC, Int_t &labelM, Int_t &labelGM);
        Bool_t PionWeighingCheck(Int_t &pdgM, Int_t &pdgGM, Int_t &labelM, Int_t &labelGM , Int_t &labelGGM, Int_t &NpureMC);
        Bool_t EtaWeighingCheck(Int_t &pdgM, Int_t &pdgGM, Int_t &labelM, Int_t &labelGM , Int_t &labelGGM, Int_t &NpureMC);
        Int_t GetNpureMCParticle(AliAODMCHeader* fMCheader);

        private:
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        AliVEvent   *fVevent;  //!event object
        AliEventCuts *fEventCuts;
        AliPIDResponse *fPIDresponse;

        TString fPeriodName;

        TClonesArray  *fTracks_tender;//Tender tracks
        TClonesArray  *fCaloClusters_tender;//Tender cluster

        AliMultSelection *fMultSelection;

        AliAODMCParticle  *fMCparticle;
        AliAODMCHeader *fMCheader;
        TClonesArray  *fMCarray;//! MC array

        Bool_t fFlagMC;
        Bool_t fFlagClsTypeEMCal;
        Bool_t fFlagClsTypeDCal;
        Bool_t fFlagEMCalCorrection;
        Bool_t fFlagEG1;
        Bool_t fFlagEG2;
        Bool_t fFlagDG1;
        Bool_t fFlagDG2;

        TH1D *fEffHadron;

        //#########################//
        //Systematic uncertainties //
        //#########################//
        Double_t TrackEtaLow, TrackEtaHigh;
        Double_t NTPCClust, NITSClust, NCrossedRow;
        Double_t DCAxy, DCAz;
        Double_t TrackMatchPhi, TrackMatchEta;
        Double_t NsigmaLow, NsigmaHigh;
        Double_t M20Low, M20High;
        Double_t EopLow,EopHigh;
        Double_t AssoMinpT;
        Int_t AssoNTPCClust;
        Double_t MassCut;

        //################//
        // Primary Vertex //
        //################//
        TH1F *fNevents;
        TH1F *fvtxZ;
        TH1F *fvtxZ_NoCut;
        TH1F *fvtxZ_NcontCut;
        TH2F *fNcont;
        TH2F *fVertexCorre;

        //############//
        // Centrality //
        //############//
        Double_t fCentrality;
        Double_t fMultiplicity;
        Double_t fCentMax;
        Double_t fCentMin;
        TH1F *fCent;
        TH1F *fMulti;
        TH2F *fCentMultiCorrl;

        //##############//
        // Cluster info //
        //##############//
        TH1F *fCaloClusterE;
        TH2F *fCaloClustEtaphi;
        TH2F *fCaloClustTiming;

        //##############//
        // Track info //
        //##############//
        TH2F *fDCAxy;
        TH2F *fDCAz;
        TH2F *fTrackphieta;
        TH1F *fTrackPt;
        TH2F *fTrackCluster_woCut;
        TH2F *fTrackCluster;
        TH2F *fTrackChi2;

        //########################//
        // Mathed track & Cluster //
        //########################//
        TH2F *fCaloTrackDiff;
        TH1F *fCaloClusterEAfterMatch;
        TH1F *fCaloClusterEincE;
        TH2F *fCaloClustEtaPhiAfterMatch;
        TH1F *fTrackPtAfterMatch;
        TH2F *fTrackphietaAfterMatch;
        TH2F *fM02;
        TH2F *fM20;
        TH2F *fdEdx;
        TH2F *fNSigmaTPC;
        TH2F *fNsigmaEta;
        TH2F *fNSigmaTPCelectron;
        TH2F *fNSigmaTPChadron;
        TH2F *fNsigmaEtaElectron;
        TH2F *fEop;
        TH2F *fEopElectron;
        TH2F *fEopHadron;
        TH2F *fNSigmaEop;
        TH2F *fElectronphieta;
        TH2F *fElectronEpT;
        //---PHE search by Invariant-mass ---//
        TH2F *fInvmassLS;
        TH2F *fInvmassULS;
        TH1F *fPHEpTLS;
        TH1F *fPHEpTULS;

        //##########################//
        // two particle correlation //
        //##########################//
        //--- for event mixing ---//
        AliEventPoolManager *fPoolMgr;
        TH2F *fMixStatCent;
        TH2F *fMixStatVtxZ;
        TH2F *fMixStatCentVtxZCorrl;
        //--- particles correlation ---//
        TH1F *fTrigpTAllHadHCorrl;
        THnSparse *fSprsAllHadHCorrl;
        THnSparse *fSprsMixAllHadHCorrl;
        TH1F *fTrigpTAllHadHCorrl_CaloMatch;
        THnSparse *fSprsAllHadHCorrl_CaloMatch;
        THnSparse *fSprsMixAllHadHCorrl_CaloMatch;
        TH1F *fTrigpTHadHCorrl;
        THnSparse *fSprsHadHCorrl;
        THnSparse *fSprsMixHadHCorrl;
        TH1F *fTrigpTInclusiveEHCorrl;
        THnSparse *fSprsInclusiveEHCorrl;
        THnSparse *fSprsMixInclusiveEHCorrl;
        TH1F *fTrigpTTagULSEHCorrl;
        THnSparse *fSprsTagULSEHCorrl;
        THnSparse *fSprsMixTagULSEHCorrl;
        THnSparse *fSprsTagULSEHCorrl_Noparter;

        //#########//
        // MC info //
        //#########//
        //--- particle information --//
        TH1F *fMCPDGcodeM;
        TH1F *fMCPDGcodeGM;
        //--- hadron efficiency ---//
        TH2F *fMCchtrack_gene;
        TH2F *fMCchtrack_reco;
        //--- effeciency correction ---//
        TH1F *fMCTrackPtelectron;
        TH1F *fMCTrackPtelectron_reco;
        TH1F *fMCTrackPtelectron_wTPCPID;
        TH1F *fMCTrackPtelectron_wmatch;
        TH1F *fMCTrackPtelectron_wEMCALPID;
        TH1F *fMCTrackPtHFE;
        TH1F *fMCTrackPtHFE_reco;
        TH1F *fMCTrackPtHFE_wTPCPID;
        TH1F *fMCTrackPtHFE_wmatch;
        TH1F *fMCTrackPtHFE_wEMCALPID;
        TH2F *fMCHFEResponceMatrix;
        // --- PID check --- //
        TH2F *fMCTPCNSigmaelectron;
        TH2F *fMCNsigmaEtaElectron;
        TH2F *fMCHFEEop;
        TH2F *fMCHFEEopwPID;
        //--- DCA for c/b separation ---//
        TH2F *fMCDCAinclusive;
        TH2F *fMCDCAconv;
        TH2F *fMCDCAdalitz;
        TH2F *fMCDCAcharm;
        TH2F *fMCDCAbeauty;
        TH2F *fMCDCAPHE;
        TH2F *fMCDCAPHEULS;
        TH2F *fMCDCAPHELS;
        //--- Non-HFE tagging efficiency ---//
        TF1 *PionWeight;
        TF1 *EtaWeight;
        TH2F *fMCPionInput;
        TH2F *fMCEtaInput;
        TH1F *fMCTrackPtPHEHijing;
        TH1F *fMCTrackPtPHEnoweighting;
        TH1F *fMCTrackPtPHEaftweighting;
        TH1F *fMCTrackPtPHEHijingwPID;
        TH1F *fMCTrackPtPHEHijingwMasscut;
        TH1F *fMCTrackPtPHEwPIDnoweighting;
        TH1F *fMCTrackPtPHEwPIDaftweighting;
        TH1F *fMCTrackPtPHEwMasscutnoweighting;
        TH1F *fMCTrackPtPHEwMasscutaftweighting;
        // --- D meson & B meson --- //
        TH2F *fMCHFhadronpT;
        // --- c->e & b->e --- //
        TH1F *fMCDdecayE;
        TH1F *fMCBdecayE;



        AliAnalysisTaskCaloHFEpPbRun2(const AliAnalysisTaskCaloHFEpPbRun2&); // not implemented
        AliAnalysisTaskCaloHFEpPbRun2& operator=(const AliAnalysisTaskCaloHFEpPbRun2&); // not implemented

        ClassDef(AliAnalysisTaskCaloHFEpPbRun2, 1);
};


class AliehDPhiBasicParticlepPbRun2 : public AliVParticle
{
  public:
    AliehDPhiBasicParticlepPbRun2(Float_t eta, Float_t phi, Float_t pt, Short_t charge)
      : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge)
    {
    }
    ~AliehDPhiBasicParticlepPbRun2() {}

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
    virtual Int_t   GetLabel()    const { AliFatal("Not implemented"); return 0; }
    // PID
    virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }
    virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }

  private:
    Float_t fEta;      // eta
    Float_t fPhi;      // phi
    Float_t fpT;       // pT
    Short_t fCharge;   // charge

    ClassDef( AliehDPhiBasicParticlepPbRun2, 1); // class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};

#endif
