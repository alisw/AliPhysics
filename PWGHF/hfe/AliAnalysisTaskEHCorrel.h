#ifndef AliAnalysisTaskEHCorrel_h
#define AliAnalysisTaskEHCorrel_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//Task for Heavy Flavour Electron-Hadron DeltaPhi Correlation in Run 2//
//                                                                    //
//  Author: Deepa Thomas (University of Texas at Austin)              //
//                                                                    //
////////////////////////////////////////////////////////////////////////

class THnSparse;
class TH2F;
class TLorentzVector;

class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliESDtrack;
class AliAODTrack;
class AliCFManager;
class AliEventPoolManager;
class AliMultSelection;

#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
#include "AliSelectNonHFE.h"

class AliAnalysisTaskEHCorrel : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskEHCorrel();
    AliAnalysisTaskEHCorrel(const char *name);
    virtual ~AliAnalysisTaskEHCorrel();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    Bool_t  PassEventSelect(AliVEvent *fVevent);
    void SetEMCalTriggerEG1(Bool_t flagTr1) { fEMCEG1=flagTr1; fEMCEG2=kFALSE;};
    void SetEMCalTriggerEG2(Bool_t flagTr2) { fEMCEG2=flagTr2; fEMCEG1=kFALSE;};
    void    CheckCentrality(AliAODEvent* fAOD, Bool_t &centralitypass);
    Bool_t  PassTrackCuts(AliAODTrack *atrack);
    void    GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff);
    Bool_t  PassEIDCuts(AliVTrack *track, AliVCluster *clust, Bool_t &Hadtrack);
    Bool_t  PassHadronCuts(AliAODTrack *HadTrack);
    void    HadronInfo(Int_t itrack);
    void    ElectronHadCorrel(Int_t itrack, AliVTrack *track, THnSparse *SparseEHCorrl);
    void    EMCalClusterInfo();
    void    SelectNonHFElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec);
    TObjArray* CloneAndReduceTrackList();
    void    MixedEvent(AliVTrack *track, THnSparse *SparseMixEHCorrl);

    void    SetCentralitySelection(Double_t centMin, Double_t centMax) {fCentralityMin = centMin; fCentralityMax = centMax;};
    void    SetTPCnsigCut(Double_t nsigMin, Double_t nsigMax) {fTPCnSigmaMin = nsigMin; fTPCnSigmaMax= nsigMax;};
    void    SetM02Cut(Double_t m02Min, Double_t m02Max) {fM02Min = m02Min; fM02Max = m02Max;};
    void    SetM20Cut(Double_t m20Min, Double_t m20Max) {fM20Min = m20Min; fM20Max = m20Max;};
    void    SetEovPCut(Double_t eovpMin, Double_t eovpMax) {fEovPMin = eovpMin; fEovPMax = eovpMax;};
    void    SetHadronCutCase(Int_t hadCutCase) {fHadCutCase = hadCutCase;};
    void    SetTriggerElePtCut(Bool_t trigElePtcut) {fTrigElePtCut = trigElePtcut;};


  private:
    AliVEvent 		    *fVevent;//!V event object
    AliAODEvent 		*fAOD;//!AOD object
    const AliVVertex    *fpVtx; //!
    AliPIDResponse      *fpidResponse; //!pid response
    AliMultSelection    *fMultSelection;//!

    Double_t            fCentrality;//!
    Double_t            fCentralityMin;//
    Double_t            fCentralityMax;//
    Double_t            fMultiplicity;//!
    Bool_t              fEMCEG1;//
    Bool_t              fEMCEG2;//
    Bool_t              fFlagClsTypeEMC;//switch to select EMC clusters
    Bool_t              fFlagClsTypeDCAL;//switch to select DCAL clusters
    Double_t            fTPCnSigma;//!
    Double_t            fTPCnSigmaMin;//
    Double_t            fTPCnSigmaMax;//
    Double_t            fM02Min;//
    Double_t            fM02Max;//
    Double_t            fM20Min;//
    Double_t            fM20Max;//
    Double_t            fEovPMin;//
    Double_t            fEovPMax;//
    Int_t               fTPCNClsHad;// Had track TPC NClusters
    Double_t            fInvmassCut;//
    Double_t            fTPCnSigmaHadMin;//
    Double_t            fTPCnSigmaHadMax;//
    Int_t               fHadCutCase;//
    AliEventPoolManager *fPoolMgr;//!
    Bool_t              fTrigElePtCut;//!

    TList       	   	*fOutputList;		//!output list
    TH1F                *fNevents;		//!no of events
    TH1F                *fVtxZ;//!
    TH1F                *fVtxX;//!
    TH1F                *fVtxY;//!
    TH1F                *fCentralityNoPass;//!
    TH1F                *fCentralityPass;//!
    TH1F                *fMultiplicityNoPass;//!
    TH1F                *fMultiplicityPass;//!
    TH1F                *fHistClustE;//!
    TH2F                *fEMCClsEtaPhi;//!

    TH1F                *fNegTrkIDPt;//!
    TH1F                *fTrkPt;//!
    TH1F                *fTrketa;//!
    TH1F                *fTrkphi;//!
    TH2F                *fdEdx;//!
    TH2F                *fTPCnsig;//!
    TH1F                *fHistPtMatch;//!
    TH2F                *fEMCTrkMatch;//!
    TH1F                *fEMCTrkPt;//!
    TH1F                *fEMCTrketa;//!
    TH1F                *fEMCTrkphi;//!
    TH2F                *fEMCTPCnsig;//!
    TH1F                *fClsEAftMatch;//!
    TH2F                *fClsEtaPhiAftMatch;//!
    TH2F                *fHistNsigEop;//!
    TH2F                *fM20EovP;//!
    TH2F                *fM02EovP;//!
    TH2F                *fHistEop;//!
    TH2F                *fM20;//!
    TH2F                *fM02;//!
    TH2F                *fHistEop_AftEID;//!
    TH1F                *fInclsElecPt;//!
    TH2F                *fHadEop;//!
    TH1F                *fHadPt_AftEID;//!

    TH2F                *fHadronPhiPt;//!
    TH1F                *fHadronPhi;//!
    TH1F                *fHadronPhiTPChalf;//!
    TH1F                *fHadronPt;//!

    TH1F                *fInvmassLS;//!
    TH1F                *fInvmassULS;//!
    TH2F                *fInvmassLSPt;//!
    TH2F                *fInvmassULSPt;//!

    TH1F                *fNoMixedEvents;//!
    TH2F                *fMixStatCent;//!
    TH2F                *fMixStatVtxZ;//!

    THnSparse           *fSprsHadHCorrl;//!
    THnSparse           *fSprsInclusiveEHCorrl;//!
    THnSparse           *fSprsLSEHCorrl;//!
    THnSparse           *fSprsULSEHCorrl;//!
    THnSparse           *fSprsMixInclusiveEHCorrl;//!
    THnSparse           *fSprsMixLSEHCorrl;//!
    THnSparse           *fSprsMixULSEHCorrl;//!


    AliAnalysisTaskEHCorrel(const AliAnalysisTaskEHCorrel&); // not implemented
    AliAnalysisTaskEHCorrel& operator=(const AliAnalysisTaskEHCorrel&); // not implemented

    ClassDef(AliAnalysisTaskEHCorrel, 2); //!example of analysis
};

class AliehDPhiBasicParticle : public AliVParticle
{
  public:
    AliehDPhiBasicParticle(Float_t eta, Float_t phi, Float_t pt, Short_t charge)
      : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge)
    {
    }
    ~AliehDPhiBasicParticle() {}

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

    ClassDef( AliehDPhiBasicParticle, 1); // class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};
#endif
