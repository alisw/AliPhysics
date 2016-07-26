#ifndef AliAnalysisTaskDielectronsPbPb_Efficiency_cxx
#define AliAnalysisTaskDielectronsPbPb_Efficiency_cxx


//========================== DIELECTRON ANALYSIS (EFFICIENCY) ==========================//
//                                                                                      //
//    Pair efficiency to correct dielectron spectra in Pb-Pb collisions at 2.76 TeV.    //
//    Additional factor to account for signal losses due to the pre-filtering.          //
//                                                                                      //
//======================================================================================//


#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTask.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "TParticle.h"
#include "TVector3.h"
#include "AliStack.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"


class AliAnalysisTaskDielectronsPbPb_Efficiency : public AliAnalysisTaskSE {

public:
    AliAnalysisTaskDielectronsPbPb_Efficiency();
    AliAnalysisTaskDielectronsPbPb_Efficiency (const char *name);
    virtual ~AliAnalysisTaskDielectronsPbPb_Efficiency();
    
    void SetCentralityRange (Double_t CentralityMin, Double_t CentralityMax)  {
        
        fCentralityMin = CentralityMin;
        fCentralityMax = CentralityMax;
    }
    
    void SetDCAparameters (Double_t DCAxy_param0, Double_t DCAxy_param1, Double_t DCAxy_param2, Double_t DCAz_max)  {
        
        fDCAxy_param0 = DCAxy_param0;
        fDCAxy_param1 = DCAxy_param1;
        fDCAxy_param2 = DCAxy_param2;
        fDCAz_max =     DCAz_max;
    }
    
    void SetTrackCuts (Int_t ITS_minNcls, Int_t TPC_minNcls, Int_t TPC_nClsdEdx, Int_t TPC_minCr, Double_t MinCrOverFindableCls, Double_t MaxGoldenChi2, Double_t MaxTPCchi2, Double_t MaxITSchi2, Double_t MaxFracSharedCls)  {
        
        fITS_minNcls =          ITS_minNcls;
        fTPC_minNcls =          TPC_minNcls;
        fTPC_nClsdEdx =         TPC_nClsdEdx;
        fTPC_minCr =            TPC_minCr;
        fMinCrOverFindableCls = MinCrOverFindableCls;
        fMaxGoldenChi2 =        MaxGoldenChi2;
        fMaxTPCchi2 =           MaxTPCchi2;
        fMaxITSchi2 =           MaxITSchi2;
        fMaxFracSharedCls =     MaxFracSharedCls;
    }
    
    void SetPIDCuts (Double_t nsigmaTOF_max, Double_t nsigmaITS_max, Double_t nsigmaTPC_min, Double_t nsigmaTPC_max)  {
        
        fnsigmaTOF_max = nsigmaTOF_max;
        fnsigmaITS_max = nsigmaITS_max;
        fnsigmaTPC_min = nsigmaTPC_min;
        fnsigmaTPC_max = nsigmaTPC_max;
    }
    
    void SetPrefilterCuts (Double_t MassLim, Double_t PhivLim)  {
        
        fMassLim = MassLim;
        fPhivLim = PhivLim;
    }

    void GetDetectorResponseMatrices (TH2F *HistoMomentum, TH2F *HistoTheta, TH2F *HistoPhi_Elec, TH2F *HistoPhi_Pos)  {
        
        fHistoDetResponseMatrix_Momentum =      HistoMomentum;
        fHistoDetResponseMatrix_Theta =         HistoTheta;
        fHistoDetResponseMatrix_Phi_Electrons = HistoPhi_Elec;
        fHistoDetResponseMatrix_Phi_Positrons = HistoPhi_Pos;
    }
    
    void GetWeightsPtDistributions (TH2F *HistoPizero, TH2F *HistoEta, TH2F *HistoEtaPrime, TH2F *HistoRho, TH2F *HistoOmega, TH2F *HistoPhi)  {
        
        fHisto_Hijing_PizeroWeight =   HistoPizero;
        fHisto_Hijing_EtaWeight =      HistoEta;
        fHisto_Hijing_EtaPrimeWeight = HistoEtaPrime;
        fHisto_Hijing_RhoWeight =      HistoRho;
        fHisto_Hijing_OmegaWeight =    HistoOmega;
        fHisto_Hijing_PhiWeight =      HistoPhi;
    }
    
    void GetCentralityBins (TH1F *HistoCentralityBins)  { fHistoCentralityBins = HistoCentralityBins; }

    
    virtual void   UserCreateOutputObjects ();
    virtual void   UserExec (Option_t *option);
    
    Bool_t   GetEvent ();
    Bool_t   IsPrimaryElectron           (TParticle *particle);
    Bool_t   IsLightFlavorParticle       (TParticle *parent);
    Bool_t   IsHeavyFlavorParticle       (TParticle *parent);
    Bool_t   IsCorrelatedPair            (TParticle *particle1,TParticle *particle2);
    TVector3 GetReconstructedMomentum    (Double_t p,Double_t theta,Double_t phi,Short_t q);
    Bool_t   IsTrackFromHijing           (AliESDtrack *track);
    Bool_t   PassedTrackQualityCuts      (AliESDtrack *track);
    Bool_t   PassedLooseTrackQualityCuts (AliESDtrack* track);
    Double_t FractionSharedClsITS        (AliESDtrack *track);
    Bool_t   PassedPIDCuts               (AliESDtrack *track);
    Bool_t   IsFromConversion            (TParticle *particle);
    Double_t GetMass                     (TVector3 P1,TVector3 P2);
    Double_t GetMass                     (AliESDtrack *track1,AliESDtrack *track2);
    Double_t GetPtee                     (AliESDtrack *track1,AliESDtrack *track2);
    Double_t GetPhiV                     (AliESDtrack *track1,AliESDtrack *track2 );
    Double_t Weight                      (TParticle *particle);
    Double_t WeightConversions           (TParticle *parent);

    virtual void   Terminate(Option_t *);
    
private:
    TList           *fOutputList;//!
    AliESDEvent     *fESDevent;//!
    AliMCEvent      *fMCEvent;//!
    AliStack        *fStack;//!
    AliPIDResponse  *fPIDResponse;//!
    AliESDtrackCuts *fESDTrackCuts;//!
    Double_t fCentralityMin;//
    Double_t fCentralityMax;//
    Double_t fDCAxy_param0;//
    Double_t fDCAxy_param1;//
    Double_t fDCAxy_param2;//
    Double_t fDCAz_max;//
    Int_t    fITS_minNcls;//
    Int_t    fTPC_minNcls;//
    Int_t    fTPC_nClsdEdx;//
    Int_t    fTPC_minCr;//
    Double_t fMinCrOverFindableCls;//
    Double_t fMaxGoldenChi2;//
    Double_t fMaxTPCchi2;//
    Double_t fMaxITSchi2;//
    Double_t fMaxFracSharedCls;//
    Double_t fnsigmaTOF_max;//
    Double_t fnsigmaITS_max;//
    Double_t fnsigmaTPC_min;//
    Double_t fnsigmaTPC_max;//
    Double_t fMassLim;//
    Double_t fPhivLim;//
    
    //Detector Response Matrices
    TH2F *fHistoDetResponseMatrix_Momentum;//
    TH2F *fHistoDetResponseMatrix_Theta;//
    TH2F *fHistoDetResponseMatrix_Phi_Electrons;//
    TH2F *fHistoDetResponseMatrix_Phi_Positrons;//
    
    //Weights Pt Distributions (Hijing To ALICE MEasurements)
    TH2F *fHisto_Hijing_PizeroWeight;//
    TH2F *fHisto_Hijing_EtaWeight;//
    TH2F *fHisto_Hijing_EtaPrimeWeight;//
    TH2F *fHisto_Hijing_RhoWeight;//
    TH2F *fHisto_Hijing_OmegaWeight;//
    TH2F *fHisto_Hijing_PhiWeight;//
    
    //Centrality Bins
    TH1F *fHistoCentralityBins;//
    
    
    //Statistics & Centrality
    TH1F *fHistoEvents;//!

    //Pair Efficiency
    TH2F *fHistoInvMass_Gen;//!
    TH2F *fHistoInvMass_Rec;//!
    TH2F *fHistoInvMass_Rec_Pref;//!
    TH2F *fHistoInvMass_Gen_noweights;//!
    TH2F *fHistoInvMass_Rec_noweights;//!
    TH2F *fHistoInvMass_Rec_Pref_noweights;//!
    
    //Conversions Residual Contribution
    TH2F *fHistoInvariantMass_Dielectrons;//!
    TH2F *fHistoInvariantMass_Conversions;//!
    TH2F *fHistoInvariantMass_Dielectrons_noweights;//!
    TH2F *fHistoInvariantMass_Conversions_noweights;//!
    
    
    AliAnalysisTaskDielectronsPbPb_Efficiency(const AliAnalysisTaskDielectronsPbPb_Efficiency&);
    AliAnalysisTaskDielectronsPbPb_Efficiency& operator=(const AliAnalysisTaskDielectronsPbPb_Efficiency&);
    
    ClassDef(AliAnalysisTaskDielectronsPbPb_Efficiency, 1);
};

#endif
