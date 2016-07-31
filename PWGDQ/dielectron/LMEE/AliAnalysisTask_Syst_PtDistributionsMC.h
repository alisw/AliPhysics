#ifndef AliAnalysisTask_Syst_PtDistributionsMC_cxx
#define AliAnalysisTask_Syst_PtDistributionsMC_cxx


//========================== PT DISTRIBUTIONS (MC) ===========================//
//                                                                            //
//    Pt distributions of electrons & positrons for different cut settings    //
//    before & after pre-filter cuts.                                         //
//                                                                            //
//============================================================================//


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


class AliAnalysisTask_Syst_PtDistributionsMC : public AliAnalysisTaskSE {

public:
    AliAnalysisTask_Syst_PtDistributionsMC();
    AliAnalysisTask_Syst_PtDistributionsMC (const char *name);
    virtual ~AliAnalysisTask_Syst_PtDistributionsMC();
    
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
    
    void SetTrackCuts (Int_t ITS_minNcls, Int_t TPC_minNcls, Int_t TPC_nClsdEdx, Int_t TPC_minCr, Double_t MinCrOverFindableCls, Double_t MaxGoldenChi2, Double_t MaxTPCchi2, Double_t MaxITSchi2, Double_t MaxFracSharedCls, const char *ITSreq)  {
        
        fITS_minNcls =          ITS_minNcls;
        fTPC_minNcls =          TPC_minNcls;
        fTPC_nClsdEdx =         TPC_nClsdEdx;
        fTPC_minCr =            TPC_minCr;
        fMinCrOverFindableCls = MinCrOverFindableCls;
        fMaxGoldenChi2 =        MaxGoldenChi2;
        fMaxTPCchi2 =           MaxTPCchi2;
        fMaxITSchi2 =           MaxITSchi2;
        fMaxFracSharedCls =     MaxFracSharedCls;
        fITSreq =               ITSreq;
    }
    
    void SetPIDCuts (Double_t nsigmaTOF_max, Double_t nsigmaITS_max, Double_t nsigmaTPC_min, Double_t nsigmaTPC_max)  {
        
        fnsigmaTOF_max = nsigmaTOF_max;
        fnsigmaITS_max = nsigmaITS_max;
        fnsigmaTPC_min = nsigmaTPC_min;
        fnsigmaTPC_max = nsigmaTPC_max;
    }
    
    void SetPrefilterCuts (Double_t MassMin, Double_t MassMax, Double_t PhivLim)  {
        
        fMassMin = MassMin;
        fMassMax = MassMax;
        fPhivLim = PhivLim;
    }
    
    void SetKinematicCuts (Double_t PtMin, Double_t PtMax, Double_t EtaLim)  {
        
        fPtMin  = PtMin;
        fPtMax  = PtMax;
        fEtaLim = EtaLim;
    }
    
    void GetCentralityBins (TH1F *HistoCentralityBins)  { fHistoCentralityBins = HistoCentralityBins; }

    void GetWeightsPtDistributions (TH2F *HistoPizero, TH2F *HistoEta, TH2F *HistoEtaPrime, TH2F *HistoRho, TH2F *HistoOmega, TH2F *HistoPhi)  {
        
        fHisto_Hijing_PizeroWeight =   HistoPizero;
        fHisto_Hijing_EtaWeight =      HistoEta;
        fHisto_Hijing_EtaPrimeWeight = HistoEtaPrime;
        fHisto_Hijing_RhoWeight =      HistoRho;
        fHisto_Hijing_OmegaWeight =    HistoOmega;
        fHisto_Hijing_PhiWeight =      HistoPhi;
    }
    
    
    virtual void   UserCreateOutputObjects ();
    virtual void   UserExec (Option_t *option);
    
    Bool_t   GetEvent ();
    Bool_t   IsTrackFromHijing           (AliESDtrack *track);
    Bool_t   PassedTrackQualityCuts      (AliESDtrack *track);
    Bool_t   PassedLooseTrackQualityCuts (AliESDtrack *track);
    Double_t FractionSharedClsITS        (AliESDtrack *track);
    Bool_t   PassedPIDcuts               (AliESDtrack *track);
    Double_t GetMass                     (AliESDtrack *track1, AliESDtrack *track2);
    Double_t GetPhiV                     (AliESDtrack *track1, AliESDtrack *track2);
    Double_t Weight                      (TParticle *particle);
    Double_t WeightConversions           (TParticle *parent);

    virtual void   Terminate(Option_t *);
    
private:
    TList           *fOutputList;//!
    AliESDEvent     *fESDevent;//!
    AliMCEvent      *fMCEvent;//!
    AliStack        *fStack;//!
    AliPIDResponse  *fPIDResponse;//!
    AliESDtrackCuts *fESDTrackCuts_Std;//!
    AliESDtrackCuts *fESDTrackCuts_Loose;//!
    
    Double_t    fCentralityMin;//
    Double_t    fCentralityMax;//
    Double_t    fDCAxy_param0;//
    Double_t    fDCAxy_param1;//
    Double_t    fDCAxy_param2;//
    Double_t    fDCAz_max;//
    Int_t       fITS_minNcls;//
    Int_t       fTPC_minNcls;//
    Int_t       fTPC_nClsdEdx;//
    Int_t       fTPC_minCr;//
    Double_t    fMinCrOverFindableCls;//
    Double_t    fMaxGoldenChi2;//
    Double_t    fMaxTPCchi2;//
    Double_t    fMaxITSchi2;//
    Double_t    fMaxFracSharedCls;//
    const char *fITSreq;//
    Double_t    fnsigmaTOF_max;//
    Double_t    fnsigmaITS_max;//
    Double_t    fnsigmaTPC_min;//
    Double_t    fnsigmaTPC_max;//
    Double_t    fMassMin;//
    Double_t    fMassMax;//
    Double_t    fPhivLim;//
    Double_t    fPtMin;//
    Double_t    fPtMax;//
    Double_t    fEtaLim;//
    
    
    //Weights Pt Distributions
    TH2F *fHisto_Hijing_PizeroWeight;//
    TH2F *fHisto_Hijing_EtaWeight;//
    TH2F *fHisto_Hijing_EtaPrimeWeight;//
    TH2F *fHisto_Hijing_RhoWeight;//
    TH2F *fHisto_Hijing_OmegaWeight;//
    TH2F *fHisto_Hijing_PhiWeight;//
    
    //Centrality Bins
    TH1F *fHistoCentralityBins;//
    
    
    //Statistics
    TH1F *fHistoEvents;//!

    
    //Pt distributions
    TH1F *fHistoPtDistribution_Electrons;//!
    TH1F *fHistoPtDistribution_Positrons;//!
    TH1F *fHistoPtDistribution_Electrons_Pref;//!
    TH1F *fHistoPtDistribution_Positrons_Pref;//!

    
    AliAnalysisTask_Syst_PtDistributionsMC(const AliAnalysisTask_Syst_PtDistributionsMC&);
    AliAnalysisTask_Syst_PtDistributionsMC& operator=(const AliAnalysisTask_Syst_PtDistributionsMC&);
    
    ClassDef(AliAnalysisTask_Syst_PtDistributionsMC, 1);
};

#endif
