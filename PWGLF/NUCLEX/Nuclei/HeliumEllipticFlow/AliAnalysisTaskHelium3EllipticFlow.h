#ifndef AliAnalysisTaskHelium3EllipticFlow_cxx
#define AliAnalysisTaskHelium3EllipticFlow_cxx

//========================== HELIUM-3 ELLIPTIC FLOW ===========================//
//                                                                             //
//    Analysis of Helium-3 elliptic flow (v2) using the event plane method.    //
//                                                                             //
//=============================================================================//

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "AliEventplane.h"
#include "AliAODVertex.h"
#include "AliEventCuts.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "TObjArray.h"
#include "TVectorD.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"

class AliAnalysisTaskHelium3EllipticFlow : public AliAnalysisTaskSE {
    
public:
    AliAnalysisTaskHelium3EllipticFlow();
    AliAnalysisTaskHelium3EllipticFlow(const char *name);
    virtual ~AliAnalysisTaskHelium3EllipticFlow();
    
    void SetEventCuts (Double_t SignMagField, Double_t VertexZmin, Double_t VertexZmax, Int_t NumberVertexContributorsMin, const char *CentralityEstimator)  {
        
        fSignMagField                = SignMagField;
        fVertexZmin                  = VertexZmin;
        fVertexZmax                  = VertexZmax;
        fNumberVertexContributorsMin = NumberVertexContributorsMin;
        fCentralityEstimator         = CentralityEstimator;
    }
    
    
    void SetTrackCuts (UInt_t FilterBitMask, Double_t PtMin, Double_t PtMax, Double_t EtaMax, Int_t NumberClustersITSMin, Int_t NumberClustersTPCMin,
                       Double_t CrossedRowsFindableClsMin, Int_t NumberClustersTPCdEdxMin, const char *ITSrequirement,
                       Double_t DCAzMax, Double_t DCAxyMax, Double_t NumberSigmaTPCmax)  {
        
        fFilterBitMask             = FilterBitMask;
        fPtMin                     = PtMin;
        fPtMax                     = PtMax;
        fEtaMax                    = EtaMax;
        fNumberClustersITSMin      = NumberClustersITSMin;
        fNumberClustersTPCMin      = NumberClustersTPCMin;
        fCrossedRowsFindableClsMin = CrossedRowsFindableClsMin;
        fNumberClustersTPCdEdxMin  = NumberClustersTPCdEdxMin;
        fITSrequirement            = ITSrequirement;
        fDCAzMax                   = DCAzMax;
        fDCAxyMax                  = DCAxyMax;
        fNumberSigmaTPCmax         = NumberSigmaTPCmax;
    }
    
    
    void SetMomentumTransformation   (TH2F *HistoPtTransformationMatrix)               { fHistoPtTransformationMatrix = HistoPtTransformationMatrix; }
    void SetCentralityRange          (Double_t CentralityMin, Double_t CentralityMax)  { fCentralityMin = CentralityMin; fCentralityMax = CentralityMax; }
    void GetCentralityEvtPlaneWeight (TH2F *HistoCentralityEvtPlaneWeight)             { fHistoCentralityEvtPlaneWeight = HistoCentralityEvtPlaneWeight; }
    void GetEventPlaneResolution     (Double_t centralityEvtPlaneWeight);

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec (Option_t *option);
    
    Double_t GetEventCentrality ();
    Double_t GetEventPlaneAngle ();
    Bool_t   GetEvent                          (Double_t &centralityEvtPlaneWeight);
    Bool_t   PassedEvtPlaneResTrackQualityCuts (AliAODTrack *track);
    Bool_t   PassedTrackQualityCuts            (AliAODTrack *track);
    Bool_t   PassedHelium3IDCuts               (AliAODTrack *track);
    Double_t GetDCAxy                          (AliAODTrack *track);
    Double_t GetDCAz                           (AliAODTrack *track);
    Bool_t   IsTrackInPlane                    (AliAODTrack *track);
    Double_t GetCorrectedPt                    (AliAODTrack *track);

    virtual void   Terminate(Option_t *);
    
private:
    AliAODEvent      *fAODevent;//!
    AliPIDResponse   *fPIDResponse;//!
    AliEventCuts      fAODeventCuts;//
    AliAnalysisUtils *fUtils;//!
    TList            *fOutputList;//!
    TList            *fQAList;//!

    Double_t        fSignMagField;//
    Double_t        fVertexZmin;//
    Double_t        fVertexZmax;//
    Int_t           fNumberVertexContributorsMin;//
    const char     *fCentralityEstimator;//
    
    UInt_t          fFilterBitMask;//
    Double_t        fPtMin;//
    Double_t        fPtMax;//
    Double_t        fEtaMax;//
    Int_t           fNumberClustersITSMin;//
    Int_t           fNumberClustersTPCMin;//
    Double_t        fCrossedRowsFindableClsMin;//
    Int_t           fNumberClustersTPCdEdxMin;//
    const char     *fITSrequirement;//
    Double_t        fDCAzMax;//
    Double_t        fDCAxyMax;//
    Double_t        fNumberSigmaTPCmax;//
    
    TH2F           *fHistoPtTransformationMatrix;//
    Double_t        fCentralityMin;//
    Double_t        fCentralityMax;//
    TH2F           *fHistoCentralityEvtPlaneWeight;//
    
    
    
    
    //Number of Events
    TH1F *fHistoEvents;//!
    
    //Centrality Distribution
    TH1F *fHistoCentralityDistribution;//!

    //Centrality & Event Plane Angle
    TH2F *fHistoCentrality_EvtPlaneAngle_NoFlattening;//!
    TH2F *fHistoCentrality_EvtPlaneAngle_Flattening;//!
    
    //Event Plane Resolution Histograms
    TH2F *fHistoCosn_PsiA_PsiB_vs_Centrality;//!
    TH2F *fHistoCosn_PsiA_PsiC_vs_Centrality;//!
    TH2F *fHistoCosn_PsiB_PsiC_vs_Centrality;//!
    
    //He3 Centrality Distributions
    TH2F *fHistoHe3_CentralityDistribution_Pt;//!
    TH2F *fHistoAntiHe3_CentralityDistribution_Pt;//!
    
    //DCAxy Distributions of Pure He3 Candidates
    TH2F *fHistoHe3_DCAxy_vs_Pt_InPlane;//!
    TH2F *fHistoHe3_DCAxy_vs_Pt_OutPlane;//!
    TH2F *fHistoAntiHe3_DCAxy_vs_Pt_InPlane;//!
    TH2F *fHistoAntiHe3_DCAxy_vs_Pt_OutPlane;//!

    //In-Plane & Out-Of-Plane He3 Candidates
    TH2F *fHistoHe3_nsigmaTPC_vs_Pt_InPlane;//!
    TH2F *fHistoHe3_nsigmaTPC_vs_Pt_OutPlane;//!
    TH2F *fHistoAntiHe3_nsigmaTPC_vs_Pt_InPlane;//!
    TH2F *fHistoAntiHe3_nsigmaTPC_vs_Pt_OutPlane;//!
    
    TH2F *fHistoHe3_nsigmaTPC_vs_Pt_PerNucleon_InPlane;//!
    TH2F *fHistoHe3_nsigmaTPC_vs_Pt_PerNucleon_OutPlane;//!
    TH2F *fHistoAntiHe3_nsigmaTPC_vs_Pt_PerNucleon_InPlane;//!
    TH2F *fHistoAntiHe3_nsigmaTPC_vs_Pt_PerNucleon_OutPlane;//!
    
    //Pt Spectra
    TH1F *fHistoHe3_Pt;//!
    TH1F *fHistoAntiHe3_Pt;//!

    
    
    
    AliAnalysisTaskHelium3EllipticFlow(const AliAnalysisTaskHelium3EllipticFlow&);
    AliAnalysisTaskHelium3EllipticFlow& operator=(const AliAnalysisTaskHelium3EllipticFlow&);
    
    ClassDef(AliAnalysisTaskHelium3EllipticFlow, 1);
};
#endif
