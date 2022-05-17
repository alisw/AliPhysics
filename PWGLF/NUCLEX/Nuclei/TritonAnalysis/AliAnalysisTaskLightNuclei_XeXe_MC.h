#ifndef   AliAnalysisTaskLightNuclei_XeXe_MC_cxx
#define   AliAnalysisTaskLightNuclei_XeXe_MC_cxx
#include "AliMCEventHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "TClonesArray.h"
#include "AliMCEvent.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "AliAODVertex.h"
#include "AliEventCuts.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"


class AliAnalysisTaskLightNuclei_XeXe_MC : public AliAnalysisTaskSE {

public:
    AliAnalysisTaskLightNuclei_XeXe_MC();
    AliAnalysisTaskLightNuclei_XeXe_MC(const char *name);
    virtual ~AliAnalysisTaskLightNuclei_XeXe_MC();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec (Option_t *option);

    Bool_t      GetInputEvent ();
    Bool_t      PassedMinimalTrackQualityCuts 		(AliAODTrack *track);
    Bool_t      PassedCandidateSelectionHe3    		(AliAODTrack *track);
    Bool_t      PassedCandidateSelectionTriton 		(AliAODTrack *track);
    Bool_t      PassedCandidateSelectionHelium3	    (AliAODTrack *track);
    Double_t    GetDCAxy                      		(AliAODTrack *track);
    Double_t    GetDCAz              		        (AliAODTrack *track);

    Bool_t      PassedTrackQualityCutsNoDCA         (AliAODTrack *track);
    Bool_t      PassedTrackQualityCuts              (AliAODTrack *track);
    Bool_t      IsHe3Candidate                      (AliAODTrack *track);
    Bool_t      IsCleanHe3Candidate                 (AliAODTrack *track);
    Double_t    Centered_nsigmaTPC                  (AliAODTrack *track);
    Double_t    Centered_nsigmaTOF                  (AliAODTrack *track);
    Bool_t      PassedTOFSelection                  (AliAODTrack *track);
    Bool_t      PassedTPCSelection                  (AliAODTrack *track);





    virtual void   Terminate(Option_t *);


    void SetCentralityRange (Double_t CentralityMin, Double_t CentralityMax)  {

        fCentralityMin = CentralityMin;
        fCentralityMax = CentralityMax;

    }

    void SetEventCuts (Double_t VertexZmin, Double_t VertexZmax, Int_t NumberVertexContributorsMin, const char *CentralityEstimator)  {

        fVertexZmin                     = VertexZmin;
        fVertexZmax                     = VertexZmax;
        fNumberVertexContributorsMin    = NumberVertexContributorsMin;
        fCentralityEstimator            = CentralityEstimator;
    }


    void SetTrackCuts (Double_t PtMin, Double_t PtMax, Double_t EtaMax, Double_t YMax, Int_t NumberClustersITSMin,
                       Int_t NumberClustersTPCMin, Int_t NumberCrossedRowsTPCMin,
                       Double_t CrossedRowsFindableClsMin, Int_t NumberClustersTPCdEdxMin, Double_t ChiSquarePerNDFMax, const char *ITSrequirement,
                       Double_t DCAzMax, Double_t DCAxyMax, Double_t nSigmaTOFmax, Double_t nSigmaTPCmax, Int_t TRDntracklets)  {

        fPtMin                          = PtMin;
        fPtMax                          = PtMax;
        fEtaMax                         = EtaMax;
        fYMax                           = YMax;
        fNumberClustersITSMin           = NumberClustersITSMin;
        fNumberClustersTPCMin           = NumberClustersTPCMin;
        fNumberCrossedRowsTPCMin        = NumberCrossedRowsTPCMin;
        fCrossedRowsFindableClsMin      = CrossedRowsFindableClsMin;
        fNumberClustersTPCdEdxMin       = NumberClustersTPCdEdxMin;
        fChiSquarePerNDFMax             = ChiSquarePerNDFMax;
        fITSrequirement                 = ITSrequirement;
        fDCAzMax                        = DCAzMax;
        fDCAxyMax                       = DCAxyMax;
        fnSigmaTOFmax                   = nSigmaTOFmax;
        fnSigmaTPCmax                   = nSigmaTPCmax;
        fTRDntracklets                  = TRDntracklets;
          }

    void RecalibrationTPCandTOF (Double_t par0_mean_TPC, Double_t par1_mean_TPC, Double_t par0_sigma_TPC, Double_t par0_mean_TOF, Double_t par1_mean_TOF, Double_t par0_sigma_TOF, Double_t par1_sigma_TOF )  {

        fpar0_mean_TPC                  = par0_mean_TPC;
        fpar1_mean_TPC                  = par1_mean_TPC;
        fpar0_sigma_TPC                 = par0_sigma_TPC;
        fpar0_mean_TOF                  = par0_mean_TOF;
        fpar1_mean_TOF                  = par1_mean_TOF;
        fpar0_sigma_TOF                 = par0_sigma_TOF;
        fpar1_sigma_TOF                 = par1_sigma_TOF;



    }


private:
    AliAODEvent       *fAODevent;//!
    AliMCEvent        *fMCevent;//!
    AliPIDResponse    *fPIDResponse;//!
    AliEventCuts       fAODeventCuts;//
    AliAnalysisUtils  *fUtils;//!
    TList             *fOutputList;//!
    TList             *fQAList;//!
    AliAODMCHeader    *fAODMCHeader;//!
    TClonesArray      *fAODArrayMCParticles;//!
    AliMCEventHandler *fMCEventHandler;//!
    TTree             *reducedTree_gen;//!
    TTree             *reducedTree_rec;//!


    Double_t    fCentralityMin;//
    Double_t    fCentralityMax;//
    Double_t    fVertexZmin;//
    Double_t    fVertexZmax;//
    Int_t       fNumberVertexContributorsMin;//
    const char *fCentralityEstimator;//
    Double_t    fPtMin;//
    Double_t    fPtMax;//
    Double_t    fEtaMax;//
    Double_t    fYMax;//
    Int_t       fNumberClustersITSMin;//
    Int_t       fNumberClustersTPCMin;//
    Int_t       fNumberCrossedRowsTPCMin;//
    Double_t    fCrossedRowsFindableClsMin;//
    Int_t       fNumberClustersTPCdEdxMin;//
    Double_t    fChiSquarePerNDFMax;//
    const char *fITSrequirement;//
    Double_t    fDCAzMax;//
    Double_t    fDCAxyMax;//
    Double_t    fnSigmaTOFmax;//
    Double_t    fnSigmaTPCmax;//
    Int_t       fTRDntracklets;//
    Double_t    fpar0_mean_TPC;//
    Double_t    fpar1_mean_TPC;//
    Double_t    fpar0_sigma_TPC;//
    Double_t    fpar0_mean_TOF;//
    Double_t    fpar1_mean_TOF;//
    Double_t    fpar0_sigma_TOF;//
    Double_t    fpar1_sigma_TOF;//

    //Histograms for He3
    TH1F        *histoNumberOfEvents;//!
    TH2F        *histoNsigmaTPCHe3_vs_pt;//!
    TH2F        *histoNsigmaTOFHe3_vs_pt;//!
    TH2F        *histoNsigmaTPCantiHe3_vs_pt;//!
    TH2F        *histoNsigmaTOFantiHe3_vs_pt;//!
    TH2F        *histoNsigmaTPCHe3_vs_pt_centered;//!
    TH2F        *histoNsigmaTPCantiHe3_vs_pt_centered;//!
    TH2F        *histoNsigmaTOFHe3_vs_pt_centered;//!
    TH2F        *histoNsigmaTOFantiHe3_vs_pt_centered;//!
    TH2F        *histoDCAxyHe3_vs_pt;//!
    TH2F        *histoDCAxyAntiHe3_vs_pt;//!
    TH2F        *histoNsigmaTOFHe3_vs_pt_trd;//!
    TH2F        *histoNsigmaTOFantiHe3_vs_pt_trd;//!
    TH2F        *histoNsigmaTPCHe3_vs_p;//!
    TH2F        *histoNsigmaTPCantiHe3_vs_p;//!
    TH2F        *histoNsigmaTOFHe3_vs_p;//!
    TH2F        *histoNsigmaTOFantiHe3_vs_p;//!
    TH2F        *histoNsigmaTPCHe3_vs_p_notof;//!
    TH2F        *histoNsigmaTPCantiHe3_vs_p_notof;//!


    // Tree variables
    Float_t centrality;//
    Float_t multPercentile_V0M;//
    Int_t   particleType;
    Float_t xVertex;//
    Float_t yVertex;//
    Float_t zVertex;//
    Int_t   ID_event;//
    Int_t   Helium3Sec;//
    Int_t   TritonSec;//
    Int_t   lp;//
    Float_t px;//
    Float_t py;//
    Float_t pz;//
    Float_t pt;//
    Float_t pt_particle;//
    Float_t pz_particle;//
    Float_t deltapt;//
    Float_t deltapz;//
    Int_t   charge;//
    Double_t dcaxy;//
    Double_t dcaz;//
    Int_t   trackType;//
    Int_t   nTPC_Clusters;//
    Int_t   nITS_Clusters;//
    Int_t   nTPC_Clusters_dEdx;//
    Int_t   nTPC_FindableClusters;//
    Int_t   nTPC_CrossedRows;//
    Int_t   HasPointOnITSLayer0;//
    Int_t   HasPointOnITSLayer1;//
    Int_t   HasPointOnITSLayer2;//
    Int_t   HasPointOnITSLayer3;//
    Int_t   HasPointOnITSLayer4;//
    Int_t   HasPointOnITSLayer5;//
    Float_t chi2_NDF;//
    Float_t chi2_ITS;//
    Float_t nSigmaITS_He3;
    Float_t nSigmaTPC_He3;//
    Float_t nSigmaTOF_He3;//
    Float_t nSigmaITS_Triton;//
    Float_t nSigmaTPC_Triton;//
    Float_t nSigmaTOF_Triton;//
    Int_t   PrimaryParticle;//
    Int_t   SecondaryMaterial;//
    Int_t   SecondaryDecay;//
    Int_t   PrimaryTrack;//
    Double_t TPC_signal;//
    Double_t ITS_signal;//
    Int_t   TRDntracklets;//
    Int_t   hasTOFhit;//



    AliAnalysisTaskLightNuclei_XeXe_MC(const AliAnalysisTaskLightNuclei_XeXe_MC&);
    AliAnalysisTaskLightNuclei_XeXe_MC& operator=(const AliAnalysisTaskLightNuclei_XeXe_MC&);
    ClassDef(AliAnalysisTaskLightNuclei_XeXe_MC, 1);
};
#endif
