
#ifndef AliAnalysisTaskTritonVsMultiplicity_XeXe_cxx
#define AliAnalysisTaskTritonVsMultiplicity_XeXe_cxx


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

class AliAnalysisTaskTritonVsMultiplicity_XeXe : public AliAnalysisTaskSE {

public:
    AliAnalysisTaskTritonVsMultiplicity_XeXe();
    AliAnalysisTaskTritonVsMultiplicity_XeXe(const char *name);
    virtual ~AliAnalysisTaskTritonVsMultiplicity_XeXe();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec (Option_t *option);

    void SetCentralityRange (Double_t CentralityMin, Double_t CentralityMax)  {

        fCentralityMin = CentralityMin;
        fCentralityMax = CentralityMax;

    }

    void SetEventCuts (Double_t VertexZmin, Double_t VertexZmax, Int_t NumberVertexContributorsMin, const char *CentralityEstimator)  {

        fVertexZmin                  = VertexZmin;
        fVertexZmax                  = VertexZmax;
        fNumberVertexContributorsMin = NumberVertexContributorsMin;
        fCentralityEstimator         = CentralityEstimator;
    }


    void SetTrackCuts (Double_t PtMin, Double_t PtMax, Double_t EtaMax, Double_t YMax, Int_t NumberClustersITSMin,
                       Int_t NumberClustersTPCMin, Int_t NumberCrossedRowsTPCMin,
                       Double_t CrossedRowsFindableClsMin, Int_t NumberClustersTPCdEdxMin, Double_t ChiSquarePerNDFMax, const char *ITSrequirement,
                       Double_t DCAzMax, Double_t DCAxyMax, Double_t nSigmaTOFmax, Double_t nSigmaTPCmax, Int_t TRDntracklets)  {

        fPtMin                     = PtMin;
        fPtMax                     = PtMax;
        fEtaMax                    = EtaMax;
        fYMax                      = YMax;
        fNumberClustersITSMin      = NumberClustersITSMin;
        fNumberClustersTPCMin      = NumberClustersTPCMin;
        fNumberCrossedRowsTPCMin   = NumberCrossedRowsTPCMin;
        fCrossedRowsFindableClsMin = CrossedRowsFindableClsMin;
        fNumberClustersTPCdEdxMin  = NumberClustersTPCdEdxMin;
        fChiSquarePerNDFMax        = ChiSquarePerNDFMax;
        fITSrequirement            = ITSrequirement;
        fDCAzMax                   = DCAzMax;
        fDCAxyMax                  = DCAxyMax;
        fnSigmaTOFmax              = nSigmaTOFmax;
        fnSigmaTPCmax              = nSigmaTPCmax;
        fTRDntracklets             = TRDntracklets;
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



    Bool_t      GetInputEvent ();
    Bool_t      PassedTrackQualityCutsNoDCA (AliAODTrack *track);
    Bool_t      PassedTrackQualityCuts      (AliAODTrack *track);
    Bool_t      IsTritonCandidate           (AliAODTrack *track);
    Bool_t      IsCleanTritonCandidate      (AliAODTrack *track);
    Double_t    Centered_nsigmaTPC          (AliAODTrack *track);
    Double_t    Centered_nsigmaTOF          (AliAODTrack *track);
    Bool_t      PassedTOFSelection          (AliAODTrack *track);
    Bool_t      PassedTPCSelection          (AliAODTrack *track);
    Double_t    GetDCAxy                    (AliAODTrack *track);
    Double_t    GetDCAz                     (AliAODTrack *track);

    virtual void   Terminate(Option_t *);

private:
    AliAODEvent      *fAODevent;//!
    AliPIDResponse   *fPIDResponse;//!
    AliEventCuts      fAODeventCuts;//
    AliAnalysisUtils *fUtils;//!
    TList            *fOutputList;//!
    //TList            *fQAList;//!

    Double_t fCentralityMin;//
    Double_t fCentralityMax;//
    Double_t fVertexZmin;//
    Double_t fVertexZmax;//
    Int_t fNumberVertexContributorsMin;//
    const char *fCentralityEstimator;//

    Double_t fPtMin;//
    Double_t fPtMax;//
    Double_t fEtaMax;//
    Double_t fYMax;//
    Int_t fNumberClustersITSMin;//
    Int_t fNumberClustersTPCMin;//
    Int_t fNumberCrossedRowsTPCMin;//
    Double_t fCrossedRowsFindableClsMin;//
    Int_t fNumberClustersTPCdEdxMin;//
    Double_t fChiSquarePerNDFMax;//
    const char *fITSrequirement;//
    Double_t fDCAzMax;//
    Double_t fDCAxyMax;//
    Double_t fnSigmaTOFmax;//
    Double_t fnSigmaTPCmax;//
    Int_t fTRDntracklets;//


    Double_t fpar0_mean_TPC;//
    Double_t fpar1_mean_TPC;//
    Double_t fpar0_sigma_TPC;//
    Double_t fpar0_mean_TOF;//
    Double_t fpar1_mean_TOF;//
    Double_t fpar0_sigma_TOF;//
    Double_t fpar1_sigma_TOF;//

    //Histograms
    TH1F *histoNumberOfEvents;//!

    TH2F *histoNsigmaTPCtriton_vs_pt;//!
    TH2F *histoNsigmaTOFtriton_vs_pt;//!

    TH2F *histoNsigmaTPCantitriton_vs_pt;//!
    TH2F *histoNsigmaTOFantitriton_vs_pt;//!

    TH2F *histoNsigmaTPCtriton_vs_pt_centered;//!
    TH2F *histoNsigmaTPCantitriton_vs_pt_centered;//!

    TH2F *histoNsigmaTOFtriton_vs_pt_centered;//!
    TH2F *histoNsigmaTOFantitriton_vs_pt_centered;//!

    TH2F *histoDCAxyTriton_vs_pt;//!
    TH2F *histoDCAxyAntiTriton_vs_pt;//!

    TH2F *histoNsigmaTOFtriton_vs_pt_trd;//!
    TH2F *histoNsigmaTOFantitriton_vs_pt_trd;//!

    TH2F *histoNsigmaTPCtriton_vs_p;//!
    TH2F *histoNsigmaTPCantitriton_vs_p;//!

    TH2F *histoNsigmaTOFtriton_vs_p;//!
    TH2F *histoNsigmaTOFantitriton_vs_p;//!

    TH2F *histoNsigmaTPCtriton_vs_p_notof;//!
    TH2F *histoNsigmaTPCantitriton_vs_p_notof;//!



    //Reduced Trees
    TTree *reducedTree_Triton;//!

    Double_t multPercentile_V0M;//
    Double_t multPercentile_V0A;//
    Double_t pt;//
    Double_t p;//
    Double_t px;//
    Double_t py;//
    Double_t pz;//
    Double_t eta;//
    Double_t y;//
    Int_t q;//
    Double_t dcaxy;//
    Double_t dcaz;//
    Int_t nTPC_Clusters;//
    Int_t nTRD_Clusters;//
    Int_t nITS_Clusters;//
    Int_t nTPC_FindableClusters;//
    Int_t nTPC_CrossedRows;//
    Int_t nTPC_Clusters_dEdx;//
    bool HasPointOnITSLayer0;//
    bool HasPointOnITSLayer1;//
    bool HasPointOnITSLayer2;//
    bool HasPointOnITSLayer3;//
    bool HasPointOnITSLayer4;//
    bool HasPointOnITSLayer5;//
    Double_t chi2_TPC;//check
    Double_t chi2_ITS;//check
    Double_t ITSsignal;
    Double_t TPCsignal;
    Double_t TOFsignal;
    Double_t nSigmaITS_Trit;//
    Double_t nSigmaTPC_Trit;//
    Double_t nSigmaTOF_Trit;//


    AliAnalysisTaskTritonVsMultiplicity_XeXe(const AliAnalysisTaskTritonVsMultiplicity_XeXe&);
    AliAnalysisTaskTritonVsMultiplicity_XeXe& operator=(const AliAnalysisTaskTritonVsMultiplicity_XeXe&);

    ClassDef(AliAnalysisTaskTritonVsMultiplicity_XeXe, 1);
};
#endif
