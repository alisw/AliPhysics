#ifndef ALIANALYSISTASKPI0V2_H
#define ALIANALYSISTASKPI0V2_H

class TH1F;
class TH2F;
class TH3F;
class TList;
class AliESDCaloCluster;
class AliVCluster;
class AliESDtrackCuts;
class AliESDEvent;
class THnSparse;
class TClonesArray;
class TString;
class TProfile;
class TProfile2D;
class AliOADBContainer;
class AliEPFlattener;
class AliEMCALGeometry;
class AliCalorimeterUtils;
class AliEventplane;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisTaskPi0V2 : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskPi0V2(const char *name);
    AliAnalysisTaskPi0V2();
    virtual ~AliAnalysisTaskPi0V2();
    
    virtual void           UserCreateOutputObjects();
    virtual void           UserExec(Option_t *option);
    virtual void           Terminate(Option_t *);

    void                   SetDebug(Int_t d)                     {fDebug          = d;}
    void                   UseV2Cluster(Bool_t e)                {fUseV2Cluster   = e;}
    void                   UseV1Cluster(Bool_t e)                {fUseV1Cluster   = e;}
    void                   UseTrack(Bool_t e)                    {fUseTrk         = e;}
    void                   SetV2ClusterName(TString n)           {fV2ClusterName  = n;} 
    void                   SetV1ClusterName(TString n)           {fV1ClusterName  = n;} 
    void                   SetTrackName(TString n)               {fTrackName      = n;}

    void                   SetEvtVzCut(Double_t v)               {fVzCut       = v;}
    void                   SetCentCut(Double_t min, Double_t max){fCentMin     = min; fCentMax = max;}
    void                   SetCentDetector(TString d)            {fCentDetector = d;}
    void                   FlattenSemiCent(Bool_t e)             {fFlattenSemiCent = e;}
    void                   UsePhosEPCali(Bool_t e)               {fUsePhosEPCali      = e;}
    void                   SetPhosEPCaliFileName(TString n)      {fPhosEPCaliFileName = n;}
    void                   SetCaloUtils(AliCalorimeterUtils* cu) {fCaloUtils = cu;}

    void                   SetClusterNCell(Double_t c)           {fNCellCut    = c;}
    void                   SetClusterE(Double_t e)               {fECut        = e;}
    void                   SetClusterEta(Double_t e)             {fEtaCut      = e;}
    void                   SetClusterDrCut(Double_t m)           {fDrCut       = m;}
    void                   SetV2M02Cut(Double_t m)               {fV2M02Cut    = m;}
    void                   SetV1M02Cut(Double_t m)               {fV1M02Cut    = m;}
    void                   CutV2ClusterPi0Asy(Double_t a)        {fPi0AsyCut   = a;}
    void                   SetNLMCut(Int_t min, Int_t max)       {fNLMCutMin   = min; fNLMCutMax = max;}
    void                   ApplySSCut(Bool_t a)                  {fApplySSCut  = a;}
    void                   SplitV1Cluster(Bool_t s)              {fSplitV1Cluster = s;}

 private:
    Int_t                  ConvertToInternalRunNum(Int_t n);
    Bool_t                 IsCentAccepted();
    void                   VZEROEventPlane(Bool_t flattenEP);
    Double_t               FlattenV0A(Double_t phi, Double_t c);
    Double_t               FlattenV0C(Double_t phi, Double_t c);
    Double_t               FlattenTPC(Double_t phi, Double_t c);

    Bool_t                 IsGoodCluster(const AliVCluster *c) const;
    Double_t               GetMaxCellEnergy(const AliVCluster *c, Short_t &id) const;
    Double_t               GetCrossEnergy(const AliVCluster *c, Short_t &idmax) const;
    Bool_t                 IsWithinFiducialVolume(Short_t id) const;
    void                   FillPionFromV2(const TLorentzVector& p1, const TLorentzVector& p2);
    Bool_t                 PassPi0SSCut(const AliVCluster *c);
    void                   FillPionFromV1(const TLorentzVector& p1, AliVCluster* c);
    void                   FillPionFromSplitV1(const TLorentzVector& p1);
    Bool_t                 IsInPi0SplitAsymmetryRange(Float_t energy, Float_t asy,  Int_t nlm) const;
    void                   GetMom(TLorentzVector& p, const AliVCluster* c, Double_t* vertex);      

    Int_t                  fDebug;

    Bool_t                 fUseV2Cluster;
    Bool_t                 fUseV1Cluster;
    Bool_t                 fUseTrk;
    TString                fV2ClusterName;
    TString                fV1ClusterName;
    TString                fTrackName;
    TClonesArray*          fV2Cluster;
    TClonesArray*          fV1Cluster;
    TClonesArray*          fTrack;

    TList*                 fOutput;
    AliAODEvent*           fAODEvent;
    AliEMCALGeometry*      fGeom;
    TString                fGeomName;
    AliOADBContainer*      fPhosEPCaliContainer;
    TString                fPhosEPCaliFileName;
    Bool_t                 fUsePhosEPCali;
    AliEventplane*         fEventPlane;
    AliCalorimeterUtils*   fCaloUtils;
    Int_t                  fRunNum;
    Int_t                  fInterRunNum;
    Double_t               fVzCut;
    Int_t                  fVzBin;
    Double_t               fCentMin;
    Double_t               fCentMax;
    TString                fCentDetector;
    Double_t               fCentrality;
    Int_t                  fCentBin;
    Bool_t                 fFlattenSemiCent;

    Double_t               fNCellCut;
    Double_t               fECut;
    Double_t               fEtaCut;
    Double_t               fV2M02Cut;
    Double_t               fV1M02Cut;
    Double_t               fDrCut;
    Bool_t                 fPi0AsyCut;
    Int_t                  fNLMCutMin;
    Int_t                  fNLMCutMax;
    Bool_t                 fApplySSCut;
    Bool_t                 fSplitV1Cluster;
    Int_t                  fBufferIndex[10][2];
    TLorentzVector         fBufferSplitV1[10][2][10];

    Double_t               fEPTPC;
    Double_t               fEPTPCReso;
    Double_t               fEPV0;
    Double_t               fEPV0A;
    Double_t               fEPV0C;
    Double_t               fEPV0AR;
    Double_t               fEPV0CR;
    Double_t               fEPV0R;
    Double_t               fEPV0AR4;
    Double_t               fEPV0AR5;
    Double_t               fEPV0AR6;
    Double_t               fEPV0AR7;
    Double_t               fEPV0CR0;
    Double_t               fEPV0CR1;
    Double_t               fEPV0CR2;
    Double_t               fEPV0CR3;
    AliEPFlattener*        fEPTPCFlat;
    AliEPFlattener*        fEPV0AFlat;
    AliEPFlattener*        fEPV0CFlat;

    TH1F*                  hEvtCount;
    TH1F*                  hCentA;
    TH1F*                  hCentB;

    TH2F*                  hEPTPC;    
    TH2F*                  hEPTPCReso;
    TH2F*                  hEPV0A;
    TH2F*                  hEPV0C;
    TH2F*                  hEPTPCFlat;
    TH2F*                  hEPV0AFlat;
    TH2F*                  hEPV0CFlat;
    TH2F*                  hEPV0;
    TH2F*                  hEPV0AR;
    TH2F*                  hEPV0CR;
    TH2F*                  hEPV0R;
    TH2F*                  hEPV0AR4;
    TH2F*                  hEPV0AR7;
    TH2F*                  hEPV0CR0;
    TH2F*                  hEPV0CR3;
    TH2F*                  hEPDiffV0A_V0CR0;
    TH2F*                  hEPDiffV0A_V0CR3;
    TH2F*                  hEPDiffV0CR0_V0CR3;
    TH2F*                  hEPDiffV0C_V0AR4;
    TH2F*                  hEPDiffV0C_V0AR7;
    TH2F*                  hEPDiffV0AR4_V0AR7;
    TProfile2D*            hEPRbrCosV0A;
    TProfile2D*            hEPRbrSinV0A;
    TProfile2D*            hEPRbrCosV0C;
    TProfile2D*            hEPRbrSinV0C;
    TProfile2D*            hEPRbrCosTPC;
    TProfile2D*            hEPRbrSinTPC;

    TH2F*                  hV2ClusterDxDzA;
    TH2F*                  hV2ClusterDxDzB;
    TH3F*                  hV2ClusterDphiV0A;
    TH3F*                  hV2ClusterDphiV0C;
    TH3F*                  hV2ClusterCos2phiV0A;
    TH3F*                  hV2ClusterCos2phiV0C;
    TH2F*                  hV1ClusterDxDzA;
    TH2F*                  hV1ClusterDxDzB;
    TH2F*                  hV1ClusterM02EA;
    TH2F*                  hV1ClusterM02EB;
    TH1F*                  hV1ClusterNlmA;
    TH1F*                  hV1ClusterNlmB;

    TH2F*                  hTrkPhiEta;
    TH1F*                  hTrkPt;
    TH3F*                  hTrkDphiEmcV0A;
    TH3F*                  hTrkDphiEmcV0C;
    TH3F*                  hTrkCos2phiEmcV0A;
    TH3F*                  hTrkCos2phiEmcV0C;
    TH3F*                  hTrkDphiOutEmcV0A;
    TH3F*                  hTrkDphiOutEmcV0C;
    TH3F*                  hTrkCos2phiOutEmcV0A;
    TH3F*                  hTrkCos2phiOutEmcV0C;

    THnSparse*             fV2ClusterV0A;
    THnSparse*             fV2ClusterV0C;
    THnSparse*             fV2ClusterTPC;
    THnSparse*             fV1ClusterV0A;
    THnSparse*             fV1ClusterV0C;
    THnSparse*             fV1ClusterTPC;

    AliAnalysisTaskPi0V2(const AliAnalysisTaskPi0V2&);
    AliAnalysisTaskPi0V2& operator=(const AliAnalysisTaskPi0V2&);

    ClassDef(AliAnalysisTaskPi0V2, 6);
};
#endif
