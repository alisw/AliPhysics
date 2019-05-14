#ifndef AliAnalysisTaskReducedTreeDS_cxx
#define AliAnalysisTaskReducedTreeDS_cxx

#include "AliAnalysisTaskSE.h"
using namespace std;

class AliAnalysisTaskReducedTreeDS : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskReducedTreeDS();
    AliAnalysisTaskReducedTreeDS(const char* name);
    virtual ~AliAnalysisTaskReducedTreeDS();
    void SetMinPtCut(Float_t min){fMinPtCut = min;}
    void SetMaxEtaCut(Float_t max){fMaxEtaCut = max;}
    void SetMaxTPCNsigmaEleCut(Float_t max){fMaxTPCNsigmaEleCut = max;}

  protected:
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    virtual void ProcessMC(Option_t *option);
    Double_t PsiPair(AliAODv0 *v0, Float_t Bz);
    Double_t PhivPair(AliAODv0 *v0, Float_t Bz);

    void GetMCInfoAOD(){
      fMCArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if(!fMCArray){
        AliError("Could not retrieve MC array!");
        return;
      }
    }

  protected:
    Float_t fMinPtCut;
    Float_t fMaxEtaCut;
    Float_t fMaxTPCNsigmaEleCut;
    TTree *fTree;
    AliPIDResponse *fPIDResponse;     //! PID response object
    AliVEvent *fEvent; 
    AliESDEvent *fESDEvent;
    AliAODEvent *fAODEvent;
    TClonesArray *fMCArray;
    Int_t fRunNumber;
    Float_t fMagneticField; //Bz in kG
    AliMultSelection *fMultSelection; 
    Float_t fCentralityV0M;
    Float_t fCentralityV0A;
    Float_t fCentralityV0C;
    Float_t fCentralityZNA;
    Float_t fCentralityZNC;
    Float_t fCentralityCL0;//SPD inner layer
    Float_t fCentralityCL1;//SPD outer layer
    Float_t fVertex[3];
    Int_t fNContributor;
    Int_t fNTPCCluster;
    Int_t fNTrackTPCout;
    Int_t fNHybridTrack08;
    Int_t fNSPDTracklet05;//|eta| < 0.5
    Int_t fNSPDTracklet10;//|eta| < 1.0
    Float_t fV0AMultiplicity;
    Float_t fV0CMultiplicity;
    Bool_t fIsPileupFromSPD;
    Bool_t fIsPileupFromSPDInMultBins;
    Bool_t fIsPileupMV;//SPD multi vertexer

    Bool_t fIskINT7;
    Bool_t fIskCentral;
    Bool_t fIskSemiCentral;

    Bool_t fIskHighMult;
    Bool_t fIskHighMultV0;
    Bool_t fIskHighMultSPD;

    vector <TVector3> fTrackMomentum;//use TVector3
    vector<Int_t>fTrackCharge;
    vector<Float_t>fTrackDCAxy;
    vector<Float_t>fTrackDCAz;

    vector<Float_t>fTrackPin;//momentum at inner wall of TPC (if available), used for PID
    vector<vector <Bool_t>> fPointOnITSLayer;//
    vector<vector<Bool_t>> fSharedPointOnITSLayer;

    vector<Float_t>fChi2TPC;
    vector<Float_t>fChi2ITS;
    vector<Int_t>fNclusterTPC;
    vector<Int_t>fNclusterITS;
    vector<Int_t>fTPCNCrossedRows;
    vector<Int_t>fTPCNFindableCluster;
    vector<Float_t>fChi2TPCConstrainedVsGlobal;

    vector<Float_t>fTPCsignal;
    vector<Float_t>fTPCNsigmaEl;
    vector<Float_t>fTPCNsigmaPi;
    vector<Float_t>fTPCNsigmaKa;
    vector<Float_t>fTPCNsigmaPr;

    vector<Float_t>fITSsignal;
    vector<Float_t>fITSNsigmaEl;
    vector<Float_t>fITSNsigmaPi;
    vector<Float_t>fITSNsigmaKa;
    vector<Float_t>fITSNsigmaPr;

    vector<Float_t>fTOFbeta;
    vector<Float_t>fTOFNsigmaEl;
    vector<Float_t>fTOFNsigmaPi;
    vector<Float_t>fTOFNsigmaKa;
    vector<Float_t>fTOFNsigmaPr;
    vector<Bool_t>fIsTOFAvailable;

    //MC track info
    vector<TVector3>fTrackMCMomentum;
    vector<TVector3>fTrackMCProdVtx;//production vertex in MC for track
    vector<Int_t>fTrackMCIndex;
    vector<Int_t>fTrackMCPdgCode;
    vector<Int_t>fTrackMCMotherIndex;
    vector<Int_t>fTrackMCMotherPdgCode;
    vector<Int_t>fTrackMCFirstMotherIndex;
    vector<Int_t>fTrackMCFirstMotherPdgCode;

    //V0 info
    vector<TVector3> fV0Momentum;//N
    vector<vector<TVector3>> fV0legMomentum;//N x 2
    vector<vector<Float_t>> fV0legPin;//N x 2
    vector<Float_t> fV0Lxy;//N 
    vector<Float_t> fV0alpha;//N 
    vector<Float_t> fV0qT;//N 
    vector<Float_t> fV0DCA;//N //DCA between daughters
    vector<Float_t> fV0PsiPair;//N //Angle cut w.r.t. to magnetic field
    vector<Float_t> fV0PhivPair;//N
    vector<Float_t> fV0PointingAngle;//N
    vector<Float_t> fV0Chi2;//N
    vector<vector<Float_t>> fV0legChi2TPCConstrainedVsGlobal;
    vector<vector<Float_t>> fV0Mass;//N x 4//K0S Lambda Anti-Lambda, Gamma
    vector<vector<Float_t>> fV0legDCAxy;//N x 2
    vector<vector<Float_t>> fV0legDCAz;//N x 2
    vector<vector<vector<Bool_t>>> fV0PointOnITSLayer;//N x 2 x 6
    vector<vector<vector<Bool_t>>> fV0SharedPointOnITSLayer;//N x 2 x 6
    vector<vector<Float_t>> fV0legTPCNsigmaEl;//N x 2//for post calibration
    vector<vector<Float_t>> fV0legTPCNsigmaPi;//N x 2//for post calibration
    vector<vector<Float_t>> fV0legTPCNsigmaKa;//N x 2//for post calibration
    vector<vector<Float_t>> fV0legTPCNsigmaPr;//N x 2//for post calibration
    vector<vector<Float_t>> fV0legITSNsigmaEl;//N x 2//for post calibration
    vector<vector<Float_t>> fV0legITSNsigmaPi;//N x 2//for post calibration
    vector<vector<Float_t>> fV0legITSNsigmaKa;//N x 2//for post calibration
    vector<vector<Float_t>> fV0legITSNsigmaPr;//N x 2//for post calibration
    vector<vector<Float_t>> fV0legTOFNsigmaEl;//N x 2//for post calibration
    vector<vector<Float_t>> fV0legTOFNsigmaPi;//N x 2//for post calibration
    vector<vector<Float_t>> fV0legTOFNsigmaKa;//N x 2//for post calibration
    vector<vector<Float_t>> fV0legTOFNsigmaPr;//N x 2//for post calibration
    vector<vector<Bool_t>>  fV0legIsTOFAvailable;//N x 2//for post calibration

    //MC V0 info //be carefull, there is no TRUE V0 object!
    vector<vector<TVector3>> fV0MClegMomentum;//N x 2
    vector<vector<TVector3>> fV0MClegProdVtx;//N x 2
    vector<vector<Int_t>> fV0MClegIndex;//N x 2
    vector<vector<Int_t>> fV0MClegPdgCode;//N x 2
    vector<vector<Int_t>> fV0MClegMotherIndex;//N x 2
    vector<vector<Int_t>> fV0MClegMotherPdgCode;//N x 2
    vector<vector<Int_t>> fV0MClegFirstMotherIndex;//N x 2
    vector<vector<Int_t>> fV0MClegFirstMotherPdgCode;//N x 2

    //MC variables for true electrons
    Float_t fMCVertex[3];//true vertex in MC
    vector<TVector3> fMCMomentum;
    vector<TVector3> fMCProdVtx;//production vertex of true electrons
    vector<Int_t> fMCIndex;
    vector<Int_t> fMCPdgCode;
    vector<Int_t> fMCMotherIndex;
    vector<Int_t> fMCMotherPdgCode;
    vector<Int_t> fMCFirstMotherIndex;
    vector<Int_t> fMCFirstMotherPdgCode;

    AliAnalysisTaskReducedTreeDS(const AliAnalysisTaskReducedTreeDS&); // not implemented
    AliAnalysisTaskReducedTreeDS& operator=(const AliAnalysisTaskReducedTreeDS&); // not implemented

    ClassDef(AliAnalysisTaskReducedTreeDS, 2);

};

#endif

