#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TString.h>
#include <TFile.h>
#include <TList.h>
#include <THashList.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TVector3.h>
#include <TClonesArray.h>

#include "AliESDtrackCuts.h"
#include "AliESDHeader.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"

#include "AliVEvent.h"
#include "AliVHeader.h"
#include "AliVTrack.h"

#include "AliMultSelection.h"
#include "AliEventplane.h"
#include "AliQnCorrectionsManager.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliQnCorrectionsQnVector.h"

#include "AliOADBContainer.h"

#include "AliPID.h"
#include "AliPIDResponse.h"

#include "AliAODMCHeader.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODInputHandler.h"
#include "AliAODTracklets.h"
#include "AliAnalysisUtils.h"

#include "AliKFVertex.h"
#include "AliDielectronPair.h"


#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskReducedTreeDS.h"

using namespace std;

//event property and single electron tree.
//Daiki Sekihata
//daiki.sekihata@cern.ch

ClassImp(AliAnalysisTaskReducedTreeDS)
//_______________________________________________________________________________________________
AliAnalysisTaskReducedTreeDS::AliAnalysisTaskReducedTreeDS():
  AliAnalysisTaskSE(),
  fMinPtCut(0.2),
  fMaxEtaCut(0.8),
  fMaxTPCNsigmaEleCut(5.),
  fTree(0x0),
  fPIDResponse(0x0),
  fFlowQnVectorMgr(0x0),
  fEvent(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fMCArray(0x0),
  fRunNumber(-1),
  fMagneticField(0),
  fMultSelection(0x0),
  fCentralityV0M(-1),
  fCentralityV0A(-1),
  fCentralityV0C(-1),
  fCentralityZNA(-1),
  fCentralityZNC(-1),
  fCentralityCL0(-1),
  fCentralityCL1(-1),
  fVertex(),
  fNContributor(-1),
  fNTPCCluster(-1),
  fNTrackTPCout(-1),
  fNHybridTrack08(-1),
  fNSPDTracklet05(-1),
  fNSPDTracklet10(-1),
  fV0AMultiplicity(-1),
  fV0CMultiplicity(-1),
  fIsPileupFromSPD(kFALSE),
  fIsPileupFromSPDInMultBins(kFALSE),
  fIsPileupMV(kFALSE),
  fIskINT7(kFALSE),
  fIskCentral(kFALSE),
  fIskSemiCentral(kFALSE),
  fIskHighMult(kFALSE),
  fIskHighMultV0(kFALSE),
  fIskHighMultSPD(kFALSE),
  fIsQnTPCAvailable(kFALSE),
  fQ2vectorTPC(),
  fQ2vectorTPCNegEta(),
  fQ2vectorTPCPosEta(),
  fIsQnV0Available(kFALSE),
  fQ2vectorV0(),
  fQ2vectorV0A(),
  fQ2vectorV0C(),
  fIsQnZDCAvailable(kFALSE),
  fQ2vectorZDCA(),
  fQ2vectorZDCC(),
  fTrackMomentum(0),
  fTrackCharge(0),
  fTrackDCAxy(0),
  fTrackDCAz(0),
  fTrackPin(0),
  fPointOnITSLayer(0),
  fSharedPointOnITSLayer(0),
  fChi2TPC(0),
  fChi2ITS(0),
  fNclusterTPC(0),
  fNclusterITS(0),
  fTPCNCrossedRows(0),
  fTPCNFindableCluster(0),
  fChi2TPCConstrainedVsGlobal(0),
  fTPCsignal(0),
  fTPCNsigmaEl(0),
  fTPCNsigmaPi(0),
  fTPCNsigmaKa(0),
  fTPCNsigmaPr(0),
  fITSsignal(0),
  fITSNsigmaEl(0),
  fITSNsigmaPi(0),
  fITSNsigmaKa(0),
  fITSNsigmaPr(0),
  fTOFbeta(0),
  fTOFNsigmaEl(0),
  fTOFNsigmaPi(0),
  fTOFNsigmaKa(0),
  fTOFNsigmaPr(0),
  fIsTOFAvailable(0),
  fTrackMCMomentum(0),
  fTrackMCProdVtx(0),
  fTrackMCIndex(0),
  fTrackMCPdgCode(0),
  fTrackMCMotherIndex(0),
  fTrackMCMotherPdgCode(0),
  fTrackMCFirstMotherIndex(0),
  fTrackMCFirstMotherPdgCode(0),
  fV0Momentum(0),
  fV0legMomentum(0),
  fV0legPin(0),
  fV0Lxy(0),
  fV0alpha(0),
  fV0qT(0),
  fV0DCA(0),
  fV0PsiPair(0),
  fV0PhivPair(0),
  fV0PointingAngle(0),
  fV0Chi2(0),
  fV0legChi2TPCConstrainedVsGlobal(0),
  fV0Mass(0),
  fV0legDCAxy(0),
  fV0legDCAz(0),
  fV0PointOnITSLayer(0),
  fV0SharedPointOnITSLayer(0),
  fV0legTPCNsigmaEl(0),
  fV0legTPCNsigmaPi(0),
  fV0legTPCNsigmaKa(0),
  fV0legTPCNsigmaPr(0),
  fV0legITSNsigmaEl(0),
  fV0legITSNsigmaPi(0),
  fV0legITSNsigmaKa(0),
  fV0legITSNsigmaPr(0),
  fV0legTOFNsigmaEl(0),
  fV0legTOFNsigmaPi(0),
  fV0legTOFNsigmaKa(0),
  fV0legTOFNsigmaPr(0),
  fV0legIsTOFAvailable(0),
  fV0MClegMomentum(0),
  fV0MClegProdVtx(0),
  fV0MClegIndex(0),
  fV0MClegPdgCode(0),
  fV0MClegMotherIndex(0),
  fV0MClegMotherPdgCode(0),
  fV0MClegFirstMotherIndex(0),
  fV0MClegFirstMotherPdgCode(0),
  fMCVertex(),
  fMCMomentum(0),
  fMCProdVtx(0),
  fMCIndex(0),
  fMCPdgCode(0),
  fMCMotherIndex(0),
  fMCMotherPdgCode(0),
  fMCFirstMotherIndex(0),
  fMCFirstMotherPdgCode(0)
{

  for(Int_t i=0;i<3;i++) fVertex[i] = 0;
  for(Int_t i=0;i<3;i++) fMCVertex[i] = 0;

  for(Int_t i=0;i<2;i++){//Qx, Qy
    fQ2vectorTPC[i] = -999;
    fQ2vectorTPCNegEta[i] = -999;
    fQ2vectorTPCPosEta[i] = -999;
    fQ2vectorV0[i] = -999;
    fQ2vectorV0A[i] = -999;
    fQ2vectorV0C[i] = -999;
    fQ2vectorZDCA[i] = -999;
    fQ2vectorZDCC[i] = -999;
  }

}
//_______________________________________________________________________________________________
AliAnalysisTaskReducedTreeDS::AliAnalysisTaskReducedTreeDS(const char *name):
  AliAnalysisTaskSE(name),
  fMinPtCut(0.2),
  fMaxEtaCut(0.8),
  fMaxTPCNsigmaEleCut(5.),
  fTree(0x0),
  fPIDResponse(0x0),
  fFlowQnVectorMgr(0x0),
  fEvent(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fMCArray(0x0),
  fRunNumber(-1),
  fMagneticField(0),
  fMultSelection(0x0),
  fCentralityV0M(-1),
  fCentralityV0A(-1),
  fCentralityV0C(-1),
  fCentralityZNA(-1),
  fCentralityZNC(-1),
  fCentralityCL0(-1),
  fCentralityCL1(-1),
  fVertex(),
  fNContributor(-1),
  fNTPCCluster(-1),
  fNTrackTPCout(-1),
  fNHybridTrack08(-1),
  fNSPDTracklet05(-1),
  fNSPDTracklet10(-1),
  fV0AMultiplicity(-1),
  fV0CMultiplicity(-1),
  fIsPileupFromSPD(kFALSE),
  fIsPileupFromSPDInMultBins(kFALSE),
  fIsPileupMV(kFALSE),
  fIskINT7(kFALSE),
  fIskCentral(kFALSE),
  fIskSemiCentral(kFALSE),
  fIskHighMult(kFALSE),
  fIskHighMultV0(kFALSE),
  fIskHighMultSPD(kFALSE),
  fIsQnTPCAvailable(kFALSE),
  fQ2vectorTPC(),
  fQ2vectorTPCNegEta(),
  fQ2vectorTPCPosEta(),
  fIsQnV0Available(kFALSE),
  fQ2vectorV0(),
  fQ2vectorV0A(),
  fQ2vectorV0C(),
  fIsQnZDCAvailable(kFALSE),
  fQ2vectorZDCA(),
  fQ2vectorZDCC(),
  fTrackMomentum(0),
  fTrackCharge(0),
  fTrackDCAxy(0),
  fTrackDCAz(0),
  fTrackPin(0),
  fPointOnITSLayer(0),
  fSharedPointOnITSLayer(0),
  fChi2TPC(0),
  fChi2ITS(0),
  fNclusterTPC(0),
  fNclusterITS(0),
  fTPCNCrossedRows(0),
  fTPCNFindableCluster(0),
  fChi2TPCConstrainedVsGlobal(0),
  fTPCsignal(0),
  fTPCNsigmaEl(0),
  fTPCNsigmaPi(0),
  fTPCNsigmaKa(0),
  fTPCNsigmaPr(0),
  fITSsignal(0),
  fITSNsigmaEl(0),
  fITSNsigmaPi(0),
  fITSNsigmaKa(0),
  fITSNsigmaPr(0),
  fTOFbeta(0),
  fTOFNsigmaEl(0),
  fTOFNsigmaPi(0),
  fTOFNsigmaKa(0),
  fTOFNsigmaPr(0),
  fIsTOFAvailable(0),
  fTrackMCMomentum(0),
  fTrackMCProdVtx(0),
  fTrackMCIndex(0),
  fTrackMCPdgCode(0),
  fTrackMCMotherIndex(0),
  fTrackMCMotherPdgCode(0),
  fTrackMCFirstMotherIndex(0),
  fTrackMCFirstMotherPdgCode(0),
  fV0Momentum(0),
  fV0legMomentum(0),
  fV0legPin(0),
  fV0Lxy(0),
  fV0alpha(0),
  fV0qT(0),
  fV0DCA(0),
  fV0PsiPair(0),
  fV0PhivPair(0),
  fV0PointingAngle(0),
  fV0Chi2(0),
  fV0legChi2TPCConstrainedVsGlobal(0),
  fV0Mass(0),
  fV0legDCAxy(0),
  fV0legDCAz(0),
  fV0PointOnITSLayer(0),
  fV0SharedPointOnITSLayer(0),
  fV0legTPCNsigmaEl(0),
  fV0legTPCNsigmaPi(0),
  fV0legTPCNsigmaKa(0),
  fV0legTPCNsigmaPr(0),
  fV0legITSNsigmaEl(0),
  fV0legITSNsigmaPi(0),
  fV0legITSNsigmaKa(0),
  fV0legITSNsigmaPr(0),
  fV0legTOFNsigmaEl(0),
  fV0legTOFNsigmaPi(0),
  fV0legTOFNsigmaKa(0),
  fV0legTOFNsigmaPr(0),
  fV0legIsTOFAvailable(0),
  fV0MClegMomentum(0),
  fV0MClegProdVtx(0),
  fV0MClegIndex(0),
  fV0MClegPdgCode(0),
  fV0MClegMotherIndex(0),
  fV0MClegMotherPdgCode(0),
  fV0MClegFirstMotherIndex(0),
  fV0MClegFirstMotherPdgCode(0),
  fMCVertex(),
  fMCMomentum(0),
  fMCProdVtx(0),
  fMCIndex(0),
  fMCPdgCode(0),
  fMCMotherIndex(0),
  fMCMotherPdgCode(0),
  fMCFirstMotherIndex(0),
  fMCFirstMotherPdgCode(0)
{

  for(Int_t i=0;i<3;i++) fVertex[i] = 0;
  for(Int_t i=0;i<3;i++) fMCVertex[i] = 0;


  for(Int_t i=0;i<2;i++){
    fQ2vectorTPC[i] = -999;
    fQ2vectorTPCNegEta[i] = -999;
    fQ2vectorTPCPosEta[i] = -999;
    fQ2vectorV0[i] = -999;
    fQ2vectorV0A[i] = -999;
    fQ2vectorV0C[i] = -999;
    fQ2vectorZDCA[i] = -999;
    fQ2vectorZDCC[i] = -999;
  }



  DefineInput(0,TChain::Class());
  DefineOutput(1, TTree::Class());  // reduced information tree

}
//_______________________________________________________________________________________________
AliAnalysisTaskReducedTreeDS::~AliAnalysisTaskReducedTreeDS()
{

}
//_______________________________________________________________________________________________
void AliAnalysisTaskReducedTreeDS::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fTree = new TTree("EventTree","Reduced event tree for dielectron");
  fTree->SetDirectory(0);//force memory-resident tree

  fTree->Branch("fRunNumber",&fRunNumber,"fRunNumber/I");
  fTree->Branch("fMagneticField",&fMagneticField,"fMagneticField/F");

  fTree->Branch("fCentralityV0M",&fCentralityV0M,"fCentralityV0M/F");
  fTree->Branch("fCentralityV0A",&fCentralityV0A,"fCentralityV0A/F");
  fTree->Branch("fCentralityV0C",&fCentralityV0C,"fCentralityV0C/F");
  fTree->Branch("fCentralityZNA",&fCentralityZNA,"fCentralityZNA/F");
  fTree->Branch("fCentralityZNC",&fCentralityZNC,"fCentralityZNC/F");
  fTree->Branch("fCentralityCL0",&fCentralityCL0,"fCentralityCL0/F");
  fTree->Branch("fCentralityCL1",&fCentralityCL1,"fCentralityCL1/F");

  fTree->Branch("fVertex",fVertex,"fVertex[3]/F");
  fTree->Branch("fNContributor",&fNContributor,"fNContributor/I");
  fTree->Branch("fNTPCCluster",&fNTPCCluster,"fNTPCCluster/I");
  fTree->Branch("fNTrackTPCout",&fNTrackTPCout,"fNTrackTPCout/I");
  fTree->Branch("fNHybridTrack08",&fNHybridTrack08,"fNHybridTrack08/I");
  fTree->Branch("fNSPDTracklet05",&fNSPDTracklet05,"fNSPDTracklet05/I");
  fTree->Branch("fNSPDTracklet10",&fNSPDTracklet10,"fNSPDTracklet10/I");
  fTree->Branch("fV0AMultiplicity",&fV0AMultiplicity,"fV0AMultiplicity/F");
  fTree->Branch("fV0CMultiplicity",&fV0CMultiplicity,"fV0CMultiplicity/F");

  fTree->Branch("fIsPileupFromSPD",&fIsPileupFromSPD,"fIsPileupFromSPD/O");
  fTree->Branch("fIsPileupFromSPDInMultBins",&fIsPileupFromSPDInMultBins,"fIsPileupFromSPDInMultBins/O");
  fTree->Branch("fIsPileupMV",&fIsPileupMV,"fIsPileupMV/O");

  fTree->Branch("fIskINT7",&fIskINT7,"fIskINT7/O");
  fTree->Branch("fIskCentral",&fIskCentral,"fIskCentral/O");
  fTree->Branch("fIskSemiCentral",&fIskSemiCentral,"fIskSemiCentral/O");
  fTree->Branch("fIskHighMult",&fIskHighMult,"fIskHighMult/O");
  fTree->Branch("fIskHighMultV0",&fIskHighMultV0,"fIskHighMultV0/O");
  fTree->Branch("fIskHighMultSPD",&fIskHighMultSPD,"fIskHighMultSPD/O");

  fTree->Branch("fIsQnTPCAvailable",&fIsQnTPCAvailable,"fIsQnTPCAvailable/O");
  fTree->Branch("fQ2vectorTPC",fQ2vectorTPC,"fQ2vectorTPC[2]/F");
  fTree->Branch("fQ2vectorTPCNegEta",fQ2vectorTPCNegEta,"fQ2vectorTPCNegEta[2]/F");
  fTree->Branch("fQ2vectorTPCPosEta",fQ2vectorTPCPosEta,"fQ2vectorTPCPosEta[2]/F");
  fTree->Branch("fIsQnV0Available",&fIsQnV0Available,"fIsQnV0Available/O");
  fTree->Branch("fQ2vectorV0",fQ2vectorV0,"fQ2vectorV0[2]/F");
  fTree->Branch("fQ2vectorV0A",fQ2vectorV0A,"fQ2vectorV0A[2]/F");
  fTree->Branch("fQ2vectorV0C",fQ2vectorV0C,"fQ2vectorV0C[2]/F");
  fTree->Branch("fIsQnZDCAvailable",&fIsQnZDCAvailable,"fIsQnZDCAvailable/O");
  fTree->Branch("fQ2vectorZDCA",fQ2vectorZDCA,"fQ2vectorZDCA[2]/F");
  fTree->Branch("fQ2vectorZDCC",fQ2vectorZDCC,"fQ2vectorZDCC[2]/F");

  fTree->Branch("fTrackMomentum",&fTrackMomentum);
  fTree->Branch("fTrackCharge",&fTrackCharge);
  fTree->Branch("fTrackDCAxy",&fTrackDCAxy);
  fTree->Branch("fTrackDCAz",&fTrackDCAz);

  fTree->Branch("fTrackPin",&fTrackPin);
  fTree->Branch("fPointOnITSLayer",&fPointOnITSLayer);
  fTree->Branch("fSharedPointOnITSLayer",&fSharedPointOnITSLayer);

  fTree->Branch("fChi2TPC",&fChi2TPC);
  fTree->Branch("fChi2ITS",&fChi2ITS);
  fTree->Branch("fNclusterTPC",&fNclusterTPC);
  fTree->Branch("fNclusterITS",&fNclusterITS);
  fTree->Branch("fTPCNCrossedRows",&fTPCNCrossedRows);
  fTree->Branch("fTPCNFindableCluster",&fTPCNFindableCluster);
  fTree->Branch("fChi2TPCConstrainedVsGlobal",&fChi2TPCConstrainedVsGlobal);

  fTree->Branch("fTPCsignal",&fTPCsignal);
  fTree->Branch("fTPCNsigmaEl",&fTPCNsigmaEl);
  fTree->Branch("fTPCNsigmaPi",&fTPCNsigmaPi);
  fTree->Branch("fTPCNsigmaKa",&fTPCNsigmaKa);
  fTree->Branch("fTPCNsigmaPr",&fTPCNsigmaPr);

  fTree->Branch("fITSsignal",&fITSsignal);
  fTree->Branch("fITSNsigmaEl",&fITSNsigmaEl);
  fTree->Branch("fITSNsigmaPi",&fITSNsigmaPi);
  fTree->Branch("fITSNsigmaKa",&fITSNsigmaKa);
  fTree->Branch("fITSNsigmaPr",&fITSNsigmaPr);

  fTree->Branch("fTOFbeta",&fTOFbeta);
  fTree->Branch("fTOFNsigmaEl",&fTOFNsigmaEl);
  fTree->Branch("fTOFNsigmaPi",&fTOFNsigmaPi);
  fTree->Branch("fTOFNsigmaKa",&fTOFNsigmaKa);
  fTree->Branch("fTOFNsigmaPr",&fTOFNsigmaPr);
  fTree->Branch("fIsTOFAvailable",&fIsTOFAvailable);

  //MC info for reconstructed tracks
  fTree->Branch("fTrackMCMomentum",&fTrackMCMomentum);
  fTree->Branch("fTrackMCProdVtx",&fTrackMCProdVtx);
  fTree->Branch("fTrackMCIndex",&fTrackMCIndex);
  fTree->Branch("fTrackMCPdgCode",&fTrackMCPdgCode);
  fTree->Branch("fTrackMCMotherIndex",&fTrackMCMotherIndex);
  fTree->Branch("fTrackMCMotherPdgCode",&fTrackMCMotherPdgCode);
  fTree->Branch("fTrackMCFirstMotherIndex",&fTrackMCFirstMotherIndex);
  fTree->Branch("fTrackMCFirstMotherPdgCode",&fTrackMCFirstMotherPdgCode);

  fTree->Branch("fV0Momentum",&fV0Momentum);
  fTree->Branch("fV0legMomentum",&fV0legMomentum);
  fTree->Branch("fV0legPin",&fV0legPin);

  fTree->Branch("fV0Lxy",&fV0Lxy);
  fTree->Branch("fV0alpha",&fV0alpha);
  fTree->Branch("fV0qT",&fV0qT);
  fTree->Branch("fV0DCA",&fV0DCA);
  fTree->Branch("fV0PsiPair",&fV0PsiPair);
  fTree->Branch("fV0PhivPair",&fV0PhivPair);

  fTree->Branch("fV0PointingAngle",&fV0PointingAngle);
  fTree->Branch("fV0Chi2",&fV0Chi2);
  fTree->Branch("fV0legChi2TPCConstrainedVsGlobal",&fV0legChi2TPCConstrainedVsGlobal);
  fTree->Branch("fV0Mass",&fV0Mass);
  fTree->Branch("fV0legDCAxy",&fV0legDCAxy);
  fTree->Branch("fV0legDCAz",&fV0legDCAz);
  fTree->Branch("fV0PointOnITSLayer",&fV0PointOnITSLayer);
  fTree->Branch("fV0SharedPointOnITSLayer",&fV0SharedPointOnITSLayer);

  fTree->Branch("fV0legTPCNsigmaEl",&fV0legTPCNsigmaEl);
  fTree->Branch("fV0legTPCNsigmaPi",&fV0legTPCNsigmaPi);
  fTree->Branch("fV0legTPCNsigmaKa",&fV0legTPCNsigmaKa);
  fTree->Branch("fV0legTPCNsigmaPr",&fV0legTPCNsigmaPr);

  fTree->Branch("fV0legITSNsigmaEl",&fV0legITSNsigmaEl);
  fTree->Branch("fV0legITSNsigmaPi",&fV0legITSNsigmaPi);
  fTree->Branch("fV0legITSNsigmaKa",&fV0legITSNsigmaKa);
  fTree->Branch("fV0legITSNsigmaPr",&fV0legITSNsigmaPr);
  fTree->Branch("fV0legTOFNsigmaEl",&fV0legTOFNsigmaEl);
  fTree->Branch("fV0legTOFNsigmaPi",&fV0legTOFNsigmaPi);
  fTree->Branch("fV0legTOFNsigmaKa",&fV0legTOFNsigmaKa);
  fTree->Branch("fV0legTOFNsigmaPr",&fV0legTOFNsigmaPr);
  fTree->Branch("fV0legIsTOFAvailable",&fV0legIsTOFAvailable);

  //MC info for reconstructed V0s
  fTree->Branch("fV0MClegMomentum",&fV0MClegMomentum);
  fTree->Branch("fV0MClegProdVtx",&fV0MClegProdVtx);
  fTree->Branch("fV0MClegIndex",&fV0MClegIndex);
  fTree->Branch("fV0MClegPdgCode",&fV0MClegPdgCode);
  fTree->Branch("fV0MClegMotherIndex",&fV0MClegMotherIndex);
  fTree->Branch("fV0MClegMotherPdgCode",&fV0MClegMotherPdgCode);
  fTree->Branch("fV0MClegFirstMotherIndex",&fV0MClegFirstMotherIndex);
  fTree->Branch("fV0MClegFirstMotherPdgCode",&fV0MClegFirstMotherPdgCode);

  //MC true info
  fTree->Branch("fMCVertex",fMCVertex,"fMCVertex[3]/F");
  fTree->Branch("fMCMomentum",&fMCMomentum);
  fTree->Branch("fMCProdVtx",&fMCProdVtx);
  fTree->Branch("fMCIndex",&fMCIndex);
  fTree->Branch("fMCPdgCode",&fMCPdgCode);
  fTree->Branch("fMCMotherIndex",&fMCMotherIndex);
  fTree->Branch("fMCMotherPdgCode",&fMCMotherPdgCode);
  fTree->Branch("fMCFirstMotherIndex",&fMCFirstMotherIndex);
  fTree->Branch("fMCFirstMotherPdgCode",&fMCFirstMotherPdgCode);

  fPIDResponse = dynamic_cast<AliPIDResponse*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetPIDResponse());
  if(!fPIDResponse){
    AliFatal("fPIDResponse does not exist!");
    return;
  }

  AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections*>(AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
  if(flowQnVectorTask != NULL){
    fFlowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
  }
  else {
    AliInfo("Flow Qn vector corrections framework is not found. Q vectors will be filled with useless values!");
  }

  PostData(1,fTree);

}
//_______________________________________________________________________________________________
void AliAnalysisTaskReducedTreeDS::UserExec(Option_t *option)
{
  // Main loop
  // Called for each event

  fEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if(!fEvent){
    AliError("event is not available.");
    return;
  }

  fESDEvent = dynamic_cast<AliESDEvent*>(fEvent);
  fAODEvent = dynamic_cast<AliAODEvent*>(fEvent);

  if(!fESDEvent && !fAODEvent){
    AliError("event type is neither ESD nor AOD. return.");
    return;
  }

  //fPIDResponse = fInputHandler->GetPIDResponse();
  //if(!fPIDResponse){
  //  AliFatal("fPIDResponse does not exist!");
  //  return;
  //}

  ClearVectorElement();

  fIskINT7        = kFALSE;
  fIskCentral     = kFALSE;
  fIskSemiCentral = kFALSE;
  fIskHighMult    = kFALSE;
  fIskHighMultV0  = kFALSE;
  fIskHighMultSPD = kFALSE;

  UInt_t SelectMask = fInputHandler->IsEventSelected();
  fIskINT7        = SelectMask & AliVEvent::kINT7;
  fIskCentral     = SelectMask & AliVEvent::kCentral;
  fIskSemiCentral = SelectMask & AliVEvent::kSemiCentral;
  fIskHighMult    = SelectMask & AliVEvent::kHighMult;
  fIskHighMultV0  = SelectMask & AliVEvent::kHighMultV0;
  fIskHighMultSPD = SelectMask & AliVEvent::kHighMultSPD;

  //AliVEventHandler* eventHandler = AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();//for ESD
  //cout << "eventHandler = " << eventHandler << endl;

  //AliAODInputHandler* aodHandler=dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  //cout << "MCEvent in AOD = " << aodHandler->MCEvent()<< endl;
  //AliMCEvent* MCevent = dynamic_cast<AliMCEvent*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->MCEvent());
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->MCEvent() != 0x0);
  AliInfo(Form("hasMC = %d",hasMC));

  if(hasMC) GetMCInfoAOD();

  fRunNumber = fEvent->GetRunNumber();
  fMagneticField = fEvent->GetMagneticField();

  //Get Centrality
  fMultSelection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
  if(!fMultSelection){
    //If you get this warning (and fCentralityV0M 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
  }
  else{
    fCentralityV0M = fMultSelection->GetMultiplicityPercentile("V0M");
    fCentralityV0A = fMultSelection->GetMultiplicityPercentile("V0A");
    fCentralityV0C = fMultSelection->GetMultiplicityPercentile("V0C");
    fCentralityZNA = fMultSelection->GetMultiplicityPercentile("ZNA");
    fCentralityZNC = fMultSelection->GetMultiplicityPercentile("ZNC");
    fCentralityCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
    fCentralityCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  }

  ExtractQnVectors();

  if(fESDEvent){
    fIsPileupFromSPD           = fESDEvent->IsPileupFromSPD();
    fIsPileupFromSPDInMultBins = fESDEvent->IsPileupFromSPDInMultBins();
  }
  else if(fAODEvent){
    fIsPileupFromSPD           = fAODEvent->IsPileupFromSPD();
    fIsPileupFromSPDInMultBins = fAODEvent->IsPileupFromSPDInMultBins();
  }

  const Int_t minContributors=5;
  const Float_t minChi2=5.;
  const Float_t minWeiZDiff=15;
  const Bool_t checkPlpFromDifferentBC=kFALSE;

  AliAnalysisUtils utils;
  utils.SetMinPlpContribMV(minContributors);
  utils.SetMaxPlpChi2MV(minChi2);
  utils.SetMinWDistMV(minWeiZDiff);
  utils.SetCheckPlpFromDifferentBCMV(checkPlpFromDifferentBC);
  fIsPileupMV = utils.IsPileUpMV(fEvent);

  AliVVZERO *V0info = (AliVVZERO*)fEvent->GetVZEROData();
  fV0AMultiplicity = V0info->GetMTotV0A();
  fV0CMultiplicity = V0info->GetMTotV0C();

  //Check SPD tracklet multilicity in |eta| < 1
  fNSPDTracklet05 = 0;
  fNSPDTracklet10 = 0;
  AliAODTracklets *tracklets = (AliAODTracklets*)fAODEvent->GetTracklets();
  const Int_t Ntl = tracklets->GetNumberOfTracklets();
  for(Int_t itl=0;itl<Ntl;itl++){
    Double_t theta = tracklets->GetTheta(itl);
    Double_t eta   = -TMath::Log(TMath::Tan(theta/2.0));
    if(TMath::Abs(eta) < 0.5) fNSPDTracklet05++;
    if(TMath::Abs(eta) < 1.0) fNSPDTracklet10++;
  }

  const AliVVertex *vVertex = fEvent->GetPrimaryVertex();
  fVertex[0] = vVertex->GetX();
  fVertex[1] = vVertex->GetY();
  fVertex[2] = vVertex->GetZ();
  fNContributor = vVertex->GetNContributors();
  if(fESDEvent)      fNTPCCluster = fESDEvent->GetNumberOfTPCClusters();
  else if(fAODEvent) fNTPCCluster = fAODEvent->GetNumberOfTPCClusters();

  const Int_t Ntrack = fEvent->GetNumberOfTracks();
  fNHybridTrack08 = 0;
  if(fESDEvent){
    AliInfo("So far, ESD is not supported. return.");
  }//end of ESD
  else if(fAODEvent){
    for(Int_t itrack=0;itrack<Ntrack;itrack++){
      AliVTrack* track = dynamic_cast<AliVTrack*>(fEvent->GetTrack(itrack));
      if(track->Pt() < 0.15) continue;
      if(TMath::Abs(track->Eta()) > 0.8) continue;

      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(track);
      if(aodtrack->IsHybridGlobalConstrainedGlobal()) fNHybridTrack08++;//hybrid track (global + complementary)
    }//end of track loop
  }//end of AOD

  fNTrackTPCout = 0;
  Float_t DCAxy = -999, DCAz = -999;
  Float_t Chi2Global = -1;
  Float_t TOFbeta = -999;
  Double32_t expt[5] = {0,0,0,0,0};

  for(Int_t itrack=0;itrack<Ntrack;itrack++){
    AliVTrack* track = dynamic_cast<AliVTrack*>(fEvent->GetTrack(itrack));
    ULong64_t status = track->GetStatus();
    if(status & AliVTrack::kTPCout) fNTrackTPCout++;

    if(track->Pt() < fMinPtCut) continue;
    if(TMath::Abs(track->Eta()) > fMaxEtaCut) continue;

    if(fESDEvent){
      AliInfo("So far, ESD is not supported. return.");
    }//end of ESD
    else if(fAODEvent){
      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(track);
      if(!aodtrack->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA)) continue;//standard cuts with very loose DCA cut //bit4

      aodtrack->GetImpactParameters(DCAxy,DCAz);
      if(TMath::Abs(DCAxy) > 1.) continue;
      if(TMath::Abs(DCAz) > 3.) continue;
      Chi2Global = aodtrack->GetChi2TPCConstrainedVsGlobal();
    }//end of AOD

    //if(!track->HasPointOnITSLayer(0)) continue;//require a hit on first SPD layer

    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron)) > fMaxTPCNsigmaEleCut) continue;//pre-select electrons to reduce data size.

    if(Chi2Global > 100) continue;
    if(track->GetNcls(0) < 4)  continue;//minimum number of ITS cluster 4
    if(track->GetNcls(1) < 70) continue;//minimum number of TPC cluster 70
    if((Double_t)(track->GetTPCchi2()) / (Double_t)(track->GetNcls(1)) > 4.) continue;//maximum chi2 per cluster TPC
    if((Double_t)(track->GetITSchi2()) / (Double_t)(track->GetNcls(0)) > 5.) continue;//maximum chi2 per cluster ITS

    vector<Float_t> vec(3,0);
    vec[0] = track->Pt();
    vec[1] = track->Eta();
    vec[2] = track->Phi();
    fTrackMomentum.push_back(vec);
    vec.clear();

    fTrackCharge.push_back(track->Charge());
    fTrackDCAxy.push_back(DCAxy);
    fTrackDCAz.push_back(DCAz);

    vector<Bool_t> PointITSLayer_tmp(6,kFALSE);
    for(Int_t il=0;il<6;il++) PointITSLayer_tmp[il] = track->HasPointOnITSLayer(il);
    fPointOnITSLayer.push_back(PointITSLayer_tmp);
    PointITSLayer_tmp.clear();

    vector<Bool_t> SharedPointITS_tmp(6,kFALSE);
    for(Int_t il=0;il<6;il++) SharedPointITS_tmp[il] = track->HasSharedPointOnITSLayer(il);
    fSharedPointOnITSLayer.push_back(SharedPointITS_tmp);
    SharedPointITS_tmp.clear();

    fTrackPin.push_back(track->GetTPCmomentum());
    fChi2TPC.push_back(track->GetTPCchi2());
    fChi2ITS.push_back(track->GetITSchi2());
    fNclusterTPC.push_back(track->GetNcls(1));
    fNclusterITS.push_back(track->GetNcls(0));
    fTPCNCrossedRows.push_back(track->GetTPCCrossedRows());
    fTPCNFindableCluster.push_back(track->GetTPCNclsF());
    fChi2TPCConstrainedVsGlobal.push_back(Chi2Global);

    //TPC PID info
    fTPCsignal.push_back(track->GetTPCsignal());
    fTPCNsigmaEl.push_back(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
    fTPCNsigmaPi.push_back(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
    fTPCNsigmaKa.push_back(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon));
    fTPCNsigmaPr.push_back(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton));

    //ITS PID info
    fITSsignal.push_back(track->GetITSsignal());
    fITSNsigmaEl.push_back(fPIDResponse->NumberOfSigmasITS(track,AliPID::kElectron));
    fITSNsigmaPi.push_back(fPIDResponse->NumberOfSigmasITS(track,AliPID::kPion));
    fITSNsigmaKa.push_back(fPIDResponse->NumberOfSigmasITS(track,AliPID::kKaon));
    fITSNsigmaPr.push_back(fPIDResponse->NumberOfSigmasITS(track,AliPID::kProton));

    //TOF info
    Bool_t isTIME   = status & AliVTrack::kTIME;
    Bool_t isTOFout = status & AliVTrack::kTOFout;
    Bool_t isTOFOK  = isTIME & isTOFout;
    track->GetIntegratedTimes(expt); 
    Double_t length = TMath::C() * expt[0] * 1e-12;// m
    Double_t time   = track->GetTOFsignal();       //ps
    time -= fPIDResponse->GetTOFResponse().GetStartTime(track->P()); // ps

    TOFbeta = -999;
    if( 
        !isTOFOK
        //time <= 0.         //time can not be measured
        //|| length < 360.e-2//too short distance
        //|| length > 800.e-2//too long distance
      ){
      TOFbeta = -999;
    }
    else{
      time *= 1e-12;//ps -> s
      Double_t velocity = length / time;
      TOFbeta = velocity / TMath::C();
    }

    fTOFbeta.push_back(TOFbeta);
    fTOFNsigmaEl.push_back(fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
    fTOFNsigmaPi.push_back(fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion));
    fTOFNsigmaKa.push_back(fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon));
    fTOFNsigmaPr.push_back(fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton));
    fIsTOFAvailable.push_back(isTOFOK);

    //MC info for reconstructed tracks
    if(hasMC){
      Int_t label = TMath::Abs(track->GetLabel());
      AliAODMCParticle *p = (AliAODMCParticle*)fMCArray->At(label);
      Int_t pdg = p->GetPdgCode();
      fTrackMCIndex.push_back(label);
      fTrackMCPdgCode.push_back(pdg);


      vector<Float_t>vec_mc(3,0);
      vec_mc[0] = p->Pt();
      vec_mc[1] = p->Eta();
      vec_mc[2] = p->Phi();
      fTrackMCMomentum.push_back(vec_mc);
      vec_mc.clear();

      vector<Float_t> mcvtx(3,0);
      mcvtx[0] = p->Xv();
      mcvtx[1] = p->Yv();
      mcvtx[2] = p->Zv();
      fTrackMCProdVtx.push_back(mcvtx);
      mcvtx.clear();

      //check mother
      Int_t mother_index = p->GetMother();
      Int_t mother_pdg   = 0;//pdgcode 0 does not exist.
      if(mother_index > -1){
        AliAODMCParticle *mp = (AliAODMCParticle*)fMCArray->At(mother_index);
        mother_pdg = mp->GetPdgCode();
        fTrackMCMotherIndex.push_back(mother_index);
        fTrackMCMotherPdgCode.push_back(mother_pdg);
        //printf("mother_index = %d , mother_pdg = %d\n",mother_index,mother_pdg);

        //check first mother
        Int_t first_mother_index     = mp->GetMother();
        Int_t first_mother_index_tmp = first_mother_index;
        Int_t first_mother_pdg   = 0;//pdgcode 0 does not exist.
        while(first_mother_index > -1){
          first_mother_index_tmp = first_mother_index;
          AliAODMCParticle *fmp = (AliAODMCParticle*)fMCArray->At(first_mother_index);
          first_mother_index = fmp->GetMother();
          first_mother_pdg   = fmp->GetPdgCode();
        }//end of mother loop
        fTrackMCFirstMotherIndex.push_back(first_mother_index_tmp);
        fTrackMCFirstMotherPdgCode.push_back(first_mother_pdg);
        //printf("first_mother_index_tmp = %d , first_mother_pdg = %d\n",first_mother_index_tmp,first_mother_pdg);
      }
      else{//no mother, i.e. no first mother, too.
        fTrackMCMotherIndex.push_back(mother_index);
        fTrackMCMotherPdgCode.push_back(mother_pdg);//pdgcode 0 does not exist.

        fTrackMCFirstMotherIndex.push_back(mother_index);
        fTrackMCFirstMotherPdgCode.push_back(mother_pdg);//pdgcode 0 does not exist.
      }
    }//end of hasMC

  }//end of track loop

  AliInfo(Form("fNSPDTracklet05 = %d , fNSPDTracklet10 = %d , fNHybridTrack08 = %d , fNTPCCluster = %d , fNTrackTPCout = %d",fNSPDTracklet05,fNSPDTracklet10,fNHybridTrack08,fNTPCCluster,fNTrackTPCout));


  AliKFVertex primaryVertexKF(*vVertex);
  Double_t secVtx[3] = {primaryVertexKF.GetX(), primaryVertexKF.GetY(), primaryVertexKF.GetZ()};
  Float_t dca = 0;
  Float_t KFchi2ndf = 999;

  //FillV0Info//select gamma conversion
  const Int_t Nv0 = fEvent->GetNumberOfV0s();  
  for(Int_t iv0=0;iv0<Nv0;iv0++){
    AliAODv0 *v0 = (AliAODv0*)fAODEvent->GetV0(iv0);

    if(!v0->GetOnFlyStatus()) continue;//select v0 reconstructed on the fly.

    dca = v0->DcaV0Daughters();
    if(dca > 0.25) continue;

    if(v0->RadiusV0() < 3.) continue;//in cm
    if(v0->RadiusV0() > 60.) continue;//in cm

    if(v0->ChargeProng(0) * v0->ChargeProng(1) > 0) continue;//reject same sign pair
    AliAODTrack *legPos = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(0));
    AliAODTrack *legNeg = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(1));
    //if(legPos->Charge() < 0 && legNeg->Charge() > 0){//swap charge sign
    if(v0->ChargeProng(0) < 0 && v0->ChargeProng(0) > 0){//swap charge sign //index0 is expect to be positive leg, index1 to be negative.//protection
      AliInfo("charge is swapped.");
      legPos = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(1));
      legNeg = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(0));
    }

    if(legPos->Pt() < fMinPtCut) continue;
    if(legNeg->Pt() < fMinPtCut) continue;

    Float_t DCAxy_leg = -999, DCAz_leg = -999;
    legPos->GetImpactParameters(DCAxy_leg,DCAz_leg);
    if(TMath::Abs(DCAxy_leg) > 1.) continue;
    if(TMath::Abs(DCAz_leg) > 3.) continue;

    DCAxy_leg = -999; DCAz_leg = -999;
    legNeg->GetImpactParameters(DCAxy_leg,DCAz_leg);
    if(TMath::Abs(DCAxy_leg) > 1.) continue;
    if(TMath::Abs(DCAz_leg) > 3.) continue;

    if((Double_t)(legPos->GetTPCchi2()) / (Double_t)(legPos->GetNcls(1)) > 4.) continue;//maximum chi2 per cluster TPC
    if((Double_t)(legNeg->GetTPCchi2()) / (Double_t)(legNeg->GetNcls(1)) > 4.) continue;//maximum chi2 per cluster TPC

    if(TMath::Abs(legPos->Eta()) > fMaxEtaCut) continue;
    if(TMath::Abs(legNeg->Eta()) > fMaxEtaCut) continue;


    if(legPos->GetNcls(1) < 70) continue;//minimum number of TPC cluster 70
    if(legNeg->GetNcls(1) < 70) continue;//minimum number of TPC cluster 70

    if(!(legPos->GetStatus() & AliVTrack::kTPCrefit)) continue;
    if(!(legNeg->GetStatus() & AliVTrack::kTPCrefit)) continue;

    Float_t ratio_pos = legPos->GetTPCNclsF() > 0 ? (Float_t)legPos->GetTPCCrossedRows() / (Float_t)legPos->GetTPCNclsF() : -1;
    Float_t ratio_neg = legNeg->GetTPCNclsF() > 0 ? (Float_t)legNeg->GetTPCCrossedRows() / (Float_t)legNeg->GetTPCNclsF() : -1;

    if(ratio_pos < 0.8) continue;
    if(ratio_neg < 0.8) continue;

    if(v0->PtArmV0() > 0.05) continue;//qT < 0.05
    if(TMath::Abs(v0->AlphaV0()) > 1.) continue;//|alpha| < 1

    AliDielectronPair* DielePair = new AliDielectronPair();
    DielePair->SetKFUsage(kTRUE); 
    DielePair->SetTracks(&(*static_cast<AliVTrack*>(legPos)),11, &(*static_cast<AliVTrack*>(legNeg)), 11);
    KFchi2ndf = DielePair->GetKFNdf() > 0 ? DielePair->GetKFChi2()/DielePair->GetKFNdf() : 999;
    delete DielePair;
    DielePair = 0x0;
    if(KFchi2ndf > 30) continue;

    if(v0->CosPointingAngle(secVtx) < 0.99) continue;

    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kElectron)) > 5.) continue;
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kElectron)) > 5.) continue;
    //if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kElectron)) > fMaxTPCNsigmaEleCut) continue;
    //if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kElectron)) > fMaxTPCNsigmaEleCut) continue;

    vector<Float_t> V0vec(3,0);
    V0vec[0] = v0->Pt();
    V0vec[1] = v0->Eta();
    V0vec[2] = v0->Phi();
    fV0Momentum.push_back(V0vec);
    V0vec.clear();

    vector<vector<Float_t>> legP_vec_tmp;//0 for legPos, 1 for legNeg
    vector<Float_t>legPos_vec(3,0);
    legPos_vec[0] = legPos->Pt();
    legPos_vec[1] = legPos->Eta();
    legPos_vec[2] = legPos->Phi();
    legP_vec_tmp.push_back(legPos_vec); 
    vector<Float_t>legNeg_vec(3,0);
    legNeg_vec[0] = legNeg->Pt();
    legNeg_vec[1] = legNeg->Eta();
    legNeg_vec[2] = legNeg->Phi();
    legP_vec_tmp.push_back(legNeg_vec); 
    fV0legMomentum.push_back(legP_vec_tmp);
    legP_vec_tmp.clear();

    fV0DCA.push_back(dca);
    fV0Lxy.push_back(v0->RadiusV0());
    
    vector<Float_t> legPin_tmp(2,0);
    legPin_tmp[0] = legPos->GetTPCmomentum();
    legPin_tmp[1] = legNeg->GetTPCmomentum();
    fV0legPin.push_back(legPin_tmp);
    legPin_tmp.clear();

    vector<Float_t> legChi2Global_tmp(2,0);
    legChi2Global_tmp[0] = legPos->GetChi2TPCConstrainedVsGlobal();
    legChi2Global_tmp[1] = legNeg->GetChi2TPCConstrainedVsGlobal();
    fV0legChi2TPCConstrainedVsGlobal.push_back(legChi2Global_tmp);
    legChi2Global_tmp.clear();

    fV0PointingAngle.push_back(v0->CosPointingAngle(secVtx));
    fV0Chi2.push_back(KFchi2ndf);

    fV0alpha.push_back(v0->AlphaV0());
    fV0qT.push_back(v0->PtArmV0());

    fV0PsiPair.push_back(PsiPair(v0,fEvent->GetMagneticField()));
    fV0PhivPair.push_back(PhivPair(v0,fEvent->GetMagneticField()));

    //info for ITS point
    vector<vector<Bool_t>> legPointOnITS_tmp;//2 x 6 elements
    for(Int_t ileg=0;ileg<2;ileg++){
      vector<Bool_t> tmp2;
      for(Int_t ilayer=0;ilayer<6;ilayer++){
        tmp2.push_back(kFALSE);
      }//end of layer loop
      legPointOnITS_tmp.push_back(tmp2);
    }//end of leg loop

    for(Int_t ilayer=0;ilayer<6;ilayer++){
      legPointOnITS_tmp[0][ilayer] = legPos->HasPointOnITSLayer(ilayer);
      legPointOnITS_tmp[1][ilayer] = legNeg->HasPointOnITSLayer(ilayer);
    }//end of layer loop
 
    fV0PointOnITSLayer.push_back(legPointOnITS_tmp);
    legPointOnITS_tmp.clear();

    //info for ITS shared point
    vector<vector<Bool_t>> legSharedPointOnITS_tmp;//2 x 6 elements
    for(Int_t ileg=0;ileg<2;ileg++){
      vector<Bool_t> tmp2;
      for(Int_t ilayer=0;ilayer<6;ilayer++){
        tmp2.push_back(kFALSE);
      }//end of layer loop
      legSharedPointOnITS_tmp.push_back(tmp2);
    }//end of leg loop

    for(Int_t ilayer=0;ilayer<6;ilayer++){
      legSharedPointOnITS_tmp[0][ilayer] = legPos->HasSharedPointOnITSLayer(ilayer);
      legSharedPointOnITS_tmp[1][ilayer] = legNeg->HasSharedPointOnITSLayer(ilayer);
    }//end of layer loop
 
    fV0SharedPointOnITSLayer.push_back(legSharedPointOnITS_tmp);
    legSharedPointOnITS_tmp.clear();

    vector<Float_t> mass_tmp(4,0);
    mass_tmp[0] = v0->InvMass2Prongs(0,1,11,11);
    mass_tmp[1] = v0->MassK0Short();
    mass_tmp[2] = v0->MassLambda();
    mass_tmp[3] = v0->MassAntiLambda();
    fV0Mass.push_back(mass_tmp);
    mass_tmp.clear();

    vector<Float_t> legDCAxy_tmp(2,999);
    vector<Float_t> legDCAz_tmp(2,999);
    legPos->GetImpactParameters(legDCAxy_tmp[0],legDCAz_tmp[0]);
    legNeg->GetImpactParameters(legDCAxy_tmp[1],legDCAz_tmp[1]);
    fV0legDCAxy.push_back(legDCAxy_tmp);
    fV0legDCAz.push_back(legDCAz_tmp);
    legDCAxy_tmp.clear();
    legDCAz_tmp.clear();

    vector<Float_t> v0_leg_TPCNsigmaEl_tmp(2,-999);//0 for legPos, 1 for legNeg
    v0_leg_TPCNsigmaEl_tmp[0] = fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kElectron);
    v0_leg_TPCNsigmaEl_tmp[1] = fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kElectron);
    fV0legTPCNsigmaEl.push_back(v0_leg_TPCNsigmaEl_tmp);
    v0_leg_TPCNsigmaEl_tmp.clear();

    vector<Float_t> v0_leg_TPCNsigmaPi_tmp(2,-999);//0 for legPos, 1 for legNeg
    v0_leg_TPCNsigmaPi_tmp[0] = fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kPion);
    v0_leg_TPCNsigmaPi_tmp[1] = fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kPion);
    fV0legTPCNsigmaPi.push_back(v0_leg_TPCNsigmaPi_tmp);
    v0_leg_TPCNsigmaPi_tmp.clear();

    vector<Float_t> v0_leg_TPCNsigmaKa_tmp(2,-999);//0 for legPos, 1 for legNeg
    v0_leg_TPCNsigmaKa_tmp[0] = fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kKaon);
    v0_leg_TPCNsigmaKa_tmp[1] = fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kKaon);
    fV0legTPCNsigmaKa.push_back(v0_leg_TPCNsigmaKa_tmp);
    v0_leg_TPCNsigmaKa_tmp.clear();

    vector<Float_t> v0_leg_TPCNsigmaPr_tmp(2,-999);//0 for legPos, 1 for legNeg
    v0_leg_TPCNsigmaPr_tmp[0] = fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kProton);
    v0_leg_TPCNsigmaPr_tmp[1] = fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kProton);
    fV0legTPCNsigmaPr.push_back(v0_leg_TPCNsigmaPr_tmp);
    v0_leg_TPCNsigmaPr_tmp.clear();

    vector<Float_t> v0_leg_ITSNsigmaEl_tmp(2,-999);//0 for legPos, 1 for legNeg
    v0_leg_ITSNsigmaEl_tmp[0] = fPIDResponse->NumberOfSigmasITS(legPos,AliPID::kElectron);
    v0_leg_ITSNsigmaEl_tmp[1] = fPIDResponse->NumberOfSigmasITS(legNeg,AliPID::kElectron);
    fV0legITSNsigmaEl.push_back(v0_leg_ITSNsigmaEl_tmp);
    v0_leg_ITSNsigmaEl_tmp.clear();

    vector<Float_t> v0_leg_ITSNsigmaPi_tmp(2,-999);//0 for legPos, 1 for legNeg
    v0_leg_ITSNsigmaPi_tmp[0] = fPIDResponse->NumberOfSigmasITS(legPos,AliPID::kPion);
    v0_leg_ITSNsigmaPi_tmp[1] = fPIDResponse->NumberOfSigmasITS(legNeg,AliPID::kPion);
    fV0legITSNsigmaPi.push_back(v0_leg_ITSNsigmaPi_tmp);
    v0_leg_ITSNsigmaPi_tmp.clear();

    vector<Float_t> v0_leg_ITSNsigmaKa_tmp(2,-999);//0 for legPos, 1 for legNeg
    v0_leg_ITSNsigmaKa_tmp[0] = fPIDResponse->NumberOfSigmasITS(legPos,AliPID::kKaon);
    v0_leg_ITSNsigmaKa_tmp[1] = fPIDResponse->NumberOfSigmasITS(legNeg,AliPID::kKaon);
    fV0legITSNsigmaKa.push_back(v0_leg_ITSNsigmaKa_tmp);
    v0_leg_ITSNsigmaKa_tmp.clear();

    vector<Float_t> v0_leg_ITSNsigmaPr_tmp(2,-999);//0 for legPos, 1 for legNeg
    v0_leg_ITSNsigmaPr_tmp[0] = fPIDResponse->NumberOfSigmasITS(legPos,AliPID::kProton);
    v0_leg_ITSNsigmaPr_tmp[1] = fPIDResponse->NumberOfSigmasITS(legNeg,AliPID::kProton);
    fV0legITSNsigmaPr.push_back(v0_leg_ITSNsigmaPr_tmp);
    v0_leg_ITSNsigmaPr_tmp.clear();

    ULong64_t status1 = legPos->GetStatus();
    ULong64_t status2 = legNeg->GetStatus();

    Bool_t isTIME1 = status1 & AliVTrack::kTIME;
    Bool_t isTIME2 = status2 & AliVTrack::kTIME;

    Bool_t isTOFout1 = status1 & AliVTrack::kTOFout;
    Bool_t isTOFout2 = status2 & AliVTrack::kTOFout;

    Bool_t isTOFOK1 = isTIME1 & isTOFout1;
    Bool_t isTOFOK2 = isTIME2 & isTOFout2;

    vector<Bool_t> v0_leg_isTOFOK_tmp(2,kFALSE);//0 for legPos, 1 for legNeg
    v0_leg_isTOFOK_tmp[0] = isTOFOK1; 
    v0_leg_isTOFOK_tmp[1] = isTOFOK2; 
    fV0legIsTOFAvailable.push_back(v0_leg_isTOFOK_tmp);
    v0_leg_isTOFOK_tmp.clear();

    vector<Float_t> v0_leg_TOFNsigmaEl_tmp(2,-999);//0 for legPos, 1 for legNeg
    vector<Float_t> v0_leg_TOFNsigmaPi_tmp(2,-999);//0 for legPos, 1 for legNeg
    vector<Float_t> v0_leg_TOFNsigmaKa_tmp(2,-999);//0 for legPos, 1 for legNeg
    vector<Float_t> v0_leg_TOFNsigmaPr_tmp(2,-999);//0 for legPos, 1 for legNeg

    v0_leg_TOFNsigmaEl_tmp[0] = fPIDResponse->NumberOfSigmasTOF(legPos,AliPID::kElectron);
    v0_leg_TOFNsigmaPi_tmp[0] = fPIDResponse->NumberOfSigmasTOF(legPos,AliPID::kPion);
    v0_leg_TOFNsigmaKa_tmp[0] = fPIDResponse->NumberOfSigmasTOF(legPos,AliPID::kKaon);
    v0_leg_TOFNsigmaPr_tmp[0] = fPIDResponse->NumberOfSigmasTOF(legPos,AliPID::kProton);
    v0_leg_TOFNsigmaEl_tmp[1] = fPIDResponse->NumberOfSigmasTOF(legNeg,AliPID::kElectron);
    v0_leg_TOFNsigmaPi_tmp[1] = fPIDResponse->NumberOfSigmasTOF(legNeg,AliPID::kPion);
    v0_leg_TOFNsigmaKa_tmp[1] = fPIDResponse->NumberOfSigmasTOF(legNeg,AliPID::kKaon);
    v0_leg_TOFNsigmaPr_tmp[1] = fPIDResponse->NumberOfSigmasTOF(legNeg,AliPID::kProton);

    fV0legTOFNsigmaEl.push_back(v0_leg_TOFNsigmaEl_tmp);
    v0_leg_TOFNsigmaEl_tmp.clear();

    fV0legTOFNsigmaPi.push_back(v0_leg_TOFNsigmaPi_tmp);
    v0_leg_TOFNsigmaPi_tmp.clear();

    fV0legTOFNsigmaKa.push_back(v0_leg_TOFNsigmaKa_tmp);
    v0_leg_TOFNsigmaKa_tmp.clear();

    fV0legTOFNsigmaPr.push_back(v0_leg_TOFNsigmaPr_tmp);
    v0_leg_TOFNsigmaPr_tmp.clear();

    //MC info
    if(hasMC){
      vector<Int_t> leglabel_tmp(2,0);
      leglabel_tmp[0] = TMath::Abs(legPos->GetLabel());
      leglabel_tmp[1] = TMath::Abs(legNeg->GetLabel());
      AliAODMCParticle *pPos = (AliAODMCParticle*)fMCArray->At(leglabel_tmp[0]);
      AliAODMCParticle *pNeg = (AliAODMCParticle*)fMCArray->At(leglabel_tmp[1]);
      fV0MClegIndex.push_back(leglabel_tmp);
      leglabel_tmp.clear();

      vector<Int_t> legpdg_tmp(2,0);
      legpdg_tmp[0] = pPos->GetPdgCode();
      legpdg_tmp[1] = pNeg->GetPdgCode();
      fV0MClegPdgCode.push_back(legpdg_tmp);
      legpdg_tmp.clear();

      vector<vector<Float_t>> legP_vec_mc_tmp;//0 for legPos, 1 for legNeg
      vector<Float_t>legPos_vec_mc(3,0);
      legPos_vec_mc[0] = pPos->Pt();
      legPos_vec_mc[1] = pPos->Eta();
      legPos_vec_mc[2] = pPos->Phi();
      legP_vec_mc_tmp.push_back(legPos_vec_mc); 
      vector<Float_t>legNeg_vec_mc(3,0);
      legNeg_vec_mc[0] = pNeg->Pt();
      legNeg_vec_mc[1] = pNeg->Eta();
      legNeg_vec_mc[2] = pNeg->Phi();
      legP_vec_mc_tmp.push_back(legNeg_vec_mc); 
      fV0MClegMomentum.push_back(legP_vec_mc_tmp);
      legP_vec_mc_tmp.clear();

      vector<vector<Float_t>> leg_vec_pv_mc_tmp;//0 for legPos, 1 for legNeg//production vertex
      vector<Float_t> legPos_vec_pv_mc(3,0);
      legPos_vec_pv_mc[0] = pPos->Xv();
      legPos_vec_pv_mc[1] = pPos->Yv();
      legPos_vec_pv_mc[2] = pPos->Zv();
      leg_vec_pv_mc_tmp.push_back(legPos_vec_pv_mc); 
      vector<Float_t> legNeg_vec_pv_mc(3,0);
      legNeg_vec_pv_mc[0] = pNeg->Xv();
      legNeg_vec_pv_mc[1] = pNeg->Yv();
      legNeg_vec_pv_mc[2] = pNeg->Zv();
      leg_vec_pv_mc_tmp.push_back(legNeg_vec_pv_mc); 
      fV0MClegProdVtx.push_back(leg_vec_pv_mc_tmp);
      leg_vec_pv_mc_tmp.clear();

      //check mother
      vector<Int_t> leg_mother_index_tmp(2,0);
      leg_mother_index_tmp[0] = pPos->GetMother();
      leg_mother_index_tmp[1] = pNeg->GetMother();

      vector<Int_t> leg_mother_pdg_tmp(2,0);
      leg_mother_pdg_tmp[0] = 0;
      leg_mother_pdg_tmp[1] = 0;

      vector<Int_t> leg_first_mother_index_tmp(2,0);
      leg_first_mother_index_tmp[0] = -1;
      leg_first_mother_index_tmp[1] = -1;

      vector<Int_t> leg_first_mother_pdg_tmp(2,0);
      leg_first_mother_pdg_tmp[0] = 0;
      leg_first_mother_pdg_tmp[1] = 0;

      if(leg_mother_index_tmp[0] > -1){//check mother for pos leg
        AliAODMCParticle *mpPos = (AliAODMCParticle*)fMCArray->At(leg_mother_index_tmp[0]);
        leg_mother_pdg_tmp[0] = mpPos->GetPdgCode();

        Int_t first_mother_index_tmp = mpPos->GetMother();
        leg_first_mother_index_tmp[0] = first_mother_index_tmp;
        while(first_mother_index_tmp > -1){
          leg_first_mother_index_tmp[0] = first_mother_index_tmp;
          AliAODMCParticle *fmpPos = (AliAODMCParticle*)fMCArray->At(first_mother_index_tmp);
          first_mother_index_tmp      = fmpPos->GetMother();
          leg_first_mother_pdg_tmp[0] = fmpPos->GetPdgCode();
        }//end of mother loop

      }//end of check mother for pos leg

      if(leg_mother_index_tmp[1] > -1){//check mother for neg leg
        AliAODMCParticle *mpNeg = (AliAODMCParticle*)fMCArray->At(leg_mother_index_tmp[1]);
        leg_mother_pdg_tmp[1] = mpNeg->GetPdgCode();

        Int_t first_mother_index_tmp = mpNeg->GetMother();
        leg_first_mother_index_tmp[1] = first_mother_index_tmp;
        while(first_mother_index_tmp > -1){
          leg_first_mother_index_tmp[1] = first_mother_index_tmp;
          AliAODMCParticle *fmpNeg = (AliAODMCParticle*)fMCArray->At(first_mother_index_tmp);
          first_mother_index_tmp      = fmpNeg->GetMother();
          leg_first_mother_pdg_tmp[1] = fmpNeg->GetPdgCode();
        }//end of mother loop

      }//end of check mother for pos leg

      fV0MClegMotherIndex.push_back(leg_mother_index_tmp);
      fV0MClegMotherPdgCode.push_back(leg_mother_pdg_tmp);
      fV0MClegFirstMotherIndex.push_back(leg_first_mother_index_tmp);
      fV0MClegFirstMotherPdgCode.push_back(leg_first_mother_pdg_tmp);

      leg_mother_index_tmp.clear();
      leg_mother_pdg_tmp.clear();
      leg_first_mother_index_tmp.clear();
      leg_first_mother_pdg_tmp.clear();

    }//end of hasMC

  }//end of V0 loop

  if(hasMC) ProcessMC(option);

  fTree->Fill();
  ClearVectorMemory();

}
//_______________________________________________________________________________________________
void AliAnalysisTaskReducedTreeDS::ProcessMC(Option_t *option)
{

  const Int_t Ntrack = fMCArray->GetEntries();
  AliInfo(Form("Ntrack in MC = %d",Ntrack));

  AliMCEvent* MCevent = dynamic_cast<AliMCEvent*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->MCEvent());//MC truth info

  const AliVVertex *MCvertex = MCevent->GetPrimaryVertex();
  fMCVertex[0] = MCvertex->GetX();//true vertex
  fMCVertex[1] = MCvertex->GetY();//true vertex
  fMCVertex[2] = MCvertex->GetZ();//true vertex

  Int_t pdg = -1;
  Float_t pT = 0, eta = 0, phi = 0;

  for(Int_t itrack=0;itrack<Ntrack;itrack++){
    AliAODMCParticle *p = (AliAODMCParticle*)fMCArray->At(itrack);
    pdg = p->GetPdgCode();
    pT  = p->Pt();
    eta = p->Eta();//pseudo-rapidity
    phi = p->Phi();

    if(pT < fMinPtCut) continue;
    if(TMath::Abs(eta) > fMaxEtaCut) continue;//select only electrons in |eta| < 0.8
    if(TMath::Abs(pdg) != 11) continue; //select only electrons

    Double_t dx = p->Xv() - fMCVertex[0];
    Double_t dy = p->Yv() - fMCVertex[1];
    //Double_t dz = p->Zv() - fMCVertex[2];
    Double_t R   = TMath::Sqrt(dx*dx + dy*dy);
    //Double_t Rho = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
    if(R > 1.) continue;//select electrons from primary vertex.

    vector<Float_t> prodvtx(3,0);
    prodvtx[0] = p->Xv();
    prodvtx[1] = p->Yv();
    prodvtx[2] = p->Zv();
    fMCProdVtx.push_back(prodvtx);
    prodvtx.clear();

    vector<Float_t> vec(3,0);
    vec[0] = pT;
    vec[1] = eta;
    vec[2] = phi;
    fMCMomentum.push_back(vec);
    vec.clear();

    fMCIndex.push_back(itrack);
    fMCPdgCode.push_back(pdg);//11 or -11

    //check mother
    Int_t mother_index = p->GetMother();
    Int_t mother_pdg   = 0;//pdgcode 0 does not exist.
    if(mother_index > -1){
      AliAODMCParticle *mp = (AliAODMCParticle*)fMCArray->At(mother_index);
      mother_pdg = mp->GetPdgCode();
      fMCMotherIndex.push_back(mother_index);
      fMCMotherPdgCode.push_back(mother_pdg);
      //printf("mother_index = %d , mother_pdg = %d\n",mother_index,mother_pdg);

      //check first mother
      Int_t first_mother_index     = mp->GetMother();
      Int_t first_mother_index_tmp = first_mother_index;
      Int_t first_mother_pdg   = 0;//pdgcode 0 does not exist.
      while(first_mother_index > -1){
        first_mother_index_tmp = first_mother_index;
        AliAODMCParticle *fmp = (AliAODMCParticle*)fMCArray->At(first_mother_index);
        first_mother_index = fmp->GetMother();
        first_mother_pdg   = fmp->GetPdgCode();
      }//end of mother loop
      fMCFirstMotherIndex.push_back(first_mother_index_tmp);
      fMCFirstMotherPdgCode.push_back(first_mother_pdg);
      //printf("first_mother_index_tmp = %d , first_mother_pdg = %d\n",first_mother_index_tmp,first_mother_pdg);
    }
    else{//no mother, i.e. no first mother, too.
      fMCMotherIndex.push_back(mother_index);
      fMCMotherPdgCode.push_back(mother_pdg);//pdgcode 0 does not exist.

      fMCFirstMotherIndex.push_back(mother_index);
      fMCFirstMotherPdgCode.push_back(mother_pdg);//pdgcode 0 does not exist.
    }


  }//end of MC track loop

}
//_______________________________________________________________________________________________
void AliAnalysisTaskReducedTreeDS::Terminate(Option_t *option)
{
  //Called once at the end of the query
  AliInfo(Form("%s is done.",GetName()));
}
//_______________________________________________________________________________________________
Double_t AliAnalysisTaskReducedTreeDS::PsiPair(AliAODv0 *v0, Float_t Bz)
{
  AliAODTrack *legPos = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(0));
  AliAODTrack *legNeg = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(1));

  AliExternalTrackParam pt; pt.CopyFromVTrack(legPos);
  AliExternalTrackParam nt; nt.CopyFromVTrack(legNeg);

  // Cartesian Momentum for each leg
  Double_t xyz[3] = {0.,0.,0.};
  v0->GetXYZ(xyz);
  const AliVVertex *vVertex = fEvent->GetPrimaryVertex();
  xyz[0] -= vVertex->GetX();
  xyz[1] -= vVertex->GetY();
  xyz[2] -= vVertex->GetZ();
  
  Double_t mn[3] = {0,0,0};
  mn[0] = v0->MomNegX();
  mn[1] = v0->MomNegY();
  mn[2] = v0->MomNegZ();
  Double_t mp[3] = {0,0,0};
  mp[0] = v0->MomPosX();
  mp[1] = v0->MomPosY();
  mp[2] = v0->MomPosZ();

  Double_t deltat = 1.;
  deltat = TMath::ATan(mp[2]/(TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1])+1.e-13)) - TMath::ATan(mn[2]/(TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1])+1.e-13));//difference of angles of the two daughter tracks with z-axis
  Double_t radiussum = TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]) + 50;//radius to which tracks shall be propagated

  Double_t momPosProp[3] = {0,0,0};
  Double_t momNegProp[3] = {0,0,0};
  Double_t psiPair = 4.;
  if(nt.PropagateTo(radiussum,Bz) == 0) psiPair = -5; //propagate tracks to the outside
  if(pt.PropagateTo(radiussum,Bz) == 0) psiPair = -5; //propagate tracks to the outside
  pt.GetPxPyPz(momPosProp);//Get momentum vectors of tracks after propagation
  nt.GetPxPyPz(momNegProp);

  Double_t pEle = TMath::Sqrt(momNegProp[0]*momNegProp[0]+momNegProp[1]*momNegProp[1]+momNegProp[2]*momNegProp[2]);//absolute momentum value of negative daughter
  Double_t pPos = TMath::Sqrt(momPosProp[0]*momPosProp[0]+momPosProp[1]*momPosProp[1]+momPosProp[2]*momPosProp[2]);//absolute momentum value of positive daughter
  Double_t scalarproduct = momPosProp[0]*momNegProp[0]+momPosProp[1]*momNegProp[1]+momPosProp[2]*momNegProp[2];//scalar product of propagated positive and negative daughters' momenta
  Double_t chipair = TMath::ACos(scalarproduct/(pEle*pPos));//Angle between propagated daughter tracks

  psiPair =  TMath::Abs(TMath::ASin(deltat/chipair));
  return psiPair;
}
//_______________________________________________________________________________________________
Double_t AliAnalysisTaskReducedTreeDS::PhivPair(AliAODv0 *v0, Float_t Bz)
{
  /// Following the idea to use opening of collinear pairs in magnetic field from e.g. PHENIX
  /// to identify conversions. Angle between ee plane and magnetic field is calculated (0 to pi).
  /// Due to tracking to the primary vertex, conversions with no intrinsic opening angle
  /// always end up as pair in "cowboy" configuration. The function as defined here then
  /// returns values close to pi.
  /// Correlated Like Sign pairs (from double conversion / dalitz + conversion) may show up
  /// at pi or at 0 depending on which leg has the higher momentum. (not checked yet)
  /// This expected ambiguity is not seen due to sorting of track arrays in this framework.
  /// To reach the same result as for ULS (~pi), the legs are flipped for LS.

  TVector3 p1 = TVector3();
  TVector3 p2 = TVector3();
  //check charge//simple protection
  //index0 is expected to be a positive leg.
  if(v0->ChargeProng(0) > 0){//expected
    p1.SetXYZ(v0->MomPosX(),v0->MomPosY(),v0->MomPosZ());
    p2.SetXYZ(v0->MomNegX(),v0->MomNegY(),v0->MomNegZ());
  }
  else{//inverted
    AliInfo("charge is swapped.");
    p2.SetXYZ(v0->MomPosX(),v0->MomPosY(),v0->MomPosZ());
    p1.SetXYZ(v0->MomNegX(),v0->MomNegY(),v0->MomNegZ());
  }
  TVector3 z = TVector3();
  if(Bz > 0) z.SetXYZ(0,0,+1.);
  else       z.SetXYZ(0,0,-1.);

  TVector3 p12 = p1+p2;
  TVector3 wexp = p12.Cross(z);//expected
  TVector3 u_tmp  = p1.Cross(p2);
  TVector3 u = TVector3(u_tmp.X()/u_tmp.Mag(), u_tmp.Y()/u_tmp.Mag(), u_tmp.Z()/u_tmp.Mag());//unit vector
  TVector3 wmeas = p12.Cross(u);//measured

  Double_t cosPhiV = wexp.Dot(wmeas)/(wexp.Mag() * wmeas.Mag());
  Double_t phiv = TMath::ACos(cosPhiV);
  return phiv;//in radian
}
//_______________________________________________________________________________________________
void AliAnalysisTaskReducedTreeDS::ClearVectorElement()
{
  AliInfo("Number of elements of vectors is cleared.");

  //clear track variables
  fTrackMomentum.clear();
  fTrackCharge.clear();
  fTrackDCAxy.clear();
  fTrackDCAz.clear();

  fTrackPin.clear();
  fPointOnITSLayer.clear();
  fSharedPointOnITSLayer.clear();

  fChi2TPC.clear();
  fChi2ITS.clear();
  fNclusterTPC.clear();
  fNclusterITS.clear();
  fTPCNCrossedRows.clear();
  fTPCNFindableCluster.clear();
  fChi2TPCConstrainedVsGlobal.clear();

  fTPCsignal.clear();
  fTPCNsigmaEl.clear();
  fTPCNsigmaPi.clear();
  fTPCNsigmaKa.clear();
  fTPCNsigmaPr.clear();

  fITSsignal.clear();
  fITSNsigmaEl.clear();
  fITSNsigmaPi.clear();
  fITSNsigmaKa.clear();
  fITSNsigmaPr.clear();

  fTOFbeta.clear();
  fTOFNsigmaEl.clear();
  fTOFNsigmaPi.clear();
  fTOFNsigmaKa.clear();
  fTOFNsigmaPr.clear();
  fIsTOFAvailable.clear();

  fTrackMCMomentum.clear();
  fTrackMCProdVtx.clear();
  fTrackMCIndex.clear();
  fTrackMCPdgCode.clear();
  fTrackMCMotherIndex.clear();
  fTrackMCMotherPdgCode.clear();
  fTrackMCFirstMotherIndex.clear();
  fTrackMCFirstMotherPdgCode.clear();


  //clear V0 variables
  fV0Momentum.clear();
  fV0legMomentum.clear();
  fV0legPin.clear();
  fV0Lxy.clear();
  fV0alpha.clear();
  fV0qT.clear();
  fV0DCA.clear();
  fV0PsiPair.clear();
  fV0PhivPair.clear();
  fV0PointingAngle.clear();
  fV0Chi2.clear();
  fV0legChi2TPCConstrainedVsGlobal.clear();
  fV0Mass.clear();
  fV0legDCAxy.clear();
  fV0legDCAz.clear();
  fV0PointOnITSLayer.clear();
  fV0SharedPointOnITSLayer.clear();
  fV0legTPCNsigmaEl.clear();
  fV0legTPCNsigmaPi.clear();
  fV0legTPCNsigmaKa.clear();
  fV0legTPCNsigmaPr.clear();
  fV0legITSNsigmaEl.clear();
  fV0legITSNsigmaPi.clear();
  fV0legITSNsigmaKa.clear();
  fV0legITSNsigmaPr.clear();
  fV0legTOFNsigmaEl.clear();
  fV0legTOFNsigmaPi.clear();
  fV0legTOFNsigmaKa.clear();
  fV0legTOFNsigmaPr.clear();
  fV0legIsTOFAvailable.clear();

  fV0MClegMomentum.clear();
  fV0MClegProdVtx.clear();
  fV0MClegIndex.clear();
  fV0MClegPdgCode.clear();
  fV0MClegMotherIndex.clear();
  fV0MClegMotherPdgCode.clear();
  fV0MClegFirstMotherIndex.clear();
  fV0MClegFirstMotherPdgCode.clear();

  //clear MC variables
  fMCMomentum.clear();
  fMCProdVtx.clear();
  fMCIndex.clear();
  fMCPdgCode.clear();
  fMCMotherIndex.clear();
  fMCMotherPdgCode.clear();
  fMCFirstMotherIndex.clear();
  fMCFirstMotherPdgCode.clear();
  AliInfo("Number of elements of vectors is cleared. DONE!");

}
//_______________________________________________________________________________________________
void AliAnalysisTaskReducedTreeDS::ClearVectorMemory()
{
  AliInfo("Memories for vector objects are swapped with null.");

  vector<vector<Float_t>>().swap(fTrackMomentum);
  vector<Int_t>().swap(fTrackCharge);
  vector<Float_t>().swap(fTrackDCAxy);
  vector<Float_t>().swap(fTrackDCAz);

  vector<Float_t>().swap(fTrackPin);
  vector<vector<Bool_t>>().swap(fPointOnITSLayer);
  vector<vector<Bool_t>>().swap(fSharedPointOnITSLayer);

  vector<Float_t>().swap(fChi2TPC);
  vector<Float_t>().swap(fChi2ITS);
  vector<Int_t>().swap(fNclusterTPC);
  vector<Int_t>().swap(fNclusterITS);
  vector<Int_t>().swap(fTPCNCrossedRows);
  vector<Int_t>().swap(fTPCNFindableCluster);
  vector<Float_t>().swap(fChi2TPCConstrainedVsGlobal);

  vector<Float_t>().swap(fTPCsignal);
  vector<Float_t>().swap(fTPCNsigmaEl);
  vector<Float_t>().swap(fTPCNsigmaPi);
  vector<Float_t>().swap(fTPCNsigmaKa);
  vector<Float_t>().swap(fTPCNsigmaPr);

  vector<Float_t>().swap(fITSsignal);
  vector<Float_t>().swap(fITSNsigmaEl);
  vector<Float_t>().swap(fITSNsigmaPi);
  vector<Float_t>().swap(fITSNsigmaKa);
  vector<Float_t>().swap(fITSNsigmaPr);

  vector<Float_t>().swap(fTOFbeta);
  vector<Float_t>().swap(fTOFNsigmaEl);
  vector<Float_t>().swap(fTOFNsigmaPi);
  vector<Float_t>().swap(fTOFNsigmaKa);
  vector<Float_t>().swap(fTOFNsigmaPr);
  vector<Bool_t>().swap(fIsTOFAvailable);

  //MC track info
  vector<vector<Float_t>>().swap(fTrackMCMomentum);
  vector<vector<Float_t>>().swap(fTrackMCProdVtx);
  vector<Int_t>().swap(fTrackMCIndex);
  vector<Int_t>().swap(fTrackMCPdgCode);
  vector<Int_t>().swap(fTrackMCMotherIndex);
  vector<Int_t>().swap(fTrackMCMotherPdgCode);
  vector<Int_t>().swap(fTrackMCFirstMotherIndex);
  vector<Int_t>().swap(fTrackMCFirstMotherPdgCode);

  //V0 info
  vector<vector<Float_t>>().swap(fV0Momentum);
  vector<vector<vector<Float_t>>>().swap(fV0legMomentum);
  vector<vector<Float_t>>().swap(fV0legPin);
  vector<Float_t>().swap(fV0Lxy);
  vector<Float_t>().swap(fV0alpha);
  vector<Float_t>().swap(fV0qT);
  vector<Float_t>().swap(fV0DCA);
  vector<Float_t>().swap(fV0PsiPair);
  vector<Float_t>().swap(fV0PhivPair);
  vector<Float_t>().swap(fV0PointingAngle);
  vector<Float_t>().swap(fV0Chi2);
  vector<vector<Float_t>>().swap(fV0legChi2TPCConstrainedVsGlobal);
  vector<vector<Float_t>>().swap(fV0Mass);
  vector<vector<Float_t>>().swap(fV0legDCAxy);
  vector<vector<Float_t>>().swap(fV0legDCAz);
  vector<vector<vector<Bool_t>>>().swap(fV0PointOnITSLayer);
  vector<vector<vector<Bool_t>>>().swap(fV0SharedPointOnITSLayer);
  vector<vector<Float_t>>().swap(fV0legTPCNsigmaEl);
  vector<vector<Float_t>>().swap(fV0legTPCNsigmaPi);
  vector<vector<Float_t>>().swap(fV0legTPCNsigmaKa);
  vector<vector<Float_t>>().swap(fV0legTPCNsigmaPr);
  vector<vector<Float_t>>().swap(fV0legITSNsigmaEl);
  vector<vector<Float_t>>().swap(fV0legITSNsigmaPi);
  vector<vector<Float_t>>().swap(fV0legITSNsigmaKa);
  vector<vector<Float_t>>().swap(fV0legITSNsigmaPr);
  vector<vector<Float_t>>().swap(fV0legTOFNsigmaEl);
  vector<vector<Float_t>>().swap(fV0legTOFNsigmaPi);
  vector<vector<Float_t>>().swap(fV0legTOFNsigmaKa);
  vector<vector<Float_t>>().swap(fV0legTOFNsigmaPr);
  vector<vector<Bool_t>>().swap(fV0legIsTOFAvailable);

  //MC V0 info //be carefull, there is no TRUE V0 object!
  vector<vector<vector<Float_t>>>().swap(fV0MClegMomentum);
  vector<vector<vector<Float_t>>>().swap(fV0MClegProdVtx);
  vector<vector<Int_t>>().swap(fV0MClegIndex);
  vector<vector<Int_t>>().swap(fV0MClegPdgCode);
  vector<vector<Int_t>>().swap(fV0MClegMotherIndex);
  vector<vector<Int_t>>().swap(fV0MClegMotherPdgCode);
  vector<vector<Int_t>>().swap(fV0MClegFirstMotherIndex);
  vector<vector<Int_t>>().swap(fV0MClegFirstMotherPdgCode);

  //MC variables for true electrons
  vector<vector<Float_t>>().swap(fMCMomentum);
  vector<vector<Float_t>>().swap(fMCProdVtx);
  vector<Int_t>().swap(fMCIndex);
  vector<Int_t>().swap(fMCPdgCode);
  vector<Int_t>().swap(fMCMotherIndex);
  vector<Int_t>().swap(fMCMotherPdgCode);
  vector<Int_t>().swap(fMCFirstMotherIndex);
  vector<Int_t>().swap(fMCFirstMotherPdgCode);
  AliInfo("Memories for vector objects are swapped with null. DONE!");
}
//_______________________________________________________________________________________________
void AliAnalysisTaskReducedTreeDS::ExtractQnVectors()
{
  AliInfo("extract Qn vectors from Qn correction framework.");
  const TString TPCEPname[3] = {"TPC","TPCNegEta","TPCPosEta"};
  const TString V0EPname[3]  = {"VZERO","VZEROA","VZEROC"};
  const TString ZDCEPname[2]  = {"ZDCA","ZDCC"};
  const TString Qnorm = "QoverM";
  const Int_t harmonics = 2;

  fIsQnTPCAvailable = kFALSE;
  fIsQnV0Available = kFALSE;
  fIsQnZDCAvailable = kFALSE;

  for(Int_t i=0;i<2;i++){//Qx, Qy
    fQ2vectorTPC[i] = -999;
    fQ2vectorTPCNegEta[i] = -999;
    fQ2vectorTPCPosEta[i] = -999;
    fQ2vectorV0[i] = -999;
    fQ2vectorV0A[i] = -999;
    fQ2vectorV0C[i] = -999;
    fQ2vectorZDCA[i] = -999;
    fQ2vectorZDCC[i] = -999;
  }

  if(fFlowQnVectorMgr == NULL){
    AliInfo("fFlowQnVectorMgr does not exist. return.");
    return;
  }

  //TPC
  TList* qnlist = fFlowQnVectorMgr->GetQnVectorList();
  //qnlist->Print();//only for debeg
  const AliQnCorrectionsQnVector *QnVectorTPCDet[3];
  for(Int_t i=0;i<3;i++){
    QnVectorTPCDet[i] = GetQnVectorFromList(qnlist,Form("%s%s",TPCEPname[i].Data(),Qnorm.Data()),"latest","plain");
  }//end of det loop

  if(!QnVectorTPCDet[0] ||!QnVectorTPCDet[1] || !QnVectorTPCDet[2]){
    AliInfo("Event is rejected because event plane is not found or bad event plane quality in TPC.");
    fIsQnTPCAvailable = kFALSE;
    fQ2vectorTPC[0]       = -999;//FullTPC
    fQ2vectorTPC[1]       = -999;//FullTPC
    fQ2vectorTPCNegEta[0] = -999;//TPCNegEta
    fQ2vectorTPCNegEta[1] = -999;//TPCNegEta
    fQ2vectorTPCPosEta[0] = -999;//TPCPosEta
    fQ2vectorTPCPosEta[1] = -999;//TPCPosEta
  }
  else{
    fIsQnTPCAvailable = kTRUE;
    fQ2vectorTPC[0]       = QnVectorTPCDet[0]->Qx(harmonics);//FullTPC
    fQ2vectorTPC[1]       = QnVectorTPCDet[0]->Qy(harmonics);//FullTPC
    fQ2vectorTPCNegEta[0] = QnVectorTPCDet[1]->Qx(harmonics);//TPCNegEta
    fQ2vectorTPCNegEta[1] = QnVectorTPCDet[1]->Qy(harmonics);//TPCNegEta
    fQ2vectorTPCPosEta[0] = QnVectorTPCDet[2]->Qx(harmonics);//TPCPosEta
    fQ2vectorTPCPosEta[1] = QnVectorTPCDet[2]->Qy(harmonics);//TPCPosEta
    for(Int_t i=0;i<3;i++) AliInfo(Form("harmonics %d | TPC sub detector name %s%s : Qx = %e, Qy = %e",harmonics,TPCEPname[i].Data(),Qnorm.Data(),QnVectorTPCDet[i]->Qx(harmonics),QnVectorTPCDet[i]->Qy(harmonics)));
  }

  //V0
  const AliQnCorrectionsQnVector *QnVectorV0Det[3];
  for(Int_t i=0;i<3;i++){
    QnVectorV0Det[i]  = GetQnVectorFromList(qnlist,Form("%s%s",V0EPname[i].Data(),Qnorm.Data()),"latest","raw");
  }//end of det loop

  if(!QnVectorV0Det[0] || !QnVectorV0Det[1] || !QnVectorV0Det[2]){
    AliInfo("Event is rejected because event plane is not found or bad event plane quality in VZERO.");
    fIsQnV0Available = kFALSE;
    fQ2vectorV0[0]  = -999;//FullV0
    fQ2vectorV0[1]  = -999;//FullV0
    fQ2vectorV0A[0] = -999;//V0A
    fQ2vectorV0A[1] = -999;//V0A
    fQ2vectorV0C[0] = -999;//V0C
    fQ2vectorV0C[1] = -999;//V0C
  }
  else{
    fIsQnV0Available = kTRUE;
    fQ2vectorV0[0]  = QnVectorV0Det[0]->Qx(harmonics);//FullV0
    fQ2vectorV0[1]  = QnVectorV0Det[0]->Qy(harmonics);//FullV0
    fQ2vectorV0A[0] = QnVectorV0Det[1]->Qx(harmonics);//V0A
    fQ2vectorV0A[1] = QnVectorV0Det[1]->Qy(harmonics);//V0A
    fQ2vectorV0C[0] = QnVectorV0Det[2]->Qx(harmonics);//V0C
    fQ2vectorV0C[1] = QnVectorV0Det[2]->Qy(harmonics);//V0C
    for(Int_t i=0;i<3;i++) AliInfo(Form("harmonics %d | V0  sub detector name %s%s : Qx = %e, Qy = %e",harmonics,V0EPname[i].Data(),Qnorm.Data(),QnVectorV0Det[i]->Qx(harmonics),QnVectorV0Det[i]->Qy(harmonics)));
  }

  //ZDC
  const AliQnCorrectionsQnVector *QnVectorZDCDet[2];
  for(Int_t i=0;i<2;i++){
    QnVectorZDCDet[i]  = GetQnVectorFromList(qnlist,Form("%s%s",ZDCEPname[i].Data(),""),"latest","raw");
  }//end of det loop

  if(!QnVectorZDCDet[0] || !QnVectorZDCDet[1]){
    AliInfo("Event is rejected because event plane is not found or bad event plane quality in VZERO.");
    fIsQnZDCAvailable = kFALSE;
    fQ2vectorZDCA[0] = -999;//ZDCA
    fQ2vectorZDCA[1] = -999;//ZDCA
    fQ2vectorZDCC[0] = -999;//ZDCC
    fQ2vectorZDCC[1] = -999;//ZDCC
  }
  else{
    fIsQnZDCAvailable = kTRUE;
    fQ2vectorZDCA[0] = QnVectorZDCDet[0]->Qx(harmonics);//ZDCA
    fQ2vectorZDCA[1] = QnVectorZDCDet[0]->Qy(harmonics);//ZDCA
    fQ2vectorZDCC[0] = QnVectorZDCDet[1]->Qx(harmonics);//ZDCC
    fQ2vectorZDCC[1] = QnVectorZDCDet[1]->Qy(harmonics);//ZDCC
    for(Int_t i=0;i<2;i++) AliInfo(Form("harmonics %d | ZDC  sub detector name %s%s : Qx = %e, Qy = %e",harmonics,ZDCEPname[i].Data(),Qnorm.Data(),QnVectorZDCDet[i]->Qx(harmonics),QnVectorZDCDet[i]->Qy(harmonics)));
  }

}
//_______________________________________________________________________________________________
const AliQnCorrectionsQnVector *AliAnalysisTaskReducedTreeDS::GetQnVectorFromList(const TList *list, const char* subdetector, const char *expcorr, const char *altcorr)
{
  AliQnCorrectionsQnVector *theQnVector = NULL;
  TList *pQvecList = dynamic_cast<TList*> (list->FindObject(subdetector));
  if(pQvecList != NULL){//sub detector is found
    if(TString(expcorr).EqualTo("latest")) theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
    else                                   theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(expcorr);

    if(theQnVector == NULL || !(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0)){ //the Qn vector for the expected correction is not found
      AliInfo(Form("expected correction (%s) is not found. use %s as an alternative step in %s.",expcorr,altcorr,subdetector));
      if(TString(altcorr).EqualTo("latest")) theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
      else                                   theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(altcorr);
    }
    //check the Qn vector quality
    if(!(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0)) theQnVector = NULL; //bad quality, discarded
  }
  return theQnVector;
}
//_______________________________________________________________________________________________

