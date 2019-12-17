#include <TObject.h>
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
#include <TMath.h>
#include <TDatabasePDG.h>

#include "AliESDtrackCuts.h"
#include "AliESDHeader.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliESDv0KineCuts.h"

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
#include "AliTimeRangeCut.h"

#include "AliAODMCHeader.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODv0KineCuts.h"

#include "AliAODMCParticle.h"
#include "AliAODInputHandler.h"
#include "AliAODTracklets.h"
#include "AliAnalysisUtils.h"

#include "AliKFVertex.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliVParticle.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskReducedTreeDS.h"

using namespace std;

//event properties and electrons tree.
//Daiki Sekihata
//daiki.sekihata@cern.ch

ClassImp(AliAnalysisTaskReducedTreeDS)
//_______________________________________________________________________________________________
AliAnalysisTaskReducedTreeDS::AliAnalysisTaskReducedTreeDS():
  AliAnalysisTaskSE(),
  fMinPtCut(0.2),
  fMaxEtaCut(0.8),
  fMinTPCNsigmaEleCut(-5.),
  fMaxTPCNsigmaEleCut(5.),
  fESDtrackCutsGlobalNoDCA(0x0),
  fESDv0KineCuts(0x0),
  fAODv0KineCuts(0x0),
  fTree(0x0),
  fPIDResponse(0x0),
  fFlowQnVectorMgr(0x0),
  fEvent(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fMCEvent(0x0),
  fHasMC(kFALSE),
  fTimeRangeCut(),
  fMCArray(0x0),
  fRunNumber(-1),
  fMagneticField(0),
  fBCNumber(-1),
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
  fNTrackTPC(-1),
  fNITSCluster(),
  fNSPDTracklet05(-1),
  fNSPDTracklet10(-1),
  fV0AMultiplicity(-1),
  fV0CMultiplicity(-1),
  fIsPileupFromSPD(kFALSE),
  fIsPileupFromSPDInMultBins(kFALSE),
  fIsPileupMV(kFALSE),
  fTPCpileupMultiplicity(),
  fTPCpileupZ(),
  fIskINT7(kFALSE),
  fIskCentral(kFALSE),
  fIskSemiCentral(kFALSE),
  fIskHighMult(kFALSE),
  fIskHighMultV0(kFALSE),
  fIskHighMultSPD(kFALSE),
  fIsBadTimeRangeTPC(kFALSE),
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
  fTPCsignalN(0),
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
  fTrackMCGeneratorIndex(0),
  fTrackMCIsPhysicalPrimary(0),
  fTrackMCIsSecondaryFromMaterial(0),
  fTrackMCIsSecondaryFromWeakDecay(0),
  fTrackMCIndex(0),
  fTrackMCPdgCode(0),
  fTrackMCMotherIndex(0),
  fTrackMCMotherPdgCode(0),
  fTrackMCFirstMotherIndex(0),
  fTrackMCFirstMotherPdgCode(0),
  fTrackMCFirstMotherMomentum(0),
  fV0OnFly(0),
  fV0legMomentum(0),
  fV0legPin(0),
  fV0Lxy(0),
  fV0alpha(0),
  fV0qT(0),
  fV0Candidate(0),
  fV0Mass(0),
  fV0legPointOnITSLayer(0),
  fV0legSharedPointOnITSLayer(0),
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
  fV0MClegGeneratorIndex(0),
  fV0MClegIndex(0),
  fV0MClegPdgCode(0),
  fV0MClegMotherIndex(0),
  fV0MClegMotherPdgCode(0),
  fV0MClegFirstMotherIndex(0),
  fV0MClegFirstMotherPdgCode(0),
  fV0MClegFirstMotherMomentum(0),
  fMCVertex(),
  fMCMomentum(0),
  fMCProdVtx(0),
  fMCGeneratorIndex(0),
  fMCGeneratorName(0),
  fMCIsPhysicalPrimary(0),
  fMCIndex(0),
  fMCPdgCode(0),
  fMCMotherIndex(0),
  fMCMotherPdgCode(0),
  fMCFirstMotherIndex(0),
  fMCFirstMotherPdgCode(0),
  fMCFirstMotherMomentum(0)
{
  for(Int_t i=0;i<3;i++) fVertex[i] = 0;
  for(Int_t i=0;i<3;i++) fMCVertex[i] = 0;
  for(Int_t i=0;i<2;i++) fNITSCluster[i] = 0;

  for(Int_t i=0;i<2;i++) fTPCpileupMultiplicity[i] = 0;
  for(Int_t i=0;i<2;i++) fTPCpileupZ[i] = 0;

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

  fESDtrackCutsGlobalNoDCA = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
  fESDtrackCutsGlobalNoDCA->SetMaxDCAToVertexXY(2.4);
  fESDtrackCutsGlobalNoDCA->SetMaxDCAToVertexZ(3.2);
  fESDtrackCutsGlobalNoDCA->SetDCAToVertex2D(kTRUE);

  fESDv0KineCuts = new AliESDv0KineCuts();
  fAODv0KineCuts = new AliAODv0KineCuts();
  fESDv0KineCuts->SetMode(AliESDv0KineCuts::kPurity,AliESDv0KineCuts::kPbPb);
  fAODv0KineCuts->SetMode(AliAODv0KineCuts::kPurity,AliAODv0KineCuts::kPbPb);

}
//_______________________________________________________________________________________________
AliAnalysisTaskReducedTreeDS::AliAnalysisTaskReducedTreeDS(const char *name):
  AliAnalysisTaskSE(name),
  fMinPtCut(0.2),
  fMaxEtaCut(0.8),
  fMinTPCNsigmaEleCut(-5.),
  fMaxTPCNsigmaEleCut(5.),
  fESDtrackCutsGlobalNoDCA(0x0),
  fESDv0KineCuts(0x0),
  fAODv0KineCuts(0x0),
  fTree(0x0),
  fPIDResponse(0x0),
  fFlowQnVectorMgr(0x0),
  fEvent(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fMCEvent(0x0),
  fHasMC(kFALSE),
  fTimeRangeCut(),
  fMCArray(0x0),
  fRunNumber(-1),
  fMagneticField(0),
  fBCNumber(-1),
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
  fNTrackTPC(-1),
  fNITSCluster(),
  fNSPDTracklet05(-1),
  fNSPDTracklet10(-1),
  fV0AMultiplicity(-1),
  fV0CMultiplicity(-1),
  fIsPileupFromSPD(kFALSE),
  fIsPileupFromSPDInMultBins(kFALSE),
  fIsPileupMV(kFALSE),
  fTPCpileupMultiplicity(),
  fTPCpileupZ(),
  fIskINT7(kFALSE),
  fIskCentral(kFALSE),
  fIskSemiCentral(kFALSE),
  fIskHighMult(kFALSE),
  fIskHighMultV0(kFALSE),
  fIskHighMultSPD(kFALSE),
  fIsBadTimeRangeTPC(kFALSE),
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
  fTPCsignalN(0),
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
  fTrackMCGeneratorIndex(0),
  fTrackMCIsPhysicalPrimary(0),
  fTrackMCIsSecondaryFromMaterial(0),
  fTrackMCIsSecondaryFromWeakDecay(0),
  fTrackMCIndex(0),
  fTrackMCPdgCode(0),
  fTrackMCMotherIndex(0),
  fTrackMCMotherPdgCode(0),
  fTrackMCFirstMotherIndex(0),
  fTrackMCFirstMotherPdgCode(0),
  fTrackMCFirstMotherMomentum(0),
  fV0OnFly(0),
  fV0legMomentum(0),
  fV0legPin(0),
  fV0Lxy(0),
  fV0alpha(0),
  fV0qT(0),
  fV0Candidate(0),
  fV0Mass(0),
  fV0legPointOnITSLayer(0),
  fV0legSharedPointOnITSLayer(0),
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
  fV0MClegGeneratorIndex(0),
  fV0MClegIndex(0),
  fV0MClegPdgCode(0),
  fV0MClegMotherIndex(0),
  fV0MClegMotherPdgCode(0),
  fV0MClegFirstMotherIndex(0),
  fV0MClegFirstMotherPdgCode(0),
  fV0MClegFirstMotherMomentum(0),
  fMCVertex(),
  fMCMomentum(0),
  fMCProdVtx(0),
  fMCGeneratorIndex(0),
  fMCGeneratorName(0),
  fMCIsPhysicalPrimary(0),
  fMCIndex(0),
  fMCPdgCode(0),
  fMCMotherIndex(0),
  fMCMotherPdgCode(0),
  fMCFirstMotherIndex(0),
  fMCFirstMotherPdgCode(0),
  fMCFirstMotherMomentum(0)
{
  for(Int_t i=0;i<3;i++) fVertex[i] = 0;
  for(Int_t i=0;i<3;i++) fMCVertex[i] = 0;
  for(Int_t i=0;i<2;i++) fNITSCluster[i] = 0;

  for(Int_t i=0;i<2;i++) fTPCpileupMultiplicity[i] = 0;
  for(Int_t i=0;i<2;i++) fTPCpileupZ[i] = 0;

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

  fESDtrackCutsGlobalNoDCA = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
  fESDtrackCutsGlobalNoDCA->SetMaxDCAToVertexXY(2.4);
  fESDtrackCutsGlobalNoDCA->SetMaxDCAToVertexZ(3.2);
  fESDtrackCutsGlobalNoDCA->SetDCAToVertex2D(kTRUE);

  fESDv0KineCuts = new AliESDv0KineCuts();
  fAODv0KineCuts = new AliAODv0KineCuts();
  fESDv0KineCuts->SetMode(AliESDv0KineCuts::kPurity,AliESDv0KineCuts::kPbPb);
  fAODv0KineCuts->SetMode(AliAODv0KineCuts::kPurity,AliAODv0KineCuts::kPbPb);

  DefineInput(0,TChain::Class());
  DefineOutput(1,TTree::Class());  // reduced information tree

}
//_______________________________________________________________________________________________
AliAnalysisTaskReducedTreeDS::~AliAnalysisTaskReducedTreeDS()
{
  delete fESDtrackCutsGlobalNoDCA;
  delete fESDv0KineCuts;
  delete fAODv0KineCuts;
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
  fTree->Branch("fBCNumber",&fBCNumber,"fBCNumber/s");//UShort_t

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
  fTree->Branch("fNTrackTPC",&fNTrackTPC,"fNTrackTPC/I");
  fTree->Branch("fNITSCluster",fNITSCluster,"fNITSCluster[2]/I");
  fTree->Branch("fNSPDTracklet05",&fNSPDTracklet05,"fNSPDTracklet05/I");
  fTree->Branch("fNSPDTracklet10",&fNSPDTracklet10,"fNSPDTracklet10/I");
  fTree->Branch("fV0AMultiplicity",&fV0AMultiplicity,"fV0AMultiplicity/F");
  fTree->Branch("fV0CMultiplicity",&fV0CMultiplicity,"fV0CMultiplicity/F");

  fTree->Branch("fIsPileupFromSPD",&fIsPileupFromSPD,"fIsPileupFromSPD/O");
  fTree->Branch("fIsPileupFromSPDInMultBins",&fIsPileupFromSPDInMultBins,"fIsPileupFromSPDInMultBins/O");
  fTree->Branch("fIsPileupMV",&fIsPileupMV,"fIsPileupMV/O");

  fTree->Branch("fTPCpileupMultiplicity",fTPCpileupMultiplicity,"fTPCpileupMultiplicity[2]/I");
  fTree->Branch("fTPCpileupZ",fTPCpileupZ,"fTPCpileupZ[2]/F");

  fTree->Branch("fIskINT7",&fIskINT7,"fIskINT7/O");
  fTree->Branch("fIskCentral",&fIskCentral,"fIskCentral/O");
  fTree->Branch("fIskSemiCentral",&fIskSemiCentral,"fIskSemiCentral/O");
  fTree->Branch("fIskHighMult",&fIskHighMult,"fIskHighMult/O");
  fTree->Branch("fIskHighMultV0",&fIskHighMultV0,"fIskHighMultV0/O");
  fTree->Branch("fIskHighMultSPD",&fIskHighMultSPD,"fIskHighMultSPD/O");
  fTree->Branch("fIsBadTimeRangeTPC",&fIsBadTimeRangeTPC,"fIsBadTimeRangeTPC/O");

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

  fTree->Branch("fTPCsignalN",&fTPCsignalN);
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
  fTree->Branch("fTrackMCGeneratorIndex",&fTrackMCGeneratorIndex);
  fTree->Branch("fTrackMCIsPhysicalPrimary",&fTrackMCIsPhysicalPrimary);
  fTree->Branch("fTrackMCIsSecondaryFromMaterial",&fTrackMCIsSecondaryFromMaterial);
  fTree->Branch("fTrackMCIsSecondaryFromWeakDecay",&fTrackMCIsSecondaryFromWeakDecay);
  fTree->Branch("fTrackMCIndex",&fTrackMCIndex);
  fTree->Branch("fTrackMCPdgCode",&fTrackMCPdgCode);
  fTree->Branch("fTrackMCMotherIndex",&fTrackMCMotherIndex);
  fTree->Branch("fTrackMCMotherPdgCode",&fTrackMCMotherPdgCode);
  fTree->Branch("fTrackMCFirstMotherIndex",&fTrackMCFirstMotherIndex);
  fTree->Branch("fTrackMCFirstMotherPdgCode",&fTrackMCFirstMotherPdgCode);
  fTree->Branch("fTrackMCFirstMotherMomentum",&fTrackMCFirstMotherMomentum);

  fTree->Branch("fV0OnFly",&fV0OnFly);
  fTree->Branch("fV0legMomentum",&fV0legMomentum);
  fTree->Branch("fV0legPin",&fV0legPin);

  fTree->Branch("fV0Lxy",&fV0Lxy);
  fTree->Branch("fV0alpha",&fV0alpha);
  fTree->Branch("fV0qT",&fV0qT);

  fTree->Branch("fV0Candidate",&fV0Candidate);
  fTree->Branch("fV0Mass",&fV0Mass);
  fTree->Branch("fV0legPointOnITSLayer",&fV0legPointOnITSLayer);
  fTree->Branch("fV0legSharedPointOnITSLayer",&fV0legSharedPointOnITSLayer);

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
  fTree->Branch("fV0MClegGeneratorIndex",&fV0MClegGeneratorIndex);
  fTree->Branch("fV0MClegIndex",&fV0MClegIndex);
  fTree->Branch("fV0MClegPdgCode",&fV0MClegPdgCode);
  fTree->Branch("fV0MClegMotherIndex",&fV0MClegMotherIndex);
  fTree->Branch("fV0MClegMotherPdgCode",&fV0MClegMotherPdgCode);
  fTree->Branch("fV0MClegFirstMotherIndex",&fV0MClegFirstMotherIndex);
  fTree->Branch("fV0MClegFirstMotherPdgCode",&fV0MClegFirstMotherPdgCode);
  fTree->Branch("fV0MClegFirstMotherMomentum",&fV0MClegFirstMotherMomentum);

  //MC true info
  fTree->Branch("fMCVertex",fMCVertex,"fMCVertex[3]/F");
  fTree->Branch("fMCMomentum",&fMCMomentum);
  fTree->Branch("fMCProdVtx",&fMCProdVtx);
  fTree->Branch("fMCGeneratorIndex",&fMCGeneratorIndex);
  fTree->Branch("fMCGeneratorName",&fMCGeneratorName);
  fTree->Branch("fMCIsPhysicalPrimary",&fMCIsPhysicalPrimary);
  fTree->Branch("fMCIndex",&fMCIndex);
  fTree->Branch("fMCPdgCode",&fMCPdgCode);
  fTree->Branch("fMCMotherIndex",&fMCMotherIndex);
  fTree->Branch("fMCMotherPdgCode",&fMCMotherPdgCode);
  fTree->Branch("fMCFirstMotherIndex",&fMCFirstMotherIndex);
  fTree->Branch("fMCFirstMotherPdgCode",&fMCFirstMotherPdgCode);
  fTree->Branch("fMCFirstMotherMomentum",&fMCFirstMotherMomentum);

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

  fIsBadTimeRangeTPC = kFALSE;
  fTimeRangeCut.InitFromEvent(InputEvent());
  fIsBadTimeRangeTPC = fTimeRangeCut.CutEvent(InputEvent());

  if(fESDEvent)      fHasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);
  else if(fAODEvent) fHasMC = (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->MCEvent() != 0x0);

  AliInfo(Form("fHasMC = %d",fHasMC));

  if(fHasMC){
    if(fAODEvent) GetMCInfoAOD();

    if(fESDEvent)      fMCEvent = dynamic_cast<AliMCEvent*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()->MCEvent());
    else if(fAODEvent) fMCEvent = dynamic_cast<AliMCEvent*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->MCEvent());
  }
  fRunNumber = fEvent->GetRunNumber();
  fMagneticField = fEvent->GetMagneticField();
  fBCNumber = fEvent->GetBunchCrossNumber();

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

  if(fESDEvent)      for(Int_t i=0;i<2;i++) fNITSCluster[i] = fESDEvent->GetMultiplicity()->GetNumberOfITSClusters(i);
  else if(fAODEvent) for(Int_t i=0;i<2;i++) fNITSCluster[i] = fAODEvent->GetMultiplicity()->GetNumberOfITSClusters(i);

  //Check SPD tracklet multilicity in |eta| < 1
  fNSPDTracklet05 = 0;
  fNSPDTracklet10 = 0;

  if(fESDEvent){
    //AliESDTracklets *tracklets = (AliESDTracklets*)fESDEvent->GetTracklets();
    const Int_t Ntl = fESDEvent->GetMultiplicity()->GetNumberOfTracklets();
    for(Int_t itl=0;itl<Ntl;itl++){
      Double_t eta   = fESDEvent->GetMultiplicity()->GetEta(itl);
      if(TMath::Abs(eta) < 0.5) fNSPDTracklet05++;
      if(TMath::Abs(eta) < 1.0) fNSPDTracklet10++;
    }
  }
  else if(fAODEvent){
    AliAODTracklets *tracklets = (AliAODTracklets*)fAODEvent->GetTracklets();
    const Int_t Ntl = tracklets->GetNumberOfTracklets();
    for(Int_t itl=0;itl<Ntl;itl++){
      Double_t theta = tracklets->GetTheta(itl);
      Double_t eta   = -TMath::Log(TMath::Tan(theta/2.0));
      if(TMath::Abs(eta) < 0.5) fNSPDTracklet05++;
      if(TMath::Abs(eta) < 1.0) fNSPDTracklet10++;
    }
  }

  const AliVVertex *vVertex = fEvent->GetPrimaryVertex();
  fVertex[0] = vVertex->GetX();
  fVertex[1] = vVertex->GetY();
  fVertex[2] = vVertex->GetZ();
  fNContributor = vVertex->GetNContributors();
  if(fESDEvent){
    fNTPCCluster = fESDEvent->GetNumberOfTPCClusters();
    fNTrackTPC   = fESDEvent->GetNumberOfTPCTracks();
  }
  else if(fAODEvent){
    fNTPCCluster = fAODEvent->GetNumberOfTPCClusters();
    fNTrackTPC   = fAODEvent->GetNumberOfTPCTracks();
  }

  Float_t DCAxy_PU = -999, DCAz_PU = -999, tgl = 0;
  fTPCpileupMultiplicity[0] = 0;
  fTPCpileupMultiplicity[1] = 0;
  fTPCpileupZ[0] = 0;
  fTPCpileupZ[1] = 0;
  vector<Float_t> vec_puZ_pos;
  vector<Float_t> vec_puZ_neg;
  vec_puZ_pos.clear();
  vec_puZ_neg.clear();

  const Int_t Ntrack = fEvent->GetNumberOfTracks();
  if(fESDEvent){
    for(Int_t itrack=0;itrack<Ntrack;itrack++){
      AliESDtrack *esdtrack = (AliESDtrack*)fEvent->GetTrack(itrack);
      if(esdtrack->IsOn(AliVTrack::kITSin)) continue;

      esdtrack->GetImpactParameters(DCAxy_PU,DCAz_PU);
      if(TMath::Abs(DCAxy_PU) < 3. && TMath::Abs(DCAz_PU) > 4.){
        tgl = esdtrack->Pz() / esdtrack->Pt();
        if(tgl > +0.1) vec_puZ_pos.push_back(esdtrack->GetZ());
        if(tgl < -0.1) vec_puZ_neg.push_back(esdtrack->GetZ());
      }
    }//end of track loop
  }//end of ESD
  else if(fAODEvent){
    for(Int_t itrack=0;itrack<Ntrack;itrack++){
      AliAODTrack *aodtrack = (AliAODTrack*)fEvent->GetTrack(itrack);
      if(aodtrack->IsOn(AliVTrack::kITSin)) continue;

      AliAODVertex *av = (AliAODVertex*)aodtrack->GetProdVertex();
      aodtrack->GetImpactParameters(DCAxy_PU,DCAz_PU);
      if(TMath::Abs(DCAxy_PU) < 3. && TMath::Abs(DCAz_PU) > 4.){
        tgl = aodtrack->Pz() / aodtrack->Pt();
        if(tgl > +0.1) vec_puZ_pos.push_back(av->GetZ());
        if(tgl < -0.1) vec_puZ_neg.push_back(av->GetZ());
      }
    }//end of track loop
  }//end of AOD
  fTPCpileupMultiplicity[0] = Int_t(vec_puZ_pos.size());
  fTPCpileupMultiplicity[1] = Int_t(vec_puZ_neg.size());
  fTPCpileupZ[0] = Median(vec_puZ_pos);
  fTPCpileupZ[1] = Median(vec_puZ_neg);

  vec_puZ_pos.clear();
  vec_puZ_neg.clear();
  vector<Float_t>().swap(vec_puZ_pos);
  vector<Float_t>().swap(vec_puZ_neg);

  fNTrackTPCout = 0;
  if(fESDEvent) fNTrackTPCout = fESDEvent->GetNTPCTrackBeforeClean();
  else if(fAODEvent) fNTrackTPCout = 0;//to be inclemented in track loop

  FillTrackInfo();
  if(fESDEvent)      FillV0InfoESD();
  else if(fAODEvent) FillV0InfoAOD();

  AliInfo(Form("fCentralityV0M = %3.2f %% , fNSPDTracklet05 = %d , fNSPDTracklet10 = %d , fNTPCCluster = %d , fNTrackTPC = %d",fCentralityV0M,fNSPDTracklet05,fNSPDTracklet10,fNTPCCluster,fNTrackTPC));

  if(fHasMC) ProcessMC(option);

  fTree->Fill();
  ClearVectorMemory();

}
//_______________________________________________________________________________________________
void AliAnalysisTaskReducedTreeDS::FillTrackInfo()
{
  Float_t DCAxy = -999, DCAz = -999;
  Float_t Chi2Global = -1;
  Float_t TOFbeta = -999;
  Double32_t expt[5] = {0,0,0,0,0};
  const Int_t Ntrack = fEvent->GetNumberOfTracks();

  AliESDVertex* vVertex = 0x0;
  if(fESDEvent) vVertex = (AliESDVertex*)fESDEvent->GetPrimaryVertexTracks();

  for(Int_t itrack=0;itrack<Ntrack;itrack++){
    AliVTrack *track = (AliVTrack*)fEvent->GetTrack(itrack);

    if(track->Pt() < fMinPtCut) continue;
    if(TMath::Abs(track->Eta()) > fMaxEtaCut) continue;
    ULong64_t status = track->GetStatus();

    if(fESDEvent){
      AliESDtrack *esdtrack = dynamic_cast<AliESDtrack*>(track);
      if(!fESDtrackCutsGlobalNoDCA->AcceptTrack(esdtrack)) continue;//standard cuts with very loose DCA cut //bit4
      Chi2Global = esdtrack->GetChi2TPCConstrainedVsGlobal(vVertex);
      esdtrack->GetImpactParameters(DCAxy,DCAz);
    }
    else if(fAODEvent){
      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(track);
      if(!aodtrack->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA)) continue;//standard cuts with very loose DCA cut //bit4
      Chi2Global = aodtrack->GetChi2TPCConstrainedVsGlobal();
      aodtrack->GetImpactParameters(DCAxy,DCAz);
      if(status & AliVTrack::kTPCout) fNTrackTPCout++;
    }

    if(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron) < fMinTPCNsigmaEleCut) continue;//pre-select electrons to reduce data size.
    if(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron) > fMaxTPCNsigmaEleCut) continue;//pre-select electrons to reduce data size.

//    if(TMath::Abs(DCAxy) > 1.) continue;
//    if(TMath::Abs(DCAz) > 3.) continue;
//    if(track->GetNcls(0) < 3)  continue;//minimum number of ITS cluster 3
//    if(track->GetNcls(1) < 70) continue;//minimum number of TPC cluster 70
//    if((Double_t)(track->GetTPCchi2()) / (Double_t)(track->GetNcls(1)) > 4.) continue;//maximum chi2 per cluster TPC
//    if((Double_t)(track->GetITSchi2()) / (Double_t)(track->GetNcls(0)) > 5.) continue;//maximum chi2 per cluster ITS

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
    fTPCsignalN.push_back(track->GetTPCsignalN());
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

    TOFbeta = -999;
    if(!isTOFOK){
      TOFbeta = -999;
    }
    else{
      expt[0] = 0; expt[1] = 0; expt[2] = 0; expt[3] = 0; expt[4] = 0;
      track->GetIntegratedTimes(expt,5); 
      Double32_t length = TMath::C() * expt[0] * 1e-12;// m
      Double32_t time   = track->GetTOFsignal();       //ps
      time -= fPIDResponse->GetTOFResponse().GetStartTime(track->P()); // ps
      time *= 1e-12;//ps -> s
      Double32_t velocity = length / time;
      TOFbeta = velocity / TMath::C();
      //AliInfo(Form("isTOFOK = %d , velocity = %e , length = %e , time = %e , TOFbeta = %e",isTOFOK, velocity, length, time, TOFbeta));
    }

    fTOFbeta.push_back(TOFbeta);
    fTOFNsigmaEl.push_back(fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
    fTOFNsigmaPi.push_back(fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion));
    fTOFNsigmaKa.push_back(fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon));
    fTOFNsigmaPr.push_back(fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton));
    fIsTOFAvailable.push_back(isTOFOK);

    //MC info for reconstructed tracks
    if(fHasMC){
      Int_t label = TMath::Abs(track->GetLabel());
      AliMCParticle *p = (AliMCParticle*)fMCEvent->GetTrack(label);
      Int_t pdg = p->PdgCode();
      Int_t genID = p->GetGeneratorIndex();
      fTrackMCGeneratorIndex.push_back(genID);
      fTrackMCIsPhysicalPrimary.push_back(p->IsPhysicalPrimary());
      fTrackMCIsSecondaryFromMaterial.push_back(p->IsSecondaryFromMaterial());
      fTrackMCIsSecondaryFromWeakDecay.push_back(p->IsSecondaryFromWeakDecay());
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
        AliMCParticle *mp = (AliMCParticle*)fMCEvent->GetTrack(mother_index);
        mother_pdg = mp->PdgCode();
        fTrackMCMotherIndex.push_back(mother_index);
        fTrackMCMotherPdgCode.push_back(mother_pdg);
        //printf("mother_index = %d , mother_pdg = %d\n",mother_index,mother_pdg);

        //check first mother
        Int_t first_mother_index = GetFirstMother(p);
        Int_t first_mother_pdg   = mp->PdgCode();
        AliMCParticle *fmp = mp;
        if(first_mother_index > -1){
          fmp = (AliMCParticle*)fMCEvent->GetTrack(first_mother_index);
          first_mother_pdg = fmp->PdgCode();
        }
        fTrackMCFirstMotherIndex.push_back(first_mother_index);
        fTrackMCFirstMotherPdgCode.push_back(first_mother_pdg);

        vector<Float_t> fmvec(3,0);
        fmvec[0] = fmp->Pt();
        fmvec[1] = fmp->Eta();
        fmvec[2] = fmp->Phi();
        if(fmvec[2] < 0) fmvec[2] += TMath::TwoPi();
        fTrackMCFirstMotherMomentum.push_back(fmvec);
        fmvec.clear();
      }
      else{//no mother, i.e. no first mother, too.
        fTrackMCMotherIndex.push_back(mother_index);
        fTrackMCMotherPdgCode.push_back(mother_pdg);//pdgcode 0 does not exist.

        fTrackMCFirstMotherIndex.push_back(mother_index);
        fTrackMCFirstMotherPdgCode.push_back(mother_pdg);//pdgcode 0 does not exist.

        vector<Float_t> fmvec(3,0);
        fmvec[0] = -999;
        fmvec[1] = -999;
        fmvec[2] = -999;
        fTrackMCFirstMotherMomentum.push_back(fmvec);
        fmvec.clear();
      }
    }//end of hasMC

  }//end of track loop

}
//_______________________________________________________________________________________________
void AliAnalysisTaskReducedTreeDS::FillV0InfoESD()
{
  const AliVVertex *vVertex = fEvent->GetPrimaryVertex();
  AliKFVertex primaryVertexKF(*vVertex);
  //Double_t secVtx[3] = {primaryVertexKF.GetX(), primaryVertexKF.GetY(), primaryVertexKF.GetZ()};
  const Double_t Me  = TDatabasePDG::Instance()->GetParticle(11)->Mass();
  const Double_t Mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  const Double_t Mp  = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  Double_t M1 = 0;
  Double_t M2 = 0;
  Double_t V0xyz[3] = {0,0,0};
  Float_t Lxy = 0;

  fESDv0KineCuts->SetEvent(InputEvent());
  fESDv0KineCuts->SetPrimaryVertex(&primaryVertexKF);

  //FillV0Info
  const Int_t Nv0 = fEvent->GetNumberOfV0s();  
  for(Int_t iv0=0;iv0<Nv0;iv0++){
    AliESDv0 *v0 = (AliESDv0*)fESDEvent->GetV0(iv0);

    //if(!v0->GetOnFlyStatus()) continue;//select v0 reconstructed on the fly.

    //if(v0->PtArmV0() > 0.3) continue;//qT < 0.3
    //if(TMath::Abs(v0->AlphaV0()) > 1.0) continue;//|alpha| < 1.0
    //if(v0->GetV0CosineOfPointingAngle(secVtx[0],secVtx[1],secVtx[2]) < 0.998) continue;
    //if(v0->GetRr() < 2. || 60. < v0->GetRr()) continue;//in cm

    //Float_t dca = v0->GetDcaV0Daughters();
    //if(dca > 0.25) continue;

    AliESDtrack* legPos = fESDEvent->GetTrack(v0->GetPindex());
    AliESDtrack* legNeg = fESDEvent->GetTrack(v0->GetNindex());

    if(legPos->Charge() * legNeg->Charge() > 0) continue;//reject same sign pair

    //if(legPos->Charge() < 0 && legNeg->Charge() >0){
    //  //AliInfo("charge is swapped.");
    //  legPos = fESDEvent->GetTrack(v0->GetNindex());
    //  legNeg = fESDEvent->GetTrack(v0->GetPindex());
    //}

    if(legPos->Pt() < fMinPtCut) continue;
    if(legNeg->Pt() < fMinPtCut) continue;
    if(TMath::Abs(legPos->Eta()) > fMaxEtaCut) continue;
    if(TMath::Abs(legNeg->Eta()) > fMaxEtaCut) continue;

    Float_t DCAxy_leg = -999, DCAz_leg = -999;
    legPos->GetImpactParameters(DCAxy_leg,DCAz_leg);
    if(TMath::Abs(DCAxy_leg) > 1.) continue;
    if(TMath::Abs(DCAz_leg)  > 3.) continue;

    DCAxy_leg = -999; DCAz_leg = -999;
    legNeg->GetImpactParameters(DCAxy_leg,DCAz_leg);
    if(TMath::Abs(DCAxy_leg) > 1.) continue;
    if(TMath::Abs(DCAz_leg)  > 3.) continue;

    if(legPos->GetKinkIndex(0) != 0) continue;
    if(legNeg->GetKinkIndex(0) != 0) continue;

    if(legPos->GetNcls(0) < 3) continue;//minimum number of ITS cluster 3
    if(legNeg->GetNcls(0) < 3) continue;//minimum number of ITS cluster 3
    if(!(legPos->GetStatus() & AliVTrack::kITSrefit)) continue;
    if(!(legNeg->GetStatus() & AliVTrack::kITSrefit)) continue;

    if((Double_t)(legPos->GetTPCchi2()) / (Double_t)(legPos->GetNcls(1)) > 4.) continue;//maximum chi2 per cluster TPC
    if((Double_t)(legNeg->GetTPCchi2()) / (Double_t)(legNeg->GetNcls(1)) > 4.) continue;//maximum chi2 per cluster TPC
    //if(legPos->GetNcls(1) < 70) continue;//minimum number of TPC cluster 70
    //if(legNeg->GetNcls(1) < 70) continue;//minimum number of TPC cluster 70
    if(legPos->GetTPCCrossedRows() < 70) continue;//minimum number of TPC crossed rows 70
    if(legNeg->GetTPCCrossedRows() < 70) continue;//minimum number of TPC crossed rows 70
    if(!(legPos->GetStatus() & AliVTrack::kTPCrefit)) continue;
    if(!(legNeg->GetStatus() & AliVTrack::kTPCrefit)) continue;

    Float_t ratio_pos = legPos->GetTPCNclsF() > 0 ? (Float_t)legPos->GetTPCCrossedRows() / (Float_t)legPos->GetTPCNclsF() : 1.0;
    Float_t ratio_neg = legNeg->GetTPCNclsF() > 0 ? (Float_t)legNeg->GetTPCCrossedRows() / (Float_t)legNeg->GetTPCNclsF() : 1.0;

    if(ratio_pos < 0.8) continue;
    if(ratio_neg < 0.8) continue;

    if(fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kElectron) < fMinTPCNsigmaEleCut) continue;
    if(fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kElectron) < fMinTPCNsigmaEleCut) continue;
    if(fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kElectron) > fMaxTPCNsigmaEleCut) continue;
    if(fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kElectron) > fMaxTPCNsigmaEleCut) continue;

    Int_t pdgV0 = 0, pdgP = 0, pdgN = 0;
    if(!fESDv0KineCuts->ProcessV0(v0,pdgV0,pdgP,pdgN)) continue;

    if(pdgV0 == 22 && TMath::Abs(pdgP) == 11 && TMath::Abs(pdgN) == 11){
      M1 = Me;
      M2 = Me;
    }
    else if(pdgV0 == 310 && TMath::Abs(pdgP) == 211 && TMath::Abs(pdgN) == 211){
      M1 = Mpi;
      M2 = Mpi;
    }
    else if(pdgV0 == 3122 && (TMath::Abs(pdgP) == 2212 || TMath::Abs(pdgP) == 211) && (TMath::Abs(pdgN) == 211 || TMath::Abs(pdgN) == 2212)){
      M1 = Mp;
      M2 = Mpi;

      if(pdgP == -211 && pdgN == 2212){//swapped
        legPos = fESDEvent->GetTrack(v0->GetNindex());//proton
        legNeg = fESDEvent->GetTrack(v0->GetPindex());//pi-
        M1 = Mpi;
        M2 = Mp;
      }
    }
    else if(pdgV0 == -3122 && (TMath::Abs(pdgP) == 2212 || TMath::Abs(pdgP) == 211) && (TMath::Abs(pdgN) == 211 || TMath::Abs(pdgN) == 2212)){
      M1 = Mpi;
      M2 = Mp;

      if(pdgP == -2212 && pdgN == 211){//swapped
        legPos = fESDEvent->GetTrack(v0->GetNindex());//pi+
        legNeg = fESDEvent->GetTrack(v0->GetPindex());//anti-proton
        M1 = Mp;
        M2 = Mpi;
      }
    }
    else continue;

    fV0OnFly.push_back(v0->GetOnFlyStatus());
    fV0Candidate.push_back(pdgV0);

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

    v0->GetXYZ(V0xyz[0],V0xyz[1],V0xyz[2]);
    Lxy = TMath::Sqrt(V0xyz[0]*V0xyz[0] + V0xyz[1]*V0xyz[1]);
    fV0Lxy.push_back(Lxy);
    //fV0Lxy.push_back(v0->GetRr());
    
    vector<Float_t> legPin_tmp(2,0);
    legPin_tmp[0] = legPos->GetTPCmomentum();
    legPin_tmp[1] = legNeg->GetTPCmomentum();
    fV0legPin.push_back(legPin_tmp);
    legPin_tmp.clear();

    fV0alpha.push_back(v0->AlphaV0());
    fV0qT.push_back(v0->PtArmV0());

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
 
    fV0legPointOnITSLayer.push_back(legPointOnITS_tmp);
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
 
    fV0legSharedPointOnITSLayer.push_back(legSharedPointOnITS_tmp);
    legSharedPointOnITS_tmp.clear();

    fV0Mass.push_back(v0->GetEffMassExplicit(M1,M2));

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
    if(fHasMC){
      vector<Int_t> leglabel_tmp(2,0);
      leglabel_tmp[0] = TMath::Abs(legPos->GetLabel());
      leglabel_tmp[1] = TMath::Abs(legNeg->GetLabel());
      AliMCParticle *pPos = (AliMCParticle*)fMCEvent->GetTrack(leglabel_tmp[0]);
      AliMCParticle *pNeg = (AliMCParticle*)fMCEvent->GetTrack(leglabel_tmp[1]);
      fV0MClegIndex.push_back(leglabel_tmp);
      leglabel_tmp.clear();

      vector<Int_t> leggenID_tmp(2,0);
      leggenID_tmp[0] = pPos->GetGeneratorIndex();
      leggenID_tmp[1] = pNeg->GetGeneratorIndex();
      fV0MClegGeneratorIndex.push_back(leggenID_tmp);
      leggenID_tmp.clear();

      vector<Int_t> legpdg_tmp(2,0);
      legpdg_tmp[0] = pPos->PdgCode();
      legpdg_tmp[1] = pNeg->PdgCode();
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
      legPos_vec_mc.clear();
      legNeg_vec_mc.clear();

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
      legPos_vec_pv_mc.clear();
      legNeg_vec_pv_mc.clear();

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

      vector<vector<Float_t>> leg_fmp_mom_tmp;//0 for legPos, 1 for legNeg//production vertex
      vector<Float_t> fmvec_Pos(3,-999);
      vector<Float_t> fmvec_Neg(3,-999);

      if(leg_mother_index_tmp[0] > -1){//check mother for pos leg
        AliMCParticle *mpPos = (AliMCParticle*)fMCEvent->GetTrack(leg_mother_index_tmp[0]);
        leg_mother_pdg_tmp[0] = mpPos->PdgCode();
        leg_first_mother_pdg_tmp[0] = mpPos->PdgCode();

        leg_first_mother_index_tmp[0] = GetFirstMother(pPos);
        AliMCParticle *fmpPos = mpPos;
        if(leg_first_mother_index_tmp[0] > -1){
          fmpPos = (AliMCParticle*)fMCEvent->GetTrack(leg_first_mother_index_tmp[0]);
          leg_first_mother_pdg_tmp[0] = fmpPos->PdgCode();
        }

        fmvec_Pos[0] = fmpPos->Pt();
        fmvec_Pos[1] = fmpPos->Eta();
        fmvec_Pos[2] = fmpPos->Phi();
        if(fmvec_Pos[2] < 0) fmvec_Pos[2] += TMath::TwoPi();
        leg_fmp_mom_tmp.push_back(fmvec_Pos);

      }//end of check mother for pos leg

      if(leg_mother_index_tmp[1] > -1){//check mother for neg leg
        AliMCParticle *mpNeg = (AliMCParticle*)fMCEvent->GetTrack(leg_mother_index_tmp[1]);
        leg_mother_pdg_tmp[1] = mpNeg->PdgCode();
        leg_first_mother_pdg_tmp[1] = mpNeg->PdgCode();

        leg_first_mother_index_tmp[1] = GetFirstMother(pNeg);
        AliMCParticle *fmpNeg = mpNeg;
        if(leg_first_mother_index_tmp[1] > -1){
          fmpNeg = (AliMCParticle*)fMCEvent->GetTrack(leg_first_mother_index_tmp[1]);
          leg_first_mother_pdg_tmp[1] = fmpNeg->PdgCode();
        }

        fmvec_Neg[0] = fmpNeg->Pt();
        fmvec_Neg[1] = fmpNeg->Eta();
        fmvec_Neg[2] = fmpNeg->Phi();
        if(fmvec_Neg[2] < 0) fmvec_Neg[2] += TMath::TwoPi();
        leg_fmp_mom_tmp.push_back(fmvec_Neg);

      }//end of check mother for neg leg

      fV0MClegMotherIndex.push_back(leg_mother_index_tmp);
      fV0MClegMotherPdgCode.push_back(leg_mother_pdg_tmp);
      fV0MClegFirstMotherIndex.push_back(leg_first_mother_index_tmp);
      fV0MClegFirstMotherPdgCode.push_back(leg_first_mother_pdg_tmp);
      fV0MClegFirstMotherMomentum.push_back(leg_fmp_mom_tmp);

      leg_mother_index_tmp.clear();
      leg_mother_pdg_tmp.clear();
      leg_first_mother_index_tmp.clear();
      leg_first_mother_pdg_tmp.clear();
      leg_fmp_mom_tmp.clear();

      fmvec_Neg.clear();
      fmvec_Pos.clear();

    }//end of hasMC

  }//end of V0 loop
}
//_______________________________________________________________________________________________
void AliAnalysisTaskReducedTreeDS::FillV0InfoAOD()
{
  const AliVVertex *vVertex = fEvent->GetPrimaryVertex();
  AliKFVertex primaryVertexKF(*vVertex);
  //Double_t secVtx[3] = {primaryVertexKF.GetX(), primaryVertexKF.GetY(), primaryVertexKF.GetZ()};

  fAODv0KineCuts->SetEvent(InputEvent());
  fAODv0KineCuts->SetPrimaryVertex(&primaryVertexKF);

  const Double_t Me  = TDatabasePDG::Instance()->GetParticle(11)->Mass();
  const Double_t Mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  const Double_t Mp  = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  Double_t M1 = 0;
  Double_t M2 = 0;

  const Int_t Nv0 = fEvent->GetNumberOfV0s();  
  for(Int_t iv0=0;iv0<Nv0;iv0++){
    AliAODv0 *v0 = (AliAODv0*)fAODEvent->GetV0(iv0);

    //if(!v0->GetOnFlyStatus()) continue;//select v0 reconstructed on the fly.

    //if(v0->PtArmV0() > 0.3) continue;//qT < 0.3
    //if(TMath::Abs(v0->AlphaV0()) > 1.0) continue;//|alpha| < 1.0
    //if(v0->CosPointingAngle(secVtx) < 0.998) continue;
    //if(v0->RadiusV0() < 2. || 60. < v0->RadiusV0()) continue;//in cm
    //Float_t dca = v0->DcaV0Daughters();
    //if(dca > 0.25) continue;

    AliAODTrack *legPos = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(0));
    AliAODTrack *legNeg = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(1));
    if(legPos->Charge() * legNeg->Charge() > 0) continue;//reject same sign pair

    //if(legPos->Charge() < 0 && legNeg->Charge() > 0){//swap charge sign //index0 is expect to be positive leg, index1 to be negative.//protection
    //  //AliInfo("charge is swapped.");
    //  legPos = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(1));
    //  legNeg = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(0));
    //}

    if(legPos->Pt() < fMinPtCut) continue;
    if(legNeg->Pt() < fMinPtCut) continue;
    if(TMath::Abs(legPos->Eta()) > fMaxEtaCut) continue;
    if(TMath::Abs(legNeg->Eta()) > fMaxEtaCut) continue;

    Float_t DCAxy_leg = -999, DCAz_leg = -999;
    legPos->GetImpactParameters(DCAxy_leg,DCAz_leg);
    if(TMath::Abs(DCAxy_leg) > 1.) continue;
    if(TMath::Abs(DCAz_leg)  > 3.) continue;

    DCAxy_leg = -999; DCAz_leg = -999;
    legNeg->GetImpactParameters(DCAxy_leg,DCAz_leg);
    if(TMath::Abs(DCAxy_leg) > 1.) continue;
    if(TMath::Abs(DCAz_leg)  > 3.) continue;

    AliAODVertex *avp = (AliAODVertex*)legPos->GetProdVertex();
    AliAODVertex *avn = (AliAODVertex*)legNeg->GetProdVertex();
    if(avp->GetType() == AliAODVertex::kKink) continue;//reject kink
    if(avn->GetType() == AliAODVertex::kKink) continue;//reject kink

    if(legPos->GetNcls(0) < 3) continue;//minimum number of ITS cluster 3
    if(legNeg->GetNcls(0) < 3) continue;//minimum number of ITS cluster 3
    if(!(legPos->GetStatus() & AliVTrack::kITSrefit)) continue;
    if(!(legNeg->GetStatus() & AliVTrack::kITSrefit)) continue;

    if((Double_t)(legPos->GetTPCchi2()) / (Double_t)(legPos->GetNcls(1)) > 4.) continue;//maximum chi2 per cluster TPC
    if((Double_t)(legNeg->GetTPCchi2()) / (Double_t)(legNeg->GetNcls(1)) > 4.) continue;//maximum chi2 per cluster TPC
    //if(legPos->GetNcls(1) < 70) continue;//minimum number of TPC cluster 70
    //if(legNeg->GetNcls(1) < 70) continue;//minimum number of TPC cluster 70
    if(legPos->GetTPCCrossedRows() < 70) continue;//minimum number of TPC crossed rows 70
    if(legNeg->GetTPCCrossedRows() < 70) continue;//minimum number of TPC crossed rows 70
    if(!(legPos->GetStatus() & AliVTrack::kTPCrefit)) continue;
    if(!(legNeg->GetStatus() & AliVTrack::kTPCrefit)) continue;

    Float_t ratio_pos = legPos->GetTPCNclsF() > 0 ? (Float_t)legPos->GetTPCCrossedRows() / (Float_t)legPos->GetTPCNclsF() : 1.0;
    Float_t ratio_neg = legNeg->GetTPCNclsF() > 0 ? (Float_t)legNeg->GetTPCCrossedRows() / (Float_t)legNeg->GetTPCNclsF() : 1.0;

    if(ratio_pos < 0.8) continue;
    if(ratio_neg < 0.8) continue;

    if(fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kElectron) < fMinTPCNsigmaEleCut) continue;
    if(fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kElectron) < fMinTPCNsigmaEleCut) continue;
    if(fPIDResponse->NumberOfSigmasTPC(legPos,AliPID::kElectron) > fMaxTPCNsigmaEleCut) continue;
    if(fPIDResponse->NumberOfSigmasTPC(legNeg,AliPID::kElectron) > fMaxTPCNsigmaEleCut) continue;

    Int_t pdgV0 = 0, pdgP = 0, pdgN = 0;
    if(!fAODv0KineCuts->ProcessV0(v0,pdgV0,pdgP,pdgN)) continue;

    if(pdgV0 == 22 && TMath::Abs(pdgP) == 11 && TMath::Abs(pdgN) == 11){
      M1 = Me;
      M2 = Me;
    }
    else if(pdgV0 == 310 && TMath::Abs(pdgP) == 211 && TMath::Abs(pdgN) == 211){
      M1 = Mpi;
      M2 = Mpi;
    }
    else if(pdgV0 == 3122 && (TMath::Abs(pdgP) == 2212 || TMath::Abs(pdgP) == 211) && (TMath::Abs(pdgN) == 211 || TMath::Abs(pdgN) == 2212)){
      M1 = Mp;
      M2 = Mpi;

      if(pdgP == -211 && pdgN == 2212){//swapped
        legPos = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(1));//pi+
        legNeg = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(0));//anti-proton
        M1 = Mpi;
        M2 = Mp;
      }
    }
    else if(pdgV0 == -3122 && (TMath::Abs(pdgP) == 2212 || TMath::Abs(pdgP) == 211) && (TMath::Abs(pdgN) == 211 || TMath::Abs(pdgN) == 2212)){
      M1 = Mpi;
      M2 = Mp;

      if(pdgP == -2212 && pdgN == 211){//swapped
        legPos = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(1));//pi-
        legNeg = dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(0));//proton
        M1 = Mp;
        M2 = Mpi;
      }
    }
    else continue;

    //AliInfo(Form("Mpos = %e GeV/c2 , Mneg = %e GeV/c2",M1,M2));

    fV0OnFly.push_back(v0->GetOnFlyStatus());
    fV0Candidate.push_back(pdgV0);

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

    fV0Lxy.push_back(v0->RadiusV0());
    
    vector<Float_t> legPin_tmp(2,0);
    legPin_tmp[0] = legPos->GetTPCmomentum();
    legPin_tmp[1] = legNeg->GetTPCmomentum();
    fV0legPin.push_back(legPin_tmp);
    legPin_tmp.clear();

    fV0alpha.push_back(v0->AlphaV0());
    fV0qT.push_back(v0->PtArmV0());

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
 
    fV0legPointOnITSLayer.push_back(legPointOnITS_tmp);
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
 
    fV0legSharedPointOnITSLayer.push_back(legSharedPointOnITS_tmp);
    legSharedPointOnITS_tmp.clear();

    fV0Mass.push_back(v0->InvMass2Prongs(0,1,TMath::Abs(pdgP),TMath::Abs(pdgN)));

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
    if(fHasMC){
      vector<Int_t> leglabel_tmp(2,0);
      leglabel_tmp[0] = TMath::Abs(legPos->GetLabel());
      leglabel_tmp[1] = TMath::Abs(legNeg->GetLabel());
      AliMCParticle *pPos = (AliMCParticle*)fMCEvent->GetTrack(leglabel_tmp[0]);
      AliMCParticle *pNeg = (AliMCParticle*)fMCEvent->GetTrack(leglabel_tmp[1]);
      fV0MClegIndex.push_back(leglabel_tmp);
      leglabel_tmp.clear();

      vector<Int_t> leggenID_tmp(2,0);
      leggenID_tmp[0] = pPos->GetGeneratorIndex();
      leggenID_tmp[1] = pNeg->GetGeneratorIndex();
      fV0MClegGeneratorIndex.push_back(leggenID_tmp);
      leggenID_tmp.clear();

      vector<Int_t> legpdg_tmp(2,0);
      legpdg_tmp[0] = pPos->PdgCode();
      legpdg_tmp[1] = pNeg->PdgCode();
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
      legPos_vec_mc.clear();
      legNeg_vec_mc.clear();

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
      legPos_vec_pv_mc.clear();
      legNeg_vec_pv_mc.clear();

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

      vector<vector<Float_t>> leg_fmp_mom_tmp;//0 for legPos, 1 for legNeg//production vertex
      vector<Float_t> fmvec_Pos(3,-999);
      vector<Float_t> fmvec_Neg(3,-999);

      if(leg_mother_index_tmp[0] > -1){//check mother for pos leg
        AliMCParticle *mpPos = (AliMCParticle*)fMCEvent->GetTrack(leg_mother_index_tmp[0]);
        leg_mother_pdg_tmp[0] = mpPos->PdgCode();
        leg_first_mother_pdg_tmp[0] = mpPos->PdgCode();

        leg_first_mother_index_tmp[0] = GetFirstMother(pPos);
        AliMCParticle *fmpPos = mpPos;
        if(leg_first_mother_index_tmp[0] > -1){
          fmpPos = (AliMCParticle*)fMCEvent->GetTrack(leg_first_mother_index_tmp[0]);
          leg_first_mother_pdg_tmp[0] = fmpPos->PdgCode();
        }

        fmvec_Pos[0] = fmpPos->Pt();
        fmvec_Pos[1] = fmpPos->Eta();
        fmvec_Pos[2] = fmpPos->Phi();
        if(fmvec_Pos[2] < 0) fmvec_Pos[2] += TMath::TwoPi();
        leg_fmp_mom_tmp.push_back(fmvec_Pos);

      }//end of check mother for pos leg

      if(leg_mother_index_tmp[1] > -1){//check mother for neg leg
        AliMCParticle *mpNeg = (AliMCParticle*)fMCEvent->GetTrack(leg_mother_index_tmp[1]);
        leg_mother_pdg_tmp[1] = mpNeg->PdgCode();
        leg_first_mother_pdg_tmp[1] = mpNeg->PdgCode();

        leg_first_mother_index_tmp[1] = GetFirstMother(pNeg);
        AliMCParticle *fmpNeg = mpNeg;
        if(leg_first_mother_index_tmp[1] > -1){
          fmpNeg = (AliMCParticle*)fMCEvent->GetTrack(leg_first_mother_index_tmp[1]);
          leg_first_mother_pdg_tmp[1] = fmpNeg->PdgCode();
        }

        fmvec_Neg[0] = fmpNeg->Pt();
        fmvec_Neg[1] = fmpNeg->Eta();
        fmvec_Neg[2] = fmpNeg->Phi();
        if(fmvec_Neg[2] < 0) fmvec_Neg[2] += TMath::TwoPi();
        leg_fmp_mom_tmp.push_back(fmvec_Neg);

      }//end of check mother for neg leg

      fV0MClegMotherIndex.push_back(leg_mother_index_tmp);
      fV0MClegMotherPdgCode.push_back(leg_mother_pdg_tmp);
      fV0MClegFirstMotherIndex.push_back(leg_first_mother_index_tmp);
      fV0MClegFirstMotherPdgCode.push_back(leg_first_mother_pdg_tmp);
      fV0MClegFirstMotherMomentum.push_back(leg_fmp_mom_tmp);

      leg_mother_index_tmp.clear();
      leg_mother_pdg_tmp.clear();
      leg_first_mother_index_tmp.clear();
      leg_first_mother_pdg_tmp.clear();
      leg_fmp_mom_tmp.clear();

      fmvec_Neg.clear();
      fmvec_Pos.clear();

    }//end of hasMC

  }//end of V0 loop
}
//_______________________________________________________________________________________________
void AliAnalysisTaskReducedTreeDS::ProcessMC(Option_t *option)
{
  const AliVVertex *MCvertex = fMCEvent->GetPrimaryVertex();
  fMCVertex[0] = MCvertex->GetX();//true vertex
  fMCVertex[1] = MCvertex->GetY();//true vertex
  fMCVertex[2] = MCvertex->GetZ();//true vertex
  //AliInfo(Form("MC vertex Z = %f (cm)",fMCVertex[2]));

  Int_t pdg = -1;
  Float_t pT = 0, eta = 0, phi = 0;
  Int_t genID = -1;

  const Int_t Ntrack_all = fMCEvent->GetNumberOfTracks();
  const Int_t Ntrack     = fMCEvent->GetNumberOfPrimaries();
  AliInfo(Form("N all track = %d , N primary track in MC = %d",Ntrack_all,Ntrack));

  if(fESDEvent){
    AliGenEventHeader* genHeader = fMCEvent->GenEventHeader();
    AliGenHijingEventHeader* hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(genHeader);

    if(hijingGenHeader == NULL){//HIJING header is not found
      AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(genHeader);
      TList* headerList = genCocktailHeader->GetHeaders();
      const Int_t Ngen = headerList->GetEntries();
      AliInfo(Form("N generators = %d",Ngen)); 
      for (Int_t igen=0; igen<Ngen; igen++) {
        AliGenEventHeader *gh = (AliGenEventHeader*)headerList->At(igen);
        AliInfo(Form("Cocktail header is found : Generator name = %s , NProduced = %d.",gh->GetName(),gh->NProduced()));
        TString genname = gh->GetName();
        fMCGeneratorName.push_back(genname);
      }
    }
    else{
      AliInfo(Form("HIJING header is found : Generator name = %s , NProduced = %d.",hijingGenHeader->GetName(),hijingGenHeader->NProduced()));
      TString genname = hijingGenHeader->GetName();
      fMCGeneratorName.push_back(genname);
    }
  }//end of ESD
  else if(fAODEvent){
    AliAODMCHeader* mcHeader = (AliAODMCHeader*)fAODEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    TList *headerList = mcHeader->GetCocktailHeaders();
    const Int_t Ngen = headerList->GetEntries();
    AliInfo(Form("N generators = %d",Ngen)); 
    for(Int_t igen=0;igen<Ngen;igen++){
      AliGenEventHeader *gh = (AliGenEventHeader*)headerList->At(igen);
      TString genname = gh->GetName();
      AliInfo(Form("Generator name = %s , NProduced = %d.",genname.Data(),gh->NProduced()));
      fMCGeneratorName.push_back(genname);
    }//end of generator loop
  }//end of AOD

  for(Int_t itrack=0;itrack<Ntrack;itrack++){
    AliMCParticle *p = (AliMCParticle*)fMCEvent->GetTrack(itrack);
    pdg = p->PdgCode();
    pT  = p->Pt();
    eta = p->Eta();//pseudo-rapidity
    phi = p->Phi();
    genID = p->GetGeneratorIndex();
    //printf("   MC i = %d , genID = %d , pdg = %d , IsPhysicalPrimary() = %d , mother index = %d , first mother index = %d , pT = %f GeV/c , eta = %f , phi = %f rad\n",itrack,genID,pdg,p->IsPhysicalPrimary(),p->GetMother(),GetFirstMother(p),pT,eta,phi);

    //AliAODMCParticle *aodp = (AliAODMCParticle*)fMCEvent->GetTrack(itrack);
    //printf("AODMC i = %d , genID = %d , pdg = %d , IsPhysicalPrimary() = %d , mother index = %d , first mother index = %d , pT = %f GeV/c , eta = %f , phi = %f rad\n",itrack,aodp->GetGeneratorIndex(),aodp->PdgCode(),aodp->IsPhysicalPrimary(),aodp->GetMother(),GetFirstMother(aodp),aodp->Pt(),aodp->Eta(),aodp->Phi());

    if(TMath::Abs(pdg) != 11) continue;//select only electrons

    if(!IsPrimaryElectron(p)) continue;

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

    fMCGeneratorIndex.push_back(genID);
    fMCIsPhysicalPrimary.push_back(p->IsPhysicalPrimary());
    fMCIndex.push_back(itrack);
    fMCPdgCode.push_back(pdg);//11 or -11

    //check mother
    Int_t mother_index = p->GetMother();
    Int_t mother_pdg   = 0;//pdgcode 0 does not exist.
    Int_t first_mother_index = p->GetMother();
    Int_t first_mother_pdg   = 0;
    if(mother_index > -1){
      AliMCParticle *mp = (AliMCParticle*)fMCEvent->GetTrack(mother_index);
      mother_pdg = mp->PdgCode();
      fMCMotherIndex.push_back(mother_index);
      fMCMotherPdgCode.push_back(mother_pdg);
      //printf("mother_index = %d , mother_pdg = %d\n",mother_index,mother_pdg);

      //check first mother
      first_mother_index = GetFirstMother(p);
      first_mother_pdg   = mp->PdgCode();
      AliMCParticle *fmp = mp;
      if(first_mother_index > -1){
        fmp = (AliMCParticle*)fMCEvent->GetTrack(first_mother_index);
        first_mother_pdg = fmp->PdgCode();
      } 
      fMCFirstMotherIndex.push_back(first_mother_index);
      fMCFirstMotherPdgCode.push_back(first_mother_pdg);

      vector<Float_t> fmvec(3,0);
      fmvec[0] = fmp->Pt();
      fmvec[1] = fmp->Eta();
      fmvec[2] = fmp->Phi();
      if(fmvec[2] < 0) fmvec[2] += TMath::TwoPi();
      fMCFirstMotherMomentum.push_back(fmvec);
      fmvec.clear();
    }
    else{//no mother, i.e. no first mother, too.
      fMCMotherIndex.push_back(mother_index);
      fMCMotherPdgCode.push_back(mother_pdg);//pdgcode 0 does not exist.

      fMCFirstMotherIndex.push_back(mother_index);
      fMCFirstMotherPdgCode.push_back(mother_pdg);//pdgcode 0 does not exist.

      vector<Float_t> fmvec(3,0);
      fmvec[0] = -999;
      fmvec[1] = -999;
      fmvec[2] = -999;
      fMCFirstMotherMomentum.push_back(fmvec);
      fmvec.clear();
    }

  }//end of MC generated track loop

}
//_______________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeDS::IsPrimaryElectron(AliVParticle *p)
{
  Int_t pdg = p->PdgCode();
  if(TMath::Abs(pdg) != 11) return kFALSE;//not electron

  //Bool_t isEW = IsEWBoson(p);//not necessary
  Bool_t isLF = IsLF(p);
  Bool_t isHF = IsSemileptonicDecayFromHF(p);

  //if(isEW || isLF || isHF) return kTRUE;
  if(isLF || isHF) return kTRUE;
  else return kFALSE; 
}
//_______________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeDS::IsEWBoson(AliVParticle *p)
{
  Int_t pdg = p->PdgCode();
  if(TMath::Abs(pdg) != 11) return kFALSE;//not electron

  Int_t mother_index = p->GetMother();
  if(mother_index < 0) return kFALSE;
  AliMCParticle *mp = (AliMCParticle*)fMCEvent->GetTrack(mother_index);
  Int_t mother_pdg = mp->PdgCode();

  if(mother_pdg == 22){
    //AliInfo("mother is photon. return kTRUE.");
    return kTRUE;
  }
  else if(mother_pdg == 23){
    //AliInfo("mother is Z boson. return kTRUE.");
    return kTRUE;
  }
  else if(mother_pdg == 24){
    //AliInfo("mother is W boson. return kTRUE.");
    return kTRUE;
  }
  return kFALSE;
}
//_______________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeDS::IsLF(AliVParticle *p)
{ 
  Int_t pdg = p->PdgCode();
  if(TMath::Abs(pdg) != 11) return kFALSE;//not electron

  //J/psi is in LF injected M.C.
  const Int_t meson[] = {111, 221, 331, 113, 223, 333, 443};//pi0, eta, eta', rho(770), omega(782), phi(1020), J/psi
  const Int_t Nm = sizeof(meson)/sizeof(meson[0]);

  Int_t mother_index = p->GetMother();
  if(mother_index < 0) return kFALSE;
  AliMCParticle *mp = (AliMCParticle*)fMCEvent->GetTrack(mother_index);
  Int_t mother_pdg = mp->PdgCode();

//  Int_t first_mother_index = GetFirstMother(p);
//  Int_t first_mother_pdg = 0;
//
//  if(mother_index != first_mother_index) return kFALSE;//reject vector meson->pi0->e
//  if(first_mother_index > -1){
//    AliMCParticle *fmp = (AliMCParticle*)fMCEvent->GetTrack(first_mother_index);
//    first_mother_pdg = fmp->PdgCode();
//  }

  //Double_t dx = mp->Xv() - fMCVertex[0];
  //Double_t dy = mp->Yv() - fMCVertex[1];
  //Double_t R   = TMath::Sqrt(dx*dx + dy*dy);
  //AliInfo(Form("R in X-Y plane = %e cm.",R));
  //if(R > 1.) return kFALSE;//select electrons from primary vertex.//1cm

  for(Int_t i=0;i<Nm;i++){
    if(TMath::Abs(mother_pdg) == meson[i]){
      //AliInfo(Form("Match with %d at Rxy = %e cm. return kTRUE.",meson[i],R));
      return kTRUE;
    }
  }//end of meson loop

  return kFALSE;
}
//_______________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeDS::IsSemileptonicDecayFromHF(AliVParticle *p)
{
  Int_t pdg = p->PdgCode();
  if(TMath::Abs(pdg) != 11) return kFALSE;//not electron

  Int_t mother_index = p->GetMother();
  if(mother_index < 0) return kFALSE;
  AliMCParticle *mp = (AliMCParticle*)fMCEvent->GetTrack(mother_index);
  Int_t mother_pdg = mp->PdgCode();

  //the last digit of pdg code is 2s+1.

  if(TMath::Abs(mother_pdg) < 100){
    //AliInfo(Form("PDG code %d is not a hadron. return kFALSE.",pdg));
    return kFALSE;
  }
 
  if(mother_pdg == 130 || mother_pdg == 310){//standard rule 2s+1 is broken for K0S and K0L.
    //AliInfo("reject K0S - K0L mixing");
    return kFALSE;
  }
 
  TString str = TString::Itoa(mother_pdg,10);
  Int_t len = str.Length();
  //Double_t dx = mp->Xv() - fMCVertex[0];
  //Double_t dy = mp->Yv() - fMCVertex[1];
  //Double_t R  = TMath::Sqrt(dx*dx + dy*dy);

  if(mother_pdg %2 == 0){//baryon
    TString str_q1 = str[len-4];
    TString str_q2 = str[len-3];
    TString str_q3 = str[len-2];

    Int_t q1 = str_q1.Atoi();
    Int_t q2 = str_q2.Atoi();
    Int_t q3 = str_q3.Atoi();

    if(q1==5 || q2==5 || q3==5){
      //AliInfo(Form("Bottom baryon %d is matched at Rxy = %e. return kTRUE.",mother_pdg,R));
      return kTRUE;
    }
    if(q1==4 || q2==4 || q3==4){
      //AliInfo(Form("Charmed baryon %d is matched at Rxy = %e. return kTRUE.",mother_pdg,R));
      return kTRUE;
    }
  }
  else{//meson
    TString str_q1 = str[len-3];
    TString str_q2 = str[len-2];
    //reject quarkonia

    Int_t q1 = str_q1.Atoi();
    Int_t q2 = str_q2.Atoi();
    if((q1==5) ^ (q2==5)){
      //AliInfo(Form("Bottom meson %d is matched at Rxy = %e. return kTRUE.",mother_pdg,R));
      return kTRUE;
    }
    if((q1==4) ^ (q2==4)){
      //AliInfo(Form("Charmed meson %d is matched at Rxy = %e. return kTRUE.",mother_pdg,R));
      return kTRUE;
    }
  }
  return kFALSE;//only for protection
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
    //AliInfo("charge is swapped.");
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
  //AliInfo("Number of elements of vectors is cleared.");

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

  fTPCsignalN.clear();
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
  fTrackMCGeneratorIndex.clear();
  fTrackMCIsPhysicalPrimary.clear();
  fTrackMCIsSecondaryFromMaterial.clear();
  fTrackMCIsSecondaryFromWeakDecay.clear();
  fTrackMCIndex.clear();
  fTrackMCPdgCode.clear();
  fTrackMCMotherIndex.clear();
  fTrackMCMotherPdgCode.clear();
  fTrackMCFirstMotherIndex.clear();
  fTrackMCFirstMotherPdgCode.clear();
  fTrackMCFirstMotherMomentum.clear();

  //clear V0 variables
  fV0OnFly.clear();
  fV0legMomentum.clear();
  fV0legPin.clear();
  fV0Lxy.clear();
  fV0alpha.clear();
  fV0qT.clear();
  fV0Candidate.clear();
  fV0Mass.clear();
  fV0legPointOnITSLayer.clear();
  fV0legSharedPointOnITSLayer.clear();
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
  fV0MClegGeneratorIndex.clear();
  fV0MClegIndex.clear();
  fV0MClegPdgCode.clear();
  fV0MClegMotherIndex.clear();
  fV0MClegMotherPdgCode.clear();
  fV0MClegFirstMotherIndex.clear();
  fV0MClegFirstMotherPdgCode.clear();
  fV0MClegFirstMotherMomentum.clear();

  //clear MC variables
  fMCMomentum.clear();
  fMCProdVtx.clear();
  fMCGeneratorIndex.clear();
  fMCGeneratorName.clear();
  fMCIsPhysicalPrimary.clear();
  fMCIndex.clear();
  fMCPdgCode.clear();
  fMCMotherIndex.clear();
  fMCMotherPdgCode.clear();
  fMCFirstMotherIndex.clear();
  fMCFirstMotherPdgCode.clear();
  fMCFirstMotherMomentum.clear();
  //AliInfo("Number of elements of vectors is cleared. DONE!");

}
//_______________________________________________________________________________________________
void AliAnalysisTaskReducedTreeDS::ClearVectorMemory()
{
  //AliInfo("Memories for vector objects are swapped with null.");

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

  vector<Int_t>().swap(fTPCsignalN);
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
  vector<Int_t>().swap(fTrackMCGeneratorIndex);
  vector<Bool_t>().swap(fTrackMCIsPhysicalPrimary);
  vector<Bool_t>().swap(fTrackMCIsSecondaryFromMaterial);
  vector<Bool_t>().swap(fTrackMCIsSecondaryFromWeakDecay);
  vector<Int_t>().swap(fTrackMCIndex);
  vector<Int_t>().swap(fTrackMCPdgCode);
  vector<Int_t>().swap(fTrackMCMotherIndex);
  vector<Int_t>().swap(fTrackMCMotherPdgCode);
  vector<Int_t>().swap(fTrackMCFirstMotherIndex);
  vector<Int_t>().swap(fTrackMCFirstMotherPdgCode);
  vector<vector<Float_t>>().swap(fTrackMCFirstMotherMomentum);

  //V0 info
  vector<Bool_t>().swap(fV0OnFly);
  vector<vector<vector<Float_t>>>().swap(fV0legMomentum);
  vector<vector<Float_t>>().swap(fV0legPin);
  vector<Float_t>().swap(fV0Lxy);
  vector<Float_t>().swap(fV0alpha);
  vector<Float_t>().swap(fV0qT);
  vector<Int_t>().swap(fV0Candidate);
  vector<Float_t>().swap(fV0Mass);
  vector<vector<vector<Bool_t>>>().swap(fV0legPointOnITSLayer);
  vector<vector<vector<Bool_t>>>().swap(fV0legSharedPointOnITSLayer);
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
  vector<vector<Int_t>>().swap(fV0MClegGeneratorIndex);
  vector<vector<Int_t>>().swap(fV0MClegIndex);
  vector<vector<Int_t>>().swap(fV0MClegPdgCode);
  vector<vector<Int_t>>().swap(fV0MClegMotherIndex);
  vector<vector<Int_t>>().swap(fV0MClegMotherPdgCode);
  vector<vector<Int_t>>().swap(fV0MClegFirstMotherIndex);
  vector<vector<Int_t>>().swap(fV0MClegFirstMotherPdgCode);
  vector<vector<vector<Float_t>>>().swap(fV0MClegFirstMotherMomentum);

  //MC variables for true electrons
  vector<vector<Float_t>>().swap(fMCMomentum);
  vector<vector<Float_t>>().swap(fMCProdVtx);
  vector<Int_t>().swap(fMCGeneratorIndex);
  vector<TString>().swap(fMCGeneratorName);
  vector<Bool_t>().swap(fMCIsPhysicalPrimary);
  vector<Int_t>().swap(fMCIndex);
  vector<Int_t>().swap(fMCPdgCode);
  vector<Int_t>().swap(fMCMotherIndex);
  vector<Int_t>().swap(fMCMotherPdgCode);
  vector<Int_t>().swap(fMCFirstMotherIndex);
  vector<Int_t>().swap(fMCFirstMotherPdgCode);
  vector<vector<Float_t>>().swap(fMCFirstMotherMomentum);

  //AliInfo("Memories for vector objects are swapped with null. DONE!");
}
//_______________________________________________________________________________________________
void AliAnalysisTaskReducedTreeDS::ExtractQnVectors()
{
  //AliInfo("extract Qn vectors from Qn correction framework.");
  const TString TPCEPname[3] = {"TPC","TPCNegEta","TPCPosEta"};
  const TString V0EPname[3]  = {"VZERO","VZEROA","VZEROC"};
  const TString ZDCEPname[2] = {"ZDCA","ZDCC"};
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
    //AliInfo("fFlowQnVectorMgr does not exist. return.");
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

