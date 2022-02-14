// $Id:
//
// Analysis task for neutral pions (into two gammas), and for direct photons by subtraction method
//
// Author: B. Sahlmueller, based on code by C. Loizides

#include "AliAnalysisTaskEMCALPi0Gamma.h"
#include <TAxis.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TProfile.h>
#include <TRegexp.h>
#include <TString.h>
#include <TVector2.h>
#include <TArray.h>
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliCDBManager.h"
#include "AliCentrality.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDCaloTrigger.h"
#include "AliESDEvent.h"
#include "AliESDUtils.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"
#include "AliEventplane.h"
#include "AliGeomManager.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliMultiplicity.h"
#include "AliMultSelection.h"
#include "AliStack.h"
#include "AliTrackerBase.h"
#include "AliTriggerAnalysis.h"
#include "AliAODConversionPhoton.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"

using std::cout;
using std::endl;
using std::max;

ClassImp(AliAnalysisTaskEMCALPi0Gamma)

//________________________________________________________________________
AliAnalysisTaskEMCALPi0Gamma::AliAnalysisTaskEMCALPi0Gamma()
: AliAnalysisTaskSE(),
fCentVar("V0M"),
fCentFrom(0),
fCentTo(100),
fVtxZMin(-10),
fVtxZMax(+10),
fUseQualFlag(1),
fClusName(),
fDoNtuple(0),
fDoAfterburner(0),
fAsymMax1(0.3),
fAsymMax2(0.7),
fAsymMax3(1),
fNminCells(1),
fMinE(0.100),
fM02(100),
fMinErat(0),
fMinEcc(-1),
fDoTrMtSmpl(0),
fDoManualRecal(0),
fCalibRun(0),
fApplyBadMapManually(0),
fGeoName("EMCAL_FIRSTYEARV1"),
fMinNClusPerTr(50),
fIsoDist(0.2),
fTrClassNames(""),
fTrCuts(0),
fPrimTrCuts(0),
fPrimTracksName(""),
fDoTrMatGeom(0),
fTrainMode(0),
fMarkCells(),
fMinL0Time(-1),
fMaxL0Time(1024),
fMcMode(0),
fEmbedMode(0),
fGeom(0),
fReco(0),
fTrigName(),
fDoPSel(kFALSE),
fIsGeoMatsSet(0),
fRotateMixed(0),
fAddedSignal(0),
fSimStudies(0),
fDataPeriod(0),
fNEvs(0),
fOutput(0),
fTrClassNamesArr(0),
fEsdEv(0),
fAodEv(0),
fRecPoints(0),
fDigits(0),
fEsdClusters(0),
fEsdCells(0),
fAodClusters(0),
fAodCells(0),
fPtRanges(0),
fSelTracks(0),
fSelPrimTracks(0),
fNtuple(0),
fHeader(0),
fPrimVert(0),
fSpdVert(0),
fTpcVert(0),
fClusters(0),
fTriggers(0),
fMcParts(0),
fHCuts(0x0),
fHVertexZ(0x0),
fHVertexZ2(0x0),
fHCent(0x0),
fHCentQual(0x0),
fHMeanClusterEnergy(0x0),
fHMeanClusterNumber(0x0),
fHCellIndexEnergy(0x0),
fHCellIndexEnergyAfterCuts(0x0),
fHClusters(0x0),
fHClustAllEtaPhi(0x0),
fHClustNoEvt(0),
fHClustAccEvt(0),
fHClustEccentricity(0),
fHClustEtaPhi(0x0),
fHClustEtaPhiAll(0x0),
fHClustEnergyPt(0x0),
fHClustEnergyPtDCal(0x0),
fHClustEnergySM(0x0),
fHClustEnergySigma(0x0),
fHClustEnergyTime(0x0),
fHClustNCellEnergyRatio(0x0),
fHClustEnergyNCell(0x0),
fHPrimTrackPt(0x0),
fHPrimTrackEta(0x0),
fHPrimTrackPhi(0x0),
fHMatchDr(0x0),
fHMatchDz(0x0),
fHMatchEp(0x0),
fHPionEtaPhi(0x0),
fHAddPionEtaPt(0x0),
fHAddPionEtaPtWgt(0x0),
fHPyPionEtaPt(0x0),
fHdr(0),
fHPionMggPt(0x0),
fHPionMggAsym(0x0),
fHPionPtAsym(0x0),
fHPionMggDgg(0x0),
fHPionInvMasses(0x0),
fHPionInvMassesSym(0x0),
fHPionInvMassesAsym(0x0),
fHPionInvMassesAdd1(0x0),
fHPionInvMassesAdd1NoWgt(0x0),
fHPionInvMassesAdd1Sym(0x0),
fHPionInvMassesAdd1Asym(0x0),
fHPionInvMassesAdd1Mult(0x0),
fHPionInvMassesAdd1MultSym(0x0),
fHPionInvMassesAdd1MultAsym(0x0),
fHPionInvMassesGamAdd1(0x0),
fHPionInvMassesGamAdd1Sym(0x0),
fHPionInvMassesGamAdd1Asym(0x0),
fHPionInvMassesGamAdd1Mult(0x0),
fHPionInvMassesGamAdd1MultSym(0x0),
fHPionInvMassesGamAdd1MultAsym(0x0),
fHPionInvMassesAdd2(0x0),
fHPionInvMassesConvElZero(0x0),
fHPionInvMassesConvElOne(0x0),
fHPionInvMassesConvElBoth(0x0),
fHPionInvMassesChargedPiZero(0x0),
fHPionInvMassesChargedPiOne(0x0),
fHPionInvMassesChargedPiBoth(0x0),
fHPionInvMassesGammaBoth(0x0),
fHPionInvMassesMix(0x0),
fHPionInvMassesMix1(0x0),
fHPionInvMassesMix2(0x0),
fHPionInvMassesDCal(0x0),
fHPionInvMassesMixDCal(0x0),
fHPionInvMassesEMCalDCal(0x0),
fHPionInvMassesMixEMCalDCal(0x0),
fHPionInvMassesEMCalCalib(0x0),
fHPionInvMassesMixEMCalCalib(0x0),
fHPionInvMassesDCalCalib(0x0),
fHPionInvMassesMixDCalCalib(0x0),
//fHJPInvMasses(0x0),
fHPrimPionInvMasses(0x0),
fHPrimPionInvMassesAsym(0x0),
//fHConversionPoint(0),
fHPionTruthPt(0x0),
fHPionTruthPtIn(0x0),
fHPionTruthPtAcc(0x0),
fHEtaTruthPt(0x0),
fHEtaTruthPtIn(0x0),
fHEtaTruthPtAcc(0x0),
fHGamTruthPt(0x0),
fHGamTruthPtIn(0x0),
fHGamTruthPtAcc(0x0),
fHPionTruthPtAdd(0x0),
fHPionTruthPtInAdd(0x0),
fHPionTruthPtAccAdd(0x0),
fHEtaTruthPtAdd(0x0),
fHEtaTruthPtInAdd(0x0),
fHEtaTruthPtAccAdd(0x0),
fHNMothers(0x0),
//fHMixRotation(0x0),
ipymin(0),
ipymax(0),
ipi0min(0),
ipi0max(0),
ietamin(0),
ietamax(0),
eventHeader(0),
pythiaHeader(0),
addedPi0Header(0),
addedEtaHeader(0),
fHPriPionInvMasses(0x0),
fHSecPionInvMasses(0x0),
fHK0PionInvMasses(0x0),
fHMatPionInvMasses(0x0),
fHMCpartfrac(0),
fHECluEMC(0x0),
fHECluEMCAddPi0(0x0),
fHECluEMCAddEta(0x0),
//  fHRecTrue(0x0),
//  fHRecTrueAddPi0(0x0),
//  fHRecTrueAddEta(0x0),
fHECluEMCnofull(0x0),
fHECluEMCnofullAdd(0x0),
fHECluEMCelectron(0x0),
fHECluEMCpion(0x0),
fHECluEMCkaon(0x0),
fHECluEMCother(0x0),
fHECluEMCpi0single(0x0),
//  fHCorrection(0x0),
//  fHPionSm(0x0),
evt(),
thisEvent(),
fHWgt(0)
{
  for(int i=0;i<nMulClass;i++){
    for(int j=0;j<nZClass;j++){
      for(int k=0;k<nPtClass;k++){
        iEvt[i][j][k] = 0;
        for(int l=0;l<nEvt;l++){
          EmcEventList[i][j][k][l].SetGlobalInfo(0,0.,0.);
        }
      }
    }
  }
  // Set bad channel map
  Char_t key[55] ;
  snprintf(key,55,"BadMap") ;
  fBadMap=new TH1D(key,"Bad Modules map",18000,0.5,18000.5) ;
  
  // Constructor.
}

//________________________________________________________________________
AliAnalysisTaskEMCALPi0Gamma::AliAnalysisTaskEMCALPi0Gamma(const char *name)
: AliAnalysisTaskSE(name),
fCentVar("V0M"),
fCentFrom(0),
fCentTo(100),
fVtxZMin(-10),
fVtxZMax(+10),
fUseQualFlag(1),
fClusName(),
fDoNtuple(0),
fDoAfterburner(0),
fAsymMax1(0.3),
fAsymMax2(0.7),
fAsymMax3(1),
fNminCells(1),
fMinE(0.100),
fM02(100),
fMinErat(0),
fMinEcc(-1),
fDoTrMtSmpl(0),
fDoManualRecal(0),
fCalibRun(0),
fApplyBadMapManually(0),
fGeoName("EMCAL_FIRSTYEARV1"),
fMinNClusPerTr(50),
fIsoDist(0.2),
fTrClassNames(""),
fTrCuts(0),
fPrimTrCuts(0),
fPrimTracksName(""),
fDoTrMatGeom(0),
fTrainMode(0),
fMarkCells(),
fMinL0Time(-1),
fMaxL0Time(1024),
fMcMode(0),
fEmbedMode(0),
fGeom(0),
fReco(0),
fTrigName(),
fDoPSel(kFALSE),
fIsGeoMatsSet(0),
fRotateMixed(0),
fAddedSignal(0),
fSimStudies(0),
fDataPeriod(0),
fNEvs(0),
fOutput(0),
fTrClassNamesArr(0),
fEsdEv(0),
fAodEv(0),
fRecPoints(0),
fDigits(0),
fEsdClusters(0),
fEsdCells(0),
fAodClusters(0),
fAodCells(0),
fPtRanges(0),
fSelTracks(0),
fSelPrimTracks(0),
fNtuple(0),
fHeader(0),
fPrimVert(0),
fSpdVert(0),
fTpcVert(0),
fClusters(0),
fTriggers(0),
fMcParts(0),
fHCuts(0x0),
fHVertexZ(0x0),
fHVertexZ2(0x0),
fHCent(0x0),
fHCentQual(0x0),
fHMeanClusterEnergy(0x0),
fHMeanClusterNumber(0x0),
fHCellIndexEnergy(0x0),
fHCellIndexEnergyAfterCuts(0x0),
fHClusters(0x0),
fHClustAllEtaPhi(0x0),
fHClustNoEvt(0),
fHClustAccEvt(0),
fHClustEccentricity(0),
fHClustEtaPhi(0x0),
fHClustEtaPhiAll(0x0),
fHClustEnergyPt(0x0),
fHClustEnergyPtDCal(0x0),
fHClustEnergySM(0x0),
fHClustEnergySigma(0x0),
fHClustEnergyTime(0x0),
fHClustNCellEnergyRatio(0x0),
fHClustEnergyNCell(0x0),
fHPrimTrackPt(0x0),
fHPrimTrackEta(0x0),
fHPrimTrackPhi(0x0),
fHMatchDr(0x0),
fHMatchDz(0x0),
fHMatchEp(0x0),
fHPionEtaPhi(0x0),
fHAddPionEtaPt(0x0),
fHAddPionEtaPtWgt(0x0),
fHPyPionEtaPt(0x0),
fHdr(0),
fHPionMggPt(0x0),
fHPionMggAsym(0x0),
fHPionPtAsym(0x0),
fHPionMggDgg(0x0),
fHPionInvMasses(0x0),
fHPionInvMassesSym(0x0),
fHPionInvMassesAsym(0x0),
fHPionInvMassesAdd1(0x0),
fHPionInvMassesAdd1NoWgt(0x0),
fHPionInvMassesAdd1Sym(0x0),
fHPionInvMassesAdd1Asym(0x0),
fHPionInvMassesAdd1Mult(0x0),
fHPionInvMassesAdd1MultSym(0x0),
fHPionInvMassesAdd1MultAsym(0x0),
fHPionInvMassesAdd2(0x0),
fHPionInvMassesConvElZero(0x0),
fHPionInvMassesConvElOne(0x0),
fHPionInvMassesConvElBoth(0x0),
fHPionInvMassesChargedPiZero(0x0),
fHPionInvMassesChargedPiOne(0x0),
fHPionInvMassesChargedPiBoth(0x0),
fHPionInvMassesGammaBoth(0x0),
fHPionInvMassesMix(0x0),
fHPionInvMassesMix1(0x0),
fHPionInvMassesMix2(0x0),
fHPionInvMassesDCal(0x0),
fHPionInvMassesMixDCal(0x0),
fHPionInvMassesEMCalDCal(0x0),
fHPionInvMassesMixEMCalDCal(0x0),
fHPionInvMassesEMCalCalib(0x0),
fHPionInvMassesMixEMCalCalib(0x0),
fHPionInvMassesDCalCalib(0x0),
fHPionInvMassesMixDCalCalib(0x0),
//fHJPInvMasses(0x0),
fHPrimPionInvMasses(0x0),
fHPrimPionInvMassesAsym(0x0),
//fHConversionPoint(0),
fHPionTruthPt(0x0),
fHPionTruthPtIn(0x0),
fHPionTruthPtAcc(0x0),
fHEtaTruthPt(0x0),
fHEtaTruthPtIn(0x0),
fHEtaTruthPtAcc(0x0),
fHGamTruthPt(0x0),
fHGamTruthPtIn(0x0),
fHGamTruthPtAcc(0x0),
fHPionTruthPtAdd(0x0),
fHPionTruthPtInAdd(0x0),
fHPionTruthPtAccAdd(0x0),
fHEtaTruthPtAdd(0x0),
fHEtaTruthPtInAdd(0x0),
fHEtaTruthPtAccAdd(0x0),
fHNMothers(0x0),
//fHMixRotation(0x0),
ipymin(0),
ipymax(0),
ipi0min(0),
ipi0max(0),
ietamin(0),
ietamax(0),
eventHeader(0),
pythiaHeader(0),
addedPi0Header(0),
addedEtaHeader(0),
fHPriPionInvMasses(0x0),
fHSecPionInvMasses(0x0),
fHK0PionInvMasses(0x0),
fHMatPionInvMasses(0x0),
fHMCpartfrac(0),
fHECluEMC(0x0),
fHECluEMCAddPi0(0x0),
fHECluEMCAddEta(0x0),
//  fHRecTrue(0x0),
//  fHRecTrueAddPi0(0x0),
//  fHRecTrueAddEta(0x0),
fHECluEMCnofull(0x0),
fHECluEMCnofullAdd(0x0),
fHECluEMCelectron(0x0),
fHECluEMCpion(0x0),
fHECluEMCkaon(0x0),
fHECluEMCother(0x0),
fHECluEMCpi0single(0x0),
//  fHCorrection(0x0),
//  fHPionSm(0x0),
evt(),
thisEvent(),
fHPionInvMassesGamAdd1(0x0),
fHPionInvMassesGamAdd1Sym(0x0),
fHPionInvMassesGamAdd1Asym(0x0),
fHPionInvMassesGamAdd1Mult(0x0),
fHPionInvMassesGamAdd1MultSym(0x0),
fHPionInvMassesGamAdd1MultAsym(0x0),
fHWgt(0)
{
  // Constructor.
  for(int i=0;i<nMulClass;i++){
    for(int j=0;j<nZClass;j++){
      for(int k=0;k<nPtClass;k++){
        iEvt[i][j][k] = 0;
        for(int l=0;l<nEvt;l++){
          EmcEventList[i][j][k][l].SetGlobalInfo(0,0.,0.);
        }
      }
    }
  }
  // Set bad channel map
  Char_t key[55] ;
  snprintf(key,55,"BadMap") ;
  fBadMap=new TH1D(key,"Bad Modules map",18000,0.5,18000.5) ;
  
  DefineOutput(1, TList::Class());
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.,SPDVertex.,TPCVertex.,EMCALCells.,Tracks,EMCALTrigger.,SPDPileupVertices,TrkPileupVertices "
  "AOD:header,vertices,emcalCells,tracks";
}

//________________________________________________________________________
AliAnalysisTaskEMCALPi0Gamma::~AliAnalysisTaskEMCALPi0Gamma()
{
  // Destructor.
  
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::UserCreateOutputObjects()
{
  // Create user objects here.
  
  cout << "AliAnalysisTaskEMCALPi0Gamma: Input settings" << endl;
  cout << " fCentVar:       " << fCentVar << endl;
  cout << " fCentFrom:      " << fCentFrom << endl;
  cout << " fCentTo:        " << fCentTo << endl;
  cout << " fVtxZMin:       " << fVtxZMin << endl;
  cout << " fVtxZMax:       " << fVtxZMax << endl;
  cout << " fUseQualFlag:   " << fUseQualFlag << endl;
  cout << " fClusName:      \"" << fClusName << "\"" << endl;
  cout << " fDoNtuple:      " << fDoNtuple << endl;
  cout << " fDoAfterburner: " << fDoAfterburner << endl;
  cout << " fAsymMax1:       " << fAsymMax1 << endl;
  cout << " fAsymMax2:       " << fAsymMax2 << endl;
  cout << " fAsymMax3:       " << fAsymMax3 << endl;
  cout << " fNminCells:     " << fNminCells << endl;
  cout << " fMinE:          " << fMinE << endl;
  cout << " fMinErat:       " << fMinErat << endl;
  cout << " fMinEcc:        " << fMinEcc << endl;
  cout << " fM02:           " << fM02 << endl;
  cout << " fGeoName:       \"" << fGeoName << "\"" << endl;
  cout << " fMinNClusPerTr: " << fMinNClusPerTr << endl;
  cout << " fIsoDist:       " << fIsoDist << endl;
  cout << " fTrClassNames:  \"" << fTrClassNames << "\"" << endl;
  cout << " fTrCuts:        " << fTrCuts << endl;
  cout << " fPrimTrCuts:    " << fPrimTrCuts << endl;
  cout << " fDoTrMatGeom:   " << fDoTrMatGeom << endl;
  cout << " fTrainMode:     " << fTrainMode << endl;
  cout << " fMarkCells:     " << fMarkCells << endl;
  cout << " fMinL0Time:     " << fMinL0Time << endl;
  cout << " fMaxL0Time:     " << fMaxL0Time << endl;
  cout << " fMcMode:        " << fMcMode << endl;
  cout << " fEmbedMode:     " << fEmbedMode << endl;
  cout << " fGeom:          " << fGeom << endl;
  cout << " fReco:          " << fReco << endl;
  cout << " fTrigName:      " << fTrigName << endl;
  cout << " fDoPSel:        " << fDoPSel << endl;
  cout << " fRotateMixed:   " << fRotateMixed << endl;
  cout << " fDataPeriod:   " << fDataPeriod << endl;
  
  if (!fGeom)
    fGeom = AliEMCALGeometry::GetInstance(fGeoName);
  else {
    if (fGeom->GetMatrixForSuperModule(0))
      fIsGeoMatsSet = kTRUE;
  }
  if (!fReco)
    fReco = new AliEMCALRecoUtils();
  fTrClassNamesArr = fTrClassNames.Tokenize(" ");
  fOutput = new TList();
  fOutput->SetOwner();
  fSelTracks = new TObjArray;
  fSelPrimTracks = new TObjArray;
  if(fMcMode){
    if (TClass::GetClass("AliStaPart"))
      fMcParts = new TClonesArray("AliStaPart");
  }
  
  if (fDoNtuple) {
    TFile *f = OpenFile(1);
    TDirectory::TContext context(f);
    if (f) {
      f->SetCompressionLevel(2);
      fNtuple = new TTree(Form("tree%.0fto%.0f",fCentFrom,fCentTo), "StandaloneTree");
      fNtuple->SetDirectory(f);
      if (fTrainMode) {
        fNtuple->SetAutoFlush(-2*1024*1024);
        fNtuple->SetAutoSave(0);
      } else {
        fNtuple->SetAutoFlush(-32*1024*1024);
        fNtuple->SetAutoSave(0);
      }
      
      fHeader = new AliStaHeader;
      fNtuple->Branch("header", &fHeader, 16*1024, 99);
      fPrimVert = new AliStaVertex;
      fNtuple->Branch("primv", &fPrimVert, 16*1024, 99);
      fSpdVert = new AliStaVertex;
      fNtuple->Branch("spdv", &fSpdVert, 16*1024, 99);
      fTpcVert = new AliStaVertex;
      fNtuple->Branch("tpcv", &fTpcVert, 16*1024, 99);
      if (TClass::GetClass("AliStaCluster"))
        //TClass::GetClass("AliStaCluster")->IgnoreTObjectStreamer();
        fClusters = new TClonesArray("AliStaCluster");
      fNtuple->Branch("clusters", &fClusters, 8*16*1024, 99);
      if (TClass::GetClass("AliStaTrigger"))
        //TClass::GetClass("AliStaTrigger")->IgnoreTObjectStreamer();
        fTriggers = new TClonesArray("AliStaTrigger");
      fNtuple->Branch("l0prim", &fTriggers, 16*1024, 99);
      if (fMcMode||fEmbedMode) {
        fNtuple->Branch("mcparts", &fMcParts, 8*16*1024, 99);
      }
    }
  }
  
  AliEMCALEMCGeometry *emc = fGeom->GetEMCGeometry();
  Double_t phimin = emc->GetArm1PhiMin()*TMath::DegToRad();
  Double_t phimax = emc->GetArm1PhiMax()*TMath::DegToRad();
  
  //  if(phimin > TMath::Pi()){
  //    phimin -= 2*TMath::Pi();
  //  }
  //
  //  if(phimax > TMath::Pi()){
  //    phimax -= 2*TMath::Pi();
  //  }
  //
  //  if(phimin > phimax){
  //    Double_t phitmp = phimin;
  //    phimin = phimax;
  //    phimax = phitmp;
  //  }
  
  cout << "phimin=" << phimin << ", phimax=" << phimax << endl;
  
  // histograms
  Bool_t th1 =   TH1::GetDefaultSumw2();
  TH1::SetDefaultSumw2(kTRUE);
  Bool_t th2 =   TH2::GetDefaultSumw2();
  TH2::SetDefaultSumw2(kTRUE);
  fHCuts = new TH1F("hEventCuts","",6,0.5,6.5);
  fHCuts->GetXaxis()->SetBinLabel(1,"All");
  fHCuts->GetXaxis()->SetBinLabel(2,"PS");
  fHCuts->GetXaxis()->SetBinLabel(3,Form("%s: %.0f-%.0f",fCentVar.Data(),fCentFrom,fCentTo));
  fHCuts->GetXaxis()->SetBinLabel(4,"Trig");
  fHCuts->GetXaxis()->SetBinLabel(5,"QFlag");
  fHCuts->GetXaxis()->SetBinLabel(6,Form("zvtx: %.0f-%.0f",fVtxZMin,fVtxZMax));
  fOutput->Add(fHCuts);
  fHVertexZ = new TH1F("hVertexZBeforeCut","",100,-25,25);
  fHVertexZ->SetXTitle("z [cm]");
  fOutput->Add(fHVertexZ);
  fHVertexZ2 = new TH1F("hVertexZAfterCut","",100,-25,25);
  fHVertexZ2->SetXTitle("z [cm]");
  fOutput->Add(fHVertexZ2);
  fHCent = new TH1F("hCentBeforeCut","",102,-1,101);
  fHCent->SetXTitle(fCentVar.Data());
  fOutput->Add(fHCent);
  fHCentQual = new TH1F("hCentAfterCut","",102,-1,101);
  fHCentQual->SetXTitle(fCentVar.Data());
  fOutput->Add(fHCentQual);
  
  // histograms for clusters
  Int_t nsm = fGeom->GetEMCGeometry()->GetNumberOfSuperModules();
  if (!fTrainMode) {
    
    // default values: LHC11a
    Double_t fRunMin = 146700.5;
    Double_t fRunMax = 146900.5;
    Int_t iRunBins = 200;
    
    // run2, PbPb 5 TeV
    if(fDataPeriod == 202){
      fRunMin = 244800.5;
      fRunMax = 247000.5;
      iRunBins = 2200;
    }
    
    fHMeanClusterEnergy = new TProfile("hMeanClusterEnergy","mean cluster energy vs. run no",
                                       iRunBins,fRunMin,fRunMax," ");
    fOutput->Add(fHMeanClusterEnergy);
    
    fHMeanClusterNumber = new TProfile("hMeanClusterNumber","number of clusters vs. run no",
                                       iRunBins,fRunMin,fRunMax," ");
    fOutput->Add(fHMeanClusterNumber);
    
    fHCellIndexEnergy = new TH2F("hCellIndexEnergy","",18001,-0.5,18000.5,100,0,20);
    fHCellIndexEnergy->SetXTitle("Cell #");
    fHCellIndexEnergy->SetYTitle("E [GeV/c]");
    fOutput->Add(fHCellIndexEnergy);
    
    fHCellIndexEnergyAfterCuts= new TH2F("hCellIndexEnergyAfterCuts","",18001,-0.5,18000.5,100,0,20);
    fHCellIndexEnergyAfterCuts->SetXTitle("Cell #");
    fHCellIndexEnergyAfterCuts->SetYTitle("E [GeV/c]");
    fOutput->Add(fHCellIndexEnergyAfterCuts);
    
    fHClusters = new TH1F("hClusters","",7,0.5,7.5);
    fHClusters->GetXaxis()->SetBinLabel(1,"All");
    fHClusters->GetXaxis()->SetBinLabel(2,"Bad Cell");
    fHClusters->GetXaxis()->SetBinLabel(3,"Min E");
    fHClusters->GetXaxis()->SetBinLabel(4,"NCells");
    fHClusters->GetXaxis()->SetBinLabel(5,"E Ratio");
    fHClusters->GetXaxis()->SetBinLabel(6,"Eccentricity");
    fHClusters->GetXaxis()->SetBinLabel(7,"M02");
    fOutput->Add(fHClusters);
    
    
    fHClustAllEtaPhi = new TH2F("hClustAllEtaPhi","",100,-1,1,630,-3.15,3.15);
    fHClustAllEtaPhi->SetXTitle("#eta");
    fHClustAllEtaPhi->SetYTitle("#varphi");
    fOutput->Add(fHClustAllEtaPhi);
    
    fHClustNoEvt = new TH1F("hClustNoEvt","",2000,0,2000);
    fHClustNoEvt->SetXTitle("# Clusters");
    fOutput->Add(fHClustNoEvt);
    fHClustAccEvt = new TH1F("hClustAccEvt","",2000,0,2000);
    fHClustAccEvt->SetXTitle("# Clusters");
    fOutput->Add(fHClustAccEvt);
    fHClustEccentricity = new TH1F("hClustEccentricity","",200,-10,10);
    fHClustEccentricity->SetXTitle("#epsilon_{C}");
    fOutput->Add(fHClustEccentricity);
    fHClustEtaPhi = new TH2F("hClustEtaPhi","",160,-0.8,0.8,100*nsm,-3.15,3.15);
    fHClustEtaPhi->SetXTitle("#eta");
    fHClustEtaPhi->SetYTitle("#varphi");
    fOutput->Add(fHClustEtaPhi);
    
    fHClustEtaPhiAll = new TH2F("hClustEtaPhiAll","",160,-0.8,0.8,100*nsm,-3.15,3.15);
    fHClustEtaPhiAll->SetXTitle("#eta");
    fHClustEtaPhiAll->SetYTitle("#varphi");
    fOutput->Add(fHClustEtaPhiAll);
    
    fHClustEnergyPt = new TH2F("hClustEnergyPt","",250,0,50,250,0,50);
    fHClustEnergyPt->SetXTitle("E [GeV]");
    fHClustEnergyPt->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHClustEnergyPt);
    fHClustEnergyPtDCal = new TH2F("hClustEnergyPtDCal","",250,0,50,250,0,50);
    fHClustEnergyPtDCal->SetXTitle("E [GeV]");
    fHClustEnergyPtDCal->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHClustEnergyPtDCal);
    fHClustEnergySM = new TH2F("hClustEnergySM","",250,0,50,21,-0.5,20.5);
    fHClustEnergySM->SetXTitle("E [GeV]");
    fHClustEnergySM->SetYTitle("SM number");
    fOutput->Add(fHClustEnergySM);
    fHClustEnergySigma = new TH2F("hClustEnergySigma","",50,0,5,250,0,50);
    fHClustEnergySigma->SetXTitle("E_{C} * #sigma_{max} [GeV*cm]");
    fHClustEnergySigma->SetYTitle("E_{C} [GeV]");
    fOutput->Add(fHClustEnergySigma);
    fHClustEnergyTime = new TH2F("hClustEnergyTime","",100,0,20,500,-50e-8,200e-8);
    fHClustEnergyTime->SetXTitle("E [GeV]");
    fHClustEnergyTime->SetYTitle("t_{Cluster} [s]");
    fOutput->Add(fHClustEnergyTime);
    fHClustNCellEnergyRatio = new TH2F("hClustNCellEnergyRatio","",27,-0.5,26.5,101,-0.05,1.05);
    fHClustNCellEnergyRatio->SetXTitle("N_{cells}");
    fHClustNCellEnergyRatio->SetYTitle("E^{max}_{cell}/E_{clus}");
    fOutput->Add(fHClustNCellEnergyRatio);
    fHClustEnergyNCell = new TH2F("hClustEnergyNCell","",200,0,100,50,0,50);
    fHClustEnergyNCell->SetXTitle("E_{clus}");
    fHClustEnergyNCell->SetYTitle("N_{cells}");
    fOutput->Add(fHClustEnergyNCell);
  }
  
  
  // histograms for primary tracks
  fHPrimTrackPt = new TH1F("hPrimTrackPt",";p_{T} [GeV/c]",500,0,50);
  fOutput->Add(fHPrimTrackPt);
  fHPrimTrackEta = new TH1F("hPrimTrackEta",";#eta",40,-2,2);
  fOutput->Add(fHPrimTrackEta);
  fHPrimTrackPhi = new TH1F("hPrimTrackPhi",";#varPhi [rad]",63,0,6.3);
  fOutput->Add(fHPrimTrackPhi);
  
  // histograms for track matching
  if (fDoTrMatGeom) {
    fHMatchDr = new TH1F("hMatchDrDist",";dR [cm]",500,0,200);
    fOutput->Add(fHMatchDr);
    fHMatchDz = new TH1F("hMatchDzDist",";dZ [cm]",500,-100,100);
    fOutput->Add(fHMatchDz);
    fHMatchEp = new TH1F("hMatchEpDist",";E/p",100,0,10);
    fOutput->Add(fHMatchEp);
  }
  
  const Int_t nbins = 160;
  const Int_t ptmax = 40;
  
  // d_r of pairs
  fHdr = new TH1F("hdr","",2000,0,2);
  fHdr->SetXTitle("pair d_r");
  fOutput->Add(fHdr);
  
  // histograms for pion candidates
  if (!fTrainMode) {
    
    const Int_t massbins = 160;
    const Double_t massmax = 0.8;
    
    fHPionEtaPhi = new TH2F("hPionEtaPhi","",100,-0.8,0.8,50*nsm,-3.15,3.15);
    fHPionEtaPhi->SetXTitle("#eta_{#gamma#gamma}");
    fHPionEtaPhi->SetYTitle("#varphi_{#gamma#gamma}");
    fOutput->Add(fHPionEtaPhi);
    
    fHPionMggPt = new TH2F("hPionMggPt","",massbins,0,massmax,nbins,0,ptmax);
    fHPionMggPt->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionMggPt->SetYTitle("p_{T}^{#gamma#gamma} [GeV/c]");
    fOutput->Add(fHPionMggPt);
    
    fHPionMggAsym = new TH2F("hPionMggAsym","",400,0,0.8,50,0,1);
    fHPionMggAsym->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionMggAsym->SetYTitle("Z_{#gamma#gamma} [GeV]");
    fOutput->Add(fHPionMggAsym);
    
    fHPionPtAsym = new TH2F("hPionPtAsym","",nbins,0,ptmax,50,0,1);
    fHPionPtAsym->SetXTitle("p_{T} [GeV/c]");
    fHPionPtAsym->SetYTitle("Z_{#gamma#gamma} [GeV]");
    fOutput->Add(fHPionPtAsym);
    
    fHPionMggDgg = new TH2F("hPionMggDgg","",100,0,1,500,0,15);
    fHPionMggDgg->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionMggDgg->SetYTitle("opening angle [grad]");
    fOutput->Add(fHPionMggDgg);
    
    // pair invariant mass < asym
    fHPionInvMasses = new TH2F("hPionInvMass","hPionInvMass",massbins,0,massmax,nbins,0,ptmax);
    fHPionInvMasses->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionInvMasses->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHPionInvMasses);
    
    // pair invariant mass > asym
    fHPionInvMassesSym = new TH2F("hPionInvMassSym","hPionInvMassSym",massbins,0,massmax,nbins,0,ptmax);
    fHPionInvMassesSym->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionInvMassesSym->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHPionInvMassesSym);
    
    // pair invariant mass > asym
    fHPionInvMassesAsym = new TH2F("hPionInvMassAsym","hPionInvMassAsym",massbins,0,massmax,nbins,0,ptmax);
    fHPionInvMassesAsym->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionInvMassesAsym->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHPionInvMassesAsym);
    
    if(fMcMode){
      // primaries, secondaries, material, k0, ...
      // pairs from primary pi0s
      fHPriPionInvMasses = new TH2F("hPriPionInvMass","hPriPionInvMass",massbins,0,massmax,nbins,0,ptmax);
      fHPriPionInvMasses->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPriPionInvMasses->SetYTitle("p_{T} [GeV/c]");
      fOutput->Add(fHPriPionInvMasses);
      
      // pairs from secondary pi0s but K0
      fHSecPionInvMasses = new TH2F("hSecPionInvMass","hSecPionInvMass",massbins,0,massmax,nbins,0,ptmax);
      fHSecPionInvMasses->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHSecPionInvMasses->SetYTitle("p_{T} [GeV/c]");
      fOutput->Add(fHSecPionInvMasses);
      
      // pairs from K0 decay pions
      fHK0PionInvMasses = new TH2F("hK0PionInvMass","hK0PionInvMass",massbins,0,massmax,nbins,0,ptmax);
      fHK0PionInvMasses->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHK0PionInvMasses->SetYTitle("p_{T} [GeV/c]");
      fOutput->Add(fHK0PionInvMasses);
      
      // pairs from "material" pi0s
      fHMatPionInvMasses = new TH2F("hMatPionInvMass","hMatPionInvMass",massbins,0,massmax,nbins,0,ptmax);
      fHMatPionInvMasses->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHMatPionInvMasses->SetYTitle("p_{T} [GeV/c]");
      fOutput->Add(fHMatPionInvMasses);
      
      // added signals
      if(fAddedSignal){
        
        // no clean gammas
        // unweighted added signals
        fHPionInvMassesAdd1NoWgt= new TH2F("hPionInvMassAdd1NoWgt","hPionInvMassAdd1NoWgt",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesAdd1NoWgt->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesAdd1NoWgt->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesAdd1NoWgt->Sumw2();
        fOutput->Add(fHPionInvMassesAdd1NoWgt);
        
        // added pair invariant mass < asym
        fHPionInvMassesAdd1 = new TH2F("hPionInvMassAdd1","hPionInvMassAdd1",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesAdd1->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesAdd1->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesAdd1->Sumw2();
        fOutput->Add(fHPionInvMassesAdd1);
        
        // added pair invariant mass > asym
        fHPionInvMassesAdd1Sym = new TH2F("hPionInvMassAdd1Sym","hPionInvMassAdd1Sym",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesAdd1Sym->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesAdd1Sym->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesAdd1Sym->Sumw2();
        fOutput->Add(fHPionInvMassesAdd1Sym);
        
        // added pair invariant mass > asym
        fHPionInvMassesAdd1Asym = new TH2F("hPionInvMassAdd1Asym","hPionInvMassAdd1Asym",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesAdd1Asym->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesAdd1Asym->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesAdd1Asym->Sumw2();
        fOutput->Add(fHPionInvMassesAdd1Asym);
        
        // "unclean" added pair invariant mass < asym
        fHPionInvMassesAdd1Mult = new TH2F("hPionInvMassAdd1Mult","hPionInvMassAdd1Mult",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesAdd1Mult->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesAdd1Mult->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesAdd1Mult->Sumw2();
        fOutput->Add(fHPionInvMassesAdd1Mult);
        
        // "unclean" added pair invariant mass > asym
        fHPionInvMassesAdd1MultSym = new TH2F("hPionInvMassAdd1MultSym","hPionInvMassAdd1MultSym",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesAdd1MultSym->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesAdd1MultSym->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesAdd1MultSym->Sumw2();
        fOutput->Add(fHPionInvMassesAdd1MultSym);
        
        // "unclean" added pair invariant mass > asym
        fHPionInvMassesAdd1MultAsym = new TH2F("hPionInvMassAdd1MultAsym","hPionInvMassAdd1MultAsym",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesAdd1MultAsym->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesAdd1MultAsym->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesAdd1MultAsym->Sumw2();
        fOutput->Add(fHPionInvMassesAdd1MultAsym);
        
        // clean two-gammas
        // added pair invariant mass < asym
        fHPionInvMassesGamAdd1 = new TH2F("hPionInvMassGamAdd1","hPionInvMassGamAdd1",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesGamAdd1->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesGamAdd1->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesGamAdd1->Sumw2();
        fOutput->Add(fHPionInvMassesGamAdd1);
        
        // added pair invariant mass > asym
        fHPionInvMassesGamAdd1Sym = new TH2F("hPionInvMassGamAdd1Sym","hPionInvMassGamAdd1Sym",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesGamAdd1Sym->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesGamAdd1Sym->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesGamAdd1Sym->Sumw2();
        fOutput->Add(fHPionInvMassesGamAdd1Sym);
        
        // added pair invariant mass > asym
        fHPionInvMassesGamAdd1Asym = new TH2F("hPionInvMassGamAdd1Asym","hPionInvMassGamAdd1Asym",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesGamAdd1Asym->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesGamAdd1Asym->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesGamAdd1Asym->Sumw2();
        fOutput->Add(fHPionInvMassesGamAdd1Asym);
        
        // "unclean" added pair invariant mass < asym
        fHPionInvMassesGamAdd1Mult = new TH2F("hPionInvMassGamAdd1Mult","hPionInvMassGamAdd1Mult",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesGamAdd1Mult->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesGamAdd1Mult->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesGamAdd1Mult->Sumw2();
        fOutput->Add(fHPionInvMassesGamAdd1Mult);
        
        // "unclean" added pair invariant mass > asym
        fHPionInvMassesGamAdd1MultSym = new TH2F("hPionInvMassGamAdd1MultSym","hPionInvMassGamAdd1MultSym",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesGamAdd1MultSym->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesGamAdd1MultSym->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesGamAdd1MultSym->Sumw2();
        fOutput->Add(fHPionInvMassesGamAdd1MultSym);
        
        // "unclean" added pair invariant mass > asym
        fHPionInvMassesGamAdd1MultAsym = new TH2F("hPionInvMassGamAdd1MultAsym","hPionInvMassGamAdd1MultAsym",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesGamAdd1MultAsym->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesGamAdd1MultAsym->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesGamAdd1MultAsym->Sumw2();
        fOutput->Add(fHPionInvMassesGamAdd1MultAsym);
        
        // added pair invariant mass stream 2, no asym
        fHPionInvMassesAdd2 = new TH2F("hPionInvMassAdd2","hPionInvMassAdd2",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesAdd2->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesAdd2->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesAdd2->Sumw2();
        fOutput->Add(fHPionInvMassesAdd2);
      }
      
      if(fSimStudies){
        
        // electrons
        // main contributor: no converted electron (control histo)
        fHPionInvMassesConvElZero = new TH2F("hPionInvMassConvElZero","hPionInvMassConvElZero",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesConvElZero->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesConvElZero->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesConvElZero->Sumw2();
        fOutput->Add(fHPionInvMassesConvElZero);

        // main contributor: converted electron (one cluster)
        fHPionInvMassesConvElOne = new TH2F("hPionInvMassConvElOne","hPionInvMassConvElOne",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesConvElOne->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesConvElOne->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesConvElOne->Sumw2();
        fOutput->Add(fHPionInvMassesConvElOne);

        // main contributor: converted electron (two clusters)
        fHPionInvMassesConvElBoth = new TH2F("hPionInvMassConvElBoth","hPionInvMassConvElBoth",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesConvElBoth->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesConvElBoth->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesConvElBoth->Sumw2();
        fOutput->Add(fHPionInvMassesConvElBoth);

        // pions
        // main contributor: no charged pion  (control histo)
        fHPionInvMassesChargedPiZero = new TH2F("hPionInvMassChargedPiZero","hPionInvMassChargedPiZero",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesChargedPiZero->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesChargedPiZero->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesChargedPiZero->Sumw2();
        fOutput->Add(fHPionInvMassesChargedPiZero);

        // main contributor: charged pion (one cluster)
        fHPionInvMassesChargedPiOne = new TH2F("hPionInvMassChargedPiOne","hPionInvMassChargedPiOne",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesChargedPiOne->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesChargedPiOne->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesChargedPiOne->Sumw2();
        fOutput->Add(fHPionInvMassesChargedPiOne);

        // main contributor: charged pion (two clusters)
        fHPionInvMassesChargedPiBoth = new TH2F("hPionInvMassChargedPiBoth","hPionInvMassChargedPiBoth",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesChargedPiBoth->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesChargedPiBoth->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesChargedPiBoth->Sumw2();
        fOutput->Add(fHPionInvMassesChargedPiBoth);

        // main contributor: photon (both clusters)
        fHPionInvMassesGammaBoth = new TH2F("hPionInvMassGammaBoth","hPionInvMassGammaBoth",massbins,0,massmax,nbins,0,ptmax);
        fHPionInvMassesGammaBoth->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
        fHPionInvMassesGammaBoth->SetYTitle("p_{T} [GeV/c]");
        fHPionInvMassesGammaBoth->Sumw2();
        fOutput->Add(fHPionInvMassesGammaBoth);
        
      }
    }
    
    // mixed events
    fHPionInvMassesMix = new TH2F("hPionInvMassMix","hPionInvMassMix",massbins,0,massmax,nbins,0,ptmax);
    fHPionInvMassesMix->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionInvMassesMix->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHPionInvMassesMix);
    
    // it does not really make sense to mix added signals, but maybe it is interesting ...
    fHPionInvMassesMix1 = new TH2F("hPionInvMassMix1","hPionInvMassMix1",massbins,0,massmax,nbins,0,ptmax);
    fHPionInvMassesMix1->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionInvMassesMix1->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHPionInvMassesMix1);
    
    fHPionInvMassesMix2 = new TH2F("hPionInvMassMix2","hPionInvMassMix2",massbins,0,massmax,nbins,0,ptmax);
    fHPionInvMassesMix2->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionInvMassesMix2->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHPionInvMassesMix2);
    
    // more histograms, for DCal
    // still might want to add more of them!
    // DCal real events
    fHPionInvMassesDCal = new TH2F("hPionInvMassDCal","hPionInvMassDCal",massbins,0,massmax,nbins,0,ptmax);
    fHPionInvMassesDCal->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionInvMassesDCal->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHPionInvMassesDCal);
    
    // DCal mixed events
    fHPionInvMassesMixDCal = new TH2F("hPionInvMassMixDCal","hPionInvMassMixDCal",massbins,0,massmax,nbins,0,ptmax);
    fHPionInvMassesMixDCal->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionInvMassesMixDCal->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHPionInvMassesMixDCal);
    
    // let's see if there is anything if we combine EMCal with DCal
    // EMCal+DCal real events
    fHPionInvMassesEMCalDCal = new TH2F("hPionInvMassEMCalDCal","hPionInvMassEMCalDCal",200,0,10,40,0,20);
    fHPionInvMassesEMCalDCal->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionInvMassesEMCalDCal->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHPionInvMassesEMCalDCal);
    
    // EMCal+DCal mixed events
    fHPionInvMassesMixEMCalDCal = new TH2F("hPionInvMassMixEMCalDCal","hPionInvMassMixEMCalDCal",100,0,10,40,0,20);
    fHPionInvMassesMixEMCalDCal->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionInvMassesMixEMCalDCal->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHPionInvMassesMixEMCalDCal);
    
    // calibration stuff
    // here we fill (E_1 + E_2)/2 instead of pT!
    if(fCalibRun){
      // EMCal
      fHPionInvMassesEMCalCalib = new TH2F("hPionInvMassesEMCalCalib","hPionInvMassesEMCalCalib",massbins,0,massmax,nbins,0,ptmax);
      fHPionInvMassesEMCalCalib->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPionInvMassesEMCalCalib->SetYTitle("(E_{1} + E_{2}/2 [GeV]");
      fOutput->Add(fHPionInvMassesEMCalCalib);
      
      fHPionInvMassesMixEMCalCalib = new TH2F("hPionInvMassesMixEMCalCalib","hPionInvMassesMixEMCalCalib",massbins,0,massmax,nbins,0,ptmax);
      fHPionInvMassesMixEMCalCalib->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPionInvMassesMixEMCalCalib->SetYTitle("(E_{1} + E_{2}/2 [GeV]");
      fOutput->Add(fHPionInvMassesMixEMCalCalib);
      
      //DCal
      fHPionInvMassesDCalCalib = new TH2F("hPionInvMassesDCalCalib","hPionInvMassesDCalCalib",massbins,0,massmax,nbins,0,ptmax);
      fHPionInvMassesDCalCalib->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPionInvMassesDCalCalib->SetYTitle("(E_{1} + E_{2}/2 [GeV]");
      fOutput->Add(fHPionInvMassesDCalCalib);
      
      fHPionInvMassesMixDCalCalib = new TH2F("hPionInvMassesMixDCalCalib","hPionInvMassesMixEMCalCalib",massbins,0,massmax,nbins,0,ptmax);
      fHPionInvMassesMixDCalCalib->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPionInvMassesMixDCalCalib->SetYTitle("(E_{1} + E_{2}/2 [GeV]");
      fOutput->Add(fHPionInvMassesMixDCalCalib);
    }
    
    // distribution of particle weights
    fHWgt = new TH1F("hWgt","hWgt",100,0,10);
    fOutput->Add(fHWgt);
    
    if(fMcMode){
      fHPrimPionInvMasses = new TH2F("hPrimPionInvMass","hPrimPionInvMass",massbins,0,massmax,nbins,0,ptmax);
      fHPrimPionInvMasses->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPrimPionInvMasses->SetYTitle("p_{T} [GeV/c]");
      fOutput->Add(fHPrimPionInvMasses);
      
      fHPrimPionInvMassesAsym = new TH2F("hPrimPionInvMassAsym","hPrimPionInvMass Asym > 0.8",massbins,0,massmax,nbins,0,ptmax);
      fHPrimPionInvMassesAsym->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPrimPionInvMassesAsym->SetYTitle("p_{T} [GeV/c]");
      fOutput->Add(fHPrimPionInvMassesAsym);
      
    }
    
    //    // histogram for conversion point
    //    fHConversionPoint = new TH2F("hConversionPoint","conversion point in xy",1000,-500,500,1000,-500,500);
    //    fHConversionPoint->SetXTitle("x");
    //    fHConversionPoint->SetYTitle("y");
    //    fOutput->Add(fHConversionPoint);
  }
  
  TH1::SetDefaultSumw2(th1);
  TH2::SetDefaultSumw2(th2);
  PostData(1, fOutput);
  
  // MC histograms
  if(fMcMode){
    
    if(fAddedSignal){
      fHAddPionEtaPt = new TH2F("hAddPionEtaPt","",150,-2.5,2.5,100,0.,10.);
      fHAddPionEtaPt->SetXTitle("#eta_{#gamma#gamma}");
      fHAddPionEtaPt->SetYTitle("p_{T}");
      fOutput->Add(fHAddPionEtaPt);
      
      fHAddPionEtaPtWgt = new TH2F("hAddPionEtaPtWgt","",150,-2.5,2.5,100,0.,10.);
      fHAddPionEtaPtWgt->SetXTitle("#eta_{#gamma#gamma}");
      fHAddPionEtaPtWgt->SetYTitle("p_{T}");
      fOutput->Add(fHAddPionEtaPtWgt);
    }
    
    fHPyPionEtaPt = new TH2F("hPyPionEtaPt","",150,-2.5,2.5,100,0.,10.);
    fHPyPionEtaPt->SetXTitle("#eta_{#gamma#gamma}");
    fHPyPionEtaPt->SetYTitle("p_{T}");
    fOutput->Add(fHPyPionEtaPt);
    
    // pi0
    fHPionTruthPt = new TH1F("hPionTruthPt","pi0 truth pT from MC",nbins,0,ptmax);
    fHPionTruthPt->SetXTitle("p_{T}");
    fOutput->Add(fHPionTruthPt);
    
    fHPionTruthPtIn = new TH1F("hPionTruthPtIn","pi0 truth pT from MC within eta range",nbins,0,ptmax);
    fHPionTruthPtIn->SetXTitle("p_{T}");
    fOutput->Add(fHPionTruthPtIn);
    
    fHPionTruthPtAcc = new TH1F("hPionTruthPtAcc","pi0 truth pT from MC within acceptance",nbins,0,ptmax);
    fHPionTruthPtAcc->SetXTitle("p_{T}");
    fOutput->Add(fHPionTruthPtAcc);
    
    // eta
    fHEtaTruthPt = new TH1F("hEtaTruthPt","eta truth pT from MC",nbins,0,ptmax);
    fHEtaTruthPt->SetXTitle("p_{T}");
    fOutput->Add(fHEtaTruthPt);
    
    fHEtaTruthPtIn = new TH1F("hEtaTruthPtIn","eta truth pT from MC within eta range",nbins,0,ptmax);
    fHEtaTruthPtIn->SetXTitle("p_{T}");
    fOutput->Add(fHEtaTruthPtIn);
    
    fHEtaTruthPtAcc = new TH1F("hEtaTruthPtAcc","eta truth pT from MC within acceptance",nbins,0,ptmax);
    fHEtaTruthPtAcc->SetXTitle("p_{T}");
    fOutput->Add(fHEtaTruthPtAcc);
    
    // direct photons
    fHGamTruthPt = new TH1F("hGamTruthPt","gamma truth pT from MC",nbins,0,ptmax);
    fHGamTruthPt->SetXTitle("p_{T}");
    fOutput->Add(fHGamTruthPt);
    
    fHGamTruthPtIn = new TH1F("hGamTruthPtIn","gamma truth pT from MC within eta range",nbins,0,ptmax);
    fHGamTruthPtIn->SetXTitle("p_{T}");
    fOutput->Add(fHGamTruthPtIn);
    
    fHGamTruthPtAcc = new TH1F("hGamTruthPtAcc","gamma truth pT from MC within acceptance",nbins,0,ptmax);
    fHGamTruthPtAcc->SetXTitle("p_{T}");
    fOutput->Add(fHGamTruthPtAcc);
    
    // added signals
    if(fAddedSignal){
      // pi0
      fHPionTruthPtAdd = new TH1F("hPionTruthPtAdd","added pi0 truth pT from MC",nbins,0,ptmax);
      fHPionTruthPtAdd->SetXTitle("p_{T}");
      fHPionTruthPtAdd->Sumw2();
      fOutput->Add(fHPionTruthPtAdd);
      
      fHPionTruthPtInAdd = new TH1F("hPionTruthPtInAdd","added pi0 truth pT from MC within eta range",nbins,0,ptmax);
      fHPionTruthPtInAdd->SetXTitle("p_{T}");
      fHPionTruthPtInAdd->Sumw2();
      fOutput->Add(fHPionTruthPtInAdd);
      
      fHPionTruthPtAccAdd = new TH1F("hPionTruthPtAccAdd","added pi0 truth pT from MC within acceptance",nbins,0,ptmax);
      fHPionTruthPtAccAdd->SetXTitle("p_{T}");
      fHPionTruthPtAccAdd->Sumw2();
      fOutput->Add(fHPionTruthPtAccAdd);
      
      // eta
      fHEtaTruthPtAdd = new TH1F("hEtaTruthPtAdd","added eta truth pT from MC",nbins,0,ptmax);
      fHEtaTruthPtAdd->SetXTitle("p_{T}");
      fHEtaTruthPtAdd->Sumw2();
      fOutput->Add(fHEtaTruthPtAdd);
      
      fHEtaTruthPtInAdd = new TH1F("hEtaTruthPtInAdd","added eta truth pT from MC within eta range",nbins,0,ptmax);
      fHEtaTruthPtInAdd->SetXTitle("p_{T}");
      fHEtaTruthPtInAdd->Sumw2();
      fOutput->Add(fHEtaTruthPtInAdd);
      
      fHEtaTruthPtAccAdd = new TH1F("hEtaTruthPtAccAdd","added eta truth pT from MC within acceptance",nbins,0,ptmax);
      fHEtaTruthPtAccAdd->SetXTitle("p_{T}");
      fHEtaTruthPtAccAdd->Sumw2();
      fOutput->Add(fHEtaTruthPtAccAdd);
    }
    // particle information
    fHMCpartfrac = new TH2F("hMCpartfrac","fraction of most energy MC particle vs. fraction of this particle to all",100,0,2,100,0.8,1.2);
    fHMCpartfrac->SetXTitle("E_{MC}/E_{Clu}");
    fHMCpartfrac->SetYTitle("E_{MC}^{highest}/E_{MC}^{all}");
    fOutput->Add(fHMCpartfrac);
    
    fHECluEMC = new TH2F("hECluEMC","energy of most contributing MC track vs. energy of cluster",200,0,10,200,0,10);
    fHECluEMC->SetXTitle("E_{MC}");
    fHECluEMC->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMC);
    
    fHECluEMCAddPi0 = new TH2F("hECluEMCAddPi0","energy of most contributing added pi0 MC track vs. energy of cluster",200,0,10,200,0,10);
    fHECluEMCAddPi0->SetXTitle("E_{MC}");
    fHECluEMCAddPi0->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMCAddPi0);
    
    fHECluEMCAddEta = new TH2F("hECluEMCAddEta","energy of most contributing added eta MC track vs. energy of cluster",200,0,10,200,0,10);
    fHECluEMCAddEta->SetXTitle("E_{MC}");
    fHECluEMCAddEta->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMCAddEta);
    
    fHECluEMCnofullAdd = new TH2F("hECluEMCnofullAdd","energy of most contributing added MC track vs. energy of cluster (more than one MC particle in cluster)",200,0,10,200,0,10);
    fHECluEMCnofullAdd->SetXTitle("E_{MC}");
    fHECluEMCnofullAdd->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMCnofullAdd);
    
    fHECluEMCnofull = new TH2F("hECluEMCnofull","energy of most contributing MC track vs. energy of cluster (more than one MC particle in cluster)",200,0,10,200,0,10);
    fHECluEMCnofull->SetXTitle("E_{MC}");
    fHECluEMCnofull->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMCnofull);
    
    fHECluEMCelectron = new TH2F("hECluEMCelectron","energy of most contributing MC track vs. energy of cluster for electrons",200,0,10,200,0,10);
    fHECluEMCelectron->SetXTitle("E_{MC}");
    fHECluEMCelectron->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMCelectron);
    
    fHECluEMCpion = new TH2F("hECluEMCpion","energy of most contributing MC track vs. energy of cluster for pions",200,0,10,200,0,10);
    fHECluEMCpion->SetXTitle("E_{MC}");
    fHECluEMCpion->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMCpion);
    
    fHECluEMCkaon = new TH2F("hECluEMCkaon","energy of most contributing MC track vs. energy of cluster for kaons",200,0,10,200,0,10);
    fHECluEMCkaon->SetXTitle("E_{MC}");
    fHECluEMCkaon->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMCkaon);
    
    fHECluEMCother = new TH2F("hECluEMCother","energy of most contributing MC track vs. energy of cluster for others",200,0,10,200,0,10);
    fHECluEMCother->SetXTitle("E_{MC}");
    fHECluEMCother->SetYTitle("E_{Clu}");
    fOutput->Add(fHECluEMCother);
    
    fHECluEMCpi0single = new TH2F("hECluEMCpi0single","E pi0 in MC vs. E of merged cluster",200,5,15,200,5,15);
    fHECluEMCpi0single->SetXTitle("p_{T}^{MC}");
    fHECluEMCpi0single->SetYTitle("p_{T}^{Clu}");
    fOutput->Add(fHECluEMCpi0single);
    
    fHNMothers = new TH2F("fHNMothers","fHNMothers",21,-0.5,20.5,60,0,30);
    fHNMothers->SetXTitle("mother iterations");
    fHNMothers->SetYTitle("Energy");
    fOutput->Add(fHNMothers);
    
  }
  
  //	if(fRotateMixed){
  //		// more histograms
  //		fHMixRotation = new TH1F("hMixRotation","rotation angle of mixed events",100,-6.28,6.28);
  //		fHMixRotation->SetXTitle("phi");
  //		fOutput->Add(fHMixRotation);
  //	}
  //  fHCorrection = new TH1F("hCorrection","correction factor for single photon",100,0,2);
  //  fHCorrection->SetXTitle("correction factor");
  //  fOutput->Add(fHCorrection);
  
  //  fHPionSm = new TH2F("hPionSm","mass shift due to energy scale",200,0,0.5,200,0,0.5);
  //  fHPionSm->SetXTitle("original m_inv");
  //  fHPionSm->SetYTitle("changed m_inv");
  //  fOutput->Add(fHPionSm);
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::UserExec(Option_t *)
{
  // Called for each event.
  
  if (!InputEvent()){
    return;
    AliError("No Input Event!");
  }
  
  // get managers and handlers
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (am->GetInputEventHandler());
  
  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (am->GetInputEventHandler());
  if (!aodH && !esdH)
  {
    Printf("ERROR: Could not get AODInputHandler");
  }
  
  if(esdH)
    fEsdEv = (AliESDEvent*)esdH->GetEvent();
  
  else if(aodH)
    fAodEv = aodH->GetEvent();
  
  else{
    AliFatal("Neither ESD nor AOD event found");
    return;
  }
  
  // set all counters to nought
  ipymin=0;
  ipymax=0;
  ipi0min=0;
  ipi0max=0;
  ietamin=0;
  ietamax=0;
  
  if(fMcMode && fAddedSignal){
    // monte carlo headers from Evi
    AliMCEvent *mcEvent = MCEvent();
    if (!mcEvent){
      cout << "no MC event" << endl;
      return;
    }
    eventHeader = mcEvent->GenEventHeader();
    
    AliGenCocktailEventHeader *cocktail = dynamic_cast<AliGenCocktailEventHeader *>(eventHeader);
    
    if(!cocktail) return ;
    TList *genHeaders = cocktail->GetHeaders();
    
    Int_t nGenerators = genHeaders->GetEntries();
    Int_t pythiaLastP=0;
    
    // it seems gen1 and gen2 are pi0 and eta added signals, respectively, when looking at lhc12i3
    for(Int_t igen = 0; igen < nGenerators; igen++)
    {
      AliGenEventHeader* eventHeader2 = (AliGenEventHeader*)genHeaders->At(igen) ;
      TString name = eventHeader2->GetName();
      
      if (name.Contains("Pythia",TString::kIgnoreCase))
      {
        pythiaLastP=eventHeader2->NProduced();
        ipymax=eventHeader2->NProduced()-1;
        //Printf("pythia partcles = %d", pythiaLastP);
        pythiaHeader = (AliGenEventHeader*)genHeaders->At(igen);
      }
      //fNMCProducedMin = pythiaLastP;
      //fNMCProducedMax= fNMCProducedMin+eventHeader2->NProduced();
      
      //Printf("Generator %d: Class Name %s, Name %s, title %s \n, nProduced %d",igen, eventHeader2->ClassName(), name.Data(), eventHeader2->GetTitle(),fNMCProducedMax-fNMCProducedMin);
      
      //    if (name.Contains("Pi0Flat",TString::kIgnoreCase)) {IsPi0Flat=1; break;}
      if (igen==1) {
        addedPi0Header = (AliGenEventHeader*)genHeaders->At(igen);
        ipi0min=ipymax+1;
        ipi0max=ipi0min+addedPi0Header->NProduced()-1;
      }
      if (igen==2) {
        addedEtaHeader = (AliGenEventHeader*)genHeaders->At(igen);
        ietamin=ipi0max+1;
        ietamax=ietamin+addedEtaHeader->NProduced()-1;
        break;
      }
    }
    
    //Printf("min = %d, max = %d",fNMCProducedMin, fNMCProducedMax);
  }
  
  
  UInt_t offtrigger = 0;
  if (fEsdEv) {
    am->LoadBranch("AliESDRun.");
    am->LoadBranch("AliESDHeader.");
    UInt_t mask1 = fEsdEv->GetESDRun()->GetDetectorsInDAQ();
    UInt_t mask2 = fEsdEv->GetESDRun()->GetDetectorsInReco();
    Bool_t desc1 = (mask1 >> 18) & 0x1;
    Bool_t desc2 = (mask2 >> 18) & 0x1;
    if (desc1==0 || desc2==0) { //AliDAQ::OfflineModuleName(18)=="EMCAL"
      AliError(Form("EMCAL not in DAQ/RECO: %u (%u)/%u (%u)",
                    mask1, fEsdEv->GetESDRun()->GetDetectorsInReco(),
                    mask2, fEsdEv->GetESDRun()->GetDetectorsInDAQ()));
      return;
    }
  }
  
  if(fAodEv){
    am->LoadBranch("header");
    offtrigger =  ((AliVAODHeader*)fAodEv->GetHeader())->GetOfflineTrigger();
  }
  
  if (!fMcMode && (offtrigger & AliVEvent::kFastOnly)) {
    AliWarning(Form("EMCAL not in fast only partition"));
    return;
  }
  
  // get EMCal geometry if necessary
  if (fDoTrMatGeom && !AliGeomManager::GetGeometry()) {
    AliWarning("Accessing geometry from OCDB, this is not very efficient!");
    AliCDBManager *cdb = AliCDBManager::Instance();
    if (!cdb->IsDefaultStorageSet())
      cdb->SetDefaultStorage("raw://");
    Int_t runno = InputEvent()->GetRunNumber();
    if (runno != cdb->GetRun())
      cdb->SetRun(runno);
    AliGeomManager::LoadGeometry();
  }
  
  // set misalignment matrices (stored in first event)
  if (!AliGeomManager::GetGeometry()&&!fIsGeoMatsSet) {
    Int_t nsm = fGeom->GetEMCGeometry()->GetNumberOfSuperModules();
    for (Int_t i=0; i<nsm; ++i) {
      const TGeoHMatrix *geom = 0;
      if (fEsdEv)
        geom = fEsdEv->GetESDRun()->GetEMCALMatrix(i);
      else{
        AliAODHeader * aodheader = dynamic_cast<AliAODHeader*>(fAodEv->GetHeader());
        if(!aodheader) AliFatal("Not a standard AOD");
        geom = aodheader->GetEMCALMatrix(i);
      }
      if (!geom)
        continue;
      geom->Print();
      fGeom->SetMisalMatrix(geom,i);
    }
    fIsGeoMatsSet = kTRUE;
  }
  
  Int_t cut = 1;
  fHCuts->Fill(cut++);
  
  TString trgclasses;
  AliESDHeader *h = dynamic_cast<AliESDHeader*>(InputEvent()->GetHeader());
  if (h) {
    trgclasses = fEsdEv->GetFiredTriggerClasses();
  } else {
    AliAODHeader *h2 = dynamic_cast<AliAODHeader*>(InputEvent()->GetHeader());
    if (h2)
      trgclasses = h2->GetFiredTriggerClasses();
  }
  
  //
  //  if(trgclasses.Contains("EMC") && trgclasses.Contains("INT")){
  //    cout << "####################################################################### EMC AND INT SET #######################################################################" << endl;
  //    cout << "trigger mask EMC: "<<  std::bitset<32>(fEsdEv->GetTriggerMask()) << " = " << fEsdEv->GetTriggerMask()<< endl;
  //    cout << trgclasses << endl;
  //  }
  //
  //  if(1){
  //  if(trgclasses.Contains("EMC")){
  //    cout << trgclasses << endl;
  //    cout << "trigger mask EMC: "<<  std::bitset<32>(fEsdEv->GetTriggerMask()) << " = " << fEsdEv->GetTriggerMask()<< endl;
  //    cout << "trigger bit 0 is " <<  std::bitset<32>(fEsdEv->GetTriggerMask()).test(0) << endl;
  //  }
  //  }
  //  if(1){
  //  if(trgclasses.Contains("INT")){
  //    cout << trgclasses << endl;
  //    cout << "trigger mask INT: "<<  std::bitset<32>(fEsdEv->GetTriggerMask()) << " = " << fEsdEv->GetTriggerMask()<< endl;
  //    cout << "trigger bit 0 is " <<  std::bitset<32>(fEsdEv->GetTriggerMask()).test(0) << endl;
  //  }
  //  }
  
  //  if (fDoPSel && offtrigger==0)
  //    return;
  
  // cut on certain events
  if(fEsdEv){
    const AliESDVertex* vtxESD = fEsdEv->GetPrimaryVertexTracks();
    if( vtxESD->GetNContributors() < 1 )
    {
      AliWarning("No vertex contributor");
    }
    if(!(fEsdEv->GetPrimaryVertex()->GetStatus())){
      AliWarning("No vertex");
      return;
    }
  }
  
  // in lhc11a, cut on even more events
  Int_t runnumber = InputEvent()->GetRunNumber();
  if ((runnumber>=144871) && (runnumber<=146860)) {
    
    AliVCaloCells *cells   = InputEvent()->GetEMCALCells();
    const Short_t nCells   = cells->GetNumberOfCells();
    
    if (InputEvent()->IsA()==AliESDEvent::Class()) AliAnalysisManager::GetAnalysisManager()->LoadBranch("EMCALCells.");
    
    AliInputEventHandler *fInputHandler=(AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!fInputHandler) return;
    
    // count cells above threshold
    Int_t nCellCount[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
    for(Int_t iCell=0; iCell<nCells; ++iCell) {
      Short_t cellId = cells->GetCellNumber(iCell);
      Double_t cellE = cells->GetCellAmplitude(cellId);
      Int_t sm       = cellId / (24*48);
      if (cellE>0.1) ++nCellCount[sm];
    }
    
    Bool_t fIsLedEvent = kFALSE;
    if (nCellCount[4] > 100) {
      fIsLedEvent = kTRUE;
    } else {
      if ((runnumber>=146858) && (runnumber<=146860)) {
        if ((fInputHandler->IsEventSelected() & AliVEvent::kMB) && (nCellCount[3]>=21))
          fIsLedEvent = kTRUE;
        else if ((fInputHandler->IsEventSelected() & AliVEvent::kEMC1) && (nCellCount[3]>=35))
          fIsLedEvent = kTRUE;
      }
    }
    if (fIsLedEvent) {
      AliWarning("LED Event");
      return;
    }
  }
  
  fHCuts->Fill(cut++);
  
  
  //  const AliCentrality *centP = InputEvent()->GetCentrality();
  //  Double_t cent = centP->GetCentralityPercentileUnchecked(fCentVar);
  //  fHCent->Fill(cent);
  //  if (cent<fCentFrom||cent>fCentTo)
  //    return;
  
  // run 2 centrality
  Double_t cent = 1;
  if(fDataPeriod == 202){
    cent = 300;
    AliMultSelection *MultSelection = 0x0;
    MultSelection = (AliMultSelection * ) InputEvent()->FindListObject("MultSelection");
    if( !MultSelection) {
      //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
      AliWarning("AliMultSelection object not found!");
    }else{
      cent = MultSelection->GetMultiplicityPercentile(fCentVar.Data());
      //    cent = MultSelection->GetMultiplicityPercentile("CL1");
    }
  }
  else{
    const AliCentrality *centP = InputEvent()->GetCentrality();
    cent = centP->GetCentralityPercentileUnchecked(fCentVar.Data());
  }
  fHCent->Fill(cent);
  
  if (cent<fCentFrom||cent>=fCentTo)
    return;
  
  fHCuts->Fill(cut++);
  
  // cut on certain events  - LHC11a and LHC13g so far
  if( ((runnumber>=144871) && (runnumber<=146860)) || ((runnumber>=197470) && (runnumber<=197692))){
    
    if(fEsdEv){
      TString trigClasses = fEsdEv->GetFiredTriggerClasses();
      // remove "fast cluster events":
      if (trigClasses.Contains("FAST")  && !trigClasses.Contains("ALL")){
        AliWarning("fast cluster event");
        return;
      }
    }
    else if(fAodEv){
      TString trigClasses = fAodEv->GetFiredTriggerClasses();
      // remove "fast cluster events":
      if (trigClasses.Contains("FAST")  && !trigClasses.Contains("ALL")){
        AliWarning("fast cluster event");
        return;
      }
    }
  }
  fHCuts->Fill(cut++);
  
  if (fUseQualFlag) {
    //    if (centP->GetQuality()>0)
    //      return;
  }
  
  fHCentQual->Fill(cent);
  fHCuts->Fill(cut++);
  
  if (fEsdEv) {
    am->LoadBranch("PrimaryVertex.");
    am->LoadBranch("SPDVertex.");
    am->LoadBranch("TPCVertex.");
  } else {
    fAodEv = dynamic_cast<AliAODEvent*>(InputEvent());
    am->LoadBranch("vertices");
    if (!fAodEv) return;
  }
  
  const AliVVertex *vertex = InputEvent()->GetPrimaryVertex();
  if (!vertex)
    return;
  
  fHVertexZ->Fill(vertex->GetZ());
  
  if(vertex->GetZ()<fVtxZMin||vertex->GetZ()>fVtxZMax)
    return;
  
  fHCuts->Fill(cut++);
  fHVertexZ2->Fill(vertex->GetZ());
  
  // count number of accepted events
  ++fNEvs;
  
  fRecPoints   = 0; // will be set if fClusName is given and AliAnalysisTaskEMCALClusterizeFast is used
  fDigits      = 0; // will be set if fClusName is given and AliAnalysisTaskEMCALClusterizeFast is used
  fEsdClusters = 0; // will be set if ESD input used and if fRecPoints are not set or if clusters are attached
  fEsdCells    = 0; // will be set if ESD input used
  fAodClusters = 0; // will be set if AOD input used and if fRecPoints are not set or if clusters are attached
  //             or if fClusName is given and AliAnalysisTaskEMCALClusterizeFast in AOD output mode
  fAodCells    = 0; // will be set if AOD input used
  
  // deal with special output from AliAnalysisTaskEMCALClusterizeFast first
  Bool_t overwrite    = 0;
  Bool_t clusattached = 0;
  Bool_t recalibrated = 0;
  if (0 && !fClusName.IsNull()) {
    AliAnalysisTaskEMCALClusterizeFast *cltask = 0;
    TObjArray *ts = am->GetTasks();
    cltask = dynamic_cast<AliAnalysisTaskEMCALClusterizeFast*>(ts->FindObject(fClusName));
    if (cltask && cltask->GetClusters()) {
      fRecPoints   = cltask->GetClusters();
      fDigits      = cltask->GetDigits();
      clusattached = cltask->GetAttachClusters();
      overwrite    = cltask->GetOverwrite();
      if (cltask->GetCalibData()!=0)
        recalibrated = kTRUE;
    }
  }
  if (1 && !fClusName.IsNull()) {
    TList *l = 0;
    if (AODEvent())
      l = AODEvent()->GetList();
    else if (fAodEv)
      l = fAodEv->GetList();
    if (l) {
      fAodClusters = dynamic_cast<TClonesArray*>(l->FindObject(fClusName));
    }
  }
  
  if (fEsdEv) { // ESD input mode
    if (1 && (!fRecPoints||clusattached)) {
      if (!clusattached && !overwrite)
        am->LoadBranch("CaloClusters");
      TList *l = fEsdEv->GetList();
      if (clusattached) {
        fEsdClusters = dynamic_cast<TClonesArray*>(l->FindObject(fClusName));
      } else {
        fEsdClusters = dynamic_cast<TClonesArray*>(l->FindObject("CaloClusters"));
      }
    }
    if (1) {
      if (!recalibrated)
        am->LoadBranch("EMCALCells.");
      fEsdCells = fEsdEv->GetEMCALCells();
    }
  } else if (fAodEv) { // AOD input mode
    if (1 && (!fAodClusters || clusattached)) {
      if (!clusattached)
        am->LoadBranch("caloClusters");
      TList *l = fAodEv->GetList();
      if (l) {
        fAodClusters = dynamic_cast<TClonesArray*>(l->FindObject("caloClusters"));
      }
    }
    if (1) {
      if (!recalibrated)
        am->LoadBranch("emcalCells");
      fAodCells = fAodEv->GetEMCALCells();
    }
  } else {
    AliFatal("Impossible to not have either pointer to ESD or AOD event");
  }
  
  
  if (1) {
    AliDebug(2,Form("fRecPoints   set: %p", fRecPoints));
    AliDebug(2,Form("fDigits      set: %p", fDigits));
    AliDebug(2,Form("fEsdClusters set: %p", fEsdClusters));
    AliDebug(2,Form("fEsdCells    set: %p", fEsdCells));
    AliDebug(2,Form("fAodClusters set: %p", fAodClusters));
    AliDebug(2,Form("fAodCells    set: %p", fAodCells));
  }
  
  if (fDoAfterburner)
    ClusterAfterburner();
  
  if (fMcMode)
    CalcMcInfo();
  
  // mixed events
  
  // some info
  Int_t vtxClass = 1;
  Double_t vtxcuts = (fVtxZMax-fVtxZMin)/nZClass;
  if(vertex->GetZ()>fVtxZMin && vertex->GetZ()<=fVtxZMin+vtxcuts) vtxClass=0;
  else if(vertex->GetZ()>fVtxZMin+vtxcuts && vertex->GetZ()<=fVtxZMin+2*vtxcuts) vtxClass=1;
  else if(vertex->GetZ()>fVtxZMin+2*vtxcuts && vertex->GetZ()<=fVtxZMin+3*vtxcuts) vtxClass=2;
  
  //    vtxClass = 0;
  
  Int_t MulClass = 0;
  
  GetMulClass(MulClass);
  
  Float_t phitrig = 0;
  Float_t thetatrig = 0;
  Double_t pt_max = 0;
  Int_t ptClass = 0;
  if (!fTrainMode) {
    
    
    pt_max = FillClusHists(phitrig, thetatrig);
    
    if(pt_max < -100)
      return;
    
    if(pt_max > 1 && pt_max<=3) ptClass = 1;
    else if(pt_max > 3 && pt_max<=6) ptClass = 2;
    else ptClass = 3;
    
    ptClass = 0;
    
    // these are where the pi0 stuff is filled
    FillPionHists();
    FillMixHists(MulClass,vtxClass,ptClass,phitrig,thetatrig);
    
    //FillOtherHists();
  }
  FillMcHists();
  if(fDoNtuple)
    FillNtuple();
  
  if(MulClass < nMulClass && vtxClass < nZClass && ptClass < nPtClass){
    AddMixEvent(MulClass, vtxClass, ptClass, iEvt[MulClass][vtxClass][ptClass], phitrig, thetatrig);
  }
  else{
    AliWarning("event outside defined classes");
  }
  
  if (fTrainMode) {
    fSelTracks->Clear();
    fSelPrimTracks->Clear();
    if (fMcParts)
      fMcParts->Clear();
    if (fTriggers)
      fTriggers->Clear();
    if (fClusters)
      fClusters->Clear();
  }
  
  thisEvent.Reset();
  PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::Terminate(Option_t *)
{
  // Terminate called at the end of analysis.
  
  if (fNtuple) {
    TFile *f = OpenFile(1);
    TDirectory::TContext context(f);
    if (f)
      fNtuple->Write();
  }
  
  AliInfo(Form("%s: Accepted %lld events          ", GetName(), fNEvs));
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::ClusterAfterburner()
{
  // Run custer reconstruction afterburner.
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  
  if (!cells)
    return;
  
  Int_t ncells = cells->GetNumberOfCells();
  if (ncells<=0)
    return;
  
  Double_t cellMeanE = 0, cellSigE = 0;
  for (Int_t i = 0; i<ncells; ++i) {
    Double_t cellE = cells->GetAmplitude(i);
    cellMeanE += cellE;
    cellSigE += cellE*cellE;
  }
  cellMeanE /= ncells;
  cellSigE /= ncells;
  cellSigE -= cellMeanE*cellMeanE;
  if (cellSigE<0)
    cellSigE = 0;
  cellSigE = TMath::Sqrt(cellSigE / ncells);
  
  Double_t subE = cellMeanE - 7*cellSigE;
  if (subE<0)
    return;
  
  for (Short_t i = 0; i<ncells; ++i) {
    Short_t id=-1;
    Int_t mclabel = -1;
    Double_t amp=0,time=0, efrac = 0;
    if (!cells->GetCell(i, id, amp, time, mclabel, efrac))
      continue;
    amp -= cellMeanE;
    if (amp<0.001)
      amp = 0;
    cells->SetCell(i, id, amp, time);
  }
  
  TObjArray *clusters = fEsdClusters;
  if (!clusters)
    clusters = fAodClusters;
  if (!clusters)
    return;
  
  Int_t nclus = clusters->GetEntries();
  for (Int_t i = 0; i<nclus; ++i) {
    AliVCluster *clus = static_cast<AliVCluster*>(clusters->At(i));
    if (!clus->IsEMCAL())
      continue;
    Int_t nc = clus->GetNCells();
    Double_t clusE = 0;
    UShort_t ids[100] = {0};
    Double_t fra[100] = {0};
    for (Int_t j = 0; j<nc; ++j) {
      Short_t id = TMath::Abs(clus->GetCellAbsId(j));
      Double_t cen = cells->GetCellAmplitude(id);
      clusE += cen;
      if (cen>0) {
        ids[nc] = id;
        ++nc;
      }
    }
    if (clusE<=0) {
      clusters->RemoveAt(i);
      continue;
    }
    
    for (Int_t j = 0; j<nc; ++j) {
      Short_t id = ids[j];
      Double_t cen = cells->GetCellAmplitude(id);
      fra[j] = cen/clusE;
    }
    clus->SetE(clusE);
    AliAODCaloCluster *aodclus = dynamic_cast<AliAODCaloCluster*>(clus);
    if (aodclus) {
      aodclus->Clear("");
      aodclus->SetNCells(nc);
      aodclus->SetCellsAmplitudeFraction(fra);
      aodclus->SetCellsAbsId(ids);
    }
  }
  clusters->Compress();
}


//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPi0Gamma::FillClusHists(Float_t& max_phi, Float_t& max_theta)
{
  
  Int_t runnumber = InputEvent()->GetRunNumber();
  
  //cout << "filling cluster histograms ";
  
  bool bprint = 0;
  
  max_phi = 0;
  max_theta = 0;

  // Fill histograms related to cluster properties.
  
  // get objects
  TObjArray *clusters = fEsdClusters;
  if (!clusters)
    clusters = fAodClusters;
  if (!clusters){
    cout << "no clusters node in event!" << endl;
    return 0;
  }
  
  // get clusters
  Int_t nclus = clusters->GetEntries();
  
  //       cout << nclus << " clusters in event ";
  
  // get vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  
  // fill cluster number amd store event properties
  fHClustNoEvt->Fill(nclus);
  thisEvent.SetGlobalInfo(0,0,0);
  
  Int_t nclusemcaldcal = 0;
  // count clusters on emcal/dcal
  for(Int_t i = 0; i<nclus; ++i) {
    AliVCluster *clus = static_cast<AliVCluster*>(clusters->At(i));
    if (!clus){
      continue;
    }
    // emcal cluster?
    if(!(clus->IsEMCAL())){
      continue;
    }
    nclusemcaldcal++;
  }
  
  // set a limit due to memory
  if(nclusemcaldcal > 1000){
    AliError("Attention! More than 1000 EMCal/DCal clusters in event!");
    return -999;
  }
  
  int nclusters = 0;
  // main cluster loop
  for(Int_t i = 0; i<nclus; ++i) {
    
    Bool_t bdcal = 0;
    // get cluster
    AliVCluster *clus = static_cast<AliVCluster*>(clusters->At(i));
    if (!clus){
      continue;
    }
    
    // fill all clusters
    TLorentzVector clusterVec;
    clus->GetMomentum(clusterVec,vertex);
    fHClustAllEtaPhi->Fill(clusterVec.Eta(),clusterVec.Phi());
    
    // emcal cluster?
    if(!(clus->IsEMCAL())){
      continue;
    }
    
    // DCal (pure geometry)
    if ( (clusterVec.Phi() < 1.2 && clusterVec.Phi() > -2.8) ){
      bdcal = 1;
    }
    
    // fill QA histograms for cells
    FillCellQAHists(clus,bdcal,0);
    
    //if(bdcal) continue;
    
    Double_t maxAxis    = 1; //clus->GetTOF(); //sigma
    Double_t clusterEcc = 1; //clus->Chi2();   //eccentricity
    fHClustEccentricity->Fill(clusterEcc);
 
    // fill clusters in the beginning
    
    fHClustEtaPhiAll->Fill(clusterVec.Eta(),clusterVec.Phi());
    
    // cluster QA, fill a histogram
    Int_t cluster = 1;
    fHClusters->Fill(cluster++);

    // look if cluster is on bad cell
    if(fApplyBadMapManually) {
      UShort_t* CellsID = clus->GetCellsAbsId();
      
      // Find Max Contributing Cell
      AliVCaloCells *vcells = fEsdCells;
      if (!vcells)
        vcells = fAodCells;
      
      Float_t cellEnergy = 0.0;
      Float_t maxEnergy = 0.0;
      Int_t maxID = -1;
      
      for(Int_t kk=0; kk < clus->GetNCells(); kk++) {
        cellEnergy = vcells->GetCellAmplitude(CellsID[kk]);
        if(cellEnergy > maxEnergy) {
          maxEnergy = cellEnergy;
          maxID = CellsID[kk];
        }
      }
      if(fBadMap->GetBinContent(fBadMap->FindBin(maxID))>0)
        continue;
    }
    fHClusters->Fill(cluster++);
    
    // apply cluster cuts now
    if (clus->E()<fMinE)
      continue;
    
    fHClusters->Fill(cluster++);
    if (clus->GetNCells()<fNminCells)
      continue;
    
    fHClusters->Fill(cluster++);
    if (GetMaxCellEnergy(clus)/clus->E()<fMinErat)
      continue;
    
    fHClusters->Fill(cluster++);
    if (clusterEcc < fMinEcc) // eccentricity cut
      continue;
    
    fHClusters->Fill(cluster++);
    if(clus->GetM02()>fM02){
      cout << clus->GetM02() << " ";
      continue;
    }
    //    if(fDoTrMtSmpl){
    //      if(clus->GetNTracksMatched()!=0){
    //        continue;
    //      }
    //    }
    fHClusters->Fill(cluster++);
    
    // fill QA histograms after badmap and cluster cuts
    FillCellQAHists(clus,bdcal,1);
    
    if(bprint)
      clusterVec.Print();
    
    // fill cluster histograms, after cuts
    // eta vs. phi
    fHClustEtaPhi->Fill(clusterVec.Eta(),clusterVec.Phi());
    if (bdcal){
      fHClustEnergyPtDCal->Fill(clusterVec.E(),clusterVec.Pt());
    }
    else{
      fHClustEnergyPt->Fill(clusterVec.E(),clusterVec.Pt());
    }
    //SM number
    Int_t modnumber = GetModuleNumber(clus);
    fHClustEnergySM->Fill(clusterVec.E(),modnumber);
    //fHClustEnergySigma->Fill(clus->E()*maxAxis,clus->E());
    fHClustEnergyTime->Fill(clusterVec.E(),clus->GetTOF());
    fHClustNCellEnergyRatio->Fill(clus->GetNCells(),GetMaxCellEnergy(clus)/clus->E());
    fHClustEnergyNCell->Fill(clus->E(),clus->GetNCells());
    nclusters++;
    
    // mainly store clusters for mixing, if data
    if(!fMcMode){
      
      Double_t En = clus->E();
      // Jasons recalibration
      if(fDoManualRecal){
        Int_t fRecalibrator = 8;
        Double_t recalScale = PrivateEnergyRecal(clus->E(), fRecalibrator);
        En = clus->E()*recalScale;// TOTAL HACK - JJ
        
        clusterVec.SetPx(clusterVec.Px()*recalScale);// TOTAL HACK - JJ
        clusterVec.SetPy(clusterVec.Py()*recalScale);// TOTAL HACK - JJ
        clusterVec.SetPz(clusterVec.Pz()*recalScale);// TOTAL HACK - JJ
        clusterVec.SetE(clusterVec.E()*recalScale);// TOTAL HACK - JJ
        //      Double_t ecorr = 1;
        //Double_t ecorr = fcorrect->Eval(2.0*clusterVec1.E());
      }
      
      fHMeanClusterEnergy->Fill(runnumber,En);
      
      TLorentzVector clusterVecCorr1(clusterVec.Px(),clusterVec.Py(),clusterVec.Pz(),En);
      
      thisEvent.hit[nclusters-1].thishit=clusterVecCorr1;
      if(bdcal){
        thisEvent.hit[nclusters-1].hittype=200;
      }
      else{
        thisEvent.hit[nclusters-1].hittype=100;
      }
      thisEvent.hit[nclusters-1].weight=1.;
      thisEvent.hit[nclusters-1].imo=1;
      thisEvent.hit[nclusters-1].smno=modnumber;
    }
    
    // go through MC information of clusters, store clusters for mixing and do cluster studies using MC info
    else{
      // MC labels
      int ilabel = -1;
      // highest contribution - this one might follow to the first mother
      ilabel = clus->GetLabel();
      // all contributors
      Int_t* mcarr =  clus->GetLabels();
      // number of contributors
      Int_t nl = clus->GetNLabels();
      
      //cout << clus->GetLabel() << endl;
      Bool_t bcl = 1;
      if(nl > 1){
        bcl = 0;
      }
      
      if(ilabel != -1){
        
        // get MC event
        AliMCEvent *mcEvent = MCEvent();
        if (!mcEvent){
          cout << "no MC event" << endl;
          return 0;
        }
        mcEvent->PreReadAll();
        // get MC particle for first label
        AliMCParticle *mcP = static_cast<AliMCParticle*>(mcEvent->GetTrack(ilabel));
        // for the case it did not work
        if (!mcP)
          continue;
        
        Int_t nTracksMC  = mcEvent->GetNumberOfTracks();
        Int_t nPTracksMC = mcEvent->GetNumberOfPrimaries();
        
        // is it generator particle or added signal?
        // one needs to look at the last aka first mother?
        
        // find original particle
        Int_t imother = 1;
        Int_t ipart = mcP->Label();
        Int_t iit = -1;
        Int_t idpi0 = -1;
        
        // loop back to the "top"
        while(imother >= 0){
          AliMCParticle *tmppart = static_cast<AliMCParticle*>(mcEvent->GetTrack(ipart));
          if(bprint){
            cout << " pid " << tmppart->Label() << ", type " << tmppart->PdgCode() << ", mid " << tmppart->GetMother() << " ==>";
          }
          imother = tmppart->GetMother();
          if(tmppart->PdgCode()==111){
            idpi0 = ipart;
            //	    cout << "mother of particle with ID " << mcP->PdgCode() << " is pi0 with " << idpi0 << endl;
          }
          if(imother >= 0){
            ipart = imother;
          }
          iit++;
        }
        fHNMothers->Fill(iit,mcP->E());
        
        // calculate weight (for added signals)
        Float_t wgt = 1.;
        
        AliMCParticle *McMo = static_cast<AliMCParticle*>(mcEvent->GetTrack(ipart));
        Double_t mcPt = McMo->Pt();
        Double_t mcEta = McMo->Eta();
        
        wgt = CalcWeight(mcPt,mcEta,1);
        
        // check if generated or added
        bool bGen = kTRUE;
        bool bAddPi0 = kFALSE;
        bool bAddEta = kFALSE;
        
        // go through MC headers to associate MC particle with correct header
        if(pythiaHeader && fAddedSignal){
          if(ipart > ipymax){
            bGen = kFALSE;
            if(ipart <= ipi0max){
              bAddPi0 = kTRUE;
              if(bprint)
                cout << " added pi ";
            }
            else if(ipart > ipi0max && ipart <= ietamax){
              bAddEta = kTRUE;
              if(bprint)
                cout << " added eta ";
            }
          }
          else{
            if(bprint)
              cout << " pythia ";
          }
        }
        if(bprint)
          cout << endl;
        
        // if not added signals, distinguish between clusters from
        // 1) primary pi0, 2) secondary pi0 (not K0), 3) pi0 from K0, 4) pi0 from material
        // still needs to be implemented for DCal, I think
        
        // loop up until pi0, then look for mother of pi0
        
        // if it is not from a pi0, flag 100
        // from primary pi0, flag 101
        // from secondary pi0, flag 102
        // from K0, flag 103
        // from material, flag 104
        // only for generator, not added signals
        Int_t tmpflag = 0;
        Bool_t bpi0 = 1;
        if(idpi0 < 0)
          bpi0 = 0;
        if(bGen){
          if(bpi0 == 0){
            tmpflag = 100;
          }
          else{
            AliMCParticle *tmppart = static_cast<AliMCParticle*>(mcEvent->GetTrack(idpi0));
            Int_t itmp = tmppart->GetMother();
            if(idpi0<nPTracksMC){
              tmpflag = 101;
            }
            else{
              tmpflag = 102;
              if( ((AliMCParticle*)mcEvent->GetTrack(tmppart->GetMother()))->PdgCode() ==  310 ||
                 ((AliMCParticle*)mcEvent->GetTrack(tmppart->GetMother()))->PdgCode() == -310  ){
                tmpflag = 103;
              }
              else if(((AliMCParticle*)mcEvent->GetTrack(tmppart->GetMother()))->PdgCode() ==  2212 || //proton
                      ((AliMCParticle*)mcEvent->GetTrack(tmppart->GetMother()))->PdgCode() == -2212 || //anti-proton
                      ((AliMCParticle*)mcEvent->GetTrack(tmppart->GetMother()))->PdgCode() ==  2112 || //neutron
                      ((AliMCParticle*)mcEvent->GetTrack(tmppart->GetMother()))->PdgCode() == -2112 || //anti-neutron
                      ((AliMCParticle*)mcEvent->GetTrack(tmppart->GetMother()))->PdgCode() ==  321  || //K+
                      ((AliMCParticle*)mcEvent->GetTrack(tmppart->GetMother()))->PdgCode() == -321  || //K-
                      ((AliMCParticle*)mcEvent->GetTrack(tmppart->GetMother()))->PdgCode() ==  211  || //pi+
                      ((AliMCParticle*)mcEvent->GetTrack(tmppart->GetMother()))->PdgCode() == -211 ){    //pi-)
                tmpflag = 104;
              }
            }
          }
        }
        // sum of energy of all MC particles in cluster
        Double_t esum = 0;
        for(int ip=0;ip<nl;ip++){
          Int_t entry = mcarr[ip];
          AliMCParticle *mcPart = static_cast<AliMCParticle*>(mcEvent->GetTrack(entry));
          esum += mcPart->E();
        }
        Double_t efrac = 0;
        
        Double_t mce = mcP->E();
        Double_t cle = clus->E();
        
        // energy fraction of "leading" particle
        if(esum!=0)
          efrac = mcP->E()/esum;
        
        //        fHRecTrue->SetXTtitle("E_{clu}");
        //        fHRecTrue->SetYTtitle("E_{tru}/E_{clu}");
        //cout << "filling histos" << endl;
        // if generator
        if(bGen){
          // if leading particle with high contribution and photon
          if(mcP->PdgCode() == 22 && efrac>=0.99){
            // lets fill histograms
            fHMCpartfrac->Fill(mce/cle,efrac);
            fHECluEMC->Fill(mce,cle);
            //            fHRecTrue->Fill(cle,mce/cle);
          }
          // if photon
          else if(mcP->PdgCode() == 22 && efrac<0.99){
            fHMCpartfrac->Fill(mce/cle,efrac);
            fHECluEMCnofull->Fill(mce,cle);
          }
          // if electron with high contribution
          else if((mcP->PdgCode() == 11 || mcP->PdgCode() == -11) && efrac>=0.99){
            fHECluEMCelectron->Fill(mce,cle);
          }
          // if pion with high contribution
          else if((mcP->PdgCode() == 211 || mcP->PdgCode() == -211) && efrac>=0.99){
            fHECluEMCpion->Fill(mce,cle);
          }
          // if kaon with high contribution
          else if((mcP->PdgCode() == 321 || mcP->PdgCode() == -321) && efrac>=0.99){
            fHECluEMCkaon->Fill(mce,cle);
          }
          // if other particle with high contribution
          else if(efrac>=0.99){
            fHECluEMCother->Fill(mce,cle);
          }
          
          // should I look for two photon clusters? ingredients: 2 photons, from same mother.
          for(int ip=0;ip<nl-1;ip++){
            Int_t entry = mcarr[ip];
            AliMCParticle *mcPart = static_cast<AliMCParticle*>(mcEvent->GetTrack(entry));
            imother = mcPart->GetMother();
            for(int jp=ip;jp<nl;jp++){
              Int_t entry2 = mcarr[jp];
              AliMCParticle *mcPart2 = static_cast<AliMCParticle*>(mcEvent->GetTrack(entry2));
              Int_t jmother = mcPart2->GetMother();
              if(imother == jmother && imother >= 0){
                AliMCParticle* mcMother = static_cast<AliMCParticle*>(mcEvent->GetTrack(imother));
                if(mcMother->PdgCode() == 111){
                  Double_t mcmop = sqrt(mcMother->P()*mcMother->P() + 0.135*0.135);
                  Double_t clue = clus->E();
                  fHECluEMCpi0single->Fill(mcmop,clue);
                }
              }
            }
          } // end two photon clusters
        } // end if generator
        if(bAddPi0){
          // if leading particle with high contribution and photon
          if(mcP->PdgCode() == 22 && efrac>0.5){
            // lets fill histograms
            fHECluEMCAddPi0->Fill(mce,cle);
            //fHRecTrueAddPi0->Fill(cle,mce/cle);
          }
        }
        if(bAddEta){
          // if leading particle with high contribution and photon
          if(mcP->PdgCode() == 22 && efrac>0.5){
            // lets fill histograms
            fHECluEMCAddEta->Fill(mce,cle);
            //fHRecTrueAddEta->Fill(cle,mce/cle);
          }
        }
        if(bAddPi0 || bAddEta){
          // if photon
          if(mcP->PdgCode() == 22 && efrac<0.99){
            fHECluEMCnofullAdd->Fill(mce,cle);
          }
        } // end else
        
        
        // Jasons recalibration
        Double_t En = clus->E();
        if(fDoManualRecal){
          Int_t fRecalibrator = 9;
          Double_t recalScale = PrivateEnergyRecal(clus->E(), fRecalibrator);
          En = clus->E()*recalScale;// TOTAL HACK - JJ
          
          //clus->GetMomentum(clusterVec,vertex);
          clusterVec.SetPx(clusterVec.Px()*recalScale);// TOTAL HACK - JJ
          clusterVec.SetPy(clusterVec.Py()*recalScale);// TOTAL HACK - JJ
          clusterVec.SetPz(clusterVec.Pz()*recalScale);// TOTAL HACK - JJ
          clusterVec.SetE(clusterVec.E()*recalScale);// TOTAL HACK - JJ
          //      Double_t ecorr = 1;
          //Double_t ecorr = fcorrect->Eval(2.0*clusterVec1.E());
        }
        TLorentzVector clusterVecCorr1(clusterVec.Px(),clusterVec.Py(),clusterVec.Pz(),En);
        
        Bool_t bkeep = 1;
        // for added signals, if particle is an electron, see if it is the higher energetic partner from conversion
        //        if(!bGen && (bAddPi0 || bAddEta)){
        //
        //          if(mcP->PdgCode() == 11 || mcP->PdgCode() == -11){
        //            imother = mcP->GetMother();
        //            AliMCParticle *tmppart2 = static_cast<AliMCParticle*>(mcEvent->GetTrack(imother));
        //            // check if mother is a photon
        //            if(tmppart2->PdgCode() == 22){
        //              Int_t d1 = tmppart2->GetDaughterFirst();
        //              Int_t d2 = tmppart2->GetDaughterLast();
        //              if (d1>0){
        //                if (d2<0){
        //                  d2=d1;
        //                }
        //                if(d2-d1 == 1){
        //                  for (Int_t ida=d1;ida<=d2;++ida) {
        //                    AliMCParticle *tmpdaughter = static_cast<AliMCParticle *>(mcEvent->GetTrack(ida));
        //                    if(tmpdaughter->E() > mcP->E())
        //                      bkeep = 0;
        //                  }
        //                }
        //              }
        //            }
        //          }
        //        }
        
        // keep only eta -> gammagamma
        if(bAddEta){
          Int_t d1 = McMo->GetDaughterFirst();
          Int_t d2 = McMo->GetDaughterLast();
          if (d1>0){
            if (d2<0){
              d2=d1;
            }
            if(d2-d1 != 1){
              bkeep = 0;
            }
            else{
              AliMCParticle *tmpdaughter1 = static_cast<AliMCParticle *>(mcEvent->GetTrack(d1));
              AliMCParticle *tmpdaughter2 = static_cast<AliMCParticle *>(mcEvent->GetTrack(d2));
              if(!(tmpdaughter1->PdgCode() == 22 && tmpdaughter2->PdgCode() == 22)){
                bkeep = 0;
              }
            }
          }
          
        }
        // store clusters for mixing
        if(1){
          thisEvent.hit[nclusters-1].thishit=clusterVecCorr1;
          thisEvent.hit[nclusters-1].imo=ipart;
          thisEvent.hit[nclusters-1].pid=mcP->PdgCode();
          thisEvent.hit[nclusters-1].smno=modnumber;
          if(!bGen){
            thisEvent.hit[nclusters-1].weight = wgt;
          }
          else{
            thisEvent.hit[nclusters-1].weight = 1.;
          }
          
          thisEvent.hit[nclusters-1].hittype=-1;
          if(bGen){
            // use hittype to flag particle history
            // if it is not from a pi0, flag 100
            // from primary pi0, flag 101
            // from secondary pi0, flag 102
            // from K0, flag 103
            // from material, flag 104
            if(bdcal)
              tmpflag += 100;
            thisEvent.hit[nclusters-1].hittype=tmpflag;
          }
          else{
            if(bAddPi0){
              thisEvent.hit[nclusters-1].hittype=1;
            }
            else if(bAddEta){
              thisEvent.hit[nclusters-1].hittype=2;
            }
          }
          thisEvent.hit[nclusters-1].bclean = bcl;
        }
      } // end ilabel
    } // end fmcmode
  } // end cluster loop
  fHMeanClusterNumber->Fill(runnumber,nclusters);
  fHClustAccEvt->Fill(nclusters);
  thisEvent.SetGlobalInfo(nclusters,max_phi,max_theta);
  if(bprint){
    thisEvent.Print();
    cout << " " << endl;
  }
  return 1;
  //  return max_pt;
}


//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::CalcMcInfo()
{
  // Get Mc truth particle information and store it
  if (!fMcMode)
    return;
  
  if (!fMcParts)
    return;
  
  fMcParts->Clear();
  
  AliEMCALEMCGeometry *emc = fGeom->GetEMCGeometry();
  Double_t etamin = emc->GetArm1EtaMin();
  Double_t etamax = emc->GetArm1EtaMax();
  Double_t phimin = emc->GetArm1PhiMin()*TMath::DegToRad();
  Double_t phimax = emc->GetArm1PhiMax()*TMath::DegToRad();
  
  //  if (fAodEv) {
  //    AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  //    am->LoadBranch(AliAODMCParticle::StdBranchName());
  //    TClonesArray *tca = dynamic_cast<TClonesArray*>(fAodEv->FindListObject(AliAODMCParticle::StdBranchName()));
  //    if (!tca)
  //      return;
  //
  //    Int_t nents = tca->GetEntries();
  //    for(int it=0; it<nents; ++it) {
  //      AliAODMCParticle *part = static_cast<AliAODMCParticle*>(tca->At(it));
  //      part->Print();
  //
  //      // pion or eta meson or direct photon
  //      if(part->GetPdgCode() == 111) {
  //      } else if(part->GetPdgCode() == 221) {
  //      } else if(part->GetPdgCode() == 22 ) {
  //      }	else
  //        continue;
  //
  //      // primary particle
  //      Double_t dR = TMath::Sqrt((part->Xv()*part->Xv())+(part->Yv()*part->Yv()));
  //      if(dR > 1.0)
  //        continue;
  //
  //      // kinematic cuts
  //      Double_t pt = part->Pt() ;
  //      if (pt<0.5)
  //        continue;
  //      Double_t eta = part->Eta();
  //      if (eta<etamin||eta>etamax)
  //        continue;
  //      Double_t phi  = part->Phi();
  //      if (phi<phimin||phi>phimax)
  //        continue;
  //
  //      ProcessDaughters(part, it, tca);
  //      PrintDaughters(part, tca, 2);
  //    }
  //    return;
  //  }
  
  // get MC event
  AliMCEvent *mcEvent = MCEvent();
  if (!mcEvent){
    cout << "no MC event" << endl;
    return;
  }
  
  // get vertex
  const AliVVertex *evtVtx = mcEvent->GetPrimaryVertex();
  if (!evtVtx)
    return;
  
  // read event
  mcEvent->PreReadAll();
  
  // get number of MC particles
  Int_t nTracks = mcEvent->GetNumberOfPrimaries();

  // loop through MC particles
  for (Int_t iTrack = 0; iTrack<nTracks; ++iTrack) {
    // get particle at index iTrack
    AliMCParticle *mcP = static_cast<AliMCParticle*>(mcEvent->GetTrack(iTrack));
    if (!mcP)
      continue;
    
    // pion or eta meson or direct photon
    if(mcP->PdgCode() == 111) {
    } else if(mcP->PdgCode() == 221) {
    } else if(mcP->PdgCode() == 22 ) {
    } else
      continue;
    
    // primary particle - should be accounted for already
    // check the radius from the vertex, if it is from decay, radius is larger than 0
    Double_t dR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX())*(mcP->Xv()-evtVtx->GetX()) +
                              (mcP->Yv()-evtVtx->GetY())*(mcP->Yv()-evtVtx->GetY()));
    
    
    if(dR > 0.1)
      continue;
    
    // use nprimarytracks (these are the first ones)
    
    // is it generator particle or added signal?
    bool bGen = kTRUE;
    bool bAddPi0 = kFALSE;
    bool bAddEta = kFALSE;
    
    if(pythiaHeader && fAddedSignal){
      if(mcP->Label() > ipymax){
        bGen = kFALSE;
        if(mcP->Label() <= ipi0max){
          bAddPi0 = kTRUE;
        }
        else if(mcP->Label() > ipi0max && mcP->Label() <= ietamax){
          bAddEta = kTRUE;
        }
      }
    }
    
    double mcPt = mcP->Pt();
    double mcEta = mcP->Eta();
    double mcY = mcP->Y();
    Float_t wgt2 = 1.;
    
    if(bGen && !(bAddPi0 || bAddEta)){
      // fill truth histogram
      if(mcP->PdgCode() == 111){
        fHPionTruthPt->Fill(mcP->Pt());
        fHPyPionEtaPt->Fill(mcP->Eta(),mcP->Pt());
      }
      if(mcP->PdgCode() == 221){
        fHEtaTruthPt->Fill(mcP->Pt());
      }
    }
    
    // everything else should be secondary
    
    if(bAddPi0 && !(bGen || bAddEta)){
      // fill truth histogram
      wgt2 = CalcWeight(mcPt,mcEta,1);
      fHWgt->Fill(wgt2);
      if(mcP->PdgCode() == 111){
        fHPionTruthPtAdd->Fill(mcP->Pt(),wgt2);
        if(fAddedSignal){
          fHAddPionEtaPt->Fill(mcP->Eta(),mcP->Pt());
          fHAddPionEtaPtWgt->Fill(mcP->Eta(),mcP->Pt(),wgt2);
        }
      }
    }
    if(bAddEta && !(bGen || bAddPi0)){
      if(mcPt > 0.01 && mcPt < 30 && mcEta > -2 && mcEta < 2){
        wgt2 = CalcWeight(mcPt,mcEta,1);
      }
      else{
        continue;
      }
      if(mcP->PdgCode() == 221){
        fHEtaTruthPtAdd->Fill(mcP->Pt(),wgt2);
      }
    }
    
    if(mcP->PdgCode() == 22){
      fHGamTruthPt->Fill(mcP->Pt());
    }
    
    // kinematic cuts
    if (mcPt<0.1)
      continue;
    if (mcY<-1.0||mcY>1.0)
      continue;
    
    //    Double_t phi  = mcP->Phi();
    //    if (phi<phimin||phi>phimax)
    //      continue;
    
    // fill truth histogram in acceptance
    // photon tracks - then continue
    if(mcP->PdgCode() == 22){
      fHGamTruthPtIn->Fill(mcP->Pt());
      Double_t phi  = mcP->Phi();
      if (phi>=phimin&&phi<=phimax){
        fHGamTruthPtAcc->Fill(mcP->Pt());
      }
      continue;
    }
    
    // for pi0 and eta, we will have to consider decay photons
    // thus, need to loop through the "daughters" and see if they are 2 photons or 1 photon and one converted photons
    // then, check if both are in acceptance
    
    if(bGen && !bAddEta && !bAddPi0){
      // fill truth histogram for input
      if(mcP->PdgCode() == 111){
        fHPionTruthPtIn->Fill(mcP->Pt());
      }
      if(mcP->PdgCode() == 221){
        fHEtaTruthPtIn->Fill(mcP->Pt());
      }
    }
    
    if(bAddPi0 && !bAddEta && !bGen){
      if(mcPt > 0.01 && mcPt < 30 && mcY > -2 && mcY < 2){
        wgt2 = CalcWeight(mcPt,mcEta,1);
      }
      else{
        continue;
      }
      // fill truth histogram for input
      if(mcP->PdgCode() == 111){
        fHPionTruthPtInAdd->Fill(mcP->Pt(),wgt2);
      }
    }
    
    if(bAddEta && !bAddPi0 && !bGen){
      if(mcP->PdgCode() == 221){
        if(mcPt > 0.01 && mcPt < 30 && mcY > -2 && mcY < 2){
          wgt2 = CalcWeight(mcPt,mcEta,1);
        }
        else{
          continue;
        }
        fHEtaTruthPtInAdd->Fill(mcP->Pt(),wgt2);
      }
      
    }
    
    bool binp = true; // does it count as input?
    
    Int_t d1 = mcP->GetDaughterFirst();
    Int_t d2 = mcP->GetDaughterLast();
    if (d1<0)
      return;
    if (d2<0)
      d2=d1;
    
    if(d2-d1 != 1){
      continue;
    }
    
    // check it particle decays into two photons (binp), decay photons are in acceptance (bacc/baccconv)
    bool bacc = true;
    for (Int_t i=d1;i<=d2;++i) {
      const AliMCParticle *dmc = static_cast<const AliMCParticle *>(mcEvent->GetTrack(i));
      Double_t eta_d = dmc->Eta();
      Double_t phi_d = dmc->Phi();
      if(!(dmc->PdgCode()==22)){
        binp = false;
      }
      if(!(dmc->PdgCode()==22 && eta_d>etamin && eta_d<etamax && phi_d>phimin && phi_d<phimax)){
        bacc = false;
      }
    }
    
    if(bGen && !bAddEta && !bAddPi0){
      // if both photons are on EMCal
      if(binp){
        if(mcP->PdgCode() == 111 && bacc){
          fHPionTruthPtAcc->Fill(mcP->Pt());
        }
        if(mcP->PdgCode() == 221 && bacc){
          fHEtaTruthPtAcc->Fill(mcP->Pt());
        }
      }
    }
    
    if(bAddPi0 && !bAddEta && !bGen){
      // if both photons are on EMCal
      if(binp){
        if(mcP->PdgCode() == 111 && bacc){
          fHPionTruthPtAccAdd->Fill(mcP->Pt(),wgt2);
        }
      }
    }
    
    
    if(bAddEta && !bAddPi0 && !bGen){
      // if both photons are on EMCal
      if(binp){
        if(mcP->PdgCode() == 221 && bacc){
          fHEtaTruthPtAccAdd->Fill(mcP->Pt(),wgt2);
        }
      }
      
    }
    ProcessDaughters(mcP, iTrack, mcEvent);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::FillNtuple()
{
  // Fill ntuple.
  
  if (!fNtuple)
    return;
  
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (fAodEv) {
    AliAODHeader * aodheader = dynamic_cast<AliAODHeader*>(fAodEv->GetHeader());
    if(!aodheader) AliFatal("Not a standard AOD");
    
    fHeader->fRun            = fAodEv->GetRunNumber();
    fHeader->fOrbit          = aodheader->GetOrbitNumber();
    fHeader->fPeriod         = aodheader->GetPeriodNumber();
    fHeader->fBx             = aodheader->GetBunchCrossNumber();
    fHeader->fL0             = aodheader->GetL0TriggerInputs();
    fHeader->fL1             = aodheader->GetL1TriggerInputs();
    fHeader->fL2             = aodheader->GetL2TriggerInputs();
    fHeader->fTrClassMask    = aodheader->GetTriggerMask();
    fHeader->fTrCluster      = aodheader->GetTriggerCluster();
    fHeader->fOffTriggers    = aodheader->GetOfflineTrigger();
    fHeader->fFiredTriggers  = aodheader->GetFiredTriggerClasses();
  } else {
    fHeader->fRun            = fEsdEv->GetRunNumber();
    fHeader->fOrbit          = fEsdEv->GetHeader()->GetOrbitNumber();
    fHeader->fPeriod         = fEsdEv->GetESDRun()->GetPeriodNumber();
    fHeader->fBx             = fEsdEv->GetHeader()->GetBunchCrossNumber();
    fHeader->fL0             = fEsdEv->GetHeader()->GetL0TriggerInputs();
    fHeader->fL1             = fEsdEv->GetHeader()->GetL1TriggerInputs();
    fHeader->fL2             = fEsdEv->GetHeader()->GetL2TriggerInputs();
    fHeader->fTrClassMask    = fEsdEv->GetHeader()->GetTriggerMask();
    fHeader->fTrCluster      = fEsdEv->GetHeader()->GetTriggerCluster();
    fHeader->fOffTriggers    = ((AliInputEventHandler*)(am->GetInputEventHandler()))->IsEventSelected();
    fHeader->fFiredTriggers  = fEsdEv->GetFiredTriggerClasses();
    Float_t v0CorrR = 0;
    fHeader->fV0 = AliESDUtils::GetCorrV0(fEsdEv,v0CorrR);
    const AliMultiplicity *mult = fEsdEv->GetMultiplicity();
    if (mult)
      fHeader->fCl1 = mult->GetNumberOfITSClusters(1);
    fHeader->fTr = AliESDtrackCuts::GetReferenceMultiplicity(fEsdEv,1);
    AliTriggerAnalysis trAn; /// Trigger Analysis
    Bool_t v0B = trAn.IsOfflineTriggerFired(fEsdEv, AliTriggerAnalysis::kV0C);
    Bool_t v0A = trAn.IsOfflineTriggerFired(fEsdEv, AliTriggerAnalysis::kV0A);
    fHeader->fV0And        = v0A && v0B;
    fHeader->fIsHT         = (fHeader->fOffTriggers & AliVEvent::kEMC1) || (fHeader->fOffTriggers & AliVEvent::kEMC7);
    am->LoadBranch("SPDPileupVertices");
    am->LoadBranch("TrkPileupVertices");
    fHeader->fIsPileup     = fEsdEv->IsPileupFromSPD(3,0.8);
    fHeader->fIsPileup2    = fEsdEv->IsPileupFromSPD(3,0.4);
    fHeader->fIsPileup4    = fEsdEv->IsPileupFromSPD(3,0.2);
    fHeader->fIsPileup8    = fEsdEv->IsPileupFromSPD(3,0.1);
    fHeader->fNSpdVertices = fEsdEv->GetNumberOfPileupVerticesSPD();
    fHeader->fNTpcVertices = fEsdEv->GetNumberOfPileupVerticesTracks();
  }
  
  AliCentrality *cent = InputEvent()->GetCentrality();
  fHeader->fV0Cent    = cent->GetCentralityPercentileUnchecked("V0M");
  fHeader->fCl1Cent   = cent->GetCentralityPercentileUnchecked("CL1");
  fHeader->fTrCent    = cent->GetCentralityPercentileUnchecked("TRK");
  fHeader->fCqual     = cent->GetQuality();
  
  AliEventplane *ep = InputEvent()->GetEventplane();
  if (ep) {
    if (ep->GetQVector())
      fHeader->fPsi     = ep->GetQVector()->Phi()/2. ;
    else
      fHeader->fPsi = -1;
    if (ep->GetQsub1()&&ep->GetQsub2())
      fHeader->fPsiRes  = ep->GetQsub1()->Phi()/2.-ep->GetQsub2()->Phi()/2.;
    else
      fHeader->fPsiRes = 0;
  }
  
  Double_t val = 0;
  TString trgclasses(fHeader->fFiredTriggers);
  for (Int_t j = 0; j<fTrClassNamesArr->GetEntries(); ++j) {
    const char *name = fTrClassNamesArr->At(j)->GetName();
    TRegexp regexp(name);
    if (trgclasses.Contains(regexp))
      val += TMath::Power(2,j);
  }
  fHeader->fTcls = (UInt_t)val;
  
  fHeader->fNSelTr     = fSelTracks->GetEntries();
  fHeader->fNSelPrimTr = fSelPrimTracks->GetEntries();
  fHeader->fNSelPrimTr1   = 0;
  fHeader->fNSelPrimTr2   = 0;
  for(int iTracks=0; iTracks < fSelPrimTracks->GetEntries(); iTracks++){
    AliESDtrack *track = static_cast<AliESDtrack*>(fSelPrimTracks->At(iTracks));
    if(track->Pt()>1)
      ++fHeader->fNSelPrimTr1;
    if(track->Pt()>2)
      ++fHeader->fNSelPrimTr2;
  }
  
  fHeader->fNCells   = 0;
  fHeader->fNCells0  = 0;
  fHeader->fNCells01 = 0;
  fHeader->fNCells03 = 0;
  fHeader->fNCells1  = 0;
  fHeader->fNCells2  = 0;
  fHeader->fNCells5  = 0;
  fHeader->fMaxCellE = 0;
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  
  if (cells) {
    Int_t ncells = cells->GetNumberOfCells();
    for(Int_t j=0; j<ncells; ++j) {
      Double_t cellen = cells->GetAmplitude(j);
      if (cellen>0.045)
        ++fHeader->fNCells0;
      if (cellen>0.1)
        ++fHeader->fNCells01;
      if (cellen>0.3)
        ++fHeader->fNCells03;
      if (cellen>1)
        ++fHeader->fNCells1;
      if (cellen>2)
        ++fHeader->fNCells2;
      if (cellen>5)
        ++fHeader->fNCells5;
      if (cellen>fHeader->fMaxCellE)
        fHeader->fMaxCellE = cellen;
    }
    fHeader->fNCells = ncells;
  }
  
  fHeader->fNClus      = 0;
  fHeader->fNClus1     = 0;
  fHeader->fNClus2     = 0;
  fHeader->fNClus5     = 0;
  fHeader->fMaxClusE   = 0;
  
  TObjArray *clusters = fEsdClusters;
  if (!clusters)
    clusters = fAodClusters;
  
  if (clusters) {
    Int_t nclus = clusters->GetEntries();
    for(Int_t j=0; j<nclus; ++j) {
      AliVCluster *clus = static_cast<AliVCluster*>(clusters->At(j));
      if (!clus->IsEMCAL())
        continue;
      Double_t clusen = clus->E();
      if (clusen>1)
        ++fHeader->fNClus1;
      if (clusen>2)
        ++fHeader->fNClus2;
      if (clusen>5)
        ++fHeader->fNClus5;
      if (clusen>fHeader->fMaxClusE)
        fHeader->fMaxClusE = clusen;
    }
    fHeader->fNClus = nclus;
  }
  
  fHeader->fMaxTrE     = 0;
  if (fTriggers) {
    Int_t ntrig = fTriggers->GetEntries();
    for (Int_t j = 0; j<ntrig; ++j) {
      AliStaTrigger *sta = static_cast<AliStaTrigger*>(fTriggers->At(j));
      if (!sta)
        continue;
      if (sta->fE>fHeader->fMaxTrE)
        fHeader->fMaxTrE = sta->fE;
    }
  }
  
  // count cells above 100 MeV on super modules
  fHeader->fNcSM0 = GetNCells(0, 0.100);
  fHeader->fNcSM1 = GetNCells(1, 0.100);
  fHeader->fNcSM2 = GetNCells(2, 0.100);
  fHeader->fNcSM3 = GetNCells(3, 0.100);
  fHeader->fNcSM4 = GetNCells(4, 0.100);
  fHeader->fNcSM5 = GetNCells(5, 0.100);
  fHeader->fNcSM6 = GetNCells(6, 0.100);
  fHeader->fNcSM7 = GetNCells(7, 0.100);
  fHeader->fNcSM8 = GetNCells(8, 0.100);
  fHeader->fNcSM9 = GetNCells(9, 0.100);
  
  if (fAodEv) {
    am->LoadBranch("vertices");
    AliAODVertex *pv = fAodEv->GetPrimaryVertex();
    FillVertex(fPrimVert, pv);
    AliAODVertex *sv = fAodEv->GetPrimaryVertexSPD();
    FillVertex(fSpdVert, sv);
  } else {
    am->LoadBranch("PrimaryVertex.");
    const AliESDVertex *pv = fEsdEv->GetPrimaryVertexTracks();
    FillVertex(fPrimVert, pv);
    am->LoadBranch("SPDVertex.");
    const AliESDVertex *sv = fEsdEv->GetPrimaryVertexSPD();
    FillVertex(fSpdVert, sv);
    am->LoadBranch("TPCVertex.");
    const AliESDVertex *tv = fEsdEv->GetPrimaryVertexTPC();
    FillVertex(fTpcVert, tv);
  }
  
  fNtuple->Fill();
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::FillPionHists()
{
  // Fill histograms related to pions.
  TLorentzVector clusterVec1;
  TLorentzVector clusterVec2;
  TLorentzVector pionVec;
  Short_t hitclass1 = 0;
  Short_t hitclass2 = 0;
  
  //    Int_t nclus = clusters->GetEntries();
  Int_t nclus = thisEvent.nHits;
  if(nclus<2){
    return;
  }
  // loop over clusters
  for (Int_t i = 0; i<nclus; ++i) {
    // get 1st cluster
    clusterVec1 = thisEvent.hit[i].thishit;
    hitclass1 = thisEvent.hit[i].hittype;
    Double_t wght = 1.;
    wght = thisEvent.hit[i].weight;
    
    // loop over 2nd clusters
    for (Int_t j = i+1; j<nclus; ++j) {
      
      // get 2nd cluster
      clusterVec2 = thisEvent.hit[j].thishit;
      hitclass2 = thisEvent.hit[j].hittype;
      
      // calculate distance between clusters
      Double_t d_phi = clusterVec1.Phi() - clusterVec2.Phi();
      Double_t d_eta = clusterVec1.Eta() - clusterVec2.Eta();
      Double_t d_r = sqrt(d_phi*d_phi + d_eta*d_eta);
      
      // calculate pair vector
      
      pionVec = clusterVec1 + clusterVec2;
      fHdr->Fill(d_r);
      
      // cut on minimum distance, should be made settable
      if(d_r < 0.01)
        continue;
      
      // asymmetry
      Double_t pionZgg = TMath::Abs(clusterVec1.E()-clusterVec2.E())/pionVec.E();
      Double_t pionOpeningAngle = clusterVec1.Angle(clusterVec2.Vect());
      
      if(fCalibRun && pionZgg < 0.2){
        //EMCal histos
        if((hitclass1 >= 100 && hitclass2 >= 100) && (hitclass1 < 200 && hitclass2 < 200)){
          fHPionInvMassesEMCalCalib->Fill(pionVec.M(),(clusterVec1.E() + clusterVec2.E())/2);
        }
        else if(hitclass1 >= 200 && hitclass2 >= 200){
          fHPionInvMassesDCalCalib->Fill(pionVec.M(),(clusterVec1.E() + clusterVec2.E())/2);
        }
      }
      
      // fill all histograms with inv masses
      if (pionZgg < fAsymMax1) {
        fHPionMggPt->Fill(pionVec.M(),pionVec.Pt());
        fHPionMggAsym->Fill(pionVec.M(),pionZgg);
        fHPionMggDgg->Fill(pionVec.M(),pionOpeningAngle);
        
        if(pionVec.M()> 0.11 && pionVec.M() < 0.17)
          fHPionPtAsym->Fill(pionVec.Pt(),pionZgg);
        
        if(fMcMode){
          if(hitclass1 == 1 && hitclass2 == 1 && fAddedSignal){
            if(thisEvent.hit[i].imo == thisEvent.hit[j].imo){
              // all, no weight
              fHPionInvMassesAdd1NoWgt->Fill(pionVec.M(),pionVec.Pt());
              
              // clean clusters only
              if(thisEvent.hit[i].bclean == 1 && thisEvent.hit[j].bclean == 1){
                // photons only
                if(thisEvent.hit[i].pid == 22 && thisEvent.hit[j].pid == 22){
                  fHPionInvMassesGamAdd1Sym->Fill(pionVec.M(),pionVec.Pt(),wght);
                }
                // one non-photon
                else{
                  fHPionInvMassesAdd1Sym->Fill(pionVec.M(),pionVec.Pt(),wght);
                }
              }
              // one non-clean cluster
              else{
                // photons only
                if(thisEvent.hit[i].pid == 22 && thisEvent.hit[j].pid == 22){
                  fHPionInvMassesGamAdd1MultSym->Fill(pionVec.M(),pionVec.Pt(),wght);
                }
                // one non-photon
                else{
                  fHPionInvMassesAdd1MultSym->Fill(pionVec.M(),pionVec.Pt(),wght);
                }
              }
            }
          }
          else if(hitclass1 == 2 && hitclass2 == 2 && fAddedSignal){
            if(thisEvent.hit[i].imo == thisEvent.hit[j].imo){
              fHPionInvMassesAdd2->Fill(pionVec.M(),pionVec.Pt(),wght);
            }
          }
          else if((hitclass1 >= 100 && hitclass2 >= 100) && (hitclass1 < 200 && hitclass2 < 200)){
            fHPionInvMassesSym->Fill(pionVec.M(),pionVec.Pt());
            if(pionVec.M()> 0.11 && pionVec.M() < 0.17)
              fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
            
            // look at the different cases now
            // use only hits from same mother particle
            if(thisEvent.hit[i].imo == thisEvent.hit[j].imo){
              if(hitclass1 == 101 && hitclass2 == 101){
                fHPriPionInvMasses->Fill(pionVec.M(),pionVec.Pt());
              }
              else if(hitclass1 == 102 && hitclass2 == 102){
                fHSecPionInvMasses->Fill(pionVec.M(),pionVec.Pt());
              }
              else if(hitclass1 == 103 && hitclass2 == 103){
                fHK0PionInvMasses->Fill(pionVec.M(),pionVec.Pt());
              }
              else if(hitclass1 == 104 && hitclass2 == 104){
                fHMatPionInvMasses->Fill(pionVec.M(),pionVec.Pt());
              }
            }
          }
          else if(hitclass1 >= 200 && hitclass2 >= 200){
            fHPionInvMassesDCal->Fill(pionVec.M(),pionVec.Pt());
            if(pionVec.M()> 0.11 && pionVec.M() < 0.17)
              fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
          }
        }
        else{
          if(hitclass1 < 200 && hitclass2 < 200){
            fHPionInvMassesSym->Fill(pionVec.M(),pionVec.Pt());
            if(pionVec.M()> 0.11 && pionVec.M() < 0.17)
              fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
          }
          else if(hitclass1 >= 200 && hitclass2 >= 200) {
            fHPionInvMassesDCal->Fill(pionVec.M(),pionVec.Pt());
            if(pionVec.M()> 0.11 && pionVec.M() < 0.17)
              fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
          }
          else{
            fHPionInvMassesEMCalDCal->Fill(pionVec.M(),pionVec.Pt());
            if(pionVec.M()> 0.11 && pionVec.M() < 0.17)
              fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
          }
        }
      }
      // asymmetric
      else if(pionZgg < fAsymMax2){
        fHPionMggPt->Fill(pionVec.M(),pionVec.Pt());
        fHPionMggAsym->Fill(pionVec.M(),pionZgg);
        fHPionMggDgg->Fill(pionVec.M(),pionOpeningAngle);
        
        if(pionVec.M()> 0.11 && pionVec.M() < 0.17)
          fHPionPtAsym->Fill(pionVec.Pt(),pionZgg);
        
        if(fMcMode){
          if(hitclass1 == 1 && hitclass2 == 1 && fAddedSignal){
            if(thisEvent.hit[i].imo == thisEvent.hit[j].imo){
              // all, no weight
              fHPionInvMassesAdd1NoWgt->Fill(pionVec.M(),pionVec.Pt());
              
              // clean clusters only
              if(thisEvent.hit[i].bclean == 1 && thisEvent.hit[j].bclean == 1){
                // photons only
                if(thisEvent.hit[i].pid == 22 && thisEvent.hit[j].pid == 22){
                  fHPionInvMassesGamAdd1->Fill(pionVec.M(),pionVec.Pt(),wght);
                }
                // one non-photon
                else{
                  fHPionInvMassesAdd1->Fill(pionVec.M(),pionVec.Pt(),wght);
                }
              }
              // one non-clean cluster
              else{
                // photons only
                if(thisEvent.hit[i].pid == 22 && thisEvent.hit[j].pid == 22){
                  fHPionInvMassesGamAdd1Mult->Fill(pionVec.M(),pionVec.Pt(),wght);
                }
                // one non-photon
                else{
                  fHPionInvMassesAdd1Mult->Fill(pionVec.M(),pionVec.Pt(),wght);
                }
              }
            }
          }
          else if(hitclass1 == 2 && hitclass2 == 2 && fAddedSignal){
            if(thisEvent.hit[i].imo == thisEvent.hit[j].imo){
              fHPionInvMassesAdd2->Fill(pionVec.M(),pionVec.Pt(),wght);
            }
          }
          else if((hitclass1 >= 100 && hitclass2 >= 100) && (hitclass1 < 200 && hitclass2 < 200)){
            fHPionInvMasses->Fill(pionVec.M(),pionVec.Pt());
            if(pionVec.M()> 0.11 && pionVec.M() < 0.17)
              fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
            
            // look at the different cases now
            // use only hits from same mother particle
            if(thisEvent.hit[i].imo == thisEvent.hit[j].imo){
              if(hitclass1 == 101 && hitclass2 == 101){
                fHPriPionInvMasses->Fill(pionVec.M(),pionVec.Pt());
              }
              else if(hitclass1 == 102 && hitclass2 == 102){
                fHSecPionInvMasses->Fill(pionVec.M(),pionVec.Pt());
              }
              else if(hitclass1 == 103 && hitclass2 == 103){
                fHK0PionInvMasses->Fill(pionVec.M(),pionVec.Pt());
              }
              else if(hitclass1 == 104 && hitclass2 == 104){
                fHMatPionInvMasses->Fill(pionVec.M(),pionVec.Pt());
              }
            }
          }
          else if(hitclass1 >= 200 && hitclass2 >= 200){
            fHPionInvMassesDCal->Fill(pionVec.M(),pionVec.Pt());
            fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
          }
        }
        else{
          if(hitclass1 < 200 && hitclass2 < 200){
            fHPionInvMasses->Fill(pionVec.M(),pionVec.Pt());
            if(pionVec.M()> 0.11 && pionVec.M() < 0.17)
              fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
          }
          else if(hitclass1 >= 200 && hitclass2 >= 200) {
            fHPionInvMassesDCal->Fill(pionVec.M(),pionVec.Pt());
            if(pionVec.M()> 0.11 && pionVec.M() < 0.17)
              fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
          }
          else{
            fHPionInvMassesEMCalDCal->Fill(pionVec.M(),pionVec.Pt());
            if(pionVec.M()> 0.11 && pionVec.M() < 0.17)
              fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
          }
        }
      }
      
      // asymmetric
      else{
        fHPionMggPt->Fill(pionVec.M(),pionVec.Pt());
        fHPionMggAsym->Fill(pionVec.M(),pionZgg);
        fHPionMggDgg->Fill(pionVec.M(),pionOpeningAngle);
        
        if(pionVec.M()> 0.11 && pionVec.M() < 0.17)
          fHPionPtAsym->Fill(pionVec.Pt(),pionZgg);
        
        if(fMcMode){
          if(hitclass1 == 1 && hitclass2 == 1 && fAddedSignal){
            if(thisEvent.hit[i].imo == thisEvent.hit[j].imo){
              // all, no weight
              fHPionInvMassesAdd1NoWgt->Fill(pionVec.M(),pionVec.Pt());
              
              // clean clusters only
              if(thisEvent.hit[i].bclean == 1 && thisEvent.hit[j].bclean == 1){
                // photons only
                if(thisEvent.hit[i].pid == 22 && thisEvent.hit[j].pid == 22){
                  fHPionInvMassesGamAdd1Asym->Fill(pionVec.M(),pionVec.Pt(),wght);
                }
                // one non-photon
                else{
                  fHPionInvMassesAdd1Asym->Fill(pionVec.M(),pionVec.Pt(),wght);
                }
              }
              // one non-clean cluster
              else{
                // photons only
                if(thisEvent.hit[i].pid == 22 && thisEvent.hit[j].pid == 22){
                  fHPionInvMassesGamAdd1MultAsym->Fill(pionVec.M(),pionVec.Pt(),wght);
                }
                // one non-photon
                else{
                  fHPionInvMassesAdd1MultAsym->Fill(pionVec.M(),pionVec.Pt(),wght);
                }
              }
            }
          }
          else if(hitclass1 == 2 && hitclass2 == 2 && fAddedSignal){
            if(thisEvent.hit[i].imo == thisEvent.hit[j].imo){
              fHPionInvMassesAdd2->Fill(pionVec.M(),pionVec.Pt(),wght);
            }
          }
          else if((hitclass1 >= 100 && hitclass2 >= 100) && (hitclass1 < 200 && hitclass2 < 200)){
            fHPionInvMassesAsym->Fill(pionVec.M(),pionVec.Pt());
            if(pionVec.M()> 0.11 && pionVec.M() < 0.17)
              fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
            
            // look at the different cases now
            // use only hits from same mother particle
            if(thisEvent.hit[i].imo == thisEvent.hit[j].imo){
              if(hitclass1 == 101 && hitclass2 == 101){
                fHPriPionInvMasses->Fill(pionVec.M(),pionVec.Pt());
              }
              else if(hitclass1 == 102 && hitclass2 == 102){
                fHSecPionInvMasses->Fill(pionVec.M(),pionVec.Pt());
              }
              else if(hitclass1 == 103 && hitclass2 == 103){
                fHK0PionInvMasses->Fill(pionVec.M(),pionVec.Pt());
              }
              else if(hitclass1 == 104 && hitclass2 == 104){
                fHMatPionInvMasses->Fill(pionVec.M(),pionVec.Pt());
              }
            }
          }
          else if(hitclass1 >= 200 && hitclass2 >= 200){
            fHPionInvMassesDCal->Fill(pionVec.M(),pionVec.Pt());
            fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
          }
        }
        else{
          if(hitclass1 < 200 && hitclass2 < 200){
            fHPionInvMassesAsym->Fill(pionVec.M(),pionVec.Pt());
            if(pionVec.M()> 0.11 && pionVec.M() < 0.17)
              fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
          }
          else if(hitclass1 >= 200 && hitclass2 >= 200) {
            fHPionInvMassesDCal->Fill(pionVec.M(),pionVec.Pt());
            if(pionVec.M()> 0.11 && pionVec.M() < 0.17)
              fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
          }
          else{
            fHPionInvMassesEMCalDCal->Fill(pionVec.M(),pionVec.Pt());
            if(pionVec.M()> 0.11 && pionVec.M() < 0.17)
              fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
          }
        }
      }
      
      // fill histograms for simulation studies (electrons, charged pions in clusters)
      if(fSimStudies && (hitclass1 >= 100 && hitclass2 >= 100)){
        // electrons
        if(abs(thisEvent.hit[i].pid) == 11 && abs(thisEvent.hit[j].pid) == 11){
          fHPionInvMassesConvElBoth->Fill(pionVec.M(),pionVec.Pt());
        }
        else if((abs(thisEvent.hit[i].pid) == 11 && abs(thisEvent.hit[j].pid) == 22) || (abs(thisEvent.hit[i].pid) == 22 && abs(thisEvent.hit[j].pid) == 11)){
          fHPionInvMassesConvElOne->Fill(pionVec.M(),pionVec.Pt());
        }
        else{
          fHPionInvMassesConvElZero->Fill(pionVec.M(),pionVec.Pt());
        }
        // charged pions
        if(abs(thisEvent.hit[i].pid) == 211 && abs(thisEvent.hit[j].pid) == 211){
          fHPionInvMassesChargedPiBoth->Fill(pionVec.M(),pionVec.Pt());
        }
        else if((abs(thisEvent.hit[i].pid) == 211 && abs(thisEvent.hit[j].pid) == 22) || (abs(thisEvent.hit[i].pid) == 22 && abs(thisEvent.hit[j].pid) == 211)){
          fHPionInvMassesChargedPiOne->Fill(pionVec.M(),pionVec.Pt());
        }
        else{
          fHPionInvMassesChargedPiZero->Fill(pionVec.M(),pionVec.Pt());
        }
        // two photons
        if(abs(thisEvent.hit[i].pid) == 22 && abs(thisEvent.hit[j].pid) == 22){
          fHPionInvMassesGammaBoth->Fill(pionVec.M(),pionVec.Pt());
        }
        

      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::FillMixHists(const Int_t MulClass, const Int_t vtxClass, const Int_t PtClass, Double_t phi0, Double_t theta0)
{
  // Fill histograms related to pions.
  
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  
  TObjArray *clusters = fEsdClusters;
  if (!clusters)
    clusters = fAodClusters;
  
  if (clusters) {
    TLorentzVector clusterVec1;
    TLorentzVector clusterVec2;
    TLorentzVector pionVec;
    Short_t hitclass1 = -1;
    Short_t hitclass2 = -1;
    
    Int_t nclus = thisEvent.nHits;
    
    // loop over clusters in current events
    for (Int_t i = 0; i<nclus; ++i) {
      
      clusterVec1 = thisEvent.hit[i].thishit;
      hitclass1 = thisEvent.hit[i].hittype;
      
      // loop over old events
      for (Int_t iOld=0;iOld<nEvt;iOld++){
        EmcEvent OldEvent = EmcEventList[MulClass][vtxClass][PtClass][iOld];
        Int_t nclusold = OldEvent.nHits;
        Double_t phirot = OldEvent.TrigPhi;
        Double_t thetarot = OldEvent.TrigTheta;
        //				if(fRotateMixed && PtClass > 0)
        //          fHMixRotation->Fill(phi0-phirot);
        
        // loop over old clusters
        for (Int_t j = 0; j<nclusold; ++j) {
          clusterVec2 = OldEvent.hit[j].thishit;
          hitclass2 = OldEvent.hit[j].hittype;
          if(fRotateMixed && PtClass > 0){
            clusterVec2.RotateZ(phi0-phirot);
            clusterVec2.SetTheta(clusterVec2.Theta() + (theta0-thetarot));
          }
          
          // calculate radius
          Double_t d_phi = clusterVec1.Phi() - clusterVec2.Phi();
          Double_t d_eta = clusterVec1.Eta() - clusterVec2.Eta();
          Double_t d_r = sqrt(d_phi*d_phi + d_eta*d_eta);
          
          
          pionVec = clusterVec1 + clusterVec2;
          
          // cut on minimum distance, should be made settable and the same as in "same" events
          if(d_r < 0.01)
            continue;
          
          Double_t pionZgg = TMath::Abs(clusterVec1.E()-clusterVec2.E())/pionVec.E();
          
          if(fCalibRun && pionZgg < 0.2){
            //EMCal histos
            if((hitclass1 >= 100 && hitclass2 >= 100) && (hitclass1 < 200 && hitclass2 < 200)){
              fHPionInvMassesMixEMCalCalib->Fill(pionVec.M(),(clusterVec1.E() + clusterVec2.E())/2);
            }
            else if(hitclass1 >= 200 && hitclass2 >= 200){
              fHPionInvMassesMixDCalCalib->Fill(pionVec.M(),(clusterVec1.E() + clusterVec2.E())/2);
            }
          }
          
          // fill all histograms (maybe split into more asymmetry classes as in same events
          if (pionZgg < fAsymMax3) {
            if(fMcMode){
              if(hitclass1 == 1 && hitclass2 == 1){
                if(thisEvent.hit[i].pid != 22 || OldEvent.hit[j].pid != 22)
                  fHPionInvMassesMix1->Fill(pionVec.M(),pionVec.Pt());
              }
              else if(hitclass1 == 2 && hitclass2 == 2){
                if(thisEvent.hit[i].pid != 22 || OldEvent.hit[j].pid != 22)
                  fHPionInvMassesMix2->Fill(pionVec.M(),pionVec.Pt());
              }
              else if(hitclass1 >= 100 && hitclass2 >= 100 && hitclass1 < 200 && hitclass2 < 200){
                fHPionInvMassesMix->Fill(pionVec.M(),pionVec.Pt());
              }
              else{
                fHPionInvMassesMixDCal->Fill(pionVec.M(),pionVec.Pt());
              }
            }
            else{ // not MC mode
              if(hitclass1 < 200 && hitclass2 < 200){
                fHPionInvMassesMix->Fill(pionVec.M(),pionVec.Pt());
              }
              else if(hitclass1 >= 200 && hitclass2 >= 200) {
                fHPionInvMassesMixDCal->Fill(pionVec.M(),pionVec.Pt());
              }
              else{
                fHPionInvMassesMixEMCalDCal->Fill(pionVec.M(),pionVec.Pt());
              }
            }
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::FillMcHists()
{
  // Fill additional MC information histograms.
  
  if (!fMcParts)
    return;
  
  // check if aod or esd mc mode and the fill histos
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::FillOtherHists()
{
  // Fill other histograms.
  
  for(int iTracks=0; iTracks < fSelPrimTracks->GetEntries(); ++iTracks){
    AliESDtrack *track = static_cast<AliESDtrack*>(fSelPrimTracks->At(iTracks));
    if(!track)
      continue;
    fHPrimTrackPt->Fill(track->Pt());
    fHPrimTrackEta->Fill(track->Eta());
    fHPrimTrackPhi->Fill(track->Phi());
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::FillTrackHists()
{
  // Fill track histograms.
  
  if (fSelPrimTracks) {
    for(int iTracks=0; iTracks < fSelPrimTracks->GetEntries(); iTracks++) {
      AliESDtrack *track = static_cast<AliESDtrack*>(fSelPrimTracks->At(iTracks));
      if(!track)
        continue;
      fHPrimTrackPt->Fill(track->Pt());
      fHPrimTrackEta->Fill(track->Eta());
      fHPrimTrackPhi->Fill(track->Phi());
    }
  }
}

//__________________________________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::FillVertex(AliStaVertex *v, const AliESDVertex *esdv)
{
  // Fill vertex from ESD vertex info.
  
  v->fVx   = esdv->GetX();
  v->fVy   = esdv->GetY();
  v->fVz   = esdv->GetZ();
  v->fVc   = esdv->GetNContributors();
  v->fDisp = esdv->GetDispersion();
  v->fZres = esdv->GetZRes();
  v->fChi2 = esdv->GetChi2();
  v->fSt   = esdv->GetStatus();
  v->fIs3D = esdv->IsFromVertexer3D();
  v->fIsZ  = esdv->IsFromVertexerZ();
}

//__________________________________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::FillVertex(AliStaVertex *v, const AliAODVertex *aodv)
{
  // Fill vertex from AOD vertex info.
  
  v->fVx   = aodv->GetX();
  v->fVy   = aodv->GetY();
  v->fVz   = aodv->GetZ();
  v->fVc   = aodv->GetNContributors();
  v->fChi2 = aodv->GetChi2();
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPi0Gamma::GetCellIsolation(Double_t cEta, Double_t cPhi, Double_t radius) const
{
  // Compute isolation based on cell content.
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  if (!cells)
    return 0;
  
  Double_t cellIsolation = 0;
  Double_t rad2 = radius*radius;
  Int_t ncells = cells->GetNumberOfCells();
  for (Int_t i = 0; i<ncells; ++i) {
    Int_t absID    = TMath::Abs(cells->GetCellNumber(i));
    Float_t eta=-1, phi=-1;
    fGeom->EtaPhiFromIndex(absID,eta,phi);
    Double_t phidiff = TVector2::Phi_mpi_pi(phi-cPhi);
    Double_t dist = (eta-cEta)*(eta-cEta)+phidiff*phidiff;
    if(dist>rad2)
      continue;
    Double_t cellE = cells->GetAmplitude(i);
    Double_t theta = 2*TMath::ATan(TMath::Exp(-eta));
    Double_t cellEt = cellE*sin(theta);
    cellIsolation += cellEt;
  }
  return cellIsolation;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPi0Gamma::GetCellIsoNxM(Double_t cEta, Double_t cPhi, Int_t N, Int_t M) const
{
  // Compute isolation based on cell content, in a NxM rectangle.
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  if (!cells)
    return 0;
  
  Double_t cellIsolation = 0;
  Int_t ncells = cells->GetNumberOfCells();
  for (Int_t i = 0; i<ncells; ++i) {
    Int_t absID    = TMath::Abs(cells->GetCellNumber(i));
    Float_t eta=-1, phi=-1;
    fGeom->EtaPhiFromIndex(absID,eta,phi);
    Double_t phidiff = TVector2::Phi_mpi_pi(phi-cPhi);
    Double_t etadiff = (eta-cEta)*(eta-cEta)+phidiff*phidiff;
    if(TMath::Abs(etadiff)/0.014>N)
      continue;
    if(TMath::Abs(phidiff)/0.014>M)
      continue;
    Double_t cellE = cells->GetAmplitude(i);
    Double_t theta = 2*TMath::ATan(TMath::Exp(-eta));
    Double_t cellEt = cellE*sin(theta);
    cellIsolation += cellEt;
  }
  return cellIsolation;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPi0Gamma::GetCellEnergy(const AliVCluster *cluster) const
{
  // Get maximum energy of attached cell.
  
  Double_t ret = 0;
  Int_t ncells = cluster->GetNCells();
  if (fEsdCells) {
    for (Int_t i=0; i<ncells; i++) {
      Double_t e = fEsdCells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
      ret += e;
    }
  } else {
    for (Int_t i=0; i<ncells; i++) {
      Double_t e = fAodCells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
      ret += e;
    }
  }
  return ret;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPi0Gamma::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
{
  // Get maximum energy of attached cell.
  
  id = -1;
  Double_t maxe = 0;
  Int_t ncells = cluster->GetNCells();
  if (fEsdCells) {
    for (Int_t i=0; i<ncells; i++) {
      Double_t e = fEsdCells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
      if (e>maxe) {
        maxe = e;
        id   = cluster->GetCellAbsId(i);
      }
    }
  } else {
    for (Int_t i=0; i<ncells; i++) {
      Double_t e = fAodCells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
      if (e>maxe)
        maxe = e;
      id   = cluster->GetCellAbsId(i);
    }
  }
  return maxe;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPi0Gamma::GetSecondMaxCellEnergy(AliVCluster *clus, Short_t &id) const
{
  // Get second maximum cell.
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  if (!cells)
    return -1;
  
  Double_t secondEmax=0, firstEmax=0;
  Double_t cellen;
  for(Int_t iCell=0;iCell<clus->GetNCells();iCell++){
    Int_t absId = clus->GetCellAbsId(iCell);
    cellen = cells->GetCellAmplitude(absId);
    if(cellen > firstEmax)
      firstEmax = cellen;
  }
  for(Int_t iCell=0;iCell<clus->GetNCells();iCell++){
    Int_t absId = clus->GetCellAbsId(iCell);
    cellen = cells->GetCellAmplitude(absId);
    if(cellen < firstEmax && cellen > secondEmax) {
      secondEmax = cellen;
      id = absId;
    }
  }
  return secondEmax;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::GetSigma(const AliVCluster *c, Double_t& sigmaMax, Double_t &sigmaMin) const
{
  // Calculate the (E) weighted variance along the longer (eigen) axis.
  
  sigmaMax = 0;          // cluster variance along its longer axis
  sigmaMin = 0;          // cluster variance along its shorter axis
  Double_t Ec  = c->E(); // cluster energy
  if(Ec<=0)
    return;
  Double_t Xc  = 0 ;     // cluster first moment along X
  Double_t Yc  = 0 ;     // cluster first moment along Y
  Double_t Sxx = 0 ;     // cluster second central moment along X (variance_X^2)
  Double_t Sxy = 0 ;     // cluster second central moment along Y (variance_Y^2)
  Double_t Syy = 0 ;     // cluster covariance^2
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  
  if (!cells)
    return;
  
  Int_t ncells = c->GetNCells();
  if (ncells==1)
    return;
  
  for(Int_t j=0; j<ncells; ++j) {
    Int_t id = TMath::Abs(c->GetCellAbsId(j));
    Double_t cellen = cells->GetCellAmplitude(id);
    TVector3 pos;
    fGeom->GetGlobal(id,pos);
    Xc  += cellen*pos.X();
    Yc  += cellen*pos.Y();
    Sxx += cellen*pos.X()*pos.X();
    Syy += cellen*pos.Y()*pos.Y();
    Sxy += cellen*pos.X()*pos.Y();
  }
  Xc  /= Ec;
  Yc  /= Ec;
  Sxx /= Ec;
  Syy /= Ec;
  Sxy /= Ec;
  Sxx -= Xc*Xc;
  Syy -= Yc*Yc;
  Sxy -= Xc*Yc;
  Sxx = TMath::Abs(Sxx);
  Syy = TMath::Abs(Syy);
  sigmaMax = (Sxx + Syy + TMath::Sqrt(TMath::Abs((Sxx-Syy)*(Sxx-Syy)+4.0*Sxy*Sxy)))/2.0;
  sigmaMax = TMath::Sqrt(TMath::Abs(sigmaMax));
  sigmaMin = TMath::Abs(Sxx + Syy - TMath::Sqrt(TMath::Abs((Sxx-Syy)*(Sxx-Syy)+4.0*Sxy*Sxy)))/2.0;
  sigmaMin = TMath::Sqrt(TMath::Abs(sigmaMin));
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::GetSigmaEtaEta(const AliVCluster *c, Double_t& sEtaEta, Double_t &sPhiPhi) const
{
  // Calculate the (E) weighted variance along the pseudorapidity.
  
  sEtaEta = 0;
  sPhiPhi = 0;
  
  Double_t Ec  = c->E(); // cluster energy
  if(Ec<=0)
    return;
  
  const Int_t ncells = c->GetNCells();
  
  Double_t EtaC    = 0;  // cluster first moment along eta
  Double_t PhiC    = 0;  // cluster first moment along phi
  Double_t Setaeta = 0;  // cluster second central moment along eta
  Double_t Sphiphi = 0;  // cluster second central moment along phi
  Double_t w[ncells];    // weight max(0,4.5*log(E_i/Ec))
  Double_t sumw = 0;
  Int_t id[ncells];
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  
  if (!cells)
    return;
  
  if (ncells==1)
    return;
  
  for(Int_t j=0; j<ncells; ++j) {
    id[j] = TMath::Abs(c->GetCellAbsId(j));
    Double_t cellen = cells->GetCellAmplitude(id[j]);
    w[j] = TMath::Max(0., 4.5+TMath::Log(cellen/Ec));
    TVector3 pos;
    fGeom->GetGlobal(id[j],pos);
    EtaC += w[j]*pos.Eta();
    PhiC += w[j]*pos.Phi();
    sumw += w[j];
  }
  EtaC /= sumw;
  PhiC /= sumw;
  
  for(Int_t j=0; j<ncells; ++j) {
    TVector3 pos;
    fGeom->GetGlobal(id[j],pos);
    Setaeta =  w[j]*(pos.Eta() - EtaC)*(pos.Eta() - EtaC);
    Sphiphi =  w[j]*(pos.Phi() - PhiC)*(pos.Phi() - PhiC);
  }
  Setaeta /= sumw;
  sEtaEta = TMath::Sqrt(Setaeta);
  Sphiphi /= sumw;
  sPhiPhi = TMath::Sqrt(Sphiphi);
}

//________________________________________________________________________
Int_t AliAnalysisTaskEMCALPi0Gamma::GetNCells(const AliVCluster *c, Double_t emin) const
{
  // Calculate number of attached cells above emin.
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  if (!cells)
    return 0;
  
  Int_t n = 0;
  Int_t ncells = c->GetNCells();
  for(Int_t j=0; j<ncells; ++j) {
    Int_t id = TMath::Abs(c->GetCellAbsId(j));
    Double_t cellen = cells->GetCellAmplitude(id);
    if (cellen>=emin)
      ++n;
  }
  return n;
}

//________________________________________________________________________
Int_t AliAnalysisTaskEMCALPi0Gamma::GetNCells(Int_t sm, Double_t emin) const
{
  // Calculate number of cells per SM above emin.
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  if (!cells)
    return 0;
  
  Int_t n = 0;
  Int_t ncells = cells->GetNumberOfCells();
  for(Int_t j=0; j<ncells; ++j) {
    Int_t id = TMath::Abs(cells->GetCellNumber(j));
    Double_t cellen = cells->GetCellAmplitude(id);
    if (cellen<emin)
      continue;
    Int_t fsm = fGeom->GetSuperModuleNumber(id);
    if (fsm != sm)
      continue;
    ++n;
  }
  return n;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPi0Gamma::GetTrackIsolation(Double_t cEta, Double_t cPhi, Double_t radius, Double_t pt) const
{
  // Compute isolation based on tracks.
  
  Double_t trkIsolation = 0;
  Double_t rad2 = radius*radius;
  Int_t ntrks = fSelPrimTracks->GetEntries();
  for(Int_t j = 0; j<ntrks; ++j) {
    AliVTrack *track = static_cast<AliVTrack*>(fSelPrimTracks->At(j));
    if (!track)
      continue;
    if (track->Pt()<pt)
      continue;
    Float_t eta = track->Eta();
    Float_t phi = track->Phi();
    Double_t phidiff = TVector2::Phi_mpi_pi(phi-cPhi);
    Double_t dist = (eta-cEta)*(eta-cEta)+phidiff*phidiff;
    if(dist>rad2)
      continue;
    trkIsolation += track->Pt();
  }
  return trkIsolation;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPi0Gamma::GetTrackIsoStrip(Double_t cEta, Double_t cPhi, Double_t dEta, Double_t dPhi, Double_t pt) const
{
  // Compute isolation based on tracks.
  
  Double_t trkIsolation = 0;
  Int_t ntrks = fSelPrimTracks->GetEntries();
  for(Int_t j = 0; j<ntrks; ++j) {
    AliVTrack *track = static_cast<AliVTrack*>(fSelPrimTracks->At(j));
    if (!track)
      continue;
    if (track->Pt()<pt)
      continue;
    Float_t eta = track->Eta();
    Float_t phi = track->Phi();
    Double_t phidiff = TVector2::Phi_mpi_pi(phi-cPhi);
    Double_t etadiff = (eta-cEta);
    if(TMath::Abs(etadiff)>dEta || TMath::Abs(phidiff)>dPhi)
      continue;
    trkIsolation += track->Pt();
  }
  return trkIsolation;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPi0Gamma::IsShared(const AliVCluster *c) const
{
  // Returns if cluster shared across super modules.
  
  if (!c)
    return 0;
  
  Int_t n = -1;
  Int_t ncells = c->GetNCells();
  for(Int_t j=0; j<ncells; ++j) {
    Int_t id = TMath::Abs(c->GetCellAbsId(j));
    Int_t got = id / (24*48);
    if (n==-1) {
      n = got;
      continue;
    }
    if (got!=n)
      return 1;
  }
  return 0;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPi0Gamma::IsIdPartOfCluster(const AliVCluster *c, Short_t id) const
{
  // Returns if id is part of cluster.
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  if (!cells)
    return 0;
  
  Int_t ncells = c->GetNCells();
  for(Int_t j=0; j<ncells; ++j) {
    Int_t cid = TMath::Abs(c->GetCellAbsId(j));
    if (cid == id)
      return 1;
  }
  return 0;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::PrintDaughters(const AliVParticle *p, const TObjArray *arr, Int_t level) const
{
  // Print recursively daughter information.
  
  if (!p || !arr)
    return;
  
  const AliAODMCParticle *amc = dynamic_cast<const AliAODMCParticle*>(p);
  if (!amc)
    return;
  for (Int_t i=0; i<level; ++i) printf(" ");
  amc->Print();
  
  Int_t n = amc->GetNDaughters();
  for (Int_t i=0; i<n; ++i) {
    Int_t d = amc->GetDaughterLabel(i);
    const AliVParticle *dmc = static_cast<const AliVParticle*>(arr->At(d));
    PrintDaughters(dmc,arr,level+1);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::PrintDaughters(const AliMCParticle *p, const AliMCEvent *arr, Int_t level) const
{
  // Print recursively daughter information.
  
  if (!p || !arr)
    return;
  
  for (Int_t i=0; i<level; ++i) printf(" ");
  Int_t d1 = p->GetDaughterFirst();
  Int_t d2 = p->GetDaughterLast();
  printf("pid=%d: %.2f %.2f %.2f (%.2f %.2f %.2f); nd=%d,%d\n",
         p->PdgCode(),p->Px(),p->Py(),p->Pz(),p->Xv(),p->Yv(),p->Zv(),d1,d2);
  if (d1<0)
    return;
  if (d2<0)
    d2=d1;
  for (Int_t i=d1;i<=d2;++i) {
    const AliMCParticle *dmc = static_cast<const AliMCParticle *>(arr->GetTrack(i));
    PrintDaughters(dmc,arr,level+1);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::ProcessDaughters(AliVParticle *p, Int_t index, const TObjArray *arr)
{
  // Process and create daughters.
  
  if (!p || !arr)
    return;
  
  AliAODMCParticle *amc = dynamic_cast<AliAODMCParticle*>(p);
  if (!amc)
    return;
  
  //amc->Print();
  
  Int_t nparts = arr->GetEntries();
  Int_t nents  = fMcParts->GetEntries();
  
  AliStaPart *newp = static_cast<AliStaPart*>(fMcParts->New(nents));
  newp->fPt  = amc->Pt();
  newp->fEta = amc->Eta();
  newp->fPhi = amc->Phi();
  if (amc->Xv() != 0 || amc->Yv() != 0 || amc->Zv() != 0) {
    TVector3 vec(amc->Xv(),amc->Yv(),amc->Zv());
    newp->fVR = vec.Perp();
    newp->fVEta = vec.Eta();
    newp->fVPhi = vec.Phi();
  }
  newp->fPid  = amc->PdgCode();
  newp->fLab  = nents;
  Int_t moi = amc->GetMother();
  if (moi>=0&&moi<nparts) {
    const AliAODMCParticle *mmc = static_cast<const AliAODMCParticle*>(arr->At(moi));
    moi = mmc->GetUniqueID();
  }
  newp->fMo = moi;
  p->SetUniqueID(nents);
  
  // TODO: Determine which detector was hit
  //newp->fDet = ???
  
  Int_t n = amc->GetNDaughters();
  for (Int_t i=0; i<n; ++i) {
    Int_t d = amc->GetDaughterLabel(i);
    if (d<=index || d>=nparts)
      continue;
    AliVParticle *dmc = static_cast<AliVParticle*>(arr->At(d));
    ProcessDaughters(dmc,d,arr);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::ProcessDaughters(AliMCParticle *p, Int_t index, const AliMCEvent *arr)
{
  // Process and create daughters.
  
  if (!p || !arr)
    return;
  
  Int_t d1 = p->GetDaughterFirst();
  Int_t d2 = p->GetDaughterLast();
  if (0) {
    printf("%d pid=%d: %.3f %.3f %.3f (%.2f %.2f %.2f); nd=%d,%d, mo=%d\n",
           index,p->PdgCode(),p->Px(),p->Py(),p->Pz(),p->Xv(),p->Yv(),p->Zv(),d1,d2, p->GetMother());
  }
  Int_t nents  = fMcParts->GetEntries();
  
  AliStaPart *newp = static_cast<AliStaPart*>(fMcParts->New(nents));
  newp->fPt  = p->Pt();
  newp->fEta = p->Eta();
  newp->fPhi = p->Phi();
  if (p->Xv() != 0 || p->Yv() != 0 || p->Zv() != 0) {
    TVector3 vec(p->Xv(),p->Yv(),p->Zv());
    newp->fVR = vec.Perp();
    newp->fVEta = vec.Eta();
    newp->fVPhi = vec.Phi();
  }
  newp->fPid  = p->PdgCode();
  newp->fLab  = nents;
  Int_t moi = p->GetMother();
  if (moi>=0) {
    const AliMCParticle *mmc = static_cast<const AliMCParticle *>(arr->GetTrack(moi));
    moi = mmc->GetUniqueID();
  }
  newp->fMo = moi;
  p->SetUniqueID(nents);
  
  Int_t nref = p->GetNumberOfTrackReferences();
  if (nref>0) {
    AliTrackReference *ref = p->GetTrackReference(nref-1);
    if (ref) {
      newp->fDet = ref->DetectorId();
    }
  }
  
  if (d1<0)
    return;
  if (d2<0)
    d2=d1;
  for (Int_t i=d1;i<=d2;++i) {
    AliMCParticle *dmc = static_cast<AliMCParticle *>(arr->GetTrack(i));
    if (dmc->P()<0.01)
      continue;
    ProcessDaughters(dmc,i,arr);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::GetMulClass(Int_t& imcl)
{
  Int_t nclus = 0;
  TObjArray *clusters = fEsdClusters;
  if (!clusters)
    clusters = fAodClusters;
  
  if (clusters)
    nclus = clusters->GetEntries();
  
  //const int MultCut[8] = {5, 15, 30, 50, 80, 120, 300, 9999};
  const int MultCut[nMulClass] = {5, 12, 20 ,50, 9999};
  
  imcl=0;
  
  for (imcl=0; imcl<nMulClass; imcl++) {
    if (nclus < MultCut[imcl]) break;
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::AddMixEvent(const Int_t MulClass, const Int_t vtxClass, const Int_t PtClass, Int_t& iEvent, const Float_t& phitrig, const Float_t& thetatrig)
{
  if(iEvent >= nEvt){
    iEvent = 0;
  }
  
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  
  TObjArray *clusters = fEsdClusters;
  if (!clusters)
    clusters = fAodClusters;
  
  Int_t nclus = 0;
  nclus = clusters->GetEntries();
  nclus = thisEvent.nHits;
  
  if(nclus > evt.nMaxHit){
    AliWarning("event has more clusters than nMaxHit!");
    //    nclus = evt.nMaxHit;
  }
  
  
  //cout << Form("%d, %d, %d, %d",MulClass,vtxClass,PtClass,iEvent) << endl;
  thisEvent.SetGlobalInfo(nclus,phitrig,thetatrig);
  EmcEventList[MulClass][vtxClass][PtClass][iEvent] = thisEvent;
  
  iEvent++;
  return;
}

//_____________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::FillCellQAHists(AliVCluster* virCluster, Bool_t isDcal, Bool_t isAfter)
{
  
  UShort_t* CellsID = virCluster->GetCellsAbsId();
  
  // Find Max Contributing Cell
  AliVCaloCells *vcells = fEsdCells;
  if (!vcells)
    vcells = fAodCells;
  
  Float_t cellEnergy = 0.0;
  Float_t maxEnergy = 0.0;
  Int_t maxID = -1;
  
  for(Int_t kk=0; kk < virCluster->GetNCells(); kk++) {
    cellEnergy = vcells->GetCellAmplitude(CellsID[kk]);
    if(cellEnergy > maxEnergy) {
      maxEnergy = cellEnergy;
      maxID = CellsID[kk];
    }
  }
  
  // Cell ID vs. energy for cells that are at the position of a cluster
  //  Int_t cellClusterPosAbsID;
  //fPHOSGeo->RelToAbsNumbering(relId, cellClusterPosAbsID); //Get AbsID of Cell in which the ClusterPosition is
  Double_t cellClusterPosE = vcells->GetCellAmplitude(maxID); //Get Energy of that cell
  
  if(!isAfter)
    fHCellIndexEnergy->Fill(maxID, maxEnergy);
  else
    fHCellIndexEnergyAfterCuts->Fill(maxID,maxEnergy);
  
  return;
}

//_____________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::SetDnDpT(Int_t i, Double_t par0, Double_t par1, Double_t par2, Double_t par3, Double_t par4)
{
  /*
   // case 1: modified hagedorn
   if(i==1){
   fPi0DnDpt = new TF1("fPi0DnDpt","[0]*pow([1]/(([1]*exp(-[3]*x)+x)),[2])",0.6,12);
   fPi0DnDpt->SetParameters(par0,par1,par2,par3);
   
   par4 += 0;
   }
   
   // get rid of the warnings
   else{
   par0+=0;
   par1+=0;
   par2+=0;
   par3+=0;
   par4+=0;
   }
   */
  return;
}

//_____________________________________________________________________
Float_t AliAnalysisTaskEMCALPi0Gamma::CalcWeight(Double_t pt, Double_t eta, Int_t i){
  Float_t weight = 1.;
  if(i==1){
    Double_t par0 = 2.52684e08;
    Double_t par1 = 0.730803;
    Double_t par2 = 5.32059;
    Double_t par3 = 0.548711;
    
    weight = (par0*pow(par1/((par1*exp(-par3*pt)+pt)),par2)) * (7.93979e-01 + eta*1.31404e-03 + eta*eta*4.96607e-01);
  }
  
  return weight;
}

//_____________________________________________________________________
Int_t AliAnalysisTaskEMCALPi0Gamma::GetModuleNumber(AliVCluster * cluster) const
{
  //Get the EMCAL/PHOS module number that corresponds to this cluster
  TLorentzVector lv;
  Double_t v[]={0.,0.,0.}; //not necessary to pass the real vertex.
  if(!cluster)
  {
    if(fDebug > 1) printf("AliCalorimeterUtils::GetModuleNumber() - NUL Cluster, please check!!!");
    return -1;
  }
  
  cluster->GetMomentum(lv,v);
  Float_t phi = lv.Phi();
  if(phi < 0) phi+=TMath::TwoPi();
  Int_t absId = -1;
  if(cluster->IsEMCAL()){
    fGeom->GetAbsCellIdFromEtaPhi(lv.Eta(),phi, absId);
    if(fDebug > 2)
      printf("AliCalorimeterUtils::GetModuleNumber() - EMCAL: cluster eta %f, phi %f, absid %d, SuperModule %d\n",
             lv.Eta(), phi*TMath::RadToDeg(),absId, fGeom->GetSuperModuleNumber(absId));
    return fGeom->GetSuperModuleNumber(absId) ;
  }//EMCAL
  else if(cluster->IsPHOS())
  {
    //    Int_t    relId[4];
    //    if ( cluster->GetNCells() > 0)
    //    {
    //      absId = cluster->GetCellAbsId(0);
    //      if(fDebug > 2)
    //        printf("AliCalorimeterUtils::GetModuleNumber() - PHOS: cluster eta %f, phi %f, e %f, absId %d\n",
    //               lv.Eta(), phi*TMath::RadToDeg(), lv.E(), absId);
    //    }
    //    else return -1;
    //
    //    if ( absId >= 0)
    //    {
    //      fPHOSGeo->AbsToRelNumbering(absId,relId);
    //      if(fDebug > 2)
    //        printf("AliCalorimeterUtils::GetModuleNumber() - PHOS: Module %d\n",relId[0]-1);
    //      return relId[0]-1;
    //    }
    //    else return -1;
    return -1;
  }//PHOS
  
  return -1;
}

// Jason's energy recalibration, use calib factor 8 for data and 9 for MC
//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPi0Gamma::PrivateEnergyRecal(Double_t energy, Int_t iCalib){
  
  double recalibfactor = 0.0;
  
  if(iCalib==0){// no recalibration!
    recalibfactor = 1.0;
  }
  else if(iCalib==1){// just a scale factor:
    recalibfactor = 0.984;
  }
  else if(iCalib==2){// Symmetric Decay Fit - corrects data to uncorrected MC.
    Double_t p[3] = {0.96968, -2.68720, -0.831607};
    recalibfactor = p[0] + exp(p[1] + p[2]*energy*2.0);
  }
  else if(iCalib==3){// Jason's fit to the LHC12f1a MC single photons - 04 Aug 2013 (call it kPi0MCv4??)
    Double_t p[7] = {1.00000e+00, 3.04925e-02, 4.69043e+00, 9.67998e-02, 2.19381e+02, 6.31604e+01, 1.00046e+00};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==4){// Jason's fit to the test beam data - 04 Aug 2013(call it kBTCv3??)
    Double_t p[7] = {9.78672e-01, 2.39745e-01, 6.41199e-01, 9.13538e-02, 1.46058e+02, 1.99469e+01, 9.72716e-01};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==5){// Based on kSDM/kTBCv3 (call it kPi0MCv4??)
    Double_t p[10] = {9.78672e-01, 2.39745e-01, 6.41199e-01, 9.13538e-02, 1.46058e+02, 1.99469e+01, 9.72716e-01, 0.96968, -2.68720, -0.831607};
    recalibfactor = ( (p[6]/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5]))))) ) / ( p[7] + exp(p[8] + p[9]*energy/2.0) );
  }
  else if(iCalib==6){// kBeamTestCorrectedv2 - in AliROOT!
    Double_t p[7] = {9.83504e-01, 2.10106e-01, 8.97274e-01, 8.29064e-02, 1.52299e+02, 3.15028e+01, 0.968};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==7){// kPi0MCv3 - in AliROOT!
    Double_t p[7] = {9.81039e-01, 1.13508e-01, 1.00173e+00, 9.67998e-02, 2.19381e+02, 6.31604e+01, 1.0};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==8){// Jason's fit to the noNL MC/data- based on kSDM and kPi0MCv5 - 28 Oct 2013 (call it... ??)
    Double_t p[10] = {1.0, 6.64778e-02, 1.57000e+00, 9.67998e-02, 2.19381e+02, 6.31604e+01, 1.01286, 0.964, -3.132, -0.435};
    //Double_t p[10] = {1.0, 6.64778e-02, 1.57000e+00, 9.67998e-02, 2.19381e+02, 6.31604e+01, 1.01286, 0.96968, -2.68720, -0.831607};//same SDM piece as iCalib==2
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5]))))) * (p[7] + exp(p[8]+p[9]*energy*2.0));
  }
  else if(iCalib==9){// Jason's fit to the LHC12f1a/b MC single photons (above 400MeV), including conversions - 28 Oct 2013 (call it kPi0MCv5??)
    Double_t p[7] = {1.0, 6.64778e-02, 1.57000e+00, 9.67998e-02, 2.19381e+02, 6.31604e+01, 1.01286};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==10){// Jason played with test beam data
    Double_t p[7] = {1.0, 0.237767, 0.651203, 0.183741, 155.427, 17.0335, 0.987054};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==11){// Jason played with test beam MC
    Double_t p[7] = {1.0, 0.0797873, 1.68322, 0.0806098, 244.586, 116.938, 1.00437};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  
  return recalibfactor;
}


//__________________________________________________________________________________________________
EmcEvent::EmcEvent()
: nHits(0),
TrigPhi(0),
TrigTheta(0)
{
  nHits = 0;
  TrigPhi = 0;
  TrigTheta = 0;
}

//__________________________________________________________________________________________________
EmcEvent::EmcEvent(const EmcEvent &obj)
: nHits(0),
TrigPhi(0),
TrigTheta(0)
{
  nHits = obj.nHits;
  TrigPhi = obj.TrigPhi;
  TrigTheta = obj.TrigTheta;
  // copy all hits
  for(int i=0;i<nHits;i++){
    hit[i].hittype = obj.hit[i].hittype;
    hit[i].imo = obj.hit[i].imo;
    hit[i].pid = obj.hit[i].pid;
    hit[i].weight = obj.hit[i].weight;
    hit[i].thishit = obj.hit[i].thishit;
    hit[i].smno = obj.hit[i].smno;
  }
}


//__________________________________________________________________________________________________
void EmcEvent::SetGlobalInfo(const Int_t& Size, const Float_t& phiTrig, const Float_t& thetaTrig)
{
  //    fCenPercent = centPer;
  //    fVtx = vtxPos;
  nHits = Size;
  TrigPhi = phiTrig;
  TrigTheta = thetaTrig;
}

void EmcEvent::Print()
{
  Printf("%d hits",nHits);
  for(int i=0;i<nHits;i++){
    hit[i].thishit.Print();
  }
}

//__________________________________________________________________________________________________
void EmcEvent::Reset()
{
  for(int i=0;i<nHits;i++){
    hit[i].hittype = 0;
    hit[i].imo = 0;
    hit[i].pid = 0;
    hit[i].weight = 1;
    hit[i].smno = 0;
    TLorentzVector lv(0,0,0,0);
    hit[i].thishit = lv;
  }
  nHits = 0;
  //    fCenPercent = -1;
  //    fVtx = -9999;
}

//__________________________________________________________________________________________________
EmcHit::EmcHit()
: thishit(),
hittype(0),
imo(0),
pid(0),
weight(1),
bclean(1)
{
}  


