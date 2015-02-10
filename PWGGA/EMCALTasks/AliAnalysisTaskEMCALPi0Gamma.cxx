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
#include <TGeoGlobalMagField.h>
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
//#include "AliEMCALRecPoint.h"
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
#include "AliMagF.h"
#include "AliMultiplicity.h"
#include "AliStack.h"
#include "AliTrackerBase.h"
#include "AliTriggerAnalysis.h"
//#include "AliConversionCuts.h"
//#include "AliV0ReaderV1.h"
//#include "AliMCAnalysisUtils.h"
//#include "AliCalorimeterUtils.h"
#include "../GammaConv/AliAODConversionPhoton.h"
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
fDoConvAna(0),
fDoAfterburner(0),
fAsymMax1(0.3),
fAsymMax2(0.7),
fAsymMax3(1),
fNminCells(1),
fMinE(0.100),
fM02(100),
fMinErat(0),
fMinEcc(0),
fDoManualRecal(0),
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
fHTclsBeforeCuts(0x0),
fHTclsAfterCuts(0x0),
fHClustNoEvt(0),
fHClustAccEvt(0),
fHClustEccentricity(0),
fHClustEtaPhi(0x0),
fHClustEnergyPt(0x0),
fHClustEnergySigma(0x0),
fHClustSigmaSigma(0x0),
fHClustNCellEnergyRatio(0x0),
fHClustEnergyNCell(0x0),
fHConvEnergyPt(0x0),
fHPrimTrackPt(0x0),
fHPrimTrackEta(0x0),
fHPrimTrackPhi(0x0),
fHMatchDr(0x0),
fHMatchDz(0x0),
fHMatchEp(0x0),
fHPionEtaPhi(0x0),
fHAddPionEtaPhi(0x0),
fHPyPionEtaPhi(0x0),
fHdr(0),
fHPionMggPt(0x0),
fHPionMggAsym(0x0),
fHPionMggDgg(0x0),
fHPionInvMasses(0x0),
fHPionInvMassesSym(0x0),
fHPionInvMassesAsym(0x0),
fHPionInvMassesAdd1(0x0),
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
fHPionInvMassesMix(0x0),
fHPionInvMassesMix1(0x0),
fHPionInvMassesMix2(0x0),
fHJPInvMasses(0x0),
fHPrimPionInvMasses(0x0),
fHPrimPionInvMassesAsym(0x0),
fHPionEtaPhiConv(0x0),
fHPionMggPtConv(0x0),
fHPionMggAsymConv(0x0),
fHPionMggDggConv(0x0),
fHPionInvMassesConv(0x0),
fHPionInvMassesConvMix(0x0),
fHPionEtaPhiConvConv(0x0),
fHPionMggPtConvConv(0x0),
fHPionMggAsymConvConv(0x0),
fHPionMggDggConvConv(0x0),
fHPionInvMassesConvConv(0x0),
fHPionInvMassesConvConvMix(0x0),
fHConversionPoint(0),
fHPionTruthPt(),
fHPionTruthPtIn(),
fHPionTruthPtAcc(),
fHPionTruthPtConvAcc(),
fHEtaTruthPt(),
fHEtaTruthPtIn(),
fHEtaTruthPtAcc(),
fHEtaTruthPtConvAcc(),
fHGamTruthPt(),
fHGamTruthPtIn(),
fHGamTruthPtAcc(),
fHPionTruthPtAdd(),
fHPionTruthPtInAdd(),
fHPionTruthPtAccAdd(),
fHPionTruthPtConvAccAdd(),
fHEtaTruthPtAdd(),
fHEtaTruthPtInAdd(),
fHEtaTruthPtAccAdd(),
fHEtaTruthPtConvAccAdd(),
fHNMothers(0x0),
fHMixRotation(),
ipymin(0),
ipymax(0),
ipi0min(0),
ipi0max(0),
ietamin(0),
ietamax(0),
  fReaderGammas(0),
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
  fHRecTrue(),
  fHRecTrueAddPi0(),
  fHRecTrueAddEta(),
  fHECluEMCnofull(),
  fHECluEMCnofullAdd(),
  fHECluEMCelectron(),
  fHECluEMCpion(),
  fHECluEMCkaon(),
  fHECluEMCother(),
  fHECluEMCpi0single(),
  fHCorrection(),
  fHPionSm(),
  evt(),
  thisEvent(),
  fHWgt(0)
{
  for(int i=0;i<nMulClass;i++){
    for(int j=0;j<nZClass;j++){
      for(int k=0;k<nPtClass;k++){
        iEvt[i][j][k] = 0;
        for(int l=0;l<nEvt;l++){
          EmcEventList[i][j][k][l].SetGlobalInfo(0,0,0.,0.);
        }
      }
    }
  }
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
fDoConvAna(0),
fDoAfterburner(0),
fAsymMax1(0.3),
fAsymMax2(0.7),
fAsymMax3(1),
fNminCells(1),
fMinE(0.100),
fM02(100),
fMinErat(0),
fMinEcc(0),
fDoManualRecal(0),
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
fHTclsBeforeCuts(0x0),
fHTclsAfterCuts(0x0),
fHClustNoEvt(0),
fHClustAccEvt(0),
fHClustEccentricity(0),
fHClustEtaPhi(0x0),
fHClustEnergyPt(0x0),
fHClustEnergySigma(0x0),
fHClustSigmaSigma(0x0),
fHClustNCellEnergyRatio(0x0),
fHClustEnergyNCell(0x0),
fHConvEnergyPt(0x0),
fHPrimTrackPt(0x0),
fHPrimTrackEta(0x0),
fHPrimTrackPhi(0x0),
fHMatchDr(0x0),
fHMatchDz(0x0),
fHMatchEp(0x0),
fHPionEtaPhi(0x0),
fHAddPionEtaPhi(0x0),
fHPyPionEtaPhi(0x0),
fHdr(0),
fHPionMggPt(0x0),
fHPionMggAsym(0x0),
fHPionMggDgg(0x0),
fHPionInvMasses(0x0),
fHPionInvMassesSym(0x0),
fHPionInvMassesAsym(0x0),
fHPionInvMassesAdd1(0x0),
fHPionInvMassesAdd1Sym(0x0),
fHPionInvMassesAdd1Asym(0x0),
fHPionInvMassesAdd1Mult(0x0),
fHPionInvMassesAdd1MultSym(0x0),
fHPionInvMassesAdd1MultAsym(0x0),
fHPionInvMassesAdd2(0x0),
fHPionInvMassesMix(0x0),
fHPionInvMassesMix1(0x0),
fHPionInvMassesMix2(0x0),
fHJPInvMasses(0x0),
fHPrimPionInvMasses(0x0),
fHPrimPionInvMassesAsym(0x0),
fHPionEtaPhiConv(0x0),
fHPionMggPtConv(0x0),
fHPionMggAsymConv(0x0),
fHPionMggDggConv(0x0),
fHPionInvMassesConv(0x0),
fHPionInvMassesConvMix(0x0),
fHPionEtaPhiConvConv(0x0),
fHPionMggPtConvConv(0x0),
fHPionMggAsymConvConv(0x0),
fHPionMggDggConvConv(0x0),
fHPionInvMassesConvConv(0x0),
fHPionInvMassesConvConvMix(0x0),
fHConversionPoint(0),
fHPionTruthPt(),
fHPionTruthPtIn(),
fHPionTruthPtAcc(),
fHPionTruthPtConvAcc(),
fHEtaTruthPt(),
fHEtaTruthPtIn(),
fHEtaTruthPtAcc(),
fHEtaTruthPtConvAcc(),
fHGamTruthPt(),
fHGamTruthPtIn(),
fHGamTruthPtAcc(),
fHPionTruthPtAdd(),
fHPionTruthPtInAdd(),
fHPionTruthPtAccAdd(),
fHPionTruthPtConvAccAdd(),
fHEtaTruthPtAdd(),
fHEtaTruthPtInAdd(),
fHEtaTruthPtAccAdd(),
fHEtaTruthPtConvAccAdd(),
fHNMothers(0x0),
fHMixRotation(),
ipymin(0),
ipymax(0),
ipi0min(0),
ipi0max(0),
ietamin(0),
ietamax(0),
  fReaderGammas(0),
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
  fHRecTrue(),
  fHRecTrueAddPi0(),
  fHRecTrueAddEta(),
  fHECluEMCnofull(),
  fHECluEMCnofullAdd(),
  fHECluEMCelectron(),
  fHECluEMCpion(),
  fHECluEMCkaon(),
  fHECluEMCother(),
  fHECluEMCpi0single(),
  fHCorrection(),
  fHPionSm(),
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
          EmcEventList[i][j][k][l].SetGlobalInfo(0,0,0.,0.);
        }
      }
    }
  }

  DefineOutput(1, TList::Class());
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.,SPDVertex.,TPCVertex.,EMCALCells.,Tracks,EMCALTrigger.,SPDPileupVertices,TrkPileupVertices "
  "AOD:header,vertices,emcalCells,tracks";
}

//________________________________________________________________________
AliAnalysisTaskEMCALPi0Gamma::~AliAnalysisTaskEMCALPi0Gamma()
{
  // Destructor.
  
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput; fOutput = 0;
  }
  delete fPtRanges; fPtRanges = 0;
  fGeom = 0; // do not delete geometry when using instance
  delete fReco; fReco = 0;
  delete fTrClassNamesArr;
  delete fSelTracks;
  delete fSelPrimTracks;
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
      //TClass::GetClass("AliStaPart")->IgnoreTObjectStreamer();
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
  
  cout << "phimin=" << phimin << ", phimax=" << phimax << endl;
  
  // histograms
  Bool_t th1 =   TH1::GetDefaultSumw2();
  TH1::SetDefaultSumw2(kTRUE);
  Bool_t th2 =   TH2::GetDefaultSumw2();
  TH2::SetDefaultSumw2(kTRUE);
  fHCuts = new TH1F("hEventCuts","",5,0.5,5.5);
  fHCuts->GetXaxis()->SetBinLabel(1,"All");
  fHCuts->GetXaxis()->SetBinLabel(2,"PS");
  fHCuts->GetXaxis()->SetBinLabel(3,Form("%s: %.0f-%.0f",fCentVar.Data(),fCentFrom,fCentTo));
  fHCuts->GetXaxis()->SetBinLabel(4,"QFlag");
  fHCuts->GetXaxis()->SetBinLabel(5,Form("zvtx: %.0f-%.0f",fVtxZMin,fVtxZMax));
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
  fHTclsBeforeCuts = new TH1F("hTclsBeforeCuts","",fTrClassNamesArr->GetEntries(),0.5,0.5+fTrClassNamesArr->GetEntries());
  fHTclsAfterCuts = new TH1F("hTclsAfterCuts","",fTrClassNamesArr->GetEntries(),0.5,0.5+fTrClassNamesArr->GetEntries());
  for (Int_t i = 0; i<fTrClassNamesArr->GetEntries(); ++i) {
    const char *name = fTrClassNamesArr->At(i)->GetName();
    fHTclsBeforeCuts->GetXaxis()->SetBinLabel(1+i,name);
    fHTclsAfterCuts->GetXaxis()->SetBinLabel(1+i,name);
  }
  fOutput->Add(fHTclsBeforeCuts);
  fOutput->Add(fHTclsAfterCuts);
  
  // histograms for clusters
  Int_t nsm = fGeom->GetEMCGeometry()->GetNumberOfSuperModules();
  if (!fTrainMode) {
    fHClustNoEvt = new TH1F("hClustNoEvt","",2000,0,2000);
    fHClustNoEvt->SetXTitle("# Clusters");
    fOutput->Add(fHClustNoEvt);
    fHClustAccEvt = new TH1F("hClustAccEvt","",2000,0,2000);
    fHClustAccEvt->SetXTitle("# Clusters");
    fOutput->Add(fHClustAccEvt);
    fHClustEccentricity = new TH1F("hClustEccentricity","",100,-0.1,1.1);
    fHClustEccentricity->SetXTitle("#epsilon_{C}");
    fOutput->Add(fHClustEccentricity);
    fHClustEtaPhi = new TH2F("hClustEtaPhi","",500,-0.8,0.8,100*nsm,phimin,phimax);
    fHClustEtaPhi->SetXTitle("#eta");
    fHClustEtaPhi->SetYTitle("#varphi");
    fOutput->Add(fHClustEtaPhi);
    fHClustEnergyPt = new TH2F("hClustEnergyPt","",250,0,50,250,0,50);
    fHClustEnergyPt->SetXTitle("E [GeV]");
    fHClustEnergyPt->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHClustEnergyPt);
    fHClustEnergySigma = new TH2F("hClustEnergySigma","",100,0,100,500,0,50);
    fHClustEnergySigma->SetXTitle("E_{C} * #sigma_{max} [GeV*cm]");
    fHClustEnergySigma->SetYTitle("E_{C} [GeV]");
    fOutput->Add(fHClustEnergySigma);
    fHClustSigmaSigma = new TH2F("hClustSigmaSigma","",500,0,50,500,0,50);
    fHClustSigmaSigma->SetXTitle("#lambda_{0} [cm]");
    fHClustSigmaSigma->SetYTitle("#sigma_{max} [cm]");
    fOutput->Add(fHClustSigmaSigma);
    fHClustNCellEnergyRatio = new TH2F("hClustNCellEnergyRatio","",27,-0.5,26.5,101,-0.05,1.05);
    fHClustNCellEnergyRatio->SetXTitle("N_{cells}");
    fHClustNCellEnergyRatio->SetYTitle("E^{max}_{cell}/E_{clus}");
    fOutput->Add(fHClustNCellEnergyRatio);
    fHClustEnergyNCell = new TH2F("hClustEnergyNCell","",200,0,100,50,0,50);
    fHClustEnergyNCell->SetXTitle("E_{clus}");
    fHClustEnergyNCell->SetYTitle("N_{cells}");
    fOutput->Add(fHClustEnergyNCell);
  }
  
  // histogram for conversion photons
  fHConvEnergyPt = new TH2F("hConvEnergyPt","",250,0,50,250,0,50);
  fHConvEnergyPt->SetXTitle("E [GeV]");
  fHConvEnergyPt->SetYTitle("p_{T} [GeV/c]");
  fOutput->Add(fHConvEnergyPt);
  
  
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
    
    const Int_t massbins = 200;
    const Double_t massmax = 1.0;

    fHPionEtaPhi = new TH2F("hPionEtaPhi","",100,-0.8,0.8,100*nsm,phimin,phimax);
    fHPionEtaPhi->SetXTitle("#eta_{#gamma#gamma}");
    fHPionEtaPhi->SetYTitle("#varphi_{#gamma#gamma}");
    fOutput->Add(fHPionEtaPhi);


    fHPionMggPt = new TH2F("hPionMggPt","",massbins,0,massmax,nbins,0,ptmax);
    fHPionMggPt->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionMggPt->SetYTitle("p_{T}^{#gamma#gamma} [GeV/c]");
    fOutput->Add(fHPionMggPt);
    fHPionMggAsym = new TH2F("hPionMggAsym","",500,0,1,100,0,1);
    fHPionMggAsym->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionMggAsym->SetYTitle("Z_{#gamma#gamma} [GeV]");
    fOutput->Add(fHPionMggAsym);
    fHPionMggDgg = new TH2F("hPionMggDgg","",100,0,1,500,0,15);
    fHPionMggDgg->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionMggDgg->SetYTitle("opening angle [grad]");
    fOutput->Add(fHPionMggDgg);
    
    if(fDoConvAna){
      fHPionEtaPhiConv = new TH2F("hPionEtaPhiConv","",100,-0.8,0.8,100*nsm,phimin,phimax);
      fHPionEtaPhiConv->SetXTitle("#eta_{#gamma#gamma}");
      fHPionEtaPhiConv->SetYTitle("#varphi_{#gamma#gamma}");
      fOutput->Add(fHPionEtaPhiConv);
      fHPionMggPtConv = new TH2F("hPionMggPtConv","",500,0,1,100,0,20.0);
      fHPionMggPtConv->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPionMggPtConv->SetYTitle("p_{T}^{#gamma#gamma} [GeV/c]");
      fOutput->Add(fHPionMggPtConv);
      fHPionMggAsymConv = new TH2F("hPionMggAsymConv","",500,0,1,100,0,1);
      fHPionMggAsymConv->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPionMggAsymConv->SetYTitle("Z_{#gamma#gamma} [GeV]");
      fOutput->Add(fHPionMggAsymConv);
      fHPionMggDggConv = new TH2F("hPionMggDggConv","",500,0,1,100,0,10);
      fHPionMggDggConv->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPionMggDggConv->SetYTitle("opening angle [grad]");
      fOutput->Add(fHPionMggDggConv);
      fHPionEtaPhiConvConv = new TH2F("hPionEtaPhiConvConv","",100,-0.8,0.8,100*nsm,phimin,phimax);
      fHPionEtaPhiConvConv->SetXTitle("#eta_{#gamma#gamma}");
      fHPionEtaPhiConvConv->SetYTitle("#varphi_{#gamma#gamma}");
      fOutput->Add(fHPionEtaPhiConvConv);
      fHPionMggPtConvConv = new TH2F("hPionMggPtConvConv","",500,0,1,100,0,20.0);
      fHPionMggPtConvConv->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPionMggPtConvConv->SetYTitle("p_{T}^{#gamma#gamma} [GeV/c]");
      fOutput->Add(fHPionMggPtConvConv);
      fHPionMggAsymConvConv = new TH2F("hPionMggAsymConvConv","",500,0,1,100,0,1);
      fHPionMggAsymConvConv->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPionMggAsymConvConv->SetYTitle("Z_{#gamma#gamma} [GeV]");
      fOutput->Add(fHPionMggAsymConvConv);
      fHPionMggDggConvConv = new TH2F("hPionMggDggConvConv","",500,0,1,100,0,10);
      fHPionMggDggConvConv->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPionMggDggConvConv->SetYTitle("opening angle [grad]");
      fOutput->Add(fHPionMggDggConvConv);
    }
    //  const Int_t nbins = 300;
    //  const Int_t ptmax = 30;
    //    Double_t xbins[nbins] = {0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,12.5,15,20,25,50};
    //    fPtRanges = new TAxis(nbins-1,xbins);
    
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

    // no clean gammas
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

    fHPionInvMassesMix = new TH2F("hPionInvMassMix","hPionInvMassMix",massbins,0,massmax,nbins,0,ptmax);
    fHPionInvMassesMix->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionInvMassesMix->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHPionInvMassesMix);
    
    fHPionInvMassesMix1 = new TH2F("hPionInvMassMix1","hPionInvMassMix1",massbins,0,massmax,nbins,0,ptmax);
    fHPionInvMassesMix1->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionInvMassesMix1->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHPionInvMassesMix1);

    fHPionInvMassesMix2 = new TH2F("hPionInvMassMix2","hPionInvMassMix2",massbins,0,massmax,nbins,0,ptmax);
    fHPionInvMassesMix2->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionInvMassesMix2->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHPionInvMassesMix2);

    // is there a JPsi?
    fHJPInvMasses = new TH2F("hJPInvMass","hJPInvMass",200,2.,4.,20,0,20);
    fHJPInvMasses->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHJPInvMasses->SetYTitle("p_{T} [GeV/c]");
    fOutput->Add(fHJPInvMasses);
    
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
    
    if(fDoConvAna){
      fHPionInvMassesConv = new TH2F("hPionInvMassConv","hPionInvMassConv",500,0,1,nbins,0,ptmax);
      fHPionInvMassesConv->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPionInvMassesConv->SetYTitle("p_{T} [GeV/c]");
      fOutput->Add(fHPionInvMassesConv);
      
      fHPionInvMassesConvMix = new TH2F("hPionInvMassConvMix","hPionInvMassConvMix",500,0,1,nbins,0,ptmax);
      fHPionInvMassesConvMix->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPionInvMassesConvMix->SetYTitle("p_{T} [GeV/c]");
      fOutput->Add(fHPionInvMassesConvMix);
      
      fHPionInvMassesConvConv = new TH2F("hPionInvMassConvConv","hPionInvMassConvConv",500,0,1,nbins,0,ptmax);
      fHPionInvMassesConvConv->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPionInvMassesConvConv->SetYTitle("p_{T} [GeV/c]");
      fOutput->Add(fHPionInvMassesConvConv);
      
      fHPionInvMassesConvConvMix = new TH2F("hPionInvMassConvConvMix","hPionInvMassConvConvMix",500,0,1,nbins,0,ptmax);
      fHPionInvMassesConvConvMix->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      fHPionInvMassesConvConvMix->SetYTitle("p_{T} [GeV/c]");
      fOutput->Add(fHPionInvMassesConvConvMix);
    }
    
    // histogram for conversion point
    fHConversionPoint = new TH2F("hConversionPoint","conversion point in xy",1000,-500,500,1000,-500,500);
    fHConversionPoint->SetXTitle("x");
    fHConversionPoint->SetYTitle("y");
    fOutput->Add(fHConversionPoint);
  }
  
  TH1::SetDefaultSumw2(th1);
  TH2::SetDefaultSumw2(th2);
  PostData(1, fOutput);
  
  //  const Int_t nbins = 20;
  //  Double_t xbins[nbins+1] = {0.,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,12.5,15,20,25,50};
  //  const Int_t nbins = 300;
  //  const Int_t ptmax = 30;
  
  // MC histograms
  if(fMcMode){
    
    fHAddPionEtaPhi = new TH2F("hAddPionEtaPhi","",150,-2.5,2.5,100,0.,6.3);
    fHAddPionEtaPhi->SetXTitle("#eta_{#gamma#gamma}");
    fHAddPionEtaPhi->SetYTitle("#varphi_{#gamma#gamma}");
    fOutput->Add(fHAddPionEtaPhi);

    fHPyPionEtaPhi = new TH2F("hPyPionEtaPhi","",150,-2.5,2.5,100,0.,6.3);
    fHPyPionEtaPhi->SetXTitle("#eta_{#gamma#gamma}");
    fHPyPionEtaPhi->SetYTitle("#varphi_{#gamma#gamma}");
    fOutput->Add(fHPyPionEtaPhi);

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
    
    fHPionTruthPtConvAcc = new TH1F("hPionTruthPtConvAcc","pi0 truth pT from MC within acceptance (eta<0.8 for converted)",nbins,0,ptmax);
    fHPionTruthPtConvAcc->SetXTitle("p_{T}");
    fOutput->Add(fHPionTruthPtConvAcc);
    
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
    
    fHEtaTruthPtConvAcc = new TH1F("hEtaTruthPtConvAcc","eta truth pT from MC within acceptance (eta<0.8 for converted)",nbins,0,ptmax);
    fHEtaTruthPtConvAcc->SetXTitle("p_{T}");
    fOutput->Add(fHEtaTruthPtConvAcc);
    
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
    
    fHPionTruthPtConvAccAdd = new TH1F("hPionTruthPtConvAccAdd","added pi0 truth pT from MC within acceptance (eta<0.8 for converted)",nbins,0,ptmax);
    fHPionTruthPtConvAccAdd->SetXTitle("p_{T}");
    fOutput->Add(fHPionTruthPtConvAccAdd);
    
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
    
    fHEtaTruthPtConvAccAdd = new TH1F("hEtaTruthPtConvAccAdd","added eta truth pT from MC within acceptance (eta<0.8 for converted)",nbins,0,ptmax);
    fHEtaTruthPtConvAccAdd->SetXTitle("p_{T}");
    fOutput->Add(fHEtaTruthPtConvAccAdd);
    
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

    fHRecTrue = new TH2F("hRecTrue","truth energy/cluster energy vs. cluster energy",200,0,10,200,0,10);
    fHRecTrue->SetXTitle("E_{clu}");
    fHRecTrue->SetYTitle("E_{tru}/E_{clu}");
    fOutput->Add(fHRecTrue);
    
    fHRecTrueAddPi0 = new TH2F("hRecTrueAddPi0","truth energy/cluster energy vs. cluster energy for added pi0 MC",200,0,10,200,0,10);
    fHRecTrueAddPi0->SetXTitle("E_{clu}");
    fHRecTrueAddPi0->SetYTitle("E_{tru}/E_{clu}");
    fOutput->Add(fHRecTrueAddPi0);

    fHRecTrueAddEta = new TH2F("hRecTrueAddEta","truth energy/cluster energy vs. cluster energy for added eta MC",200,0,10,200,0,10);
    fHRecTrueAddEta->SetXTitle("E_{clu}");
    fHRecTrueAddEta->SetYTitle("E_{tru}/E_{clu}");
    fOutput->Add(fHRecTrueAddEta);
    
    // also for clusters with more than one contribution
    //    fHMCpartfracnofull = new TH2F("hMCpartfracnofull","fraction of most energy MC particle vs. fraction of this particle to all",100,0,2,100,0.8,1.2);
    //    fHMCpartfracnofull->SetXTitle("E_{MC}/E_{Clu}");
    //    fHMCpartfracnofull->SetYTitle("E_{MC}^{highest}/E_{MC}^{all}");
    //    fOutput->Add(fHMCpartfracnofull);
    
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
  
	if(fRotateMixed){
		// more histograms
		fHMixRotation = new TH1F("hMixRotation","rotation angle of mixed events",100,-6.28,6.28);
		fHMixRotation->SetXTitle("phi");
		fOutput->Add(fHMixRotation);
	}
  fHCorrection = new TH1F("hCorrection","correction factor for single photon",100,0,2);
  fHCorrection->SetXTitle("correction factor");
  fOutput->Add(fHCorrection);
  
  fHPionSm = new TH2F("hPionSm","mass shift due to energy scale",200,0,0.5,200,0,0.5);
  fHPionSm->SetXTitle("original m_inv");
  fHPionSm->SetYTitle("changed m_inv");
  fOutput->Add(fHPionSm);
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::UserExec(Option_t *)
{
  // Called for each event.
  
  if (!InputEvent())
    return;
  
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  
  
	AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (am->GetInputEventHandler());
	//if (!esdH)
	//{
	//	Printf("ERROR: Could not get ESDInputHandler");
	//}
  
	AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (am->GetInputEventHandler());
	if (!aodH && !esdH)
	{
		Printf("ERROR: Could not get AODInputHandler");
	}
  
	if(esdH)
    fEsdEv = esdH->GetEvent();
  
  else if(aodH) fAodEv = aodH->GetEvent();
  
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
    //Int_t IsPi0Flat=0;
    //printf("N generators %d \n", nGenerators);
    //Int_t fNMCProducedMin = 0;
    //Int_t fNMCProducedMax = 0;
    
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
  
  //  // monte carlo headers from Philipp
  //  TList *list = NULL;
  //  if(esdH)
  //    list = fEsdEv->GetList();
  //  else if(aodH)
  //    list = fAodEv->GetList();
  //
  //  //hijingGenHeader = NULL;
  //  //PythiaGenHeader = NULL;
  //  if(fMcMode && list){
  //      //AliGenEventHeader* header = (AliGenEventHeader*)list->FindObject(AliGenEventHeader::StdBranchName());
  //      AliAODMCHeader* header = (AliAODMCHeader*)list->FindObject(AliAODMCHeader::StdBranchName());
  //
  //    if(!header){cout << "no gen event header" << endl;
  //      return;
  //    }
  //      TList* headerList = header->GetCocktailHeaders();
  //
  //      for(Int_t i = 0; i < headerList->GetEntries(); i++)
  //      {
  //          hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(headerList->At(i));
  //          if(hijingGenHeader) break;
  //      }
  //      if(!hijingGenHeader){
  //        for(Int_t i = 0; i < headerList->GetEntries(); i++)
  //        {
  //          PythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(headerList->At(i));
  //          if(PythiaGenHeader) break;
  //        }
  //      }
  //  }
  
  
  
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
  /*
   fEsdEv = dynamic_cast<AliESDEvent*>(InputEvent());
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
   offtrigger = ((AliInputEventHandler*)(am->GetInputEventHandler()))->IsEventSelected();
   }
   else {
   fAodEv = dynamic_cast<AliAODEvent*>(InputEvent());
   if (!fAodEv) {
   AliFatal("Neither ESD nor AOD event found");
   return;
   }
   am->LoadBranch("header");
   offtrigger =  fAodEv->GetHeader()->GetOfflineTrigger();
   }
   */
  if (!fMcMode && (offtrigger & AliVEvent::kFastOnly)) {
    AliWarning(Form("EMCAL not in fast only partition"));
    return;
  }
  
  
  if (fDoTrMatGeom && !AliGeomManager::GetGeometry()) { // get geometry
    AliWarning("Accessing geometry from OCDB, this is not very efficient!");
    AliCDBManager *cdb = AliCDBManager::Instance();
    if (!cdb->IsDefaultStorageSet())
      cdb->SetDefaultStorage("raw://");
    Int_t runno = InputEvent()->GetRunNumber();
    if (runno != cdb->GetRun())
      cdb->SetRun(runno);
    AliGeomManager::LoadGeometry();
  }
  
  if (!AliGeomManager::GetGeometry()&&!fIsGeoMatsSet) { // set misalignment matrices (stored in first event)
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
  
  if (!TGeoGlobalMagField::Instance()->GetField()) { // construct field map
    if (fEsdEv) {
      const AliESDRun *erun = fEsdEv->GetESDRun();
      AliMagF *field = AliMagF::CreateFieldMap(erun->GetCurrentL3(),
                                               erun->GetCurrentDip(),
                                               AliMagF::kConvLHC,
                                               kFALSE,
                                               erun->GetBeamEnergy(),
                                               erun->GetBeamType());
      TGeoGlobalMagField::Instance()->SetField(field);
    }
    else {
      Double_t pol = -1; //polarity
      Double_t be = -1;  //beam energy
      AliMagF::BeamType_t btype = AliMagF::kBeamTypepp;
      Int_t runno = fAodEv->GetRunNumber();
      if (runno>=136851 && runno<138275) {
        pol = -1;
        be = 2760;
        btype = AliMagF::kBeamTypeAA;
      } else if (runno>=138275 && runno<=139517) {
        pol = +1;
        be = 2760;
        btype = AliMagF::kBeamTypeAA;
      } else {
        AliError(Form("Do not know the bfield parameters for run %d! Using defaults!!!", runno));
      }
      TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", pol, pol, AliMagF::k5kG, btype, be));
    }
  }
  
  /*
   AliMCEvent* fMC;
   if(fMcMode)
   {
   AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
   if (!mcH)
   {
   printf("ERROR: Could not get MCInputHandler");
   //return;
   }
   else
   fMC = mcH->MCEvent();
   
   }
   */
  
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
  for (Int_t i = 0; i<fTrClassNamesArr->GetEntries(); ++i) {
    const char *name = fTrClassNamesArr->At(i)->GetName();
    TRegexp regexp(name);
    if (trgclasses.Contains(regexp))
      fHTclsBeforeCuts->Fill(1+i);
  }
  
  if (fDoPSel && offtrigger==0)
    return;
  
  fHCuts->Fill(cut++);
  
  
  const AliCentrality *centP = InputEvent()->GetCentrality();
  Double_t cent = centP->GetCentralityPercentileUnchecked(fCentVar);
  fHCent->Fill(cent);
  if (cent<fCentFrom||cent>fCentTo)
    return;
  
  fHCuts->Fill(cut++);
  
  if (fUseQualFlag) {
    if (centP->GetQuality()>0)
      return;
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
  
  for (Int_t i = 0; i<fTrClassNamesArr->GetEntries(); ++i) {
    const char *name = fTrClassNamesArr->At(i)->GetName();
    TRegexp regexp(name);
    if (trgclasses.Contains(regexp))
      fHTclsAfterCuts->Fill(1+i);
  }
  
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
  if (1 && !fClusName.IsNull()) {
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
  
  // need to get converted photons from something, maybe apply quality cuts, and then use them for invariant mass calculation
  // also need space in event mixing to store all those candidates
  //fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
  //if(!fV0Reader){printf("Error: No V0 Reader");} // GetV0Reader
  //Int_t eventQuality = ((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetEventQuality();
  
  //if(fV0Reader)
    //fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut // GetGoodGammas instead?
  //cout << "number of conversion photons: " << fV0Reader->GetNReconstructedGammas() << endl;
  
  //else
  fReaderGammas = NULL;
  
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
  Int_t vtxClass = 2;
  Double_t vtxcuts = (fVtxZMax-fVtxZMin)/nZClass;
  if(vertex->GetZ()>fVtxZMin && vertex->GetZ()<=fVtxZMin+vtxcuts) vtxClass=0;
  else if(vertex->GetZ()>fVtxZMin+vtxcuts && vertex->GetZ()<=fVtxZMin+2*vtxcuts) vtxClass=1;
  else if(vertex->GetZ()>fVtxZMin+2*vtxcuts && vertex->GetZ()<=fVtxZMin+3*vtxcuts) vtxClass=2;
  else if(vertex->GetZ()>fVtxZMin+3*vtxcuts && vertex->GetZ()<=fVtxZMin+4*vtxcuts) vtxClass=3;
  else if(vertex->GetZ()>fVtxZMin+4*vtxcuts && vertex->GetZ()<fVtxZMax) vtxClass=4;
  
  //    vtxClass = 0;
  
  Int_t MulClass = 0;
  
  GetMulClass(MulClass);
  
	Float_t phitrig = 0;
  Float_t thetatrig = 0;
  Double_t pt_max = 0;
  Int_t ptClass = 0;
  if (!fTrainMode) {
    pt_max = FillClusHists(phitrig, thetatrig);
    //FillConvHists();
    
    if(pt_max > 1 && pt_max<=3) ptClass = 1;
    else if(pt_max > 3 && pt_max<=6) ptClass = 2;
    else ptClass = 3;

    ptClass = 0;
    
    // these are where the pi0 stuff is filled
    FillPionHists();
    if(fDoConvAna){
      FillPionConv();
      FillConvConv();
    }
    FillMixHists(MulClass,vtxClass,ptClass,phitrig,thetatrig);
    if(fDoConvAna){
      FillMixConv(MulClass,vtxClass,ptClass);
      FillMixConvConv(MulClass,vtxClass,ptClass);
    }
    
    FillOtherHists();
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
  
  //cout << "filling cluster histograms ";
  
  bool bprint = 0;
  
  max_phi = 0;
  max_theta = 0;
	//Double_t max_pt = 0;
	// Fill histograms related to cluster properties.
  
  TObjArray *clusters = fEsdClusters;
  if (!clusters)
    clusters = fAodClusters;
  if (!clusters){
    cout << "no clusters node in event!" << endl;
    return 0;
  }
  
  Int_t nclus = clusters->GetEntries();
  
  //       cout << nclus << " clusters in event ";
  
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  
  fHClustNoEvt->Fill(nclus);
  thisEvent.SetGlobalInfo(0,0,0,0);
  
  int nclusters = 0;
  
  for(Int_t i = 0; i<nclus; ++i) {
    AliVCluster *clus = static_cast<AliVCluster*>(clusters->At(i));
    if (!clus)
      continue;
    // emcal cluster?
    if (!clus->IsEMCAL()){
      //      cout << "at " << i << " cluster not EMCAL" << endl;
      continue;
    }
    //    Double_t ecorr = 1;
    //Double_t ecorr = fcorrect->Eval(2.0*clusterVec.E());
    //    fHCorrection->Fill(ecorr);
    //    TLorentzVector clusterVecCorr(clusterVec.Px()*ecorr,clusterVec.Py()*ecorr,clusterVec.Pz()*ecorr,clusterVec.E()*ecorr);
    //clusterVecCorr.SetPxPyPzE(clusterVec.Px()*ecorr,clusterVec.Py()*ecorr,clusterVec.Pz()*ecorr,clusterVec.E()*ecorr);
    Double_t maxAxis    = clus->GetTOF(); //sigma
    Double_t clusterEcc = clus->Chi2();   //eccentricity
    fHClustEccentricity->Fill(clusterEcc);
    // here we only fill after fulfilling cuts)
    /*
     if(clusterVecCorr.Pt()>max_pt){
     max_phi = (float)clusterVecCorr.Phi();
     max_theta = (float)clusterVecCorr.Theta();
     max_pt = clusterVecCorr.Pt();
     }
     */
    // fill clusters into this event

    
    // apply cluster cuts first
    if (clus->E()<fMinE)
      continue;
    if (clus->GetNCells()<fNminCells)
      continue;
    if (GetMaxCellEnergy(clus)/clus->E()<fMinErat)
      continue;
    if (clus->Chi2()<fMinEcc) // eccentricity cut
      continue;
    if(clus->GetM02()>fM02)
      continue;

    TLorentzVector clusterVec;
    clus->GetMomentum(clusterVec,vertex);

    // fill cluster histograms
    fHClustEtaPhi->Fill(clusterVec.Eta(),clusterVec.Phi());
    fHClustEnergyPt->Fill(clusterVec.E(),clusterVec.Pt());
    fHClustEnergySigma->Fill(clus->E()*maxAxis,clus->E());
    fHClustSigmaSigma->Fill(max(clus->GetM02(),clus->GetM20()),clus->E()*maxAxis);
    fHClustNCellEnergyRatio->Fill(clus->GetNCells(),GetMaxCellEnergy(clus)/clus->E());
    fHClustEnergyNCell->Fill(clus->E(),clus->GetNCells());
    nclusters++;

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
      
      TLorentzVector clusterVecCorr1(clusterVec.Px(),clusterVec.Py(),clusterVec.Pz(),En);

      thisEvent.hit[i].thishit=clusterVecCorr1;
      thisEvent.hit[i].hittype=100;
      thisEvent.hit[i].weight=1.;
      thisEvent.hit[i].imo=1;
    }
    // go through MC information of clusters
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
        while(imother >= 0){
          AliMCParticle *tmppart = static_cast<AliMCParticle*>(mcEvent->GetTrack(ipart));
          if(bprint){
            cout << " pid " << tmppart->Label() << ", type " << tmppart->PdgCode() << ", mid " << tmppart->GetMother() << " ==>";
          }
          imother = tmppart->GetMother();
          if(imother >= 0){
            ipart = imother;
          }
          iit++;
        }
        fHNMothers->Fill(iit,mcP->E());
        Double_t wgt = 1.;
        
        AliMCParticle *McMo = static_cast<AliMCParticle*>(mcEvent->GetTrack(ipart));
        Double_t mcPt = McMo->Pt();
        Double_t mcEta = McMo->Eta();
        
	wgt = CalcWeight(mcPt,mcEta,1);
        
        // check if generated or added
        bool bGen = kTRUE;
        bool bAddPi0 = kFALSE;
        bool bAddEta = kFALSE;
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

        // loop up until pi0, then look for mother of pi0

        Int_t idmother = 1;
        Int_t ithispart = mcP->Label();
        Int_t iPID = 0;
        Bool_t bpi0 = 1;
        Short_t tmpflag = 0;
        while(idmother != 111){
          AliMCParticle *tmppart = static_cast<AliMCParticle*>(mcEvent->GetTrack(ithispart));
          Int_t itmp = tmppart->GetMother();
	  if(itmp < 1){
	      bpi0 = 0;
	      break;
	  }
          AliMCParticle *tmppartmo = static_cast<AliMCParticle*>(mcEvent->GetTrack(itmp));
          idmother = tmppartmo->PdgCode();
          ithispart = itmp;
        }
        // if it is not from a pi0, flag 100
        // from primary pi0, flag 101
        // from secondary pi0, flag 102
        // from K0, flag 103
        // from material, flag 104
        // only for generator, not added signals
	if(bGen){
	  if(bpi0 == 0){
	    tmpflag = 100;
	  }
	  else{
	    AliMCParticle *tmppart = static_cast<AliMCParticle*>(mcEvent->GetTrack(ithispart));
	    Int_t itmp = tmppart->GetMother();
	    if(itmp<nPTracksMC){
	      tmpflag = 101;
	    }
	    else{
	      if( ((AliMCParticle*)mcEvent->GetTrack(tmppart->GetMother()))->PdgCode() ==  310 ||
		  ((AliMCParticle*)mcEvent->GetTrack(tmppart->GetMother()))->PdgCode() == -310  ){
		tmpflag = 103;
	      }
	      else if(((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() ==  2212 || //proton
		      ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() == -2212 || //anti-proton
		      ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() ==  2112 || //neutron
		      ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() == -2112 || //anti-neutron
		      ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() ==  321  || //K+
		      ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() == -321  || //K-
		      ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() ==  211  || //pi+
		      ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() == -211 ){    //pi-)
		tmpflag = 104;
	      }
	      else{
		tmpflag = 102;
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
            fHRecTrue->Fill(cle,mce/cle);
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
              if(imother == jmother){
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
            fHRecTrueAddPi0->Fill(cle,mce/cle);
          }
        }
        if(bAddEta){
          // if leading particle with high contribution and photon
          if(mcP->PdgCode() == 22 && efrac>0.5){
            // lets fill histograms
            fHECluEMCAddEta->Fill(mce,cle);
            fHRecTrueAddEta->Fill(cle,mce/cle);
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
//              Int_t d1 = tmppart2->GetFirstDaughter();
//              Int_t d2 = tmppart2->GetLastDaughter();
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
          Int_t d1 = McMo->GetFirstDaughter();
          Int_t d2 = McMo->GetLastDaughter();
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
        if(1){
          thisEvent.hit[i].thishit=clusterVecCorr1;
          thisEvent.hit[i].imo=ipart;
          thisEvent.hit[i].pid=mcP->PdgCode();
          if(!bGen){
            thisEvent.hit[i].weight = wgt;
          }
          else{
            thisEvent.hit[i].weight = 1.;
          }
          
          thisEvent.hit[i].hittype=-1;
          if(bGen){
            // use hittype to flag particle history
            // if it is not from a pi0, flag 100
            // from primary pi0, flag 101
            // from secondary pi0, flag 102
            // from K0, flag 103
            // from material, flag 104
            
            thisEvent.hit[i].hittype=tmpflag;
          }
          else{
            if(bAddPi0){
              thisEvent.hit[i].hittype=1;
            }
            else if(bAddEta){
              thisEvent.hit[i].hittype=2;
            }
          }
          thisEvent.hit[i].bclean = bcl;
        }
      } // end ilabel
    } // end fmcmode
  } // end cluster loop
  fHClustAccEvt->Fill(nclusters);
  thisEvent.SetGlobalInfo(nclusters,0,max_phi,max_theta);
  //cout << "this event size: " << thisEvent.evsize() << endl;
  return 1;
  //  return max_pt;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::FillConvHists()
{
  // Fill histograms related to cluster properties.
  TLorentzVector clusterVec;
  
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  
  if(fReaderGammas){
    for(Int_t j = 0; j < fReaderGammas->GetEntriesFast(); j++){
      AliAODConversionPhoton* PhotonCandidate = NULL;//(AliAODConversionPhoton*) fReaderGammas->At(j);
      Double_t x = PhotonCandidate->GetConversionX();
      Double_t y = PhotonCandidate->GetConversionY();
      fHConversionPoint->Fill(x,y);
      if(!PhotonCandidate) continue;
      clusterVec.SetPxPyPzE(PhotonCandidate->GetPx(),PhotonCandidate->GetPy(),PhotonCandidate->GetPz(),PhotonCandidate->GetPhotonP());
      fHConvEnergyPt->Fill(clusterVec.E(),clusterVec.Pt());
      //cout << "photon converted at (" << x << "," << y << ")" << endl;
    }
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::CalcMcInfo()
{
  // Get Mc truth particle information.
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
  
  AliMCEvent *mcEvent = MCEvent();
  if (!mcEvent){
    cout << "no MC event" << endl;
    return;
  }
  
  const AliVVertex *evtVtx = mcEvent->GetPrimaryVertex();
  if (!evtVtx)
    return;
  
  mcEvent->PreReadAll();
  
  Int_t nTracks = mcEvent->GetNumberOfPrimaries();
  for (Int_t iTrack = 0; iTrack<nTracks; ++iTrack) {
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
    Double_t dR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX())*(mcP->Xv()-evtVtx->GetX()) +
                              (mcP->Yv()-evtVtx->GetY())*(mcP->Yv()-evtVtx->GetY()));
    
    
    if(dR > 0.1)
      continue;
    
    // use nprimarytracks (these are the first ones)
    
    // is it generator particle or added signal?
    bool bGen = kTRUE;
    bool bAddPi0 = kFALSE;
    bool bAddEta = kFALSE;
    if(pythiaHeader){
      if(mcP->Label() > ipymax){
        bGen = kFALSE;
        if(mcP->Label() <= ipi0max)
          bAddPi0 = kTRUE;
        else if(mcP->Label() > ipi0max && mcP->Label() <= ietamax)
          bAddEta = kTRUE;
      }
    }
    
    double mcPt = mcP->Pt();
    double mcEta = mcP->Eta();
    double wgt2 = 1.;

    if(bGen && !(bAddPi0 || bAddEta)){
      // fill truth histogram
      if(mcP->PdgCode() == 111){
        fHPionTruthPt->Fill(mcP->Pt());
        fHPyPionEtaPhi->Fill(mcP->Eta(),mcP->Phi());
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
        fHAddPionEtaPhi->Fill(mcP->Eta(),mcP->Phi(),wgt2);
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
    Double_t pt = mcP->Pt() ;
    if (pt<0.3)
      continue;
    Double_t eta = mcP->Eta();
    if (eta<-1.0||eta>1.0)
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
      if(mcPt > 0.01 && mcPt < 30 && mcEta > -2 && mcEta < 2){
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
      if(mcPt > 0.01 && mcPt < 30 && mcEta > -2 && mcEta < 2){
	wgt2 = CalcWeight(mcPt,mcEta,1);
      }
      else{
	continue;
      }
        fHEtaTruthPtInAdd->Fill(mcP->Pt(),wgt2);
      }
      
    }

    bool binp = true; // does it count as input?

    Int_t d1 = mcP->GetFirstDaughter();
    Int_t d2 = mcP->GetLastDaughter();
    if (d1<0)
      return;
    if (d2<0)
      d2=d1;
    
    if(d2-d1 != 1){
      continue;
    }
    
    // check it particle decays into two photons (binp), decay photons are in acceptance (bacc/baccconv)
    bool bacc = true;
    bool baccconv = true;
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
      if(i==d1 &&!(dmc->PdgCode()==22 && eta_d>-0.8 && eta_d<0.8)){
        baccconv = false;
      }
      if(i==d2 &&!(dmc->PdgCode()==22 && eta_d>etamin && eta_d<etamax && phi_d>phimin && phi_d<phimax)){
        baccconv = false;
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
      
        if(mcP->PdgCode() == 111 && baccconv){
          fHPionTruthPtConvAcc->Fill(mcP->Pt());
        }
        if(mcP->PdgCode() == 221 && baccconv){
          fHEtaTruthPtConvAcc->Fill(mcP->Pt());
        }
      }
    }
    
    
    if(bAddPi0 && !bAddEta && !bGen){
      // if both photons are on EMCal
      if(binp){
        if(mcP->PdgCode() == 111 && bacc){
          fHPionTruthPtAccAdd->Fill(mcP->Pt(),wgt2);
        }
      
        if(mcP->PdgCode() == 111 && baccconv){
          fHPionTruthPtConvAccAdd->Fill(mcP->Pt(),wgt2);
        }
      }
    }
    
    
    if(bAddEta && !bAddPi0 && !bGen){
      // if both photons are on EMCal
      if(binp){
        if(mcP->PdgCode() == 221 && bacc){
          fHEtaTruthPtAccAdd->Fill(mcP->Pt(),wgt2);
        }
        if(mcP->PdgCode() == 221 && baccconv){
          fHEtaTruthPtConvAccAdd->Fill(mcP->Pt(),wgt2);
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
  
//  Double_t vertex[3] = {0};
//  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
//  
//  TObjArray *clusters = fEsdClusters;
//  if (!clusters)
//    clusters = fAodClusters;
//  
//  if (clusters) {
    TLorentzVector clusterVec1;
    TLorentzVector clusterVec2;
    TLorentzVector pionVec;
    Short_t hitclass1 = 0;
    Short_t hitclass2 = 0;
    //TLorentzVector pionVecCorr;
    
//    Int_t nclus = clusters->GetEntries();
    Int_t nclus = thisEvent.nHits;
    if(nclus<2){
      return;
    }
    for (Int_t i = 0; i<nclus; ++i) {
//      AliVCluster *clus1 = static_cast<AliVCluster*>(clusters->At(i));
//      if (!clus1)
//        continue;
//      if (!clus1->IsEMCAL())
//        continue;
//      if (clus1->E()<fMinE)
//        continue;
//      if (clus1->GetNCells()<fNminCells)
//        continue;
//      if (GetMaxCellEnergy(clus1)/clus1->E()<fMinErat)
//        continue;
//      if (clus1->Chi2()<fMinEcc) // eccentricity cut
//        continue;
//      clus1->GetMomentum(clusterVec1,vertex);
      //Double_t ecorr1 = fcorrect->Eval(2.0*clusterVec1.E());
      //      Double_t ecorr1 = 1;
      //      TLorentzVector clusterVecCorr1(clusterVec1.Px()*ecorr1,clusterVec1.Py()*ecorr1,clusterVec1.Pz()*ecorr1,clusterVec1.E()*ecorr1);
      //      clusterVecCorr1.SetPxPyPzE(clusterVec1.Px()*ecorr1,clusterVec1.Py()*ecorr1,clusterVec1.Pz()*ecorr1,clusterVec1.E()*ecorr1);
      //Bool_t bkeep1 = 1;
      clusterVec1 = thisEvent.hit[i].thishit;
      hitclass1 = thisEvent.hit[i].hittype;
      Double_t wght = 1.;
      wght = thisEvent.hit[i].weight;

      //thisEvent.hit[i].Print();
      
      for (Int_t j = i+1; j<nclus; ++j) {
//        AliVCluster *clus2 = static_cast<AliVCluster*>(clusters->At(j));
//        if (!clus2)
//          continue;
//        if (!clus2->IsEMCAL())
//          continue;
//        if (clus2->E()<fMinE)
//          continue;
//        if (clus2->GetNCells()<fNminCells)
//          continue;
//        if (GetMaxCellEnergy(clus2)/clus2->E()<fMinErat)
//          continue;
//        if (clus2->Chi2()<fMinEcc) // eccentricity cut
//          continue;
//        clus2->GetMomentum(clusterVec2,vertex);
        //Bool_t bkeep2 = 1;
        clusterVec2 = thisEvent.hit[j].thishit;
        hitclass2 = thisEvent.hit[j].hittype;

        
        //Double_t ecorr2 = fcorrect->Eval(2.0*clusterVec2.E());
        //        Double_t ecorr2 = 1;
        //        TLorentzVector clusterVecCorr2(clusterVec2.Px()*ecorr2,clusterVec2.Py()*ecorr2,clusterVec2.Pz()*ecorr2,clusterVec2.E()*ecorr2);
        //clusterVecCorr2.SetPxPyPzE(clusterVec2.Px()*ecorr2,clusterVec2.Py()*ecorr2,clusterVec2.Pz()*ecorr2,clusterVec2.E()*ecorr1);
        // calculate radius
        Double_t d_phi = clusterVec1.Phi() - clusterVec2.Phi();
        Double_t d_eta = clusterVec1.Eta() - clusterVec2.Eta();
        Double_t d_r = sqrt(d_phi*d_phi + d_eta*d_eta);
        
        pionVec = clusterVec1 + clusterVec2;
        //        pionVecCorr = clusterVecCorr1 + clusterVecCorr2;
        //        if(pionVec.Pt()>1)
        //          fHPionSm->Fill(pionVec.M(),pionVecCorr.M());
        fHdr->Fill(d_r);
        
//        if(d_r < 0.005 && pionVec.Pt()<8)
//          continue;

        if(d_r < 0.01)
          continue;
        
        Double_t pionZgg = TMath::Abs(clusterVec1.E()-clusterVec2.E())/pionVec.E();
        Double_t pionOpeningAngle = clusterVec1.Angle(clusterVec2.Vect());
        if (pionZgg < fAsymMax1) {
          fHPionMggPt->Fill(pionVec.M(),pionVec.Pt());
          fHPionMggAsym->Fill(pionVec.M(),pionZgg);
          fHPionMggDgg->Fill(pionVec.M(),pionOpeningAngle);

          if(fMcMode){
            if(hitclass1 == 1 && hitclass2 == 1){
              if(thisEvent.hit[i].imo == thisEvent.hit[j].imo){
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
            else if(hitclass1 == 2 && hitclass2 == 2){
              if(thisEvent.hit[i].imo == thisEvent.hit[j].imo){
                fHPionInvMassesAdd2->Fill(pionVec.M(),pionVec.Pt(),wght);
              }
            }
            else if(hitclass1 >= 100 && hitclass2 >= 100){
              fHPionInvMassesSym->Fill(pionVec.M(),pionVec.Pt());
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
          }
          else{
            fHPionInvMassesSym->Fill(pionVec.M(),pionVec.Pt());
            fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
          }
        }
        // asymmetric
        else if(pionZgg < fAsymMax2){
          fHPionMggPt->Fill(pionVec.M(),pionVec.Pt());
          fHPionMggAsym->Fill(pionVec.M(),pionZgg);
          fHPionMggDgg->Fill(pionVec.M(),pionOpeningAngle);
          
          if(fMcMode){
            if(hitclass1 == 1 && hitclass2 == 1){
              if(thisEvent.hit[i].imo == thisEvent.hit[j].imo){
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
            else if(hitclass1 == 2 && hitclass2 == 2){
              if(thisEvent.hit[i].imo == thisEvent.hit[j].imo){
                fHPionInvMassesAdd2->Fill(pionVec.M(),pionVec.Pt(),wght);
              }
            }
            else if(hitclass1 >= 100 && hitclass2 >= 100){
              fHPionInvMasses->Fill(pionVec.M(),pionVec.Pt());
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
          }
          else{
            fHPionInvMasses->Fill(pionVec.M(),pionVec.Pt());
            fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
          }
        }
        
        // asymmetric
        else{
          fHPionMggPt->Fill(pionVec.M(),pionVec.Pt());
          fHPionMggAsym->Fill(pionVec.M(),pionZgg);
          fHPionMggDgg->Fill(pionVec.M(),pionOpeningAngle);
          
          if(fMcMode){
            if(hitclass1 == 1 && hitclass2 == 1){
              if(thisEvent.hit[i].imo == thisEvent.hit[j].imo){
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
            else if(hitclass1 == 2 && hitclass2 == 2){
              if(thisEvent.hit[i].imo == thisEvent.hit[j].imo){
                fHPionInvMassesAdd2->Fill(pionVec.M(),pionVec.Pt(),wght);
              }
            }
            else if(hitclass1 >= 100 && hitclass2 >= 100){
              fHPionInvMassesAsym->Fill(pionVec.M(),pionVec.Pt());
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
          }
          else{
            fHPionInvMassesAsym->Fill(pionVec.M(),pionVec.Pt());
            fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi());
          }
        }
        
        // look for matched tracks (JPsi?)
//        if(clus1->GetNTracksMatched()>=1 && clus2->GetNTracksMatched()>=1){
//          fHJPInvMasses->Fill(pionVec.M(),pionVec.Pt());
//        }
        
        // try to find only cluster pairs from primary pi0's
        // should I look for two photon clusters? ingredients: 2 photons, from same mother.
        
//        if(fMcMode){
          
//          AliMCEvent *mcEvent = MCEvent();
//          
//          mcEvent->PreReadAll();
//          
//          // MC labels
//          int ilabel1 = -1;
//          int ilabel2 = -1;
//          // highest contribution
//          ilabel1 = clus1->GetLabel();
//          ilabel2 = clus2->GetLabel();
//          
//          if(ilabel1 == -1 || ilabel2 == -1){
//            //AliWarning("label not set");
//            continue;
//          }
//          // get the two photons
//          AliMCParticle *g1 = static_cast<AliMCParticle*>(mcEvent->GetTrack(ilabel1));
//          AliMCParticle *g2 = static_cast<AliMCParticle*>(mcEvent->GetTrack(ilabel2));
//          if(g1->PdgCode()!=22 || g2->PdgCode()!=22){
//            //AliWarning("not two photons");
//            continue;
//          }
//          
//          Int_t imother1 = g1->GetMother();
//          Int_t imother2 = g2->GetMother();
//          
//          if(imother1 != imother2){
//            //AliWarning("different mother ID");
//            continue;
//          }
//          
//          AliMCParticle* mcMother = static_cast<AliMCParticle*>(mcEvent->GetTrack(imother1));
//          if(mcMother->PdgCode() != 111){
//            //AliWarning("mother not a pi0");
//            continue;
//          }
//          
//          if(!mcEvent->IsPhysicalPrimary(imother1)){
//            //AliWarning("not physical primary");
//            continue;
//          }
//          if (pionZgg < fAsymMax)
//            fHPrimPionInvMasses->Fill(pionVec.M(),pionVec.Pt());
//          else
//            fHPrimPionInvMassesAsym->Fill(pionVec.M(),pionVec.Pt());
          
//        }
        
      }
    }
 // }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::FillPionConv()
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
    
    Int_t nclus = clusters->GetEntries();
    for (Int_t i = 0; i<nclus; ++i) {
      AliVCluster *clus1 = static_cast<AliVCluster*>(clusters->At(i));
      if (!clus1)
        continue;
      if (!clus1->IsEMCAL())
        continue;
      if (clus1->E()<fMinE)
        continue;
      if (clus1->GetNCells()<fNminCells)
        continue;
      if (GetMaxCellEnergy(clus1)/clus1->E()<fMinErat)
        continue;
      if (clus1->Chi2()<fMinEcc) // eccentricity cut
        continue;
      clus1->GetMomentum(clusterVec1,vertex);
      //      Double_t ecorr = 1;
      //Double_t ecorr = fcorrect->Eval(2.0*clusterVec1.E());
      //      TLorentzVector clusterVecCorr1;
      //      clusterVecCorr1.SetPxPyPzE(clusterVec1.Px()*ecorr,clusterVec1.Py()*ecorr,clusterVec1.Pz()*ecorr,clusterVec1.E()*ecorr);
      
      if(fReaderGammas){
        for(Int_t j = 0; j < fReaderGammas->GetEntriesFast(); j++){
          AliAODConversionPhoton* PhotonCandidate = NULL;//(AliAODConversionPhoton*) fReaderGammas->At(j);
          if(!PhotonCandidate) continue;
          clusterVec2.SetPxPyPzE(PhotonCandidate->GetPx(),PhotonCandidate->GetPy(),PhotonCandidate->GetPz(),PhotonCandidate->GetPhotonP());
          if (clusterVec2.E()<fMinE)
            continue;
          pionVec = clusterVec1 + clusterVec2;
          Double_t pionZgg = TMath::Abs(clusterVec1.E()-clusterVec2.E())/pionVec.E();
          Double_t pionOpeningAngle = clusterVec1.Angle(clusterVec2.Vect());
          if (pionZgg < fAsymMax2) {
            fHPionEtaPhiConv->Fill(pionVec.Eta(),pionVec.Phi());
            fHPionMggPtConv->Fill(pionVec.M(),pionVec.Pt());
            fHPionMggAsymConv->Fill(pionVec.M(),pionZgg);
            fHPionMggDggConv->Fill(pionVec.M(),pionOpeningAngle);
            //            Int_t bin = fPtRanges->FindBin(pionVec.Pt());
            fHPionInvMassesConv->Fill(pionVec.M(),pionVec.Pt());
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::FillConvConv()
{
  // Fill histograms related to pions.
  
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  
  TLorentzVector clusterVec1;
  TLorentzVector clusterVec2;
  TLorentzVector pionVec;
  
  if(fReaderGammas){
    Int_t nclus = fReaderGammas->GetEntriesFast();
    for (Int_t i = 0; i<nclus; ++i) {
      AliAODConversionPhoton* PhotonCandidate1 = NULL;//(AliAODConversionPhoton*) fReaderGammas->At(i);
      if(!PhotonCandidate1) continue;
      clusterVec1.SetPxPyPzE(PhotonCandidate1->GetPx(),PhotonCandidate1->GetPy(),PhotonCandidate1->GetPz(),PhotonCandidate1->GetPhotonP());
      if (clusterVec1.E()<fMinE)
        continue;
      
      for(Int_t j = 0; j < fReaderGammas->GetEntriesFast(); j++){
        AliAODConversionPhoton* PhotonCandidate2 = NULL;//(AliAODConversionPhoton*) fReaderGammas->At(j);
        if(!PhotonCandidate2) continue;
        clusterVec2.SetPxPyPzE(PhotonCandidate2->GetPx(),PhotonCandidate2->GetPy(),PhotonCandidate2->GetPz(),PhotonCandidate2->GetPhotonP());
        if (clusterVec2.E()<fMinE)
          continue;
        pionVec = clusterVec1 + clusterVec2;
        Double_t pionZgg = TMath::Abs(clusterVec1.E()-clusterVec2.E())/pionVec.E();
        Double_t pionOpeningAngle = clusterVec1.Angle(clusterVec2.Vect());
        if (pionZgg < fAsymMax2) {
          fHPionEtaPhiConvConv->Fill(pionVec.Eta(),pionVec.Phi());
          fHPionMggPtConvConv->Fill(pionVec.M(),pionVec.Pt());
          fHPionMggAsymConvConv->Fill(pionVec.M(),pionZgg);
          fHPionMggDggConvConv->Fill(pionVec.M(),pionOpeningAngle);
          //            Int_t bin = fPtRanges->FindBin(pionVec.Pt());
          fHPionInvMassesConvConv->Fill(pionVec.M(),pionVec.Pt());
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
    
    Int_t nclus = clusters->GetEntries();
    for (Int_t i = 0; i<nclus; ++i) {
//      AliVCluster *clus1 = static_cast<AliVCluster*>(clusters->At(i));
//      if (!clus1)
//        continue;
//      if (!clus1->IsEMCAL())
//        continue;
//      if (clus1->E()<fMinE)
//        continue;
//      if (clus1->GetNCells()<fNminCells)
//        continue;
//      if (GetMaxCellEnergy(clus1)/clus1->E()<fMinErat)
//        continue;
//      if (clus1->Chi2()<fMinEcc) // eccentricity cut
//        continue;
//      clus1->GetMomentum(clusterVec1,vertex);

      //      Double_t ecorr = 1;
      //Double_t ecorr = fcorrect->Eval(2.0*clusterVec1.E());
      //      TLorentzVector clusterVecCorr1;
      //      clusterVecCorr1.SetPxPyPzE(clusterVec1.Px()*ecorr,clusterVec1.Py()*ecorr,clusterVec1.Pz()*ecorr,clusterVec1.E()*ecorr);
      
      clusterVec1 = thisEvent.hit[i].thishit;
      hitclass1 = thisEvent.hit[i].hittype;

      // loop over old events
      for (Int_t iOld=0;iOld<nEvt;iOld++){
        EmcEvent OldEvent = EmcEventList[MulClass][vtxClass][PtClass][iOld];
        Int_t nclusold = OldEvent.nHits;
        Double_t phirot = OldEvent.TrigPhi;
        Double_t thetarot = OldEvent.TrigTheta;
				if(fRotateMixed && PtClass > 0)
          fHMixRotation->Fill(phi0-phirot);
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
          
//          if(d_r < 0.005 && pionVec.Pt()<8)
//            continue;
          if(d_r < 0.01)
            continue;
          
          Double_t pionZgg = TMath::Abs(clusterVec1.E()-clusterVec2.E())/pionVec.E();
          //                    Double_t pionOpeningAngle = clusterVec1.Angle(clusterVec2.Vect());
          if (pionZgg < fAsymMax3) {
            /*
             fHPionEtaPhiMix->Fill(pionVec.Eta(),pionVec.Phi());
             fHPionMggPtMix->Fill(pionVec.M(),pionVec.Pt());
             fHPionMggAsymMix->Fill(pionVec.M(),pionZgg);
             fHPionMggDggMix->Fill(pionVec.M(),pionOpeningAngle);
             */
            //           Int_t bin = fPtRanges->FindBin(pionVec.Pt());
            if(fMcMode){
              if(hitclass1 == 1 && hitclass2 == 1){
                if(thisEvent.hit[i].pid != 22 || OldEvent.hit[j].pid != 22)
                  fHPionInvMassesMix1->Fill(pionVec.M(),pionVec.Pt());
              }
              else if(hitclass1 == 2 && hitclass2 == 2){
                if(thisEvent.hit[i].pid != 22 || OldEvent.hit[j].pid != 22)
                  fHPionInvMassesMix2->Fill(pionVec.M(),pionVec.Pt());
              }
              else{
                fHPionInvMassesMix->Fill(pionVec.M(),pionVec.Pt());
              }
            }
            else{ // not MC mode
              fHPionInvMassesMix->Fill(pionVec.M(),pionVec.Pt());
            }
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::FillMixConv(const Int_t MulClass, const Int_t vtxClass, const Int_t PtClass)
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

    Int_t nclus = thisEvent.nHits;
    
    // use clusters from "this" event and pair with v0 from old events
    // Int_t nclus = clusters->GetEntries();

    for (Int_t i = 0; i<nclus; ++i) {
//      AliVCluster *clus1 = static_cast<AliVCluster*>(clusters->At(i));
//      if (!clus1)
//        continue;
//      if (!clus1->IsEMCAL())
//        continue;
//      if (clus1->E()<fMinE)
//        continue;
//      if (clus1->GetNCells()<fNminCells)
//        continue;
//      if (GetMaxCellEnergy(clus1)/clus1->E()<fMinErat)
//        continue;
//      if (clus1->Chi2()<fMinEcc) // eccentricity cut
//        continue;
//      clus1->GetMomentum(clusterVec1,vertex);
      //      Double_t ecorr = 1;
      //Double_t ecorr = fcorrect->Eval(2.0*clusterVec1.E());
      //      TLorentzVector clusterVecCorr1;
      //      clusterVecCorr1.SetPxPyPzE(clusterVec1.Px()*ecorr,clusterVec1.Py()*ecorr,clusterVec1.Pz()*ecorr,clusterVec1.E()*ecorr);
      
      
      // loop over old events
      for (Int_t iOld=0;iOld<nEvt;iOld++){
        EmcEvent OldEvent = EmcEventList[MulClass][vtxClass][PtClass][iOld];
        Int_t nclusold = OldEvent.nV0Hits;
        for (Int_t j = 0; j<nclusold; ++j) {
          clusterVec2 = OldEvent.hitv0[j].thishit;
          pionVec = clusterVec1 + clusterVec2;
          Double_t pionZgg = TMath::Abs(clusterVec1.E()-clusterVec2.E())/pionVec.E();
          //                    Double_t pionOpeningAngle = clusterVec1.Angle(clusterVec2.Vect());
          if (pionZgg < fAsymMax2) {
            /*
             fHPionEtaPhiMix->Fill(pionVec.Eta(),pionVec.Phi());
             fHPionMggPtMix->Fill(pionVec.M(),pionVec.Pt());
             fHPionMggAsymMix->Fill(pionVec.M(),pionZgg);
             fHPionMggDggMix->Fill(pionVec.M(),pionOpeningAngle);
             */
            //            Int_t bin = fPtRanges->FindBin(pionVec.Pt());
            fHPionInvMassesConvMix->Fill(pionVec.M(),pionVec.Pt());
          }
        }
      }
    }
    // also loop over v0 from "this" event and pair with old clusters?
    if(fReaderGammas){
      for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
        AliAODConversionPhoton* PhotonCandidate = NULL;//(AliAODConversionPhoton*) fReaderGammas->At(i);
        if(!PhotonCandidate) continue;
        clusterVec1.SetPxPyPzE(PhotonCandidate->GetPx(),PhotonCandidate->GetPy(),PhotonCandidate->GetPz(),PhotonCandidate->GetPhotonP());
        if (clusterVec1.E()<fMinE)
          continue;
        
        // loop over old events
        for (Int_t iOld=0;iOld<nEvt;iOld++){
          EmcEvent OldEvent = EmcEventList[MulClass][vtxClass][PtClass][iOld];
          Int_t nclusold = OldEvent.nHits;
          for (Int_t j = 0; j<nclusold; ++j) {
            clusterVec2 = OldEvent.hit[j].thishit;
            pionVec = clusterVec1 + clusterVec2;
            Double_t pionZgg = TMath::Abs(clusterVec1.E()-clusterVec2.E())/pionVec.E();
            //                    Double_t pionOpeningAngle = clusterVec1.Angle(clusterVec2.Vect());
            if (pionZgg < fAsymMax2) {
              //              Int_t bin = fPtRanges->FindBin(pionVec.Pt());
              fHPionInvMassesConvMix->Fill(pionVec.M(),pionVec.Pt());
            }
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0Gamma::FillMixConvConv(const Int_t MulClass, const Int_t vtxClass, const Int_t PtClass)
{
  // Fill histograms related to pions.
  
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  
  if(fReaderGammas){
    
    TLorentzVector clusterVec1;
    TLorentzVector clusterVec2;
    TLorentzVector pionVec;
    
    // use v0 from "this" event and pair with v0 from old events
    Int_t nclus = fReaderGammas->GetEntriesFast();
    for (Int_t i = 0; i<nclus; ++i) {
      AliAODConversionPhoton* PhotonCandidate1 = NULL;//(AliAODConversionPhoton*) fReaderGammas->At(i);
      if(!PhotonCandidate1) continue;
      clusterVec1.SetPxPyPzE(PhotonCandidate1->GetPx(),PhotonCandidate1->GetPy(),PhotonCandidate1->GetPz(),PhotonCandidate1->GetPhotonP());
      if (clusterVec1.E()<fMinE)
        continue;
      
      // loop over old events
      for (Int_t iOld=0;iOld<nEvt;iOld++){
        EmcEvent OldEvent = EmcEventList[MulClass][vtxClass][PtClass][iOld];
        Int_t nclusold = OldEvent.nV0Hits;
        for (Int_t j = 0; j<nclusold; ++j) {
          clusterVec2 = OldEvent.hitv0[j].thishit;
          pionVec = clusterVec1 + clusterVec2;
          Double_t pionZgg = TMath::Abs(clusterVec1.E()-clusterVec2.E())/pionVec.E();
          //                    Double_t pionOpeningAngle = clusterVec1.Angle(clusterVec2.Vect());
          if (pionZgg < fAsymMax2) {
            /*
             fHPionEtaPhiMix->Fill(pionVec.Eta(),pionVec.Phi());
             fHPionMggPtMix->Fill(pionVec.M(),pionVec.Pt());
             fHPionMggAsymMix->Fill(pionVec.M(),pionZgg);
             fHPionMggDggMix->Fill(pionVec.M(),pionOpeningAngle);
             */
            //            Int_t bin = fPtRanges->FindBin(pionVec.Pt());
            fHPionInvMassesConvConvMix->Fill(pionVec.M(),pionVec.Pt());
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
    Int_t d = amc->GetDaughter(i);
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
  Int_t d1 = p->GetFirstDaughter();
  Int_t d2 = p->GetLastDaughter();
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
    Int_t d = amc->GetDaughter(i);
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
  
  Int_t d1 = p->GetFirstDaughter();
  Int_t d2 = p->GetLastDaughter();
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
  const int MultCut[5] = {3, 7, 14, 30, 9999};
  
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


  Int_t nv0 = 0;
  if(fReaderGammas)
    nv0 = fReaderGammas->GetEntriesFast();

// add new event
  //	evt.SetGlobalInfo(MulClass, vtxClass, nclus, nv0, phitrig, thetatrig);
//	evt.SetGlobalInfo(nclus, nv0, phitrig, thetatrig);
  
//  if (clusters) {
//    int ncl = 0;
//    for (Int_t i = 0; i<nclus; ++i) {
//      TLorentzVector clusterVec1;
//      AliVCluster *clus1 = static_cast<AliVCluster*>(clusters->At(i));
//      if (!clus1)
//        continue;
//      if (!clus1->IsEMCAL())
//        continue;
//      if (clus1->E()<fMinE)
//        continue;
//      if (clus1->GetNCells()<fNminCells)
//        continue;
//      if (GetMaxCellEnergy(clus1)/clus1->E()<fMinErat)
//        continue;
//      if (clus1->Chi2()<fMinEcc) // eccentricity cut
//        continue;
//      clus1->GetMomentum(clusterVec1,vertex);
//      //      Double_t ecorr = 1;
//      //Double_t ecorr = fcorrect->Eval(2.0*clusterVec1.E());
//      //      TLorentzVector clusterVecCorr1;
//      //      clusterVecCorr1.SetPxPyPzE(clusterVec1.Px()*ecorr,clusterVec1.Py()*ecorr,clusterVec1.Pz()*ecorr,clusterVec1.E()*ecorr);
//      
//      evt.hit[i].thishit=clusterVec1;
//      ncl++;
//    }
//    EmcEventList[MulClass][vtxClass][PtClass][iEvent] = evt;
//    iEvent++;
//    //        cout << "added event with " << ncl << " clusters, in mulclass " << MulClass << ", zclass " << vtxClass << endl;
//  }
  
  //cout << Form("%d, %d, %d, %d",MulClass,vtxClass,PtClass,iEvent) << endl;
  thisEvent.SetGlobalInfo(nclus,nv0,phitrig,thetatrig);
  EmcEventList[MulClass][vtxClass][PtClass][iEvent] = thisEvent;

//  for(int i=0;i<nclus;i++){
//    EmcEventList[MulClass][vtxClass][PtClass][iEvent].hit[i].hittype = thisEvent.hit[i].hittype;
//    EmcEventList[MulClass][vtxClass][PtClass][iEvent].hit[i].imo = thisEvent.hit[i].imo;
//    EmcEventList[MulClass][vtxClass][PtClass][iEvent].hit[i].pid = thisEvent.hit[i].pid;
//    EmcEventList[MulClass][vtxClass][PtClass][iEvent].hit[i].weight = thisEvent.hit[i].weight;
//    EmcEventList[MulClass][vtxClass][PtClass][iEvent].hit[i].thishit = thisEvent.hit[i].thishit;
//  }

  iEvent++;
  if(fReaderGammas){
    for(Int_t j = 0; j < fReaderGammas->GetEntriesFast(); j++){
      AliAODConversionPhoton* PhotonCandidate = NULL;//(AliAODConversionPhoton*) fReaderGammas->At(j);
      if(!PhotonCandidate) continue;
      TLorentzVector clusterVec2;
      clusterVec2.SetPxPyPzE(PhotonCandidate->GetPx(),PhotonCandidate->GetPy(),PhotonCandidate->GetPz(),PhotonCandidate->GetPhotonP());
      if (clusterVec2.E()<fMinE)
        continue;
      evt.hitv0[j].thishit=clusterVec2;
    }
  }
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
Double_t AliAnalysisTaskEMCALPi0Gamma::CalcWeight(Double_t pt, Double_t eta, Int_t i){
  Double_t weight = 1.;
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

/*
 //__________________________________________________________________________________________________
 void AliStaCluster::GetMom(TLorentzVector& p, Double_t *vertex)
 {
 // Calculate momentum.
 
 TVector3 pos;
 pos.SetPtEtaPhi(fR,fEta,fPhi);
 
 if(vertex){ //calculate direction relative to vertex
 pos -= vertex;
 }
 
 Double_t r = pos.Mag();
 p.SetPxPyPzE(fE*pos.x()/r, fE*pos.y()/r, fE*pos.z()/r, fE);
 }
 
 //__________________________________________________________________________________________________
 void AliStaCluster::GetMom(TLorentzVector& p, AliStaVertex *vertex)
 {
 // Calculate momentum.
 
 Double_t v[3] = {0,0,0};
 if (vertex) {
 v[0] = vertex->fVx;
 v[1] = vertex->fVy;
 v[2] = vertex->fVz;
 }
 GetMom(p, v);
 }
 */

//__________________________________________________________________________________________________
EmcEvent::EmcEvent()
: nHits(0),
  nV0Hits(0),
  TrigPhi(0),
  TrigTheta(0)
{
  nHits = 0;
  nV0Hits = 0;
  TrigPhi = 0;
  TrigTheta = 0;
}

//__________________________________________________________________________________________________
EmcEvent::EmcEvent(const EmcEvent &obj)
: nHits(0),
  nV0Hits(0),
  TrigPhi(0),
  TrigTheta(0)
{
  nHits = obj.nHits;
  nV0Hits = obj.nV0Hits;
  TrigPhi = obj.TrigPhi;
  TrigTheta = obj.TrigTheta;
  // copy all hits
  for(int i=0;i<nHits;i++){
    hit[i].hittype = obj.hit[i].hittype;
    hit[i].imo = obj.hit[i].imo;
    hit[i].pid = obj.hit[i].pid;
    hit[i].weight = obj.hit[i].weight;
    hit[i].thishit = obj.hit[i].thishit;
  }
}


//__________________________________________________________________________________________________
void EmcEvent::SetGlobalInfo(const Int_t& Size, const Int_t& V0Size, const Float_t& phiTrig, const Float_t& thetaTrig)
{
  //    fCenPercent = centPer;
  //    fVtx = vtxPos;
  nHits = Size;
  nV0Hits = V0Size;
	TrigPhi = phiTrig;
  TrigTheta = thetaTrig;
}

//__________________________________________________________________________________________________
void EmcEvent::Reset()
{
  for(int i=0;i<nHits;i++){
    hit[i].hittype = 0;
    hit[i].imo = 0;
    hit[i].pid = 0;
    hit[i].weight = 1;
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

//__________________________________________________________________________________________________
V0Hit::V0Hit()
: thishit()
{
}


