// $Id$

#include "AliAnalysisTaskEMCALPi0PbPb.h"
#include <TAxis.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TGeoGlobalMagField.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TProfile.h>
#include <TRegexp.h>
#include <TString.h>
#include <TVector2.h>
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliCDBManager.h"
#include "AliCentrality.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDCaloTrigger.h"
#include "AliESDEvent.h"
#include "AliESDUtils.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
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

ClassImp(AliAnalysisTaskEMCALPi0PbPb)

//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::AliAnalysisTaskEMCALPi0PbPb() 
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
    fAsymMax(1),
    fNminCells(1),
    fMinE(0.100),
    fMinErat(0),
    fMinEcc(0),
    fGeoName("EMCAL_FIRSTYEARV1"),
    fMinNClusPerTr(50),
    fIsoDist(0.2),
    fTrClassNames(""),
    fTrCuts(0),
    fPrimTrCuts(0),
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
    fDoPSel(kTRUE),
    fIsGeoMatsSet(0),
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
    fHColuRow(0x0),
    fHColuRowE(0x0),
    fHCellMult(0x0),
    fHCellE(0x0),
    fHCellH(0x0),
    fHCellM(0x0),
    fHCellM2(0x0),
    fHCellFreqNoCut(0x0),
    fHCellFreqCut100M(0x0),
    fHCellFreqCut300M(0x0),
    fHCellFreqE(0x0),
    fHCellCheckE(0x0),
    fHClustEccentricity(0),
    fHClustEtaPhi(0x0),
    fHClustEnergyPt(0x0),
    fHClustEnergySigma(0x0),
    fHClustSigmaSigma(0x0),
    fHClustNCellEnergyRatio(0x0),
    fHClustEnergyNCell(0x0),
    fHPrimTrackPt(0x0),
    fHPrimTrackEta(0x0),
    fHPrimTrackPhi(0x0),
    fHMatchDr(0x0),
    fHMatchDz(0x0),
    fHMatchEp(0x0),
    fHPionEtaPhi(0x0),
    fHPionMggPt(0x0),
    fHPionMggAsym(0x0),
    fHPionMggDgg(0x0)
{
  // Constructor.
}

//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::AliAnalysisTaskEMCALPi0PbPb(const char *name) 
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
    fAsymMax(1),
    fNminCells(1),
    fMinE(0.100),
    fMinErat(0),
    fMinEcc(0),
    fGeoName("EMCAL_FIRSTYEARV1"),
    fMinNClusPerTr(50),
    fIsoDist(0.2),
    fTrClassNames(""),
    fTrCuts(0),
    fPrimTrCuts(0),
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
    fDoPSel(kTRUE),
    fIsGeoMatsSet(0),
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
    fHColuRow(0x0),
    fHColuRowE(0x0),
    fHCellMult(0x0),
    fHCellE(0x0),
    fHCellH(0x0),
    fHCellM(0x0),
    fHCellM2(0x0),
    fHCellFreqNoCut(0x0),
    fHCellFreqCut100M(0x0),
    fHCellFreqCut300M(0x0),
    fHCellFreqE(0x0),
    fHCellCheckE(0x0),
    fHClustEccentricity(0),
    fHClustEtaPhi(0x0),
    fHClustEnergyPt(0x0),
    fHClustEnergySigma(0x0),
    fHClustSigmaSigma(0x0),
    fHClustNCellEnergyRatio(0x0),
    fHClustEnergyNCell(0x0),
    fHPrimTrackPt(0x0),
    fHPrimTrackEta(0x0),
    fHPrimTrackPhi(0x0),
    fHMatchDr(0x0),
    fHMatchDz(0x0),
    fHMatchEp(0x0),
    fHPionEtaPhi(0x0),
    fHPionMggPt(0x0),
    fHPionMggAsym(0x0),
    fHPionMggDgg(0x0)
{
  // Constructor.

  DefineOutput(1, TList::Class());
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.,SPDVertex.,TPCVertex.,EMCALCells.,Tracks,EMCALTrigger.,SPDPileupVertices,TrkPileupVertices "
               "AOD:header,vertices,emcalCells,tracks";
}

//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::~AliAnalysisTaskEMCALPi0PbPb()
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
  delete [] fHColuRow;
  delete [] fHColuRowE;
  delete [] fHCellMult;
  delete [] fHCellFreqNoCut;
  delete [] fHCellFreqCut100M;
  delete [] fHCellFreqCut300M;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::UserCreateOutputObjects()
{
  // Create user objects here.

  cout << "AliAnalysisTaskEMCALPi0PbPb: Input settings" << endl;
  cout << " fCentVar:       " << fCentVar << endl;
  cout << " fCentFrom:      " << fCentFrom << endl;
  cout << " fCentTo:        " << fCentTo << endl;
  cout << " fVtxZMin:       " << fVtxZMin << endl;
  cout << " fVtxZMax:       " << fVtxZMax << endl;
  cout << " fUseQualFlag:   " << fUseQualFlag << endl;
  cout << " fClusName:      \"" << fClusName << "\"" << endl;
  cout << " fDoNtuple:      " << fDoNtuple << endl;
  cout << " fDoAfterburner: " << fDoAfterburner << endl;
  cout << " fAsymMax:       " << fAsymMax << endl;
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
  fOutput->SetOwner(1);
  fSelTracks = new TObjArray;
  fSelPrimTracks = new TObjArray;

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
        TClass::GetClass("AliStaCluster")->IgnoreTObjectStreamer();
      fClusters = new TClonesArray("AliStaCluster");
      fNtuple->Branch("clusters", &fClusters, 8*16*1024, 99);
      if (TClass::GetClass("AliStaTrigger"))
        TClass::GetClass("AliStaTrigger")->IgnoreTObjectStreamer();
      fTriggers = new TClonesArray("AliStaTrigger");
      fNtuple->Branch("l0prim", &fTriggers, 16*1024, 99);
      if (fMcMode||fEmbedMode) {
        if (TClass::GetClass("AliStaPart"))
          TClass::GetClass("AliStaPart")->IgnoreTObjectStreamer();
        fMcParts = new TClonesArray("AliStaPart");
        fNtuple->Branch("mcparts", &fMcParts, 8*16*1024, 99);
      }
    }  
  }

  AliEMCALEMCGeometry *emc = fGeom->GetEMCGeometry();
  Double_t phimin = emc->GetArm1PhiMin()*TMath::DegToRad();
  Double_t phimax = emc->GetArm1PhiMax()*TMath::DegToRad();

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

  // histograms for cells
  Int_t nsm = fGeom->GetEMCGeometry()->GetNumberOfSuperModules();
  fHColuRow   = new TH2*[nsm];
  fHColuRowE  = new TH2*[nsm];
  fHCellMult  = new TH1*[nsm];
  for (Int_t i = 0; i<nsm; ++i) {
    fHColuRow[i] = new TH2F(Form("hColRow_Mod%d",i),"",48,0,48,24,0.,24);
    fHColuRow[i]->SetTitle(Form("Module %d: Occupancy", i));
    fHColuRow[i]->SetXTitle("col (i#eta)");
    fHColuRow[i]->SetYTitle("row (i#phi)");
    fHColuRowE[i] = new TH2F(Form("hColRowE_Mod%d", i),"",48,0,48,24,0,24);
    fHColuRowE[i]->SetTitle(Form("Module %d: Cell energy",i));
    fHColuRowE[i]->SetXTitle("col (i#eta)");
    fHColuRowE[i]->SetYTitle("row (i#phi)");
    fHCellMult[i] = new TH1F(Form("hCellMult_Mod%d",i),"",1000,0,1000); 
    fHCellMult[i]->SetTitle(Form("Module %d: Cell multiplicity",i));
    fHCellMult[i]->SetXTitle("# of cells");
    fOutput->Add(fHColuRow[i]);
    fOutput->Add(fHColuRowE[i]);
    fOutput->Add(fHCellMult[i]);
  }  
  fHCellE = new TH1F("hCellE","",250,0.,25.);
  fHCellE->SetXTitle("E_{cell} [GeV]");
  fOutput->Add(fHCellE);
  fHCellH = new TH1F ("hCellHighestE","",250,0.,25.);
  fHCellH->SetXTitle("E^{max}_{cell} [GeV]");
  fOutput->Add(fHCellH);
  fHCellM = new TH1F ("hCellMeanEperHitCell","",250,0.,2.5);
  fHCellM->SetXTitle("#LT E_{cell}#GT [GeV]");
  fOutput->Add(fHCellM);
  fHCellM2 = new TH1F ("hCellMeanEperAllCells","",250,0.,1);
  fHCellM2->SetXTitle("1/N_{cells} #Sigma E_{cell} [GeV]");
  fOutput->Add(fHCellM2);

  fHCellFreqNoCut   = new TH1*[nsm];
  fHCellFreqCut100M = new TH1*[nsm];
  fHCellFreqCut300M = new TH1*[nsm];
  fHCellFreqE       = new TH1*[nsm];
  for (Int_t i = 0; i<nsm; ++i){
    Double_t lbin = i*24*48-0.5;
    Double_t hbin = lbin+24*48;
    fHCellFreqNoCut[i]   = new TH1F(Form("hCellFreqNoCut_SM%d",i),    
                                    Form("Frequency SM%d (no cut);id;#",i),  1152, lbin, hbin);
    fHCellFreqCut100M[i] = new TH1F(Form("hCellFreqCut100M_SM%d",i), 
                                    Form("Frequency SM%d (>0.1GeV);id;#",i), 1152, lbin, hbin);
    fHCellFreqCut300M[i] = new TH1F(Form("hCellFreqCut300M_SM%d",i), 
                                    Form("Frequency SM%d (>0.3GeV);id;#",i), 1152, lbin, hbin);
    fHCellFreqE[i]       = new TH1F(Form("hCellFreqE_SM%d",i), 
                                    Form("Frequency SM%d (E weighted);id;#",i), 1152, lbin, hbin);
    fOutput->Add(fHCellFreqNoCut[i]);
    fOutput->Add(fHCellFreqCut100M[i]);
    fOutput->Add(fHCellFreqCut300M[i]);
    fOutput->Add(fHCellFreqE[i]);
  }
  if (!fMarkCells.IsNull()) {
    fHCellCheckE = new TH1*[24*48*nsm];
    memset(fHCellCheckE,0,24*48*nsm*sizeof(TH1*));
    TObjArray *cells = fMarkCells.Tokenize(" ");
    Int_t n = cells->GetEntries();
    Int_t *tcs = new Int_t[n];
    for (Int_t i=0;i<n;++i) {
      TString name(cells->At(i)->GetName());
      tcs[i]=name.Atoi();
    }
    for (Int_t i = 0; i<n; ++i) {
      Int_t c=tcs[i];
      if (c<24*48*nsm) {
        fHCellCheckE[i] = new TH1F(Form("hCellE_id%d",c), Form("Cell %d;E [GeV/c];#",c), 1000, 0, 10);
        fOutput->Add(fHCellCheckE[i]);
      }
    }
    delete cells;
    delete [] tcs;
  }

  // histograms for clusters
  if (!fTrainMode) {
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

  // histograms for primary tracks
  fHPrimTrackPt = new TH1F("hPrimTrackPt",";p_{T} [GeV/c]",500,0,50);
  fOutput->Add(fHPrimTrackPt);
  fHPrimTrackEta = new TH1F("hPrimTrackEta",";#eta",40,-2,2);
  fOutput->Add(fHPrimTrackEta);
  fHPrimTrackPhi = new TH1F("hPrimTrackPhi",";#varPhi [rad]",63,0,6.3);
  fOutput->Add(fHPrimTrackPhi);
  // histograms for track matching
  fHMatchDr = new TH1F("hMatchDrDist",";dR [cm]",500,0,200);
  fOutput->Add(fHMatchDr);
  fHMatchDz = new TH1F("hMatchDzDist",";dZ [cm]",500,-100,100);
  fOutput->Add(fHMatchDz);
  fHMatchEp = new TH1F("hMatchEpDist",";E/p",100,0,10);
  fOutput->Add(fHMatchEp);

  // histograms for pion candidates
  if (!fTrainMode) {
    fHPionEtaPhi = new TH2F("hPionEtaPhi","",100,-0.8,0.8,100*nsm,phimin,phimax);
    fHPionEtaPhi->SetXTitle("#eta_{#gamma#gamma}");
    fHPionEtaPhi->SetYTitle("#varphi_{#gamma#gamma}");
    fOutput->Add(fHPionEtaPhi);
    fHPionMggPt = new TH2F("hPionMggPt","",1000,0,2,100,0,20.0);
    fHPionMggPt->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionMggPt->SetYTitle("p_{T}^{#gamma#gamma} [GeV/c]");
    fOutput->Add(fHPionMggPt);
    fHPionMggAsym = new TH2F("hPionMggAsym","",1000,0,2,100,0,1);
    fHPionMggAsym->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionMggAsym->SetYTitle("Z_{#gamma#gamma} [GeV]");
    fOutput->Add(fHPionMggAsym);
    fHPionMggDgg = new TH2F("hPionMggDgg","",1000,0,2,100,0,10);
    fHPionMggDgg->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    fHPionMggDgg->SetYTitle("opening angle [grad]");
    fOutput->Add(fHPionMggDgg);
    const Int_t nbins = 20; 
    Double_t xbins[nbins] = {0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,12.5,15,20,25,50};
    fPtRanges = new TAxis(nbins-1,xbins);
    for (Int_t i = 0; i<=nbins; ++i) {
      fHPionInvMasses[i] = new TH1F(Form("hPionInvMass%d",i),"",1000,0,2);
      fHPionInvMasses[i]->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
      if (i==0)
        fHPionInvMasses[i]->SetTitle(Form("0 < p_{T}^{#gamma#gamma} <%.1f",xbins[0]));
      else if (i==nbins)
        fHPionInvMasses[i]->SetTitle(Form("p_{T}^{#gamma#gamma} > 50"));
      else 
        fHPionInvMasses[i]->SetTitle(Form("%.1f < p_{T}^{#gamma#gamma} <%.1f",xbins[i-1],xbins[i]));
      fOutput->Add(fHPionInvMasses[i]);
    }
  }

  TH1::SetDefaultSumw2(th1);
  TH2::SetDefaultSumw2(th2);
  PostData(1, fOutput); 
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::UserExec(Option_t *) 
{
  // Called for each event.

  if (!InputEvent())
    return;

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
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
  } else {
    fAodEv = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAodEv) {
      AliFatal("Neither ESD nor AOD event found");
      return;
    }
    am->LoadBranch("header");
    offtrigger =  fAodEv->GetHeader()->GetOfflineTrigger();
  }
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
      else 
        geom = fAodEv->GetHeader()->GetEMCALMatrix(i);
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
    } else {
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

  CalcMcInfo();
  CalcCaloTriggers();
  CalcPrimTracks();
  CalcTracks();
  CalcClusterProps();

  FillCellHists();
  if (!fTrainMode) {
    FillClusHists();
    FillPionHists();
    FillOtherHists();
  }
  FillMcHists();
  FillNtuple();

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

  PostData(1, fOutput);
}      

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::Terminate(Option_t *) 
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
void AliAnalysisTaskEMCALPi0PbPb::CalcCaloTriggers()
{
  // Calculate triggers

  if (fAodEv)
    return; // information not available in AOD

  if (!fTriggers)
    return;

  fTriggers->Clear();

  if (fTrigName.Length()<=0)
    return;

  TClonesArray *arr = dynamic_cast<TClonesArray*>(fEsdEv->FindListObject(fTrigName));
  if (!arr) {
    AliError(Form("Could not get array with name %s", fTrigName.Data()));
    return;
  }

  Int_t nNumberOfCaloClusters = arr->GetEntries();
  for(Int_t j = 0, ntrigs = 0; j < nNumberOfCaloClusters; ++j) {
    AliVCluster *cl = dynamic_cast<AliVCluster*>(arr->At(j));
    if (!cl)
      continue;
    if (!cl->IsEMCAL())
      continue;
    if (cl->E()<1)
      continue;
    AliStaTrigger *trignew = static_cast<AliStaTrigger*>(fTriggers->New(ntrigs++));
    Float_t pos[3] = {0,0,0};
    cl->GetPosition(pos);  
    TVector3 vpos(pos); 
    trignew->fE       = cl->E();
    trignew->fEta     = vpos.Eta();
    trignew->fPhi     = vpos.Phi();
    Short_t  id    = -1;
    GetMaxCellEnergy(cl, id);
    trignew->fIdMax    = id;
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::CalcClusterProps()
{
  // Calculate cluster properties

  if (!fClusters)
    return;

  fClusters->Clear();

  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;
  if (!cells)
    return;

  TObjArray *clusters = fEsdClusters;
  if (!clusters)
    clusters = fAodClusters;
  if (!clusters)
    return;

  Int_t ncells = cells->GetNumberOfCells();
  Int_t nclus  = clusters->GetEntries();
  Int_t ntrks  = fSelTracks->GetEntries();
  Int_t btracks[6][ntrks];
  memset(btracks,0,sizeof(Int_t)*6*ntrks);

  std::map<Short_t,Short_t> map;
  for (Short_t pos=0;pos<ncells;++pos) {
    Short_t id = cells->GetCellNumber(pos);
    map[id]=pos;
  }

  TObjArray filtMcParts;
  if (fMcParts) {
    Int_t nmc = fMcParts->GetEntries();
    for (Int_t i=0; i<nmc; ++i) {
      AliStaPart *pa = static_cast<AliStaPart*>(fMcParts->At(i));
      if (pa->OnEmcal())
        filtMcParts.Add(pa);
    }
  }

  for(Int_t i=0, ncl=0; i<nclus; ++i) {
    AliVCluster *clus = static_cast<AliVCluster*>(clusters->At(i));

    if (!clus)
      continue;
    if (!clus->IsEMCAL())
      continue;
    if (clus->E()<fMinE)
      continue;

    Float_t clsPos[3] = {0};
    clus->GetPosition(clsPos);
    TVector3 clsVec(clsPos);
    Double_t vertex[3] = {0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    TLorentzVector clusterVec;
    clus->GetMomentum(clusterVec,vertex);
    Double_t clsEta = clusterVec.Eta();

    AliStaCluster *cl = static_cast<AliStaCluster*>(fClusters->New(ncl++));
    cl->fE        = clus->E();
    cl->fR        = clsVec.Perp();
    cl->fEta      = clsVec.Eta();
    cl->fPhi      = clsVec.Phi();
    cl->fN        = clus->GetNCells();
    cl->fN1       = GetNCells(clus,0.100);
    cl->fN3       = GetNCells(clus,0.300);
    Short_t id    = -1;
    Double_t emax = GetMaxCellEnergy(clus, id);
    cl->fIdMax    = id;
    cl->fSM       = fGeom->GetSuperModuleNumber(id);
    cl->fEmax     = emax;
    Short_t id2   = -1;
    cl->fE2max    = GetSecondMaxCellEnergy(clus,id2);
    cl->fTmax     = cells->GetCellTime(id);
    if (clus->GetDistanceToBadChannel()<10000)
      cl->fDbc    = clus->GetDistanceToBadChannel();
    if (!TMath::IsNaN(clus->GetDispersion()))
      cl->fDisp   = clus->GetDispersion();
    if (!TMath::IsNaN(clus->GetM20()))
      cl->fM20    = clus->GetM20();
    if (!TMath::IsNaN(clus->GetM02()))
      cl->fM02    = clus->GetM02();
    Double_t maxAxis = -1, minAxis = -1;
    GetSigma(clus,maxAxis,minAxis);
    clus->SetTOF(maxAxis);     // store sigma in TOF for later plotting
    cl->fSig      = maxAxis;
    Double_t sEtaEta = -1;
    Double_t sPhiPhi = -1;
    GetSigmaEtaEta(clus, sEtaEta, sPhiPhi);
    cl->fSigEtaEta = sEtaEta;
    cl->fSigPhiPhi = sPhiPhi;
    Double_t clusterEcc = -1;
    if (maxAxis > 0)
      clusterEcc = TMath::Sqrt(1.0 - minAxis*minAxis/(maxAxis*maxAxis));
    clus->SetChi2(clusterEcc); // store ecc in chi2 for later plotting
    cl->fEcc        = clusterEcc;
    cl->fTrIso      = GetTrackIsolation(clusterVec.Eta(),clusterVec.Phi(),fIsoDist);
    cl->fTrIso1     = GetTrackIsolation(clusterVec.Eta(),clusterVec.Phi(),fIsoDist, 1);
    cl->fTrIso2     = GetTrackIsolation(clusterVec.Eta(),clusterVec.Phi(),fIsoDist, 2);
    cl->fTrIsoD1    = GetTrackIsolation(clusterVec.Eta(),clusterVec.Phi(),fIsoDist-0.1);
    cl->fTrIso1D1   = GetTrackIsolation(clusterVec.Eta(),clusterVec.Phi(),fIsoDist-0.1, 1);
    cl->fTrIso2D1   = GetTrackIsolation(clusterVec.Eta(),clusterVec.Phi(),fIsoDist-0.1, 2);
    cl->fTrIsoD3    = GetTrackIsolation(clusterVec.Eta(),clusterVec.Phi(),fIsoDist+0.1);
    cl->fTrIso1D3   = GetTrackIsolation(clusterVec.Eta(),clusterVec.Phi(),fIsoDist+0.1, 1);
    cl->fTrIso2D3   = GetTrackIsolation(clusterVec.Eta(),clusterVec.Phi(),fIsoDist+0.1, 2);
    cl->fTrIsoD4    = GetTrackIsolation(clusterVec.Eta(),clusterVec.Phi(),fIsoDist+0.2);
    cl->fTrIso1D4   = GetTrackIsolation(clusterVec.Eta(),clusterVec.Phi(),fIsoDist+0.2, 1);
    cl->fTrIso2D4   = GetTrackIsolation(clusterVec.Eta(),clusterVec.Phi(),fIsoDist+0.2, 2);
    cl->fTrIsoStrip = GetTrackIsoStrip(clusterVec.Eta(), clusterVec.Phi());
    cl->fCeCore     = GetCellIsolation(clsVec.Eta(),clsVec.Phi(),0.05);
    cl->fCeIso      = GetCellIsolation(clsVec.Eta(),clsVec.Phi(),fIsoDist);
    cl->fCeIso1     = GetCellIsolation(clsVec.Eta(),clsVec.Phi(),0.10);
    cl->fCeIso3     = GetCellIsolation(clsVec.Eta(),clsVec.Phi(),0.30);
    cl->fCeIso4     = GetCellIsolation(clsVec.Eta(),clsVec.Phi(),0.40);
    cl->fCeIso4x4   = GetCellIsoNxM(clsVec.Eta(),clsVec.Phi(), 4, 4);
    cl->fCeIso5x5   = GetCellIsoNxM(clsVec.Eta(),clsVec.Phi(), 5, 5);
    cl->fCeIso3x22  = GetCellIsoNxM(clsVec.Eta(),clsVec.Phi(), 3, 22);
    cl->fIsShared   = IsShared(clus);
    cl->fTrigId     = -1;
    cl->fTrigE      = 0;
    if (fTriggers) {
      Int_t ntrig = fTriggers->GetEntries();
      for (Int_t j = 0; j<ntrig; ++j) {
        AliStaTrigger *sta = static_cast<AliStaTrigger*>(fTriggers->At(j));
        if (!sta)
          continue;
        Short_t idmax = sta->fIdMax;
        Bool_t inc = IsIdPartOfCluster(clus, idmax);
        if (inc) {
          cl->fTrigId     = j;
          cl->fTrigE      = sta->fE;
          break;
        }
      }
    }

    // track matching
    Double_t mind2 = 1e10;
    for(Int_t j = 0; j<ntrks; ++j) {
      AliVTrack *track = static_cast<AliVTrack*>(fSelTracks->At(j));
      if (!track)
        continue;

      if (TMath::Abs(clsEta-track->Eta())>0.5)
        continue;

      TVector3 vec(clsPos);
      Int_t index =  (Int_t)(vec.Phi()*TMath::RadToDeg()/20);
      if (btracks[index-4][j]) {
        continue;
      }

      Float_t tmpR=-1, tmpZ=-1;
      Double_t dedx = 0;
      if (!fDoTrMatGeom) {
        AliExternalTrackParam *tParam = 0;
        if (fEsdEv) {
          AliESDtrack *esdTrack = static_cast<AliESDtrack*>(track);
          tParam = new AliExternalTrackParam(*esdTrack->GetTPCInnerParam());
	  dedx = esdTrack->GetTPCsignal();
        } else 
          tParam = new AliExternalTrackParam(track);

        Double_t bfield[3] = {0};
        track->GetBxByBz(bfield);
        Double_t alpha = (index+0.5)*20*TMath::DegToRad();
        vec.RotateZ(-alpha);   //Rotate the cluster to the local extrapolation coordinate system
        tParam->Rotate(alpha); //Rotate the track to the same local extrapolation system
        Bool_t ret = tParam->PropagateToBxByBz(vec.X(), bfield);
        if (!ret) {
          btracks[index-4][j]=1;
          delete tParam;
          continue;
        }
        Double_t trkPos[3] = {0};
        tParam->GetXYZ(trkPos); //Get the extrapolated global position
        tmpR = TMath::Sqrt( TMath::Power(clsPos[0]-trkPos[0],2) + 
                            TMath::Power(clsPos[1]-trkPos[1],2) +
                            TMath::Power(clsPos[2]-trkPos[2],2) );
        tmpZ = clsPos[2]-trkPos[2];
        delete tParam;
      } else {
        if (TMath::Abs(clsEta-track->Eta())>fIsoDist)
          continue;
        AliExternalTrackParam tParam(track);
        if (!fReco->ExtrapolateTrackToCluster(&tParam, clus, tmpR, tmpZ))
          continue;
      }

      Double_t d2 = tmpR;
      if (mind2>d2) {
        mind2=d2;
        cl->fTrDz   = tmpZ;
        cl->fTrDr   = TMath::Sqrt(tmpR*tmpR-tmpZ*tmpZ);
        cl->fTrEp   = clus->E()/track->P();
	cl->fTrDedx = dedx;
        cl->fIsTrackM = 1;
      }
    }
    
    if (cl->fIsTrackM) {
      fHMatchDr->Fill(cl->fTrDr);
      fHMatchDz->Fill(cl->fTrDz);
      fHMatchEp->Fill(cl->fTrEp);
    }

    //mc matching
    if (fMcParts) {
      Int_t nmc = filtMcParts.GetEntries();
      Double_t diffR2 = 1e9;
      AliStaPart *msta = 0;
      for (Int_t j=0; j<nmc; ++j) {
        AliStaPart *pa = static_cast<AliStaPart*>(filtMcParts.At(j));
        Double_t t1=clsVec.Eta()-pa->fVEta;
        Double_t t2=TVector2::Phi_mpi_pi(clsVec.Phi()-pa->fVPhi);
        Double_t tmp = t1*t1+t2*t2;
        if (tmp<diffR2) {
          diffR2 = tmp;
          msta   = pa;
        }
      }
      if (diffR2<10 && msta!=0) {
        cl->fMcLabel = msta->fLab;
      }
    }

    cl->fEmbE = 0;
    if (fDigits && fEmbedMode) {
      for(Int_t j=0; j<cl->fN; ++j) {
        Short_t cid = TMath::Abs(clus->GetCellAbsId(j));
        Short_t pos = -1;
        std::map<Short_t,Short_t>::iterator it = map.find(cid);
        if (it!=map.end())
          pos = it->second;
        if (pos<0)
          continue;
        AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigits->At(pos));
        if (!digit)
          continue;
        if (digit->GetId() != cid) {
          AliError(Form("Ids should be equal: %d %d", cid, digit->GetId()));
          continue;
        }
        if (digit->GetType()<-1) {
          cl->fEmbE += digit->GetChi2();
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::CalcPrimTracks()
{
  // Calculate track properties for primary tracks.

  fSelPrimTracks->Clear();

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  TClonesArray *tracks = 0;
  if (fEsdEv) {
    am->LoadBranch("Tracks");
    TList *l = fEsdEv->GetList();
    tracks = dynamic_cast<TClonesArray*>(l->FindObject("Tracks"));
  } else {
    am->LoadBranch("tracks");
    TList *l = fAodEv->GetList();
    tracks = dynamic_cast<TClonesArray*>(l->FindObject("tracks"));
  }

  if (!tracks)
    return;

  AliEMCALEMCGeometry *emc = fGeom->GetEMCGeometry();
  Double_t phimin = emc->GetArm1PhiMin()*TMath::DegToRad()-fIsoDist*1.25;
  Double_t phimax = emc->GetArm1PhiMax()*TMath::DegToRad()+fIsoDist*1.25;

  if (fEsdEv) {
    fSelPrimTracks->SetOwner(kTRUE);
    am->LoadBranch("PrimaryVertex.");
    const AliESDVertex *vtx = fEsdEv->GetPrimaryVertexSPD();
    am->LoadBranch("SPDVertex.");
    const AliESDVertex *vtxSPD = fEsdEv->GetPrimaryVertexSPD();
    am->LoadBranch("Tracks");
    const Int_t Ntracks = tracks->GetEntries();
    for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
      AliESDtrack *track = static_cast<AliESDtrack*>(tracks->At(iTracks));
      if (!track) {
        AliWarning(Form("Could not receive track %d\n", iTracks));
        continue;
      }
      if (fTrCuts && !fTrCuts->IsSelected(track))
        continue;
      Double_t eta = track->Eta();
      if (eta<-1||eta>1)
        continue;
      if (track->Phi()<phimin||track->Phi()>phimax)
        continue;

      AliESDtrack copyt(*track);
      Double_t bfield[3];
      copyt.GetBxByBz(bfield);
      AliExternalTrackParam tParam;
      Bool_t relate = copyt.RelateToVertexBxByBz(vtxSPD,bfield,kVeryBig,&tParam);
      if (!relate)
        continue;
      copyt.Set(tParam.GetX(),tParam.GetAlpha(),tParam.GetParameter(),tParam.GetCovariance());

      Double_t p[3]      = { 0. };
      copyt.GetPxPyPz(p);
      Double_t pos[3]    = { 0. };      
      copyt.GetXYZ(pos);
      Double_t covTr[21] = { 0. };
      copyt.GetCovarianceXYZPxPyPz(covTr);
      Double_t pid[10]   = { 0. };  
      copyt.GetESDpid(pid);
      AliAODTrack *aTrack = new AliAODTrack(copyt.GetID(),
                                            copyt.GetLabel(),
                                            p,
                                            kTRUE,
                                            pos,
                                            kFALSE,
                                            covTr, 
                                            (Short_t)copyt.GetSign(),
                                            copyt.GetITSClusterMap(), 
                                            pid,
                                            0,/*fPrimaryVertex,*/
                                            kTRUE, // check if this is right
                                            vtx->UsesTrack(copyt.GetID()));
      aTrack->SetTPCClusterMap(copyt.GetTPCClusterMap());
      aTrack->SetTPCSharedMap (copyt.GetTPCSharedMap());
      Float_t ndf = copyt.GetTPCNcls() + 1 - 5 ;
      if(ndf>0)
        aTrack->SetChi2perNDF(copyt.GetTPCchi2()/ndf);
      else
        aTrack->SetChi2perNDF(-1);
      aTrack->SetFlags(copyt.GetStatus());
      aTrack->SetTPCPointsF(copyt.GetTPCNclsF());
      fSelPrimTracks->Add(aTrack);
    }
  } else {
    Int_t ntracks = tracks->GetEntries();
    for (Int_t i=0; i<ntracks; ++i) {
      AliAODTrack *track = static_cast<AliAODTrack*>(tracks->At(i));
      if (!track)
        continue;
      Double_t eta = track->Eta();
      if (eta<-1||eta>1)
        continue;
      if (track->Phi()<phimin||track->Phi()>phimax)
        continue;
      if(track->GetTPCNcls()<fMinNClusPerTr)
        continue;
      //todo: Learn how to set/filter AODs for prim/sec tracks
      fSelPrimTracks->Add(track);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::CalcTracks()
{
  // Calculate track properties (including secondaries).

  fSelTracks->Clear();

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  TClonesArray *tracks = 0;
  if (fEsdEv) {
    am->LoadBranch("Tracks");
    TList *l = fEsdEv->GetList();
    tracks = dynamic_cast<TClonesArray*>(l->FindObject("Tracks"));
  } else {
    am->LoadBranch("tracks");
    TList *l = fAodEv->GetList();
    tracks = dynamic_cast<TClonesArray*>(l->FindObject("tracks"));
  }

  if (!tracks)
    return;

  if (fEsdEv) {
    const Int_t Ntracks = tracks->GetEntries();
    for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
      AliESDtrack *track = static_cast<AliESDtrack*>(tracks->At(iTracks));
      if (!track) {
        AliWarning(Form("Could not receive track %d\n", iTracks));
        continue;
      }
      if (fTrCuts && !fTrCuts->IsSelected(track))
        continue;
      Double_t eta = track->Eta();
      if (eta<-1||eta>1)
        continue;
      fSelTracks->Add(track);
    }
  } else {
    Int_t ntracks = tracks->GetEntries();
    for (Int_t i=0; i<ntracks; ++i) {
      AliAODTrack *track = static_cast<AliAODTrack*>(tracks->At(i));
      if (!track)
        continue;
      Double_t eta = track->Eta();
      if (eta<-1||eta>1)
        continue;
      if(track->GetTPCNcls()<fMinNClusPerTr)
        continue;

      if (0 && (track->Pt()>=0.6) && (track->PxAtDCA()==-999)) { // compute position on EMCAL 
        AliExternalTrackParam tParam(track);
        if (AliTrackerBase::PropagateTrackToBxByBz(&tParam, 438, 0.139, 1, kTRUE)) {
          Double_t trkPos[3];
          tParam.GetXYZ(trkPos);
          track->SetPxPyPzAtDCA(trkPos[0],trkPos[1],trkPos[2]);
        }
      }
      fSelTracks->Add(track);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::ClusterAfterburner()
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

  for (Int_t i = 0; i<ncells; ++i) {
    Short_t id=-1;
    Double_t amp=0,time=0;
    if (!cells->GetCell(i, id, amp, time))
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
void AliAnalysisTaskEMCALPi0PbPb::FillCellHists()
{
  // Fill histograms related to cell properties.
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;

  if (!cells)
    return;

  Int_t cellModCount[12] = {0};
  Double_t cellMaxE = 0; 
  Double_t cellMeanE = 0; 
  Int_t ncells = cells->GetNumberOfCells();
  if (ncells==0)
    return;

  Int_t nsm = fGeom->GetEMCGeometry()->GetNumberOfSuperModules();

  for (Int_t i = 0; i<ncells; ++i) {
    Int_t absID    = TMath::Abs(cells->GetCellNumber(i));
    Double_t cellE = cells->GetAmplitude(i);
    fHCellE->Fill(cellE);
    if (cellE>cellMaxE) 
      cellMaxE = cellE;
    cellMeanE += cellE;
    
    Int_t iSM=-1, iTower=-1, nIphi=-1, nIeta=-1;
    Bool_t ret = fGeom->GetCellIndex(absID, iSM, iTower, nIphi, nIeta);
    if (!ret) {
      AliError(Form("Could not get cell index for %d", absID));
      continue;
    }
    ++cellModCount[iSM];
    Int_t iPhi=-1, iEta=-1;
    fGeom->GetCellPhiEtaIndexInSModule(iSM, iTower, nIphi, nIeta, iPhi, iEta);   
    fHColuRow[iSM]->Fill(iEta,iPhi,1);
    fHColuRowE[iSM]->Fill(iEta,iPhi,cellE);
    fHCellFreqNoCut[iSM]->Fill(absID);
    if (cellE > 0.1) fHCellFreqCut100M[iSM]->Fill(absID);
    if (cellE > 0.3) fHCellFreqCut300M[iSM]->Fill(absID);
    if (fHCellCheckE && fHCellCheckE[absID])
      fHCellCheckE[absID]->Fill(cellE);
    fHCellFreqE[iSM]->Fill(absID, cellE);
  }    
  fHCellH->Fill(cellMaxE);
  cellMeanE /= ncells;
  fHCellM->Fill(cellMeanE);
  fHCellM2->Fill(cellMeanE*ncells/24/48/nsm);
  for (Int_t i=0; i<nsm; ++i) 
    fHCellMult[i]->Fill(cellModCount[i]);
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::FillClusHists()
{
  // Fill histograms related to cluster properties.

  TObjArray *clusters = fEsdClusters;
  if (!clusters)
    clusters = fAodClusters;
  if (!clusters)
    return;

  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  Int_t nclus = clusters->GetEntries();
  for(Int_t i = 0; i<nclus; ++i) {
    AliVCluster *clus = static_cast<AliVCluster*>(clusters->At(i));
    if (!clus)
      continue;
    if (!clus->IsEMCAL()) 
      continue;
    TLorentzVector clusterVec;
    clus->GetMomentum(clusterVec,vertex);
    Double_t maxAxis    = clus->GetTOF(); //sigma
    Double_t clusterEcc = clus->Chi2();   //eccentricity
    fHClustEccentricity->Fill(clusterEcc); 
    fHClustEtaPhi->Fill(clusterVec.Eta(),clusterVec.Phi());
    fHClustEnergyPt->Fill(clusterVec.E(),clusterVec.Pt());
    fHClustEnergySigma->Fill(clus->E()*maxAxis,clus->E());
    fHClustSigmaSigma->Fill(max(clus->GetM02(),clus->GetM20()),clus->E()*maxAxis);
    fHClustNCellEnergyRatio->Fill(clus->GetNCells(),GetMaxCellEnergy(clus)/clus->E());
    fHClustEnergyNCell->Fill(clus->E(),clus->GetNCells());  
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::CalcMcInfo()
{
  // Get Mc truth particle information.

  if (!fMcParts)
    return;

  fMcParts->Clear();

  AliEMCALEMCGeometry *emc = fGeom->GetEMCGeometry();
  Double_t etamin = -0.7;
  Double_t etamax = +0.7;
  Double_t phimin = emc->GetArm1PhiMin()*TMath::DegToRad();
  Double_t phimax = emc->GetArm1PhiMax()*TMath::DegToRad();

  if (fAodEv) {
    AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
    am->LoadBranch(AliAODMCParticle::StdBranchName());
    TClonesArray *tca = dynamic_cast<TClonesArray*>(fAodEv->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!tca) 
      return;

    Int_t nents = tca->GetEntries();
    for(int it=0; it<nents; ++it) {
      AliAODMCParticle *part = static_cast<AliAODMCParticle*>(tca->At(it));
      part->Print();

      // pion or eta meson or direct photon
      if(part->GetPdgCode() == 111) {
      } else if(part->GetPdgCode() == 221) {
      } else if(part->GetPdgCode() == 22 ) {
      }	else
        continue;

      // primary particle
      Double_t dR = TMath::Sqrt((part->Xv()*part->Xv())+(part->Yv()*part->Yv()));
      if(dR > 1.0)
        continue;

      // kinematic cuts
      Double_t pt = part->Pt() ;
      if (pt<0.5)
        continue;
      Double_t eta = part->Eta();
      if (eta<etamin||eta>etamax)
        continue;
      Double_t phi  = part->Phi();
      if (phi<phimin||phi>phimax)
        continue;

      ProcessDaughters(part, it, tca);
    } 
    return;
  }

  AliMCEvent *mcEvent = MCEvent();
  if (!mcEvent)
    return;

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

    // primary particle
    Double_t dR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX())*(mcP->Xv()-evtVtx->GetX()) + 
                              (mcP->Yv()-evtVtx->GetY())*(mcP->Yv()-evtVtx->GetY()));
    if(dR > 1.0)
      continue;

    // kinematic cuts
    Double_t pt = mcP->Pt() ;
    if (pt<0.5)
      continue;
    Double_t eta = mcP->Eta();
    if (eta<etamin||eta>etamax)
      continue;
    Double_t phi  = mcP->Phi();
    if (phi<phimin||phi>phimax)
      continue;

    ProcessDaughters(mcP, iTrack, mcEvent);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::FillNtuple()
{
  // Fill ntuple.

  if (!fNtuple)
    return;

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (fAodEv) {
    fHeader->fRun            = fAodEv->GetRunNumber();
    fHeader->fOrbit          = fAodEv->GetHeader()->GetOrbitNumber(); 
    fHeader->fPeriod         = fAodEv->GetHeader()->GetPeriodNumber();
    fHeader->fBx             = fAodEv->GetHeader()->GetBunchCrossNumber();
    fHeader->fL0             = fAodEv->GetHeader()->GetL0TriggerInputs();
    fHeader->fL1             = fAodEv->GetHeader()->GetL1TriggerInputs();
    fHeader->fL2             = fAodEv->GetHeader()->GetL2TriggerInputs();
    fHeader->fTrClassMask    = fAodEv->GetHeader()->GetTriggerMask();
    fHeader->fTrCluster      = fAodEv->GetHeader()->GetTriggerCluster();
    fHeader->fOffTriggers    = fAodEv->GetHeader()->GetOfflineTrigger();
    fHeader->fFiredTriggers  = fAodEv->GetHeader()->GetFiredTriggerClasses();
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
void AliAnalysisTaskEMCALPi0PbPb::FillPionHists()
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
      for (Int_t j = i+1; j<nclus; ++j) {
        AliVCluster *clus2 = static_cast<AliVCluster*>(clusters->At(j));
        if (!clus2)
          continue;
        if (!clus2->IsEMCAL()) 
          continue;
        if (clus2->E()<fMinE)
          continue;
        if (clus2->GetNCells()<fNminCells)
          continue;
        if (GetMaxCellEnergy(clus2)/clus2->E()<fMinErat)
          continue;
        if (clus2->Chi2()<fMinEcc) // eccentricity cut
          continue;
        clus2->GetMomentum(clusterVec2,vertex);
        pionVec = clusterVec1 + clusterVec2;
        Double_t pionZgg = TMath::Abs(clusterVec1.E()-clusterVec2.E())/pionVec.E();
        Double_t pionOpeningAngle = clusterVec1.Angle(clusterVec2.Vect());
        if (pionZgg < fAsymMax) {
          fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi()); 
          fHPionMggPt->Fill(pionVec.M(),pionVec.Pt()); 
          fHPionMggAsym->Fill(pionVec.M(),pionZgg); 
          fHPionMggDgg->Fill(pionVec.M(),pionOpeningAngle); 
          Int_t bin = fPtRanges->FindBin(pionVec.Pt());
          fHPionInvMasses[bin]->Fill(pionVec.M());
        }
      }
    }
  } 
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::FillMcHists()
{
  // Fill additional MC information histograms.

  if (!fMcParts)
    return;

  // check if aod or esd mc mode and the fill histos
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::FillOtherHists()
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
void AliAnalysisTaskEMCALPi0PbPb::FillTrackHists()
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
void AliAnalysisTaskEMCALPi0PbPb::FillVertex(AliStaVertex *v, const AliESDVertex *esdv)
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
void AliAnalysisTaskEMCALPi0PbPb::FillVertex(AliStaVertex *v, const AliAODVertex *aodv)
{
  // Fill vertex from AOD vertex info.

  v->fVx   = aodv->GetX();
  v->fVy   = aodv->GetY();
  v->fVz   = aodv->GetZ();
  v->fVc   = aodv->GetNContributors();
  v->fChi2 = aodv->GetChi2();
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPi0PbPb::GetCellIsolation(Double_t cEta, Double_t cPhi, Double_t radius) const
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
Double_t AliAnalysisTaskEMCALPi0PbPb::GetCellIsoNxM(Double_t cEta, Double_t cPhi, Int_t N, Int_t M) const
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
Double_t AliAnalysisTaskEMCALPi0PbPb::GetCellEnergy(const AliVCluster *cluster) const
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
Double_t AliAnalysisTaskEMCALPi0PbPb::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
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
Double_t AliAnalysisTaskEMCALPi0PbPb::GetSecondMaxCellEnergy(AliVCluster *clus, Short_t &id) const
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
void AliAnalysisTaskEMCALPi0PbPb::GetSigma(const AliVCluster *c, Double_t& sigmaMax, Double_t &sigmaMin) const
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
void AliAnalysisTaskEMCALPi0PbPb::GetSigmaEtaEta(const AliVCluster *c, Double_t& sEtaEta, Double_t &sPhiPhi) const
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
Int_t AliAnalysisTaskEMCALPi0PbPb::GetNCells(const AliVCluster *c, Double_t emin) const
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
Int_t AliAnalysisTaskEMCALPi0PbPb::GetNCells(Int_t sm, Double_t emin) const
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
Double_t AliAnalysisTaskEMCALPi0PbPb::GetTrackIsolation(Double_t cEta, Double_t cPhi, Double_t radius, Double_t pt) const
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
Double_t AliAnalysisTaskEMCALPi0PbPb::GetTrackIsoStrip(Double_t cEta, Double_t cPhi, Double_t dEta, Double_t dPhi, Double_t pt) const
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
Bool_t AliAnalysisTaskEMCALPi0PbPb::IsShared(const AliVCluster *c) const
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
Bool_t AliAnalysisTaskEMCALPi0PbPb::IsIdPartOfCluster(const AliVCluster *c, Short_t id) const
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
void AliAnalysisTaskEMCALPi0PbPb::PrintDaughters(const AliVParticle *p, const TObjArray *arr, Int_t level) const
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
void AliAnalysisTaskEMCALPi0PbPb::PrintDaughters(const AliMCParticle *p, const AliMCEvent *arr, Int_t level) const
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
void AliAnalysisTaskEMCALPi0PbPb::PrintTrackRefs(AliMCParticle *p) const
{
  // Print track ref array.

  if (!p)
    return;

  Int_t n = p->GetNumberOfTrackReferences();
  for (Int_t i=0; i<n; ++i) {
    AliTrackReference *ref = p->GetTrackReference(i);
    if (!ref)
      continue;
    ref->SetUserId(ref->DetectorId());
    ref->Print();
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::ProcessDaughters(AliVParticle *p, Int_t index, const TObjArray *arr)
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
void AliAnalysisTaskEMCALPi0PbPb::ProcessDaughters(AliMCParticle *p, Int_t index, const AliMCEvent *arr)
{
  // Process and create daughters.

  if (!p || !arr)
    return;

  Int_t d1 = p->GetFirstDaughter();
  Int_t d2 = p->GetLastDaughter();
  if (0) {
    printf("%d pid=%d: %.3f %.3f %.3f (%.2f %.2f %.2f); nd=%d,%d, mo=%d\n",
           index,p->PdgCode(),p->Px(),p->Py(),p->Pz(),p->Xv(),p->Yv(),p->Zv(),d1,d2, p->GetMother());
    PrintTrackRefs(p);
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

//__________________________________________________________________________________________________
void AliStaCluster::GetMom(TLorentzVector& p, Double_t *vertex) 
{
  // Calculate momentum.

  TVector3 pos;
  pos.SetPtEtaPhi(fR,fEta,fPhi);

  if(vertex){ //calculate direction relative to  vertex
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
