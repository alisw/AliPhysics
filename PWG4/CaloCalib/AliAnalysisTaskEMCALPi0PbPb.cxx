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
#include <TVector2.h>
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliCDBManager.h"
#include "AliCentrality.h"
#include "AliEMCALGeoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliGeomManager.h"
#include "AliLog.h"
#include "AliMagF.h"
#include "AliTrackerBase.h"

ClassImp(AliAnalysisTaskEMCALPi0PbPb)

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
    fMinNClustPerTrack(50),
    fMinPtPerTrack(1.0), 
    fIsoDist(0.2),
    fTrCuts(0),
    fNEvs(0),
    fGeom(0),
    fReco(0),
    fOutput(0),
    fEsdEv(0),
    fAodEv(0),
    fRecPoints(0),
    fEsdClusters(0),
    fEsdCells(0),
    fAodClusters(0),
    fAodCells(0),
    fPtRanges(0),
    fNtuple(0),
    fSelTracks(0),
    fHCuts(0x0),
    fHVertexZ(0x0),
    fHVertexZ2(0x0),
    fHCent(0x0),
    fHCentQual(0x0),
    fHColuRow(0x0),
    fHColuRowE(0x0),
    fHCellMult(0x0),
    fHCellE(0x0),
    fHCellH(0x0),
    fHCellM(0x0),
    fHCellM2(0x0),
    fHCellFreqNoCut(0x0),
    fHCellFrequCut100M(0x0),
    fHCellFrequCut300M(0x0),
    fHCellCheckE(0x0),
    fHClustEccentricity(0),
    fHClustEtaPhi(0x0),
    fHClustEnergyPt(0x0),
    fHClustEnergySigma(0x0),
    fHClustSigmaSigma(0x0),
    fHClustNCellEnergyRatio(0x0),
    fHPionEtaPhi(0x0),
    fHPionMggPt(0x0),
    fHPionMggAsym(0x0),
    fHPionMggDgg(0x0)
{
  // Constructor.

  if (!name)
    return;
  SetName(name);
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.,SPDVertex.,TPCVertex.,EMCALCells.,Tracks "
               "AOD:header,vertices,emcalCells,tracks";
}

//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::~AliAnalysisTaskEMCALPi0PbPb()
{
  // Destructor.

  delete fOutput; fOutput = 0;
  delete fPtRanges; fPtRanges = 0;
  delete fGeom; fGeom = 0;
  delete fReco; fReco = 0;
  delete fSelTracks;
  delete [] fHColuRow;
  delete [] fHColuRowE;
  delete [] fHCellMult;
  delete [] fHCellFreqNoCut;
  delete [] fHCellFrequCut100M;
  delete [] fHCellFrequCut300M;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::UserCreateOutputObjects()
{
  // Create user objects here.

  fGeom = new AliEMCALGeoUtils(fGeoName,"EMCAL");
  fReco = new AliEMCALRecoUtils();
  fOutput = new TList();
  fOutput->SetOwner();
  fSelTracks = new TObjArray;

  if (fDoNtuple) {
    TFile *f = OpenFile(1);
    if (f) {
      f->SetCompressionLevel(2);
      fNtuple = new TNtuple(Form("nt%.0fto%.0f",fCentFrom,fCentTo),"nt",
                            "run:evt:l0:cent:pt:eta:phi:e:emax:n:n1:nsm:db:disp:mn:ms:ecc:sig:tkdz:tkdr:tkep:tkiso:ceiso");
      fNtuple->SetDirectory(f);
      fNtuple->SetAutoFlush(-1024*1024*1024);
      fNtuple->SetAutoSave(-1024*1024*1024);
    }  
  }

  // histograms
  TH1::SetDefaultSumw2(kTRUE);
  TH2::SetDefaultSumw2(kTRUE);
  fHCuts = new TH1F("hEventCuts","",4,0.5,4.5);
  fHCuts->GetXaxis()->SetBinLabel(1,"All (PS)");
  fHCuts->GetXaxis()->SetBinLabel(2,Form("%s: %.0f-%.0f",fCentVar.Data(),fCentFrom,fCentTo));
  fHCuts->GetXaxis()->SetBinLabel(3,"QFlag");
  fHCuts->GetXaxis()->SetBinLabel(4,Form("zvtx: %.0f-%.0f",fVtxZMin,fVtxZMax));
  fOutput->Add(fHCuts);
  fHVertexZ = new TH1F("hVertexZBeforeCut","",100,-25,25);
  fHVertexZ->SetXTitle("z [cm]");
  fOutput->Add(fHVertexZ);
  fHVertexZ2 = new TH1F("hVertexZAfterCut","",100,-25,25);
  fHVertexZ2->SetXTitle("z [cm]");
  fOutput->Add(fHVertexZ2);
  fHCent = new TH1F("hCentBeforeCut","",101,-1,100);
  fHCent->SetXTitle(fCentVar.Data());
  fOutput->Add(fHCent);
  fHCentQual = new TH1F("hCentAfterCut","",101,-1,100);
  fHCentQual->SetXTitle(fCentVar.Data());
  fOutput->Add(fHCentQual);

  // histograms for cells
  Int_t nsm = fGeom->GetEMCGeometry()->GetNumberOfSuperModules();
  fHColuRow   = new TH2*[nsm];
  fHColuRowE  = new TH2*[nsm];
  fHCellMult  = new TH1*[nsm];
  for (Int_t i = 0; i<nsm; ++i) {
    fHColuRow[i] = new TH2F(Form("hColRow_Mod%d",i),"",49,0,49,25,0,25);
    fHColuRow[i]->SetTitle(Form("Module %d: Occupancy", i));
    fHColuRow[i]->SetXTitle("col (i#eta)");
    fHColuRow[i]->SetYTitle("row (i#phi)");
    fHColuRowE[i] = new TH2F(Form("hColRowE_Mod%d", i),"",49,0,49,25,0,25);
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
  fHCellE = new TH1F("hCellE","",150,0.,15.);
  fHCellE->SetXTitle("E_{cell} [GeV]");
  fOutput->Add(fHCellE);
  fHCellH = new TH1F ("fHCellHighestE","",150,0.,15.);
  fHCellH->SetXTitle("E^{max}_{cell} [GeV]");
  fOutput->Add(fHCellH);
  fHCellM = new TH1F ("fHCellMeanEperHitCell","",250,0.,2.5);
  fHCellM->SetXTitle("#LT E_{cell}#GT [GeV]");
  fOutput->Add(fHCellM);
  fHCellM2 = new TH1F ("fHCellMeanEperAllCells","",250,0.,1);
  fHCellM2->SetXTitle("1/N_{cells} #Sigma E_{cell} [GeV]");
  fOutput->Add(fHCellM2);

  fHCellFreqNoCut    = new TH1*[nsm];
  fHCellFrequCut100M = new TH1*[nsm];
  fHCellFrequCut300M = new TH1*[nsm];
  for (Int_t i = 0; i<nsm; ++i){
    Double_t lbin = i*24*48-0.5;
    Double_t hbin = lbin+24*48;
    fHCellFreqNoCut[i]    = new TH1F(Form("fHCellFreqNoCut_SM%d",i),    
                                     Form("Frequency SM%d (no cut);id;#",i),  1152, lbin, hbin);
    fHCellFrequCut100M[i] = new TH1F(Form("fHCellFreqCut100M_SM%d",i), 
                                     Form("Frequency SM%d (>0.1GeV);id;#",i), 1152, lbin, hbin);
    fHCellFrequCut300M[i] = new TH1F(Form("fHCellFreqCut300M_SM%d",i), 
                                     Form("Frequency SM%d (>0.3GeV);id;#",i), 1152, lbin, hbin);
    fOutput->Add(fHCellFreqNoCut[i]);
    fOutput->Add(fHCellFrequCut100M[i]);
    fOutput->Add(fHCellFrequCut300M[i]);
  }
  if (1) {
    fHCellCheckE = new TH1*[24*48*nsm];
    memset(fHCellCheckE,0,sizeof(fHCellCheckE));
    Int_t tcs[1] = {4102};
    for (UInt_t i = 0; i<sizeof(tcs)/sizeof(Int_t); ++i){
      Int_t c=tcs[i];
      if (c<24*48*nsm) {
        fHCellCheckE[i] = new TH1F(Form("fHCellE_id%d",c), Form("Cell %d;E [GeV/c];#",c), 500, 0, 8);
        fOutput->Add(fHCellCheckE[i]);
      }
    }
  }

  // histograms for clusters
  fHClustEccentricity = new TH1F("hClustEccentricity","",100,-0.1,1.1);
  fHClustEccentricity->SetXTitle("#epsilon_{C}");
  fOutput->Add(fHClustEccentricity);
  fHClustEtaPhi = new TH2F("hClustEtaPhi","",500,-0.8,0.8,500,1.2,2.2);
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

  // histograms for pion candidates
  fHPionEtaPhi = new TH2F("hPionEtaPhi","",100,-0.8,0.8,100,1.2,2.2);
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
  if (fEsdEv) {
    am->LoadBranch("AliESDRun.");
    am->LoadBranch("AliESDHeader.");
  } else {
    fAodEv = dynamic_cast<AliAODEvent*>(InputEvent());
    am->LoadBranch("header");
  }

  if (fHCuts->GetEntries()==0) {
    if (!AliGeomManager::GetGeometry()) { // get geometry 
      AliWarning("Accessing geometry from OCDB, this is not very efficient!");
      AliCDBManager *cdb = AliCDBManager::Instance();
      if (!cdb->IsDefaultStorageSet())
        cdb->SetDefaultStorage("raw://");
      Int_t runno = InputEvent()->GetRunNumber();
      if (runno != cdb->GetRun())
        cdb->SetRun(runno);
      AliGeomManager::LoadGeometry();
    }

    if (fEsdEv) {  // set misalignment matrices (stored in first event)
      for (Int_t i=0; i<fGeom->GetEMCGeometry()->GetNumberOfSuperModules(); ++i)
        fGeom->SetMisalMatrix(fEsdEv->GetESDRun()->GetEMCALMatrix(i),i);
    } else {
      for (Int_t i=0; i<fGeom->GetEMCGeometry()->GetNumberOfSuperModules(); ++i)
        fGeom->SetMisalMatrix(fAodEv->GetHeader()->GetEMCALMatrix(i),i);
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
  }

  Int_t cut = 1;
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

  // count number of accepted events
  ++fNEvs;

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

  fRecPoints   = 0; // will be set if fClusName is given and AliAnalysisTaskEMCALClusterizeFast is used
  fEsdClusters = 0; // will be set if ESD input used and if fRecPoints are not set of if clusters are attached
  fEsdCells    = 0; // will be set if ESD input used
  fAodClusters = 0; // will be set if AOD input used and if fRecPoints are not set of if clusters are attached
                    //             or if fClusName is given and AliAnalysisTaskEMCALClusterizeFast in AOD output mode
  fAodCells    = 0; // will be set if AOD input used

  // deal with special output from AliAnalysisTaskEMCALClusterizeFast first
  Bool_t clusattached = 0;
  Bool_t recalibrated = 0;
  if (1 && !fClusName.IsNull()) {
    AliAnalysisTaskEMCALClusterizeFast *cltask = 0;
    TObjArray *ts = am->GetTasks();
    cltask = dynamic_cast<AliAnalysisTaskEMCALClusterizeFast*>(ts->FindObject(fClusName));
    if (cltask && cltask->GetClusters()) {
      fRecPoints = const_cast<TObjArray*>(cltask->GetClusters());
      clusattached = cltask->GetAttachClusters();
      if (cltask->GetCalibData()!=0)
        recalibrated = kTRUE;
    }
  }
  if (1 && AODEvent() && !fClusName.IsNull()) {
    TList *l = AODEvent()->GetList();
    TClonesArray *clus = 0;
    if (l) {
      clus = dynamic_cast<TClonesArray*>(l->FindObject(fClusName));
      fAodClusters = clus;
    }
  }

  if (fEsdEv) { // ESD input mode
    if (1 && (!fRecPoints||clusattached)) {
      if (!clusattached)
        am->LoadBranch("CaloClusters");
      TList *l = fEsdEv->GetList();
      TClonesArray *clus = 0;
      if (l) {
        clus = dynamic_cast<TClonesArray*>(l->FindObject("CaloClusters"));
        fEsdClusters = clus;
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
      TClonesArray *clus = 0;
      if (l) {
        clus = dynamic_cast<TClonesArray*>(l->FindObject("caloClusters"));
        fAodClusters = clus;
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
    AliDebug(2,Form("fEsdClusters set: %p", fEsdClusters));
    AliDebug(2,Form("fEsdCells    set: %p", fEsdCells));
    AliDebug(2,Form("fAodClusters set: %p", fAodClusters));
    AliDebug(2,Form("fAodCells    set: %p", fAodCells));
  }

  if (fDoAfterburner)
    ClusterAfterburner();

  CalcTracks();
  CalcClusterProps();

  FillCellHists();
  FillClusHists();
  FillPionHists();
  FillOtherHists();

  PostData(1, fOutput);
}      

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::Terminate(Option_t *) 
{
  // Terminate called at the end of analysis.

  if (fNtuple) {
    TFile *f = OpenFile(1);
    if (f) 
      fNtuple->Write();
  }

  AliInfo(Form("\n%s: Accepted %lld events", GetName(), fNEvs));
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::CalcTracks()
{
  // Calculate track properties.

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

  AliEMCALEMCGeometry *emc = fGeom->GetEMCGeometry();
  Double_t phimin = emc->GetArm1PhiMin()*TMath::DegToRad()-0.25;
  Double_t phimax = emc->GetArm1PhiMax()*TMath::DegToRad()+0.25;

  if (fEsdEv) {
    fSelTracks->SetOwner(kTRUE);
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
      Double_t pt  = track->Pt();
      if (pt<fMinPtPerTrack) 
        continue;
      if(track->GetTPCNcls()<fMinNClustPerTrack)
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
      fSelTracks->Add(aTrack);
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
      Double_t pt  = track->Pt();
      if (pt<fMinPtPerTrack) 
        continue;
      if(track->GetTPCNcls()<fMinNClustPerTrack)
        continue;

      if (0 && (pt>=0.6) && (track->PxAtDCA()==-999)) { // compute position on EMCAL 
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
void AliAnalysisTaskEMCALPi0PbPb::CalcClusterProps()
{
  // Calculate cluster properties

  TObjArray *clusters = fEsdClusters;
  if (!clusters)
    clusters = fAodClusters;
  if (!clusters)
    return;

  Int_t nclus = clusters->GetEntries();
  Int_t ntrks = fSelTracks->GetEntries();

  for(Int_t i = 0; i<nclus; ++i) {
    fClusProps[i].Reset();

    AliVCluster *clus = static_cast<AliVCluster*>(clusters->At(i));
    if (!clus)
      continue;
    if (!clus->IsEMCAL())
      continue;
    if (clus->E()<fMinE)
      continue;

    Float_t  clsPos[3] = {0};
    clus->GetPosition(clsPos);
    TVector3 clsVec(clsPos);
    Double_t vertex[3] = {0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    TLorentzVector clusterVec;
    clus->GetMomentum(clusterVec,vertex);
    Double_t clsEta = clusterVec.Eta();

    Double_t mind2 = 1e10;
    for(Int_t j = 0; j<ntrks; ++j) {
      AliVTrack *track = static_cast<AliVTrack*>(fSelTracks->At(j));
      if (!track)
        continue;
      if (track->Pt()<0.6)
        continue;
      if (TMath::Abs(clsEta-track->Eta())>fIsoDist)
        continue;

      AliExternalTrackParam tParam(track);
      Float_t tmpR=-1, tmpZ=-1;
#if 1
      if (1) {
        TVector3 vec(clsPos);
        Double_t alpha =  ((int)(vec.Phi()*TMath::RadToDeg()/20)+0.5)*20*TMath::DegToRad();
        vec.RotateZ(-alpha);  //Rotate the cluster to the local extrapolation coordinate system
        tParam.Rotate(alpha); //Rotate the track to the same local extrapolation system
        if (!AliTrackerBase::PropagateTrackToBxByBz(&tParam, vec.X(), 0.139, 1, kFALSE)) 
          continue;
        Double_t trkPos[3];
        tParam.GetXYZ(trkPos); //Get the extrapolated global position
        tmpR = TMath::Sqrt( TMath::Power(clsPos[0]-trkPos[0],2)+TMath::Power(clsPos[1]-trkPos[1],2)+TMath::Power(clsPos[2]-trkPos[2],2) );
        tmpZ = clsPos[2]-trkPos[2];
      }
#else
      if (!fReco->ExtrapolateTrackToCluster(&tParam, clus, tmpR, tmpZ))
        continue;
#endif
      Double_t d2 = tmpR;
      if (mind2>d2) {
        mind2=d2;
        fClusProps[i].fTrIndex = j;
        fClusProps[i].fTrDz    = tmpZ;
        fClusProps[i].fTrDr    = TMath::Sqrt(tmpR*tmpR-tmpZ*tmpZ);
        fClusProps[i].fTrDist  = d2;
        fClusProps[i].fTrEp    = clus->E()/track->P();
      }
    }

    if (0 && (fClusProps[i].fTrIndex>=0)) {
      cout << i << " " << fClusProps[i].fTrIndex << ": Dr " << fClusProps[i].fTrDr << " " << " Dz " << fClusProps[i].fTrDz << endl;
    }
    fClusProps[i].fTrIso      = GetTrackIsolation(clusterVec.Eta(),clusterVec.Phi(),fIsoDist);
    fClusProps[i].fTrLowPtIso = 0;
    fClusProps[i].fCellIso    = GetCellIsolation(clsVec.Eta(),clsVec.Phi(),fIsoDist);
  }
}

//________________________________________________________________________
void  AliAnalysisTaskEMCALPi0PbPb::ClusterAfterburner()
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
    if (cellE > 0.1) fHCellFrequCut100M[iSM]->Fill(absID);
    if (cellE > 0.3) fHCellFrequCut300M[iSM]->Fill(absID);
    if (fHCellCheckE && fHCellCheckE[absID])
      fHCellCheckE[absID]->Fill(cellE);
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
    Double_t maxAxis = 0, minAxis = 0;
    GetSigma(clus,maxAxis,minAxis);
    Double_t clusterEcc = 0;
    if (maxAxis > 0)
      clusterEcc = TMath::Sqrt(1.0 - minAxis*minAxis/(maxAxis*maxAxis));
    clus->SetChi2(clusterEcc); // store ecc in chi2
    fHClustEccentricity->Fill(clusterEcc); 
    fHClustEtaPhi->Fill(clusterVec.Eta(),clusterVec.Phi());
    fHClustEnergyPt->Fill(clusterVec.E(),clusterVec.Pt());
    fHClustEnergySigma->Fill(clus->E()*maxAxis,clus->E());
    fHClustSigmaSigma->Fill(max(clus->GetM02(),clus->GetM20()),clus->E()*maxAxis);
    fHClustNCellEnergyRatio->Fill(clus->GetNCells(),GetMaxCellEnergy(clus)/clus->E());
    if (fNtuple) {
      if (clus->E()<fMinE)
        continue;
      Float_t vals[23];
      vals[0]  = InputEvent()->GetRunNumber();
      vals[1]  = (((UInt_t)InputEvent()->GetOrbitNumber()  << 12) | (UInt_t)InputEvent()->GetBunchCrossNumber()); 
      if (vals[1]<=0) 
        vals[1] = fNEvs;
      AliESDHeader *h = dynamic_cast<AliESDHeader*>(InputEvent()->GetHeader());
      if (h)
        vals[2]  = h->GetL0TriggerInputs();
      else {
        AliAODHeader *h2 = dynamic_cast<AliAODHeader*>(InputEvent()->GetHeader());
        if (h2)
          vals[2]  = h2->GetL0TriggerInputs();
        else 
          vals[2] = 0;
      }
      vals[3]  = InputEvent()->GetCentrality()->GetCentralityPercentileUnchecked(fCentVar);
      vals[4]  = clusterVec.Pt();
      vals[5]  = clusterVec.Eta();
      vals[6]  = clusterVec.Phi();
      vals[7]  = clusterVec.E();
      vals[8]  = GetMaxCellEnergy(clus);
      vals[9]  = clus->GetNCells();
      vals[10] = GetNCells(clus,0.100);
      vals[11] = fGeom->GetSuperModuleNumber(clus->GetCellAbsId(0));
      vals[12] = clus->GetDistanceToBadChannel();
      vals[13] = clus->GetDispersion();
      vals[14] = clus->GetM20();
      vals[15] = clus->GetM02();
      vals[16] = clusterEcc;
      vals[17] = maxAxis;
      vals[18] = fClusProps[i].fTrDz; 
      vals[19] = fClusProps[i].fTrDr;
      vals[20] = fClusProps[i].fTrEp;
      vals[21] = fClusProps[i].fTrIso;
      vals[22] = fClusProps[i].fCellIso;
      fNtuple->Fill(vals);
    }
  }
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
void AliAnalysisTaskEMCALPi0PbPb::FillOtherHists()
{
  // Fill histograms related to cell properties.
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
    Double_t cellE = cells->GetAmplitude(i);
    Float_t eta, phi;
    fGeom->EtaPhiFromIndex(absID,eta,phi);
    Double_t phidiff = TVector2::Phi_0_2pi(phi-cPhi);
    Double_t dist = (eta-cEta)*(eta-cEta)+phidiff*phidiff;
    if(dist>rad2)
      continue;
    cellIsolation += cellE;
  }
  return cellIsolation;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPi0PbPb::GetMaxCellEnergy(AliVCluster *cluster) const
{
  // Get maximum energy of attached cell.

  Double_t maxe = 0;
  Int_t ncells = cluster->GetNCells();
  if (fEsdCells) {
    for (Int_t i=0; i<ncells; i++) {
      Double_t e = fEsdCells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
      if (e>maxe) {
        maxe = e;
      }
    }
  } else {
    for (Int_t i=0; i<ncells; i++) {
      Double_t e = fAodCells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
      if (e>maxe)
        maxe = e;
    }
  }
  return maxe;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::GetSigma(AliVCluster *c, Double_t& sigmaMax, Double_t &sigmaMin) const
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

  TVector3 pos;
  for(Int_t j=0; j<ncells; ++j) {
    Int_t id = TMath::Abs(c->GetCellAbsId(j));
    fGeom->GetGlobal(id,pos);
    Double_t cellen = cells->GetCellAmplitude(id);
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
Int_t AliAnalysisTaskEMCALPi0PbPb::GetNCells(AliVCluster *c, Double_t emin) const
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
Double_t AliAnalysisTaskEMCALPi0PbPb::GetTrackIsolation(Double_t cEta, Double_t cPhi, Double_t radius) const
{
  // Compute isolation based on tracks.
  
  Double_t trkIsolation = 0;
  Double_t rad2 = radius*radius;
  Int_t ntrks = fSelTracks->GetEntries();
  for(Int_t j = 0; j<ntrks; ++j) {
    AliVTrack *track = static_cast<AliVTrack*>(fSelTracks->At(j));
    if (!track)
      continue;
    Float_t eta = track->Eta();
    Float_t phi = track->Phi();
    Double_t phidiff = TVector2::Phi_0_2pi(phi-cPhi);
    Double_t dist = (eta-cEta)*(eta-cEta)+phidiff*phidiff;
    if(dist>rad2)
      continue;
    trkIsolation += track->Pt();
  } 
  return trkIsolation;
}
