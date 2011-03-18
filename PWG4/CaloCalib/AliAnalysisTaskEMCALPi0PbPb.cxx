// $Id$

#include "AliAnalysisTaskEMCALPi0PbPb.h"
#include <TAxis.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliCentrality.h"
#include "AliEMCALGeoUtils.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliLog.h"

ClassImp(AliAnalysisTaskEMCALPi0PbPb)

//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::AliAnalysisTaskEMCALPi0PbPb(const char *name) 
  : AliAnalysisTaskSE(name),
    fAsymMax(1),
    fCentVar("V0M"),
    fCentFrom(0),
    fCentTo(100),
    fVtxZMin(-10),
    fVtxZMax(+10),
    fUseQualFlag(1),
    fClusName(),
    fDoNtuple(0),
    fDoAfterburner(0),
    fNminCells(1),
    fGeom(0),
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
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.,SPDVertex.,TPCVertex.,EMCALCells. "
               "AOD:header,vertices,emcalCells";
}

//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::~AliAnalysisTaskEMCALPi0PbPb()
{
  // Destructor.

  delete fOutput; fOutput = 0;
  delete fPtRanges; fPtRanges = 0;
  delete fGeom; fGeom = 0;
  delete [] fHColuRow;
  delete [] fHColuRowE;
  delete [] fHCellMult;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::UserCreateOutputObjects()
{
  // Create user objects here.

  fOutput = new TList();
  fOutput->SetOwner();

  if (fDoNtuple) {
    TFile *f = OpenFile(1);
    if (f)
      fNtuple = new TNtuple(Form("nt%.0fto%.0f",fCentFrom,fCentTo),"nt",
                            "run:evt:cent:pt:e:emax:n:db:disp:mn:ms:chi:cpv:ecc:sig:eta:phi");
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
  fHColuRow   = new TH2F*[4];
  fHColuRowE  = new TH2F*[4];
  fHCellMult  = new TH1F*[4];
  for (Int_t i = 0; i<4; ++i) {
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
  fHCellM2->SetXTitle("1/N_{cells} #Sum E_{cell} [GeV]");
  fOutput->Add(fHCellM2);

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
  const Int_t nbins = 19; 
  Double_t xbins[nbins] = {0.5,1,1.5,2,2.5,3,3.5,4,4.5,6,7,8,9,10,12.5,15,20,25,50};
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

  if (!fGeom) { // set misalignment matrices (stored in first event)
    fGeom = new AliEMCALGeoUtils("EMCAL_FIRSTYEARV1","EMCAL");
    if (fEsdEv) {
      for (Int_t i=0; i<4; ++i)
        fGeom->SetMisalMatrix(fEsdEv->GetESDRun()->GetEMCALMatrix(i),i);
    } else {
      for (Int_t i=0; i<4; ++i)
        fGeom->SetMisalMatrix(fAodEv->GetHeader()->GetEMCALMatrix(i),i);
    }
    fGeom->GetEMCGeometry();
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

  if (fEsdEv) {
    am->LoadBranch("PrimaryVertex.");
    am->LoadBranch("SPDVertex.");
    am->LoadBranch("TPCVertex.");
  } else {
    fAodEv = dynamic_cast<AliAODEvent*>(InputEvent());
    am->LoadBranch("vertices");
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

  FillCellHists();
  FillClusHists();
  FillPionHists();

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
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::FillCellHists()
{
  // Fill histograms related to cell properties.
  
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;

  if (cells) {
    Int_t cellModCount[4] = {0,0,0,0};
    Double_t cellMaxE = 0; 
    Double_t cellMeanE = 0; 
    Int_t ncells = cells->GetNumberOfCells();
    if (ncells==0)
      return;

    for (Int_t i = 0; i<ncells; ++i ) {
      Int_t absID    = cells->GetCellNumber(i);
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
    }    
    fHCellH->Fill(cellMaxE);
    cellMeanE /= ncells;
    fHCellM->Fill(cellMeanE);
    fHCellM2->Fill(cellMeanE*ncells/24/48/4); //hard-coded but there is a way to figure out from geometry
    for (Int_t i=0; i<4; ++i) 
      fHCellMult[i]->Fill(cellModCount[i]);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::FillClusHists()
{
  // Fill histograms related to cluster properties.
  Double_t clusterEcc = 0;


  Double_t vertex[3] = {0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  TObjArray *clusters = fEsdClusters;
  if (!clusters)
    clusters = fAodClusters;

  if (clusters) {
    TLorentzVector clusterVec;
    Int_t nclus = clusters->GetEntries();
    for(Int_t i = 0; i<nclus; ++i) {
      AliVCluster *clus = static_cast<AliVCluster*>(clusters->At(i));
      if (!clus)
        continue;
      if (!clus->IsEMCAL()) 
        continue;
      clus->GetMomentum(clusterVec,vertex);
      Double_t maxAxis = 0, minAxis = 0;
      GetSigma(clus,maxAxis,minAxis);
      if (maxAxis > 0)
        clusterEcc = TMath::Sqrt(1.0 - minAxis*minAxis/(maxAxis*maxAxis));

      fHClustEccentricity->Fill(clusterEcc); 
      fHClustEtaPhi->Fill(clusterVec.Eta(),clusterVec.Phi());
      fHClustEnergyPt->Fill(clusterVec.E(),clusterVec.Pt());
      fHClustEnergySigma->Fill(clus->E()*maxAxis,clus->E());
      fHClustSigmaSigma->Fill(max(clus->GetM02(),clus->GetM20()),clus->E()*maxAxis);
      fHClustNCellEnergyRatio->Fill(clus->GetNCells(),GetMaxCellEnergy(clus)/clus->E());
      if (fNtuple) {
        Float_t vals[18];
        vals[0]  = InputEvent()->GetRunNumber();
        vals[1]  = (((UInt_t)InputEvent()->GetOrbitNumber()  << 12) | (UInt_t)InputEvent()->GetBunchCrossNumber()); 
        vals[2]  = InputEvent()->GetCentrality()->GetCentralityPercentileUnchecked(fCentVar);
        vals[3]  = clusterVec.Pt();
        vals[4]  = clusterVec.E();
        vals[5]  = GetMaxCellEnergy(clus);
        vals[6]  = clus->GetNCells();
        vals[7]  = clus->GetDistanceToBadChannel();
        vals[8]  = clus->GetDispersion();
        vals[9]  = clus->GetM20();
        vals[10] = clus->GetM02();
        vals[11] = clus->Chi2();
        vals[12] = clus->GetEmcCpvDistance();
        vals[13] = clusterEcc;
        vals[14] = GetMaxCellEnergy(clus)/clus->E();
        vals[15] = maxAxis;
        vals[16] = clusterVec.Eta();
        vals[17] = clusterVec.Phi();
        fNtuple->Fill(vals);
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::FillPionHists()
{
  // Fill histograms related to pions.

  Double_t vertex[3] = {0,0,0};
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
      if (clus1->E()<0.010)
        continue;
      if (clus1->GetNCells()<fNminCells)
        continue;
      clus1->GetMomentum(clusterVec1,vertex);
      for (Int_t j = i+1; j<nclus; ++j) {
        AliVCluster *clus2 = static_cast<AliVCluster*>(clusters->At(j));
        if (!clus2)
          continue;
        if (!clus2->IsEMCAL()) 
          continue;
        if (clus2->E()<0.010)
          continue;
        if (clus2->GetNCells()<fNminCells)
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
void  AliAnalysisTaskEMCALPi0PbPb::ClusterAfterburner()
{
  // 

  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;

  if (!cells)
    return;

  Int_t ncells = cells->GetNumberOfCells();
  if (ncells<=0)
    return;

  Double_t cellMeanE = 0, cellSigE = 0;
  for (Int_t i = 0; i<ncells; ++i ) {
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
    Int_t nc = clus->GetNCells();
    Double_t oldE = clus->E();
    Double_t clusE = 0;
    Short_t nc2 = 0;
    UShort_t ids[100] = {0};
    Double_t fra[100] = {0};
    for (Int_t j = 0; j<nc; ++j) {
      Short_t id = clus->GetCellAbsId(j);
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
    for (Int_t j = 0; j<nc2; ++j) {
      Short_t id = ids[j];
      Double_t cen = cells->GetCellAmplitude(id);
      fra[j] = cen/clusE;
    }
    clus->SetE(clusE);
    AliAODCaloCluster *aodclus = dynamic_cast<AliAODCaloCluster*>(clus);
    if (aodclus) {
      aodclus->Clear("");
      aodclus->SetNCells(nc2);
      aodclus->SetCellsAmplitudeFraction(fra);
      aodclus->SetCellsAbsId(ids);
    }
    cout << i << " " << clusE << " " << oldE << endl;
  }

  clusters->Compress();
}

// virtual void     AddAt(TObject *obj, Int_t idx);
//   virtual TObject *RemoveAt(Int_t idx);

//________________________________________________________________________
Double_t  AliAnalysisTaskEMCALPi0PbPb::GetMaxCellEnergy(AliVCluster *cluster)
{
  // Get maximum energy of attached cell.

  Double_t maxe = 0;
  Int_t ncells = cluster->GetNCells();
  if (fEsdCells) {
    for (Int_t i=0; i<ncells; i++) {
      Double_t e = fEsdCells->GetCellAmplitude(cluster->GetCellAbsId(i));
      if (e>maxe)
        maxe = e;
    }
  } else {
    for (Int_t i=0; i<ncells; i++) {
      Double_t e = fAodCells->GetCellAmplitude(cluster->GetCellAbsId(i));
      if (e>maxe)
        maxe = e;
    }
  }
  return maxe;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::GetSigma(AliVCluster *cluster, Double_t& sigmaMax, Double_t &sigmaMin)
{
  // Calculate the (E) weighted variance along the longer (eigen) axis.

  sigmaMax = 0;               // cluster variance along its longer axis
  sigmaMin = 0;               // cluster variance along its shorter axis
  Double_t Ec = cluster->E(); // cluster energy
  Double_t Xc = 0 ;           // cluster first moment along X
  Double_t Yc = 0 ;           // cluster first moment along Y
  Double_t Sxx = 0 ;          // cluster second central moment along X (variance_X^2)
  Double_t Sxy = 0 ;          // cluster second central moment along Y (variance_Y^2)
  Double_t Syy = 0 ;          // cluster covariance^2

  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;

  if (!cells)
    return;

  Int_t ncells = cluster->GetNCells();
  if (ncells==1)
    return;

  TVector3 pos;
  for(Int_t j=0; j<ncells; ++j) {
    fGeom->GetGlobal(cluster->GetCellAbsId(j),pos);
    Int_t id = cluster->GetCellAbsId(j);
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
  sigmaMax = (Sxx + Syy + TMath::Sqrt((Sxx-Syy)*(Sxx-Syy)+4.0*Sxy*Sxy))/2.0;
  sigmaMax = TMath::Sqrt(sigmaMax); 
  sigmaMin = TMath::Abs(Sxx + Syy - TMath::Sqrt((Sxx-Syy)*(Sxx-Syy)+4.0*Sxy*Sxy))/2.0;
  sigmaMin = TMath::Sqrt(sigmaMin); 
}

