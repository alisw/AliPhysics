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
    fHClustEtaPhi(0x0),
    fHClustEnergyPt(0x0),
    fHClustEnergySigma(0x0),
    fHClustSigmaSigma(0x0),
    fHClustNTowEnergyRatio(0x0),
    fHPionEtaPhi(0x0),
    fHPionMggPt(0x0),
    fHPionMggAsym(0x0)
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
      fNtuple = new TNtuple("nt","nt","cent:m:pt:e1:e2:e1m:e2m:n1:n2:db1:db2:disp1:disp2:mm1:mm2:ms1:ms2");
  }

  // histograms
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
  fHCellE = new TH1F("hCellE","",100,0.,10.);
  fHCellE->SetXTitle("E_{cell} [GeV]");
  fOutput->Add(fHCellE);
  fHCellH = new TH1F ("fHCellHighestE","",100,0.,10.);
  fHCellH->SetXTitle("E^{max}_{cell} [GeV] ");
  fOutput->Add(fHCellH);

  // histograms for clusters
  fHClustEtaPhi = new TH2F("hClustEtaPhi","",100,-0.8,0.8,100,1.2,2.2);
  fHClustEtaPhi->SetXTitle("#eta");
  fHClustEtaPhi->SetYTitle("#varphi");
  fOutput->Add(fHClustEtaPhi);
  fHClustEnergyPt = new TH2F("hClustEnergyPt","",250,0,50,250,0,50);
  fHClustEnergyPt->SetXTitle("E [GeV]");
  fHClustEnergyPt->SetYTitle("p_{T}^{} [GeV/c]");
  fOutput->Add(fHClustEnergyPt);
  fHClustEnergySigma = new TH2F("hClustEnergySigma","",100,0,100,100,0,30);
  fHClustEnergySigma->SetXTitle("E_{C}^{}*#sigma_{max}^{} [GeV*cm]");
  fHClustEnergySigma->SetYTitle("E_{C}^{} [GeV]");
  fOutput->Add(fHClustEnergySigma);
  fHClustSigmaSigma = new TH2F("hClustSigmaSigma","",100,0,100,100,0,20);
  fHClustSigmaSigma->SetXTitle("#lambda_{0}^{} [cm]");
  fHClustSigmaSigma->SetYTitle("#sigma_{max}^{} [cm]");
  fOutput->Add(fHClustSigmaSigma);
  fHClustNTowEnergyRatio = new TH2F("hClustNTowEnergyRatio","",15,-0.5,14.5,101,-0.05,1.05);
  fHClustNTowEnergyRatio->SetXTitle("N towers");
  fHClustNTowEnergyRatio->SetYTitle("Energy ratio");
  fOutput->Add(fHClustNTowEnergyRatio);

  // histograms for pion candidates
  fHPionEtaPhi = new TH2F("hPionEtaPhi","",100,-0.8,0.8,100,1.2,2.2);
  fHPionEtaPhi->SetXTitle("#eta^{#gamma#gamma}");
  fHPionEtaPhi->SetYTitle("#varphi^{#gamma#gamma}");
  fOutput->Add(fHPionEtaPhi);
  fHPionMggPt = new TH2F("hPionMggPt","",100,0,2,100,0,20.0);
  fHPionMggPt->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
  fHPionMggPt->SetYTitle("p_{T}^{#gamma#gamma} [GeV/c]");
  fOutput->Add(fHPionMggPt);
  fHPionMggAsym = new TH2F("hPionMggAsym","",100,0,2,100,0,1);
  fHPionMggAsym->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
  fHPionMggAsym->SetYTitle("Z_{#gamma#gamma} [GeV]");
  fOutput->Add(fHPionMggAsym);
  const Int_t nbins = 19; 
  Double_t xbins[nbins] = {0.5,1,1.5,2,2.5,3,3.5,4,4.5,6,7,8,9,10,12.5,15,20,25,50};
  fPtRanges = new TAxis(nbins-1,xbins);
  for (Int_t i = 0; i<=nbins; ++i) {
    fHPionInvMasses[i] = new TH1F(Form("HPionInvMass%d",i),"",100,0,2);
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
    Int_t ncells = cells->GetNumberOfCells();

    for (Int_t i = 0; i<ncells; ++i ) {
      Int_t absID    = cells->GetCellNumber(i);
      Double_t cellE = cells->GetAmplitude(i);
      fHCellE->Fill(cellE);
      if (cellE>cellMaxE) 
        cellMaxE = cellE;

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
    for (Int_t i=0; i<4; ++i) 
      fHCellMult[i]->Fill(cellModCount[i]);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::FillClusHists()
{
  // Fill histograms related to cluster properties.

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
      
      fHClustEtaPhi->Fill(clusterVec.Eta(),clusterVec.Phi());
      fHClustEnergyPt->Fill(clusterVec.E(),clusterVec.Pt());
      fHClustEnergySigma->Fill(clus->E()*GetSigmaMax(clus),clus->E());
      fHClustSigmaSigma->Fill(max(clus->GetM02(),clus->GetM20()),clus->E()*GetSigmaMax(clus));
      fHClustNTowEnergyRatio->Fill(clus->GetNCells(),GetMaxCellEnergy(clus));
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
      clus1->GetMomentum(clusterVec1,vertex);
      for (Int_t j = i+1; j<nclus; ++j) {
        AliVCluster *clus2 = static_cast<AliVCluster*>(clusters->At(j));
        if (!clus2)
          continue;
        if (!clus2->IsEMCAL()) 
          continue;
        clus2->GetMomentum(clusterVec2,vertex);
        pionVec = clusterVec1 + clusterVec2;
        Double_t pionZgg = TMath::Abs(clusterVec1.E()-clusterVec2.E())/pionVec.E();
        if (pionZgg < fAsymMax) {
          fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi()); 
          fHPionMggPt->Fill(pionVec.M(),pionVec.Pt()); 
          fHPionMggAsym->Fill(pionVec.M(),pionZgg); 
          Int_t bin = fPtRanges->FindBin(pionVec.Pt());
          fHPionInvMasses[bin]->Fill(pionVec.M());
        }

        if (fNtuple) {
          Double_t mass = pionVec.M();
          if (mass>0.08 && mass<0.2) {
            Float_t vals[17];
            vals[0]  = InputEvent()->GetCentrality()->GetCentralityPercentileUnchecked(fCentVar);
            vals[1]  = mass;
            vals[2]  = pionVec.Pt();
            vals[3]  = clusterVec1.E();
            vals[4]  = clusterVec2.E();
            vals[5]  = GetMaxCellEnergy(clus1);
            vals[6]  = GetMaxCellEnergy(clus2);
            vals[7]  = clus1->GetNCells();
            vals[8]  = clus2->GetNCells();
            vals[9]  = clus1->GetDistanceToBadChannel();
            vals[10] = clus2->GetDistanceToBadChannel();
            vals[11] = clus1->GetDispersion();
            vals[12] = clus2->GetDispersion();
            vals[13] = clus1->GetM20();
            vals[14] = clus2->GetM20();
            vals[15] = clus1->GetM02();
            vals[16] = clus2->GetM02();
            fNtuple->Fill(vals);
          }
        }
      }
    }
  } 
}

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
Double_t  AliAnalysisTaskEMCALPi0PbPb::GetSigmaMax(AliVCluster *cluster)
{
  // Calculate the (E) weighted variance along the longer (eigen) axis.

  Double_t sigmaMax = 0;      // cluster variance along its longer axis
  Double_t Ec = cluster->E(); // cluster energy
  Double_t Xc = 0 ;           // cluster first moment along X
  Double_t Yc = 0 ;           // cluster first moment along Y
  Double_t Sxx = 0 ;          // cluster second central moment along X (variance_X^2)
  Double_t Sxy = 0 ;          // cluster second central moment along Y (variance_Y^2)
  Double_t Syy = 0 ;          // cluster covariance^2

  Double_t testenergy = 0;
  AliVCaloCells *cells = fEsdCells;
  if (!cells)
    cells = fAodCells;

  if (cells) {
    TVector3 pos;
    Int_t ncells = cluster->GetNCells();
    if (ncells==1)
      return 0;
    for(Int_t j=0; j<ncells; ++j) {
      fGeom->GetGlobal(cluster->GetCellAbsId(j),pos);
      Int_t id = cluster->GetCellAbsId(j);
      Double_t cellen = cells->GetCellAmplitude(id);
      testenergy += cellen;
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

    printf("%f %f\n",testenergy,cluster->E());
    cout << testenergy << " " << cluster->E() << endl;
  } 
  return sigmaMax;
}

