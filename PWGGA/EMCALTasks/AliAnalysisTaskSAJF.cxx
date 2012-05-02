// $Id$
//
// Jet finder analysis task (S.Aiola).
//
//

#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TParticle.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliESDtrack.h"
#include "AliEmcalJet.h"
#include "AliAODTrack.h"
#include "AliEmcalJet.h"
#include "AliVEventHandler.h"
#include "AliPicoTrack.h"

#include "AliAnalysisTaskSAJF.h"

ClassImp(AliAnalysisTaskSAJF)

//________________________________________________________________________
AliAnalysisTaskSAJF::AliAnalysisTaskSAJF() : 
  AliAnalysisTaskSE("AliAnalysisTaskSAJF"),
  fAnaType(kEMCAL),
  fOutput(0),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fKtJetsName("KtJets"),
  fTrgClusName("ClustersL1GAMMAFEE"),
  fTracks(0),
  fCaloClusters(0),
  fJets(0),
  fKtJets(0),
  fTrgClusters(0),
  fCent(0),
  fHistCentrality(0),
  fHistJetPhiEta(0),
  fHistRhoClusVSleadJetE(0),
  fHistRhoTracksVSleadJetE(0),
  Ptbins(400),
  Ptlow(0),
  Ptup(200),
  Ebins(400),
  Elow(0),
  Eup(200)
{
  // Default constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistJetsE[i] = 0;
    fHistJetsNE[i] = 0;
    fHistJetsNEF[i] = 0;
    fHistJetsZ[i] = 0;
    fHistLeadingJetE[i] = 0;
    fHistTracksPtLJ[i] = 0;
    fHistClusELJ[i] = 0;
    fHistTracksPtBkg[i] = 0;
    fHistClusEBkg[i] = 0;
    fHistBkgClusPhiCorr[i] = 0;
    fHistBkgTracksPhiCorr[i] = 0;
    fHistBkgClusMeanRho[i] = 0;
    fHistBkgTracksMeanRho[i] = 0;
    fHistBkg2JetPhiCorr[i] = 0;
    fHistMedianEnergyKt[i] = 0;
  }

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class()); 
}

//________________________________________________________________________
AliAnalysisTaskSAJF::AliAnalysisTaskSAJF(const char *name) : 
  AliAnalysisTaskSE(name),
  fAnaType(kEMCAL),
  fOutput(0),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fKtJetsName("KtJets"),
  fTrgClusName("ClustersL1GAMMAFEE"),
  fTracks(0),
  fCaloClusters(0),
  fJets(0),
  fKtJets(0),
  fTrgClusters(0),
  fCent(0),
  fHistCentrality(0),
  fHistJetPhiEta(0),
  fHistRhoClusVSleadJetE(0),
  fHistRhoTracksVSleadJetE(0),
  Ptbins(400),
  Ptlow(0),
  Ptup(200),
  Ebins(400),
  Elow(0),
  Eup(200)
{
  // Standard constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistJetsE[i] = 0;
    fHistJetsNE[i] = 0;
    fHistJetsNEF[i] = 0;
    fHistJetsZ[i] = 0;
    fHistLeadingJetE[i] = 0;
    fHistTracksPtLJ[i] = 0;
    fHistClusELJ[i] = 0;
    fHistTracksPtBkg[i] = 0;
    fHistClusEBkg[i] = 0;
    fHistBkgClusPhiCorr[i] = 0;
    fHistBkgTracksPhiCorr[i] = 0;
    fHistBkgClusMeanRho[i] = 0;
    fHistBkgTracksMeanRho[i] = 0;
    fHistBkg2JetPhiCorr[i] = 0;
    fHistMedianEnergyKt[i] = 0;
  }

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class()); 
}

//________________________________________________________________________
AliAnalysisTaskSAJF::~AliAnalysisTaskSAJF()
{
  // Destructor
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::UserCreateOutputObjects()
{
  // Create histograms
  
  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!
  
  fHistCentrality = new TH1F("fHistCentrality","Event centrality distribution", Ebins, 0, 100);
  fHistCentrality->GetXaxis()->SetTitle("Centrality (%)");
  fHistCentrality->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistCentrality);

  fHistJetPhiEta = new TH2F("fHistJetPhiEta","Phi-Eta distribution of jets", 20, -2, 2, 32, 0, 6.4);
  fHistJetPhiEta->GetXaxis()->SetTitle("Eta");
  fHistJetPhiEta->GetYaxis()->SetTitle("Phi");
  fOutput->Add(fHistJetPhiEta);

  fHistRhoClusVSleadJetE = new TH2F("fHistRhoClusVSleadJetE","fHistRhoClusVSleadJetE", Ebins, Elow, Eup, Ebins, Elow, Eup);
  fHistRhoClusVSleadJetE->GetXaxis()->SetTitle("#rho [GeV]");
  fHistRhoClusVSleadJetE->GetYaxis()->SetTitle("Leading jet energy [GeV]");
  fOutput->Add(fHistRhoClusVSleadJetE);

  fHistRhoTracksVSleadJetE = new TH2F("fHistRhoTracksVSleadJetE","fHistRhoTracksVSleadJetE", Ptbins, Ptlow, Ptup, Ebins, Elow, Eup);
  fHistRhoTracksVSleadJetE->GetXaxis()->SetTitle("#rho [GeV/c]");
  fHistRhoTracksVSleadJetE->GetYaxis()->SetTitle("Leading jet energy");
  fOutput->Add(fHistRhoTracksVSleadJetE);

  TString histname;

  for (Int_t i = 0; i < 4; i++) {
    histname = "fHistJetsE_";
    histname += i;
    fHistJetsE[i] = new TH1F(histname.Data(), histname.Data(), Ebins, Elow, Eup);
    fHistJetsE[i]->GetXaxis()->SetTitle("E [GeV]");
    fHistJetsE[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistJetsE[i]);
    
    histname = "fHistJetsNE_";
    histname += i;
    fHistJetsNE[i] = new TH1F(histname.Data(), histname.Data(), Ebins, Elow, Eup);
    fHistJetsNE[i]->GetXaxis()->SetTitle("E [GeV]");
    fHistJetsNE[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistJetsNE[i]);
    
    histname = "fHistJetsNEF_";
    histname += i;
    fHistJetsNEF[i] = new TH1F(histname.Data(), histname.Data(), Ebins, 0, 1.2);
    fHistJetsNEF[i]->GetXaxis()->SetTitle("NEF");
    fHistJetsNEF[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistJetsNEF[i]);

    histname = "fHistJetsZ_";
    histname += i;
    fHistJetsZ[i] = new TH1F(histname.Data(), histname.Data(), Ebins, 0, 1.2);
    fHistJetsZ[i]->GetXaxis()->SetTitle("Z");
    fHistJetsZ[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistJetsZ[i]);

    histname = "fHistLeadingJetE_";
    histname += i;
    fHistLeadingJetE[i] = new TH1F(histname.Data(), histname.Data(), Ebins, Elow, Eup);
    fHistLeadingJetE[i]->GetXaxis()->SetTitle("E [GeV]");
    fHistLeadingJetE[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistLeadingJetE[i]);
    
    histname = "fHistTracksPtLJ_";
    histname += i;
    fHistTracksPtLJ[i] = new TH1F(histname.Data(), histname.Data(), Ptbins, Ptlow, Ptup);
    fHistTracksPtLJ[i]->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    fHistTracksPtLJ[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistTracksPtLJ[i]);
    
    histname = "fHistClusELJ_";
    histname += i;
    fHistClusELJ[i] = new TH1F(histname.Data(), histname.Data(), Ebins, Elow, Eup);
    fHistClusELJ[i]->GetXaxis()->SetTitle("E [GeV]");
    fHistClusELJ[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistClusELJ[i]);

    histname = "fHistTracksPtBkg_";
    histname += i;
    fHistTracksPtBkg[i] = new TH1F(histname.Data(), histname.Data(), Ptbins, Ptlow, Ptup);
    fHistTracksPtBkg[i]->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    fHistTracksPtBkg[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistTracksPtBkg[i]);
    
    histname = "fHistClusEBkg_";
    histname += i;
    fHistClusEBkg[i] = new TH1F(histname.Data(), histname.Data(), Ebins, Elow, Eup);
    fHistClusEBkg[i]->GetXaxis()->SetTitle("E [GeV]");
    fHistClusEBkg[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistClusEBkg[i]);

    histname = "fHistBkgClusPhiCorr_";
    histname += i;
    fHistBkgClusPhiCorr[i] = new TH1F(histname.Data(), histname.Data(), 128, -1.6, 4.8);
    fHistBkgClusPhiCorr[i]->GetXaxis()->SetTitle("#Delta#phi");
    fHistBkgClusPhiCorr[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistBkgClusPhiCorr[i]);

    histname = "fHistBkgTracksPhiCorr_";
    histname += i;
    fHistBkgTracksPhiCorr[i] = new TH1F(histname.Data(), histname.Data(), 128, -1.6, 4.8);
    fHistBkgTracksPhiCorr[i]->GetXaxis()->SetTitle("#Delta#phi");
    fHistBkgTracksPhiCorr[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistBkgTracksPhiCorr[i]);

    histname = "fHistBkgClusMeanRho_";
    histname += i;
    fHistBkgClusMeanRho[i] = new TH1F(histname.Data(), histname.Data(), Ebins, Elow, Eup);
    fHistBkgClusMeanRho[i]->GetXaxis()->SetTitle("E [GeV]");
    fHistBkgClusMeanRho[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistBkgClusMeanRho[i]);

    histname = "fHistBkgTracksMeanRho_";
    histname += i;
    fHistBkgTracksMeanRho[i] = new TH1F(histname.Data(), histname.Data(), Ptbins, Ptlow, Ptup);
    fHistBkgTracksMeanRho[i]->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    fHistBkgTracksMeanRho[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistBkgTracksMeanRho[i]);

    histname = "fHistBkg2JetPhiCorr_";
    histname += i;
    fHistBkg2JetPhiCorr[i] = new TH1F(histname.Data(), histname.Data(), 128, -1.6, 4.8);
    fHistBkg2JetPhiCorr[i]->GetXaxis()->SetTitle("#Delta#phi");
    fHistBkg2JetPhiCorr[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistBkg2JetPhiCorr[i]);

    histname = "fHistMedianEnergyKt_";
    histname += i;
    fHistMedianEnergyKt[i] = new TH1F(histname.Data(), histname.Data(), Ebins, Elow, Eup);
    fHistMedianEnergyKt[i]->GetXaxis()->SetTitle("E [GeV]");
    fHistMedianEnergyKt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistMedianEnergyKt[i]);
  }

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

void AliAnalysisTaskSAJF::RetrieveEventObjects()
{
  fCaloClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fCaloName));
  if (!fCaloClusters) {
    AliWarning(Form("Could not retrieve clusters!")); 
  }

  fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName));
  if (!fTracks) {
    AliWarning(Form("Could not retrieve tracks!")); 
  }

  fJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsName));
  if (!fJets) {
    AliWarning(Form("Could not retrieve jets!")); 
  }

  if (strcmp(fKtJetsName,"")) {
    fKtJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fKtJetsName));
    if (!fKtJets) {
      AliWarning(Form("Could not retrieve Kt jets!")); 
    }
  }
    
  if (strcmp(fTrgClusName,"")) {
    fTrgClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTrgClusName));
    if (!fTrgClusters) {
      AliWarning(Form("Could not retrieve trigger clusters!")); 
    }
  }

  fCent = InputEvent()->GetCentrality();
  if (!fCent) {
    AliWarning(Form("Could not retrieve centrality information!")); 
  }
}

AliVTrack* AliAnalysisTaskSAJF::GetTrack(const Int_t i) const
{
  if (fTracks)
    return dynamic_cast<AliVTrack*>(fTracks->At(i));
  else
    return 0;
}

Int_t AliAnalysisTaskSAJF::GetNumberOfTracks() const
{
  if (fTracks)
    return fTracks->GetEntriesFast();
  else
    return 0;
}

AliVCluster* AliAnalysisTaskSAJF::GetCaloCluster(const Int_t i) const
{ 
  if (fCaloClusters)
    return dynamic_cast<AliVCluster*>(fCaloClusters->At(i));
  else
    return 0;
}

Int_t AliAnalysisTaskSAJF::GetNumberOfCaloClusters() const
{ 
  if (fCaloClusters)
    return fCaloClusters->GetEntriesFast();
  else
    return 0;
}

AliEmcalJet* AliAnalysisTaskSAJF::GetJet(const Int_t i) const
{
  if (fJets)
    return dynamic_cast<AliEmcalJet*>(fJets->At(i));
  else
    return 0;
}

Int_t AliAnalysisTaskSAJF::GetNumberOfJets() const
{
  if (fJets)
    return fJets->GetEntriesFast();
  else
    return 0;
}

AliEmcalJet* AliAnalysisTaskSAJF::GetKtJet(const Int_t i) const
{
  if (fKtJets)
    return dynamic_cast<AliEmcalJet*>(fKtJets->At(i));
  else
    return 0;
}

Int_t AliAnalysisTaskSAJF::GetNumberOfKtJets() const
{
  if (fKtJets)
    return fKtJets->GetEntriesFast();
  else
    return 0;
}

AliVCluster* AliAnalysisTaskSAJF::GetTrgCluster(const Int_t i) const
{
  if (fTrgClusters)
    return dynamic_cast<AliVCluster*>(fTrgClusters->At(i));
  else
    return 0;
}

Int_t AliAnalysisTaskSAJF::GetNumberOfTrgClusters() const
{
  if (fTrgClusters)
    return fTrgClusters->GetEntriesFast();
  else
    return 0;
}

void AliAnalysisTaskSAJF::FillHistograms()
{
  Float_t cent = 100;
  
  if (fCent)
    cent = fCent->GetCentralityPercentile("V0M");
  else
    AliWarning("Centrality not available!");

  fHistCentrality->Fill(cent);

  Int_t centbin=-1;
  if(cent>=0 && cent<10) centbin=0;
  else if(cent>=10 && cent<30) centbin=1;
  else if(cent>=30 && cent<50) centbin=2;
  else if(cent>=50 && cent<=100) centbin=3;

  // Kt Jet loop
  Int_t nktjets =  GetNumberOfKtJets();

  if (nktjets > 0) {
    
    TArrayF ktJets(nktjets);
    for (Int_t ij = 0; ij < nktjets; ij++) {
      
      AliEmcalJet* jet = GetKtJet(ij);
      
      if (!jet) {
	printf("ERROR: Could not receive jet %d\n", ij);
	continue;
      }  
      
      if (jet->E() <= 0)
	continue;
      
      if (!AcceptJet(jet))
	continue;
      
      Int_t i = nktjets - 1;
      while (jet->E() < ktJets[i] && i > 0)
	i--;
      memmove(ktJets.GetArray() + nktjets - ij - 1, ktJets.GetArray() + nktjets - ij, (ij + i - nktjets + 1) * sizeof(Float_t));
      ktJets[i] = jet->E();
      
    } //kt jet loop 

    Float_t ktJetsMedian;
    if (nktjets % 2)
      ktJetsMedian = ktJets[nktjets / 2 - 1];
    else
      ktJetsMedian = (ktJets[nktjets / 2 - 1] + ktJets[nktjets / 2]) / 2;
    
    fHistMedianEnergyKt[centbin]->Fill(ktJetsMedian);
  }

  // Jet loop
  Int_t njets =  GetNumberOfJets();
  //cout << njets << " jets" << endl;
  Float_t maxJetEnergy = 0;
  Int_t maxJetIndex = -1;
  Float_t max2JetEnergy = 0;
  Int_t max2JetIndex = -1;
  for (Int_t ij = 0; ij < njets; ij++) {

    AliEmcalJet* jet = GetJet(ij);

    if (!jet) {
      printf("ERROR: Could not receive jet %d\n", ij);
      continue;
    }  
    
    if (jet->E() <= 0)
      continue;

    if (!AcceptJet(jet))
      continue;

    fHistJetPhiEta->Fill(jet->Eta(), jet->Phi());

    fHistJetsE[centbin]->Fill(jet->E());
    fHistJetsNEF[centbin]->Fill(jet->NEF());
    fHistJetsNE[centbin]->Fill(jet->E() * jet->NEF());

    //cout << "****** jet id = " << ij << " ******" << endl;

    for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
      Int_t trackid = jet->TrackAt(it);
      AliVTrack *track = GetTrack(trackid);
      //cout << "track id = " << trackid << endl;
      if (track)
	fHistJetsZ[centbin]->Fill(track->Pt() / jet->E());
    }
    for (Int_t ic = 0; ic < jet->GetNumberOfClusters(); ic++) {
      Int_t clusterid = jet->ClusterAt(ic);
      //cout << "cluster id = " << clusterid << endl;
      AliVCluster *cluster = GetCaloCluster(clusterid);
      if (cluster)
	fHistJetsZ[centbin]->Fill(cluster->E() / jet->E());
    }

    if (maxJetEnergy < jet->E()) {
      max2JetEnergy = maxJetEnergy;
      max2JetIndex = maxJetIndex;
      maxJetEnergy = jet->E();
      maxJetIndex = ij;
    }
    else if (max2JetEnergy < jet->E()) {
      max2JetEnergy = jet->E();
      max2JetIndex = ij;
    }
  } //jet loop 

  
  if (!(maxJetEnergy > 0) || maxJetIndex < 0) 
    return;

  fHistLeadingJetE[centbin]->Fill(maxJetEnergy);
  
  AliEmcalJet* jet = GetJet(maxJetIndex);
  if (!jet)
    return;

  AliEmcalJet* jet2 = GetJet(max2JetIndex);

  jet->SortConstituents();
  if (jet2) jet2->SortConstituents();
  
  // Track loop 
  Float_t rhoTracks = 0;
  Int_t ntracks = GetNumberOfTracks();
  for(Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliVTrack* track = GetTrack(iTracks);         
    if(!track) {
      AliError(Form("ERROR: Could not retrieve track %d",iTracks)); 
      continue; 
    }
    
    if (!AcceptTrack(track)) continue;
    
    if (IsJetTrack(jet, iTracks)) {
      fHistTracksPtLJ[centbin]->Fill(track->Pt());
    }
    else if (!jet2 || !IsJetTrack(jet2, iTracks)) {
      fHistTracksPtBkg[centbin]->Fill(track->Pt());
      rhoTracks += track->Pt();
      
      Float_t dphijet = jet->Phi() - track->Phi();
      if (dphijet < -1.6) dphijet += TMath::Pi() * 2;
      if (dphijet > 4.8) dphijet -= TMath::Pi() * 2;
      fHistBkg2JetPhiCorr[centbin]->Fill(dphijet);

      for(Int_t it2 = iTracks+1; it2 < ntracks; it2++) {
	AliVTrack* track2 = GetTrack(it2);         
	if(!track2) {
	  AliError(Form("ERROR: Could not retrieve track %d", it2)); 
	  continue; 
	}
	
	if (!AcceptTrack(track2)) continue;
	
	if (IsJetTrack(jet, it2)) continue;

	if (jet2 && IsJetTrack(jet2, it2)) continue;
	
	Float_t dphi = track->Phi() - track2->Phi();
	if (dphi < -1.6) dphi += TMath::Pi() * 2;
	if (dphi > 4.8) dphi -= TMath::Pi() * 2;
	fHistBkgTracksPhiCorr[centbin]->Fill(dphi);
      } // second track loop
    }
  } // track loop

  rhoTracks /= GetArea();
  fHistBkgTracksMeanRho[centbin]->Fill(rhoTracks);

  fHistRhoTracksVSleadJetE->Fill(rhoTracks, maxJetEnergy);

  if (fAnaType == kFullAcceptance) return;
  
  // Cluster loop
  Float_t rhoClus = 0;
  Int_t nclusters =  GetNumberOfCaloClusters();
  for (Int_t iClusters = 0; iClusters < nclusters; iClusters++) {
    AliVCluster* cluster = GetCaloCluster(iClusters);
    if (!cluster) {
      printf("ERROR: Could not receive cluster %d\n", iClusters);
      continue;
    }  
    
    if (!(cluster->IsEMCAL())) continue;
    
    if (IsJetCluster(jet, iClusters)) {
      fHistClusELJ[centbin]->Fill(cluster->E());
    }
    else if (!jet2 || !IsJetCluster(jet2, iClusters)) {
      fHistClusEBkg[centbin]->Fill(cluster->E());
      rhoClus += cluster->E();

      Float_t pos1[3];
      cluster->GetPosition(pos1);
      TVector3 clusVec1(pos1);

      Float_t dphijet = jet->Phi() - clusVec1.Phi();
      if (dphijet < -1.6) dphijet += TMath::Pi() * 2;
      if (dphijet > 4.8) dphijet -= TMath::Pi() * 2;
      fHistBkg2JetPhiCorr[centbin]->Fill(dphijet);

      for(Int_t ic2 = iClusters+1; ic2 < nclusters; ic2++) {
	AliVCluster* cluster2 = GetCaloCluster(ic2);
	if (!cluster2) {
	  printf("ERROR: Could not receive cluster %d\n", ic2);
	  continue;
	}  
	
	if (!(cluster2->IsEMCAL())) continue;
	
	if (IsJetCluster(jet, ic2)) continue;
	
	if (jet2 && IsJetCluster(jet2, ic2)) continue;

	Float_t pos2[3];
	cluster2->GetPosition(pos2);
	TVector3 clusVec2(pos2);

	Float_t dphi = clusVec1.Phi() - clusVec2.Phi();
	if (dphi < -1.6) dphi += TMath::Pi() * 2;
	if (dphi > 4.8) dphi -= TMath::Pi() * 2;
	fHistBkgClusPhiCorr[centbin]->Fill(dphi);
      }
    }
  } //cluster loop 
  
  rhoClus /= GetArea();
  fHistBkgClusMeanRho[centbin]->Fill(rhoClus);

  fHistRhoClusVSleadJetE->Fill(rhoClus, maxJetEnergy);
}

//________________________________________________________________________
Float_t AliAnalysisTaskSAJF::GetArea() const
{
  if (fAnaType == kFullAcceptance) {
    return 2 * TMath::Pi() * 2;
  }
  else if (fAnaType == kEMCAL) {
    return 1.4 * 100 * TMath::DegToRad();
  }
  else if (fAnaType == kEMCAL) {
    return 1.4 * 100 * TMath::DegToRad();
  }
  else {
    AliWarning("Analysis type not recognized! Assuming kFullAcceptance...");
    return 2 * TMath::Pi() * 2;
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted) const
{
  for (Int_t i = 0; i < jet->GetNumberOfTracks(); i++) {
    Int_t ijettrack = jet->TrackAt(i);
    if (sorted && ijettrack > itrack)
      return kFALSE;
    if (ijettrack == itrack)
      return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted) const
{
  for (Int_t i = 0; i < jet->GetNumberOfClusters(); i++) {
    Int_t ijetclus = jet->ClusterAt(i);
    if (sorted && ijetclus > iclus)
      return kFALSE;
    if (ijetclus == iclus)
      return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::AcceptJet(AliEmcalJet* jet) const
{
  if (fAnaType == kFullAcceptance) {
    return kTRUE;
  }
  else if (fAnaType == kEMCAL) {
    return (Bool_t)(TMath::Abs(jet->Eta()) < 0.7 && jet->Phi() * TMath::RadToDeg() > 80 && jet->Phi() * TMath::RadToDeg() < 180);
  }
  else if (fAnaType == kEMCALFiducial) {
    return (Bool_t)(TMath::Abs(jet->Eta()) < 0.7 && jet->Phi() * TMath::RadToDeg() > 80 && jet->Phi() * TMath::RadToDeg() < 180);
  }
  else {
    AliWarning("Analysis type not recognized! Assuming kFullAcceptance...");
    return kTRUE;
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::AcceptTrack(AliVTrack* track) const
{
  if (fAnaType == kFullAcceptance) {
    return kTRUE;
  }
  else if (fAnaType == kEMCAL) {
    return (Bool_t)(TMath::Abs(track->Eta()) < 0.7 && track->Phi() * TMath::RadToDeg() > 80 && track->Phi() * TMath::RadToDeg() < 180);
  }
  else if (fAnaType == kEMCALFiducial) {
    return (Bool_t)(TMath::Abs(track->Eta()) < 0.7 && track->Phi() * TMath::RadToDeg() > 80 && track->Phi() * TMath::RadToDeg() < 180);
  }
  else {
    AliWarning("Analysis type not recognized! Assuming kFullAcceptance...");
    return kTRUE;
  }
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::UserExec(Option_t *) 
{
  // Main loop, called for each event.
  // Add jets to event if not yet there

  RetrieveEventObjects();

  FillHistograms();
    
  // information for this iteration of the UserExec in the container
  PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
