#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPicoTrack.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSoftDrop.h"

ClassImp(AliAnalysisTaskSoftDrop)

//________________________________________________________________________
AliAnalysisTaskSoftDrop::AliAnalysisTaskSoftDrop() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskSoftDrop", kTRUE),
  fHistTracksPt(0),
  fHistClustersPt(0),
  fHistLeadingJetPt(0),
  fHistJetsPhiEta(0),
  fHistJetsPtArea(0),
  fHistJetsPtLeadHad(0),
  fHistJetsCorrPtArea(0),
  fHistPtDEtaDPhiTrackClus(0),
  fHistPtDEtaDPhiClusTrack(0),
  fHistNTracks(0),
  fHistClustDx(0),
  fHistClustDz(0),
  fNAccJets(0),
  fhZg(0),
  fJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0)

{
  // Default constructor.

  fHistTracksPt       = new TH1*[fNcentBins];
  fHistNTracks        = new TH1*[fNcentBins];
  fHistClustersPt     = new TH1*[fNcentBins];
  fHistLeadingJetPt   = new TH1*[fNcentBins];
  fHistJetsPhiEta     = new TH2*[fNcentBins];
  fHistJetsPtArea     = new TH2*[fNcentBins];
  fHistJetsPtLeadHad  = new TH2*[fNcentBins];
  fHistJetsCorrPtArea = new TH2*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fHistTracksPt[i] = 0;
    fHistNTracks[i] = 0;
    fHistClustersPt[i] = 0;
    fHistLeadingJetPt[i] = 0;
    fHistJetsPhiEta[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistJetsPtLeadHad[i] = 0;
    fHistJetsCorrPtArea[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskSoftDrop::AliAnalysisTaskSoftDrop(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistTracksPt(0),
  fHistClustersPt(0),
  fHistLeadingJetPt(0),
  fHistJetsPhiEta(0),
  fHistJetsPtArea(0),
  fHistJetsPtLeadHad(0),
  fHistJetsCorrPtArea(0),
  fHistPtDEtaDPhiTrackClus(0),
  fHistPtDEtaDPhiClusTrack(0),
  fHistNTracks(0),
  fHistClustDx(0),
  fHistClustDz(0),
  fNAccJets(0),
  fhZg(0),
  fJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0)
{
  // Standard constructor.

  fHistTracksPt       = new TH1*[fNcentBins];
  fHistNTracks        = new TH1*[fNcentBins];
  fHistClustersPt     = new TH1*[fNcentBins];
  fHistLeadingJetPt   = new TH1*[fNcentBins];
  fHistJetsPhiEta     = new TH2*[fNcentBins];
  fHistJetsPtArea     = new TH2*[fNcentBins];
  fHistJetsPtLeadHad  = new TH2*[fNcentBins];
  fHistJetsCorrPtArea = new TH2*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fHistTracksPt[i] = 0;
    fHistNTracks[i] = 0;
    fHistClustersPt[i] = 0;
    fHistLeadingJetPt[i] = 0;
    fHistJetsPhiEta[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistJetsPtLeadHad[i] = 0;
    fHistJetsCorrPtArea[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskSoftDrop::~AliAnalysisTaskSoftDrop()
{
  // Destructor.
}

AliAnalysisTaskSoftDrop* AliAnalysisTaskSoftDrop::AddTaskSoftDrop(
  const char *ntracks,
  const char *nclusters,
  const char *njets,
  const char *nrho,
  Int_t       nCentBins,
  Double_t    jetradius,
  Double_t    jetptcut,
  Double_t    jetareacut,
  const char *type,
  Int_t       leadhadtype,
  const char *taskname
)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetSample", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetSample", "This task requires an input event handler");
    return NULL;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString trackName(ntracks);
  TString clusName(nclusters);

  TString name(taskname);
  if (strcmp(njets,"")) {
    name += "_";
    name += njets;
  }
  if (strcmp(nrho,"")) {
    name += "_";
    name += nrho;
  }
  if (strcmp(type,"")) {
    name += "_";
    name += type;
  }

  Printf("name: %s",name.Data());

  AliAnalysisTaskSoftDrop* jetTask = new AliAnalysisTaskSoftDrop(name);
  //jetTask->SetCentRange(0.,100.);
  //jetTask->SetNCentBins(nCentBins);

  AliParticleContainer* partCont = 0;
  if (trackName == "mcparticles") {
    partCont = jetTask->AddMCParticleContainer(ntracks);
  }
  else if (trackName == "tracks" || trackName == "Tracks") {
    partCont = jetTask->AddTrackContainer(ntracks);
  }
  else if (!trackName.IsNull()) {
    partCont = new AliParticleContainer(trackName);
  }

  AliClusterContainer *clusterCont = jetTask->AddClusterContainer(nclusters);

  TString strType(type);
  AliJetContainer *jetCont = jetTask->AddJetContainer(njets,strType,jetradius);
  if(jetCont) {
    jetCont->SetRhoName(nrho);
    jetCont->ConnectParticleContainer(partCont);
    jetCont->ConnectClusterContainer(clusterCont);
    //jetCont->SetZLeadingCut(0.98,0.98);
    jetCont->SetPercAreaCut(jetareacut);
    jetCont->SetJetPtCut(jetptcut);
    jetCont->SetLeadingHadronType(leadhadtype);
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(jetTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  mgr->ConnectOutput (jetTask, 1, coutput1 );

  return jetTask;
}

//________________________________________________________________________
void AliAnalysisTaskSoftDrop::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fJetsCont           = GetJetContainer(0);
  if(fJetsCont) { //get particles and clusters connected to jets
    fTracksCont       = fJetsCont->GetParticleContainer();
    fCaloClustersCont = fJetsCont->GetClusterContainer();
  } else {        //no jets, just analysis tracks and clusters
    fTracksCont       = GetParticleContainer(0);
    fCaloClustersCont = GetClusterContainer(0);
  }
  if(fTracksCont) fTracksCont->SetClassName("AliVTrack");
  if(fCaloClustersCont) fCaloClustersCont->SetClassName("AliVCluster");

  TString histname;

  for (Int_t i = 0; i < fNcentBins; i++) {
    if (fParticleCollArray.GetEntriesFast()>0) {
      histname = "fHistTracksPt_";
      histname += i;
      fHistTracksPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistTracksPt[i]->GetXaxis()->SetTitle("p_{T,track} (GeV/c)");
      fHistTracksPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistTracksPt[i]);
      
      histname = "fHistNTracks_";
      histname += i;
      fHistNTracks[i] = new TH1F(histname.Data(), histname.Data(), 200, 0., 199.);
      fOutput->Add(fHistNTracks[i]);
    }

    if (fClusterCollArray.GetEntriesFast()>0) {
      histname = "fHistClustersPt_";
      histname += i;
      fHistClustersPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistClustersPt[i]->GetXaxis()->SetTitle("p_{T,clus} (GeV/c)");
      fHistClustersPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistClustersPt[i]);
    }

    if (fJetCollArray.GetEntriesFast()>0) {
      histname = "fHistLeadingJetPt_";
      histname += i;
      fHistLeadingJetPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
      fHistLeadingJetPt[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
      fHistLeadingJetPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistLeadingJetPt[i]);
      
      histname = "fHistJetsPhiEta_";
      histname += i;
      fHistJetsPhiEta[i] = new TH2F(histname.Data(), histname.Data(), 50, -1, 1, 101, 0, TMath::Pi()*2 + TMath::Pi()/200);
      fHistJetsPhiEta[i]->GetXaxis()->SetTitle("#eta");
      fHistJetsPhiEta[i]->GetYaxis()->SetTitle("#phi");
      fOutput->Add(fHistJetsPhiEta[i]);
      
      histname = "fHistJetsPtArea_";
      histname += i;
      fHistJetsPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, 50, 0, 1);
      fHistJetsPtArea[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
      fHistJetsPtArea[i]->GetYaxis()->SetTitle("area");
      fOutput->Add(fHistJetsPtArea[i]);

      histname = "fHistJetsPtLeadHad_";
      histname += i;
      fHistJetsPtLeadHad[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistJetsPtLeadHad[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
      fHistJetsPtLeadHad[i]->GetYaxis()->SetTitle("p_{T,lead} (GeV/c)");
      fHistJetsPtLeadHad[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistJetsPtLeadHad[i]);
    
      if (!(GetJetContainer()->GetRhoName().IsNull())) {
	histname = "fHistJetsCorrPtArea_";
	histname += i;
	fHistJetsCorrPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins*2, -fMaxBinPt, fMaxBinPt, 50, 0, 1);
	fHistJetsCorrPtArea[i]->GetXaxis()->SetTitle("p_{T}^{corr} [GeV/c]");
	fHistJetsCorrPtArea[i]->GetYaxis()->SetTitle("area");
	fOutput->Add(fHistJetsCorrPtArea[i]);
      }
    }
  }

  histname = "fHistPtDEtaDPhiTrackClus";
  fHistPtDEtaDPhiTrackClus = new TH3F(histname.Data(),Form("%s;#it{p}_{T}^{track};#Delta#eta;#Delta#varphi",histname.Data()),100,0.,100.,100,-0.1,0.1,100,-0.1,0.1);
  fOutput->Add(fHistPtDEtaDPhiTrackClus);

  histname = "fHistPtDEtaDPhiClusTrack";
  fHistPtDEtaDPhiClusTrack = new TH3F(histname.Data(),Form("%s;#it{p}_{T}^{clus};#Delta#eta;#Delta#varphi",histname.Data()),100,0.,100.,100,-0.1,0.1,100,-0.1,0.1);
  fOutput->Add(fHistPtDEtaDPhiClusTrack);

  fHistClustDx = new TH1F("fHistClustDx","fHistClustDx;Dx",1000,0.,1.);
  fOutput->Add(fHistClustDx);

  fHistClustDz = new TH1F("fHistClustDz","fHistClustDz;Dz",1000,0.,1.);
  fOutput->Add(fHistClustDz);

  fNAccJets = new TH1F("fNAccJets","fNAccJets;N/ev",11,-0.5, 9.5);
  fOutput->Add(fNAccJets);
  
  fhZg = new TH1F("fhZg", "#it{Z}_{g}; #it{Z}_{g}; Entries", 200, 0., 0.5);
  fOutput->Add(fhZg);
  
  fhCorrPtZg = new TH2F("fhCorrPtZg", "#it{Z}_{g}; p_{T}^{corr} [GeV/c]; #it{Z}_{g}", 16, 0, 160, 20, 0., 0.5);
  fOutput->Add(fhCorrPtZg);

  fhCorrPtZg2 = new TH2F("fhCorrPtZg2", "#it{Z}_{g}; p_{T}^{corr} [GeV/c]; #it{Z}_{g}", 16, 0, 160, 20, 0., 0.5);
  fOutput->Add(fhCorrPtZg2);

  fhCorrPtZgD = new TH2F("fhCorrPtZgD", "#it{Z}_{g}; p_{T}^{corr} [GeV/c]; #it{Z}_{g}", 16, 0, 160, 20, 0., 0.5);
  fOutput->Add(fhCorrPtZgD);

  fhCorrPtRg = new TH2F("fhCorrPtRg", "R_{g}; p_{T}^{corr} [GeV/c]; R_{g}", 16, 0, 160, 40, 0., 0.5);
  fOutput->Add(fhCorrPtRg);

  fhCorrPtRgD = new TH2F("fhCorrPtRgD", "R_{g}; p_{T}^{corr} [GeV/c]; R_{g}", 16, 0, 160, 40, 0., 0.5);
  fOutput->Add(fhCorrPtRgD);

  fhCorrPtZgRg = new TH3F("fhCorrPtZgRg", "fhCorrPtZgRg", 8, 0, 160, 10, 0., 0.5, 10, 0, 0.5);
  fOutput->Add(fhCorrPtZgRg);

  fhCorrPtZgSDstep = new TH3F("fhCorrPtZgSDstep", "fhCorrPtZgSDstep", 16, 0, 160, 20, 0., 0.5, 20, 0, 20);
  fOutput->Add(fhCorrPtZgSDstep);

  fhCorrPtRgSDstep = new TH3F("fhCorrPtRgSDstep", "fhCorrPtRgSDstep", 16, 0, 160, 20, 0., 0.5, 20, 0, 20);
  fOutput->Add(fhCorrPtRgSDstep);

  fhCorrPtPtfrac = new TH2F("fhCorrPtPtfrac", "#deltap_{T}; p_{T}^{corr} [GeV/c]; #deltap_{T}", 16, 0, 160, 80, 0., 1.0);
  fOutput->Add(fhCorrPtPtfrac);

  fhCorrPtDropCount = new TH2F("fhCorrPtDropCount", "fhCorrPtDropCount; p_{T}^{corr} [GeV/c]; Counts", 16, 0, 160, 50, 0., 50);
  fOutput->Add(fhCorrPtDropCount);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSoftDrop::FillHistograms()
{
  // Fill histograms.

  if (fTracksCont) {
     fHistNTracks[fCentBin]->Fill(fTracksCont->GetNAcceptedParticles());
     fTracksCont->ResetCurrentID();
    AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
    while(track) {
      fHistTracksPt[fCentBin]->Fill(track->Pt()); 
      track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
      
    }
  }

  if (fCaloClustersCont) {
    AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster(); 
    while(cluster) {
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);
      fHistClustersPt[fCentBin]->Fill(nPart.Pt());
      Double_t dx = cluster->GetTrackDx();
      Double_t dz = cluster->GetTrackDz();
      fHistClustDx->Fill(dx);
      fHistClustDz->Fill(dz);
      cluster = fCaloClustersCont->GetNextAcceptCluster();
    }
  }

  if (fJetsCont) {
    Int_t count = 0;
    for (auto jet : fJetsCont->accepted() ) {
      count++;
      fHistJetsPtArea[fCentBin]->Fill(jet->Pt(), jet->Area());
      fHistJetsPhiEta[fCentBin]->Fill(jet->Eta(), jet->Phi());

      Float_t ptLeading = fJetsCont->GetLeadingHadronPt(jet);
      fHistJetsPtLeadHad[fCentBin]->Fill(jet->Pt(), ptLeading);

      if (fHistJetsCorrPtArea[fCentBin]) {
	Float_t corrPt = jet->Pt() - fJetsCont->GetRhoVal() * jet->Area();
	fHistJetsCorrPtArea[fCentBin]->Fill(corrPt, jet->Area());
      }

      Double_t jetpt_ungrmd = jet->Pt() / ( jet->GetShapeProperties()->GetSoftDropPtfrac() );

      std::vector<fastjet::PseudoJet> particles;
      UShort_t ntracks = jet->GetNumberOfTracks();
      for (int j = 0; j < ntracks; j++) {
        particles.push_back( fastjet::PseudoJet( jet->Track(j)->Px(), jet->Track(j)->Py(), jet->Track(j)->Pz(), jet->Track(j)->E() ) );
      }
      fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 0.4, fastjet::E_scheme);
      fastjet::ClusterSequence cs(particles, jet_def);
      std::vector<fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

      if (jets.size() > 0) {
        fSDM = 0;
        SoftDropDeepDeclustering( jets[0], jets[0].pt() );
        fhCorrPtZg2->Fill( jets[0].pt(), SoftDropDeclustering(jets[0], 0.5, 1.5) );
      }

      fhZg->Fill(jet->GetShapeProperties()->GetSoftDropZg());
      fhCorrPtZg->Fill(jetpt_ungrmd - fJetsCont->GetRhoVal() * jet->Area(), jet->GetShapeProperties()->GetSoftDropZg() );
      fhCorrPtRg->Fill(jetpt_ungrmd - fJetsCont->GetRhoVal() * jet->Area(), jet->GetShapeProperties()->GetSoftDropdR() );
      fhCorrPtZgRg->Fill(jetpt_ungrmd - fJetsCont->GetRhoVal() * jet->Area(), jet->GetShapeProperties()->GetSoftDropZg(), jet->GetShapeProperties()->GetSoftDropdR() );
      fhCorrPtPtfrac->Fill(jetpt_ungrmd - fJetsCont->GetRhoVal() * jet->Area(), jet->GetShapeProperties()->GetSoftDropPtfrac() );
      fhCorrPtDropCount->Fill(jetpt_ungrmd - fJetsCont->GetRhoVal() * jet->Area(), jet->GetShapeProperties()->GetSoftDropDropCount() );
    }
    fNAccJets->Fill(count);
    auto jet = (AliEmcalJet*)fJetsCont->GetLeadingJet();
    if(jet) fHistLeadingJetPt[fCentBin]->Fill(jet->Pt());
  }

  CheckClusTrackMatching();

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskSoftDrop::CheckClusTrackMatching()
{
  
  if(!fTracksCont || !fCaloClustersCont)
    return;

  Double_t deta = 999;
  Double_t dphi = 999;

  //Get closest cluster to track
  AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle()); 
  while(track) {
    //Get matched cluster
    Int_t emc1 = track->GetEMCALcluster();
    if(fCaloClustersCont && emc1>=0) {
      AliVCluster *clusMatch = fCaloClustersCont->GetCluster(emc1);
      if(clusMatch) {
	AliPicoTrack::GetEtaPhiDiff(track, clusMatch, dphi, deta);
	fHistPtDEtaDPhiTrackClus->Fill(track->Pt(),deta,dphi);
      }
    }
    track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
  }
  
  //Get closest track to cluster
  AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster(); 
  while(cluster) {
    TLorentzVector nPart;
    cluster->GetMomentum(nPart, fVertex);
    fHistClustersPt[fCentBin]->Fill(nPart.Pt());
    
    //Get matched track
    AliVTrack *mt = NULL;      
    AliAODCaloCluster *acl = dynamic_cast<AliAODCaloCluster*>(cluster);
    if(acl) {
      if(acl->GetNTracksMatched()>1)
	mt = static_cast<AliVTrack*>(acl->GetTrackMatched(0));
    }
    else {
      AliESDCaloCluster *ecl = dynamic_cast<AliESDCaloCluster*>(cluster);
      Int_t im = ecl->GetTrackMatchedIndex();
      if(fTracksCont && im>=0) {
	mt = static_cast<AliVTrack*>(fTracksCont->GetParticle(im));
      }
    }
    if(mt) {
      AliPicoTrack::GetEtaPhiDiff(mt, cluster, dphi, deta);
      fHistPtDEtaDPhiClusTrack->Fill(nPart.Pt(),deta,dphi);
      
      /* //debugging
	 if(mt->IsEMCAL()) {
	 Int_t emc1 = mt->GetEMCALcluster();
	 Printf("current id: %d  emc1: %d",fCaloClustersCont->GetCurrentID(),emc1);
	 AliVCluster *clm = fCaloClustersCont->GetCluster(emc1);
	 AliPicoTrack::GetEtaPhiDiff(mt, clm, dphi, deta);
	 Printf("deta: %f dphi: %f",deta,dphi);
	 }
      */
    }
    cluster = fCaloClustersCont->GetNextAcceptCluster();
  }
}

void AliAnalysisTaskSoftDrop::SoftDropDeepDeclustering(fastjet::PseudoJet jet, const Float_t inpt) {

  fastjet::PseudoJet jet1;
  fastjet::PseudoJet jet2;

  if ( jet.has_parents(jet1, jet2) ) {

    Float_t pt1 = jet1.pt();
    Float_t pt2 = jet2.pt();

    Float_t dr = TMath::Sqrt( jet1.plain_distance(jet2) );

    Float_t z;
    if (pt1 < pt2) z = pt1/(pt1+pt2);
    else z = pt2/(pt1+pt2);

    if (z > 0.1) {
      fSDM++;
      fhCorrPtZgD->Fill(inpt, z);
      fhCorrPtRgD->Fill(inpt, dr);
      fhCorrPtZgSDstep->Fill(inpt, z,  fSDM);
      fhCorrPtRgSDstep->Fill(inpt, dr, fSDM);
    }

    if (pt1 > pt2) SoftDropDeepDeclustering(jet1, inpt);
    else SoftDropDeepDeclustering(jet2, inpt);

  }

}

Float_t AliAnalysisTaskSoftDrop::SoftDropDeclustering(fastjet::PseudoJet jet, const Float_t zcut, const Float_t beta) {

  fastjet::PseudoJet jet1;
  fastjet::PseudoJet jet2;

  if ( !(jet.has_parents(jet1, jet2)) ) {
    return 0.0;
  }
  else {
    Float_t pt1 = jet1.pt();
    Float_t pt2 = jet2.pt();

    Float_t dr = TMath::Sqrt( jet1.plain_distance(jet2) );
    Float_t angular_term = TMath::Power(dr/0.4, beta);

    Float_t z;
    if (pt1 < pt2) z = pt1/(pt1+pt2);
    else z = pt2/(pt1+pt2);
    
    if ( z > (zcut*angular_term) ) return z;
    else {
      if (pt1 > pt2) return SoftDropDeclustering(jet1, zcut, beta);
      else return SoftDropDeclustering(jet2, zcut, beta);
    }
  }

}

//________________________________________________________________________
void AliAnalysisTaskSoftDrop::ExecOnce() {

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskSoftDrop::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//________________________________________________________________________
void AliAnalysisTaskSoftDrop::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
