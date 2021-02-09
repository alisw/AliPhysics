/*************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//
// Soft Drop observables response taks based on the AliJetResponseMaker task 
//
// Author: Kirill Lapidus, Yale University, kirill.lapdidus@cern.ch

#include <THnSparse.h>

#include "AliAnalysisTaskSoftDropResponse.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSoftDrop.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"

ClassImp(AliAnalysisTaskSoftDropResponse)

//________________________________________________________________________
AliAnalysisTaskSoftDropResponse::AliAnalysisTaskSoftDropResponse() : 
  AliJetResponseMaker("AliAnalysisTaskSoftDropResponse"),
  fZ2gAxis(0)
{
  // Default constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskSoftDropResponse::AliAnalysisTaskSoftDropResponse(const char *name) : 
  AliJetResponseMaker(name),
  fZ2gAxis(0)
{
  // Standard constructor.

  SetMakeGeneralHistograms(kTRUE);
}


//________________________________________________________________________
void AliAnalysisTaskSoftDropResponse::AllocateTHnSparse()
{
  // Allocate THnSparse histograms.

  TString title[25]= {""};
  Int_t nbins[25]  = {0};
  Double_t min[25] = {0.};
  Double_t max[25] = {0.};
  Int_t dim = 0;

  title[dim] = "#phi";
  nbins[dim] = fNbins/4;
  min[dim] = 0;
  max[dim] = 2*TMath::Pi()*(1 + 1./(nbins[dim]-1));
  dim++;

  title[dim] = "#eta";
  nbins[dim] = fNbins/4;
  min[dim] = -1;
  max[dim] = 1;
  dim++;

  title[dim] = "p_{T}";
  nbins[dim] = fNbins;
  min[dim] = 0;
  max[dim] = 250;
  dim++;

  title[dim] = "A_{jet}";
  nbins[dim] = fNbins/4;
  min[dim] = 0;
  max[dim] = 1.5;
  dim++;

  Int_t dim1 = dim, dim2 = dim;

  if (fIsJet1Rho) {
    title[dim1] = "p_{T}^{corr}";
    nbins[dim1] = fNbins*2;
    min[dim1] = -250;
    max[dim1] = 250;
    dim1++;
  }

  if (fIsEmbedded) {
    title[dim1] = "p_{T}^{MC}";
    nbins[dim1] = fNbins;
    min[dim1] = 0;
    max[dim1] = 250;
    dim1++;
  }

  fHistJets1 = new THnSparseD("fHistJets1","fHistJets1",dim1,nbins,min,max);
  for (Int_t i = 0; i < dim1; i++) fHistJets1->GetAxis(i)->SetTitle(title[i]);
  fOutput->Add(fHistJets1);

  if (fIsJet2Rho) {
    title[dim2] = "p_{T}^{corr}";
    nbins[dim2] = fNbins*2;
    min[dim2] = -250;
    max[dim2] = 250;
    dim2++;
  }

  fHistJets2 = new THnSparseD("fHistJets2","fHistJets2",dim2,nbins,min,max);
  for (Int_t i = 0; i < dim2; i++) fHistJets2->GetAxis(i)->SetTitle(title[i]);
  fOutput->Add(fHistJets2);

  // Matching

  dim = 0;

  title[dim] = "p_{T,1}";
  nbins[dim] = fNbins;
  min[dim] = 0;
  max[dim] = 250;
  dim++;

  title[dim] = "p_{T,2}";
  nbins[dim] = fNbins;
  min[dim] = 0;
  max[dim] = 250;
  dim++;

  title[dim] = "A_{jet,1}";
  nbins[dim] = fNbins/4;
  min[dim] = 0;
  max[dim] = 1.5;
  dim++;

  title[dim] = "A_{jet,2}";
  nbins[dim] = fNbins/4;
  min[dim] = 0;
  max[dim] = 1.5;
  dim++;

  title[dim] = "distance";
  nbins[dim] = fNbins/2;
  min[dim] = 0;
  max[dim] = 1.2;
  dim++;

  if (fDeltaPtAxis) {
    title[dim] = "#deltaA_{jet}";
    nbins[dim] = fNbins/2;
    min[dim] = -1;
    max[dim] = 1;
    dim++;

    title[dim] = "#deltap_{T}";
    nbins[dim] = fNbins*2;
    min[dim] = -250;
    max[dim] = 250;
    dim++;
  }
  if (fIsJet1Rho) {
    title[dim] = "p_{T,1}^{corr}";
    nbins[dim] = fNbins*2;
    min[dim] = -250;
    max[dim] = 250;
    dim++;
  }
  if (fIsJet2Rho) {
    title[dim] = "p_{T,2}^{corr}";
    nbins[dim] = fNbins*2;
    min[dim] = -250;
    max[dim] = 250;
    dim++;
  }
  if (fDeltaPtAxis && (fIsJet1Rho || fIsJet2Rho)) {
    title[dim] = "#deltap_{T}^{corr}";
    nbins[dim] = fNbins*2;
    min[dim] = -250;
    max[dim] = 250;
    dim++;
  }
  if (fDeltaEtaDeltaPhiAxis) {
    title[dim] = "#delta#eta";
    nbins[dim] = fNbins/2;
    min[dim] = -1;
    max[dim] = 1;
    dim++;

    title[dim] = "#delta#phi";
    nbins[dim] = fNbins/2;
    min[dim] = -TMath::Pi()/2;
    max[dim] = TMath::Pi()*3/2;
    dim++;
  }
  if (fIsEmbedded) {
    title[dim] = "p_{T,1}^{MC}";
    nbins[dim] = fNbins;
    min[dim] = 0;
    max[dim] = 250;
    dim++;

    if (fDeltaPtAxis) {
      title[dim] = "#deltap_{T}^{MC}";
      nbins[dim] = fNbins*2;
      min[dim] = -250;
      max[dim] = 250;
      dim++;
    }
  }

  if (fZgAxis) {
    title[dim] = "Z_{g,1}";
    nbins[dim] = 20;
    min[dim] = 0.0;
    max[dim] = 1.0;
    dim++;
    title[dim] = "Z_{g,2}";
    nbins[dim] = 20;
    min[dim] = 0.0;
    max[dim] = 1.0;
    dim++;
  }

  if (fZ2gAxis) {
    title[dim] = "Z2_{g,1}";
    nbins[dim] = 20;
    min[dim] = 0.0;
    max[dim] = 1.0;
    dim++;
    title[dim] = "Z2_{g,2}";
    nbins[dim] = 20;
    min[dim] = 0.0;
    max[dim] = 1.0;
    dim++;
  }

  if (fdRAxis) {
    title[dim] = "R_{g,1}";
    nbins[dim] = 20;
    min[dim] = 0.0;
    max[dim] = 0.5;
    dim++;
    title[dim] = "R_{g,2}";
    nbins[dim] = 20;
    min[dim] = 0.0;
    max[dim] = 0.5;
    dim++;
  }

  if (fPtgAxis) {
    title[dim] = "p_{T,g,1}";
    nbins[dim] = 16;
    min[dim] = 0.0;
    max[dim] = 160.0;
    dim++;
    title[dim] = "p_{T,g,2}";
    nbins[dim] = 16;
    min[dim] = 0.0;
    max[dim] = 160.0;
    dim++;
  }

  if (fDBCAxis) {
    title[dim] = "DBC_{1}";
    nbins[dim] = 20;
    min[dim] = 0.0;
    max[dim] = 20.0;
    dim++;
    title[dim] = "DBC_{2}";
    nbins[dim] = 20;
    min[dim] = 0.0;
    max[dim] = 20.0;
    dim++;
  }
  
  title[dim] = "nsdsteps_{jet,1}";
  nbins[dim] = 10;
  min[dim]   = 1;
  max[dim]   = 10;
  dim++;
  
  title[dim] = "nsdsteps_{jet,2}";
  nbins[dim] = 10;
  min[dim]   = 1;
  max[dim]   = 10;
  dim++;
  
  title[dim] = "nhsplits_{jet,1}";
  nbins[dim] = 10;
  min[dim]   = 1;
  max[dim]   = 10;
  dim++;
  
  title[dim] = "nhsplits_{jet,2}";
  nbins[dim] = 10;
  min[dim]   = 1;
  max[dim]   = 10;
  dim++;
  
  title[dim] = "Z^{2}_{g,1}";
  nbins[dim] = 20;
  min[dim] = 0.0;
  max[dim] = 0.5;
  dim++;
  
  title[dim] = "Z^{2}_{g,2}";
  nbins[dim] = 20;
  min[dim] = 0.0;
  max[dim] = 0.5;
  dim++;
  
  title[dim] = "R^{2}_{g,1}";
  nbins[dim] = 20;
  min[dim] = 0.0;
  max[dim] = 0.5;
  dim++;
  
  title[dim] = "R^{2}_{g,2}";
  nbins[dim] = 20;
  min[dim] = 0.0;
  max[dim] = 0.5;
  dim++;

  fHistMatching = new THnSparseD("fHistMatching","fHistMatching",dim,nbins,min,max);

  for (Int_t i = 0; i < dim; i++) fHistMatching->GetAxis(i)->SetTitle(title[i]);

  fOutput->Add(fHistMatching);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSoftDropResponse::FillHistograms()
{
  // Fill histograms.

  AliJetContainer *jets1 = static_cast<AliJetContainer*>(fJetCollArray.At(0));
  AliJetContainer *jets2 = static_cast<AliJetContainer*>(fJetCollArray.At(1));

  if (!jets1 || !jets1->GetArray() || !jets2 || !jets2->GetArray()) return kFALSE;

  AliEmcalJet* jet1 = 0;  
  AliEmcalJet* jet2 = 0;

  jets2->ResetCurrentID();
  while ((jet2 = jets2->GetNextJet())) {

    AliDebug(2,Form("Processing jet (2) %d", jets2->GetCurrentID()));

    if (jet2->Pt() < jets2->GetJetPtCut()) continue;

    UInt_t rejectionReason = 0;
    if (jets2->AcceptJet(jet2, rejectionReason))
      FillJetHisto(jet2, 2);
    else
      fHistRejectionReason2->Fill(jets2->GetRejectionReasonBitPosition(rejectionReason), jet2->Pt());

    jet1 = jet2->MatchedJet();

    if (!jet1) continue;
    rejectionReason = 0;
    if (!jets1->AcceptJet(jet1, rejectionReason)) continue;
    if (jet1->MCPt() < fMinJetMCPt) continue;

    Double_t d=-1, ce1=-1, ce2=-1;
    if (jet2->GetMatchingType() == kGeometrical) {
      if (jets2->GetIsParticleLevel() && !jets1->GetIsParticleLevel()) // the other way around is not supported
        GetMCLabelMatchingLevel(jet1, jet2, ce1, ce2);
      else if (jets1->GetIsParticleLevel() == jets2->GetIsParticleLevel())
        GetSameCollectionsMatchingLevel(jet1, jet2, ce1, ce2);

      d = jet2->ClosestJetDistance();
    }
    else if (jet2->GetMatchingType() == kMCLabel || jet2->GetMatchingType() == kSameCollections) {
      GetGeometricalMatchingLevel(jet1, jet2, d);

      ce1 = jet1->ClosestJetDistance();
      ce2 = jet2->ClosestJetDistance();
    }

    FillMatchingHistos(jet1, jet2, d, ce1, ce2);
  }

  jets1->ResetCurrentID();
  while ((jet1 = jets1->GetNextJet())) {
    UInt_t rejectionReason = 0;
    if (!jets1->AcceptJet(jet1, rejectionReason)) {
      fHistRejectionReason1->Fill(jets1->GetRejectionReasonBitPosition(rejectionReason), jet1->Pt());
      continue;
    }
    if (jet1->MCPt() < fMinJetMCPt) continue;
    AliDebug(2,Form("Processing jet (1) %d", jets1->GetCurrentID()));

    FillJetHisto(jet1, 1);
  }
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskSoftDropResponse::FillMatchingHistos(AliEmcalJet* jet1, AliEmcalJet* jet2, Double_t d, Double_t CE1, Double_t CE2)
{
  AliJetContainer* jets1 = GetJetContainer(0);
  AliJetContainer* jets2 = GetJetContainer(1);

  Double_t corrpt1 = jet1->Pt() - jets1->GetRhoVal() * jet1->Area();
  Double_t corrpt2 = jet2->Pt() - jets2->GetRhoVal() * jet2->Area();

  //Float_t z2g1 = -1.0;
  //Float_t z2g2 = -1.0;
  //CalculateZg(jet1, 0.5, 1.5, z2g1);
  //CalculateZg(jet2, 0.5, 1.5, z2g2);
 
  FillZgRgVectors(jet1);
  
  Int_t nsdsteps_jet1 = 0;
  Int_t nhsplits_jet1 = fZg_values.size();
  Float_t zg2_1 = 0.0;
  Float_t rg2_1 = 0.0;
  if (nhsplits_jet1 > 0) nsdsteps_jet1 = fSDsteps_values[0];
  if (nhsplits_jet1 > 1) {
    zg2_1 = fZg_values[1];
    rg2_1 = fRg_values[1];
  }
  
  FillZgRgVectors(jet2);
  
  Int_t nsdsteps_jet2 = 0;
  Int_t nhsplits_jet2 = fZg_values.size();
  Float_t zg2_2 = 0.0;
  Float_t rg2_2 = 0.0;
  if (nhsplits_jet2 > 0) nsdsteps_jet2 = fSDsteps_values[0];
  if (nhsplits_jet2 > 1) {
    zg2_2 = fZg_values[1];
    rg2_2 = fRg_values[1];
  }
  
  if (fHistoType==1) {
    Double_t contents[25]={0};

    for (Int_t i = 0; i < fHistMatching->GetNdimensions(); i++) {
      TString title(fHistMatching->GetAxis(i)->GetTitle());
      if (title=="p_{T,1}")
        contents[i] = jet1->Pt();
      else if (title=="p_{T,2}")
        contents[i] = jet2->Pt();
      else if (title=="A_{jet,1}")
        contents[i] = jet1->Area();
      else if (title=="A_{jet,2}")
        contents[i] = jet2->Area();
      else if (title=="distance")
        contents[i] = d;
      else if (title=="#deltaA_{jet}")
        contents[i] = jet1->Area()-jet2->Area();
      else if (title=="#deltap_{T}")
        contents[i] = jet1->Pt()-jet2->Pt();
      else if (title=="#delta#eta")
        contents[i] = jet1->Eta()-jet2->Eta();
      else if (title=="#delta#phi")
        contents[i] = jet1->Phi()-jet2->Phi();
      else if (title=="p_{T,1}^{corr}")
        contents[i] = corrpt1;
      else if (title=="p_{T,2}^{corr}")
        contents[i] = corrpt2;
      else if (title=="#deltap_{T}^{corr}")
        contents[i] = corrpt1-corrpt2;
      else if (title=="p_{T,1}^{MC}")
        contents[i] = jet1->MCPt();
      else if (title=="#deltap_{T}^{MC}")
        contents[i] = jet1->MCPt()-jet2->Pt();
      else if (title=="Z_{g,1}")
        contents[i] = jet1->GetShapeProperties()->GetSoftDropZg();
      else if (title=="Z_{g,2}")
        contents[i] = jet2->GetShapeProperties()->GetSoftDropZg();
      else if (title=="Z2_{g,1}")
        contents[i] = 0.0;
      else if (title=="Z2_{g,2}")
        contents[i] = 0.0;
      else if (title=="R_{g,1}")
        contents[i] = jet1->GetShapeProperties()->GetSoftDropdR();
      else if (title=="R_{g,2}")
        contents[i] = jet2->GetShapeProperties()->GetSoftDropdR();
      else if (title=="p_{T,g,1}")
        contents[i] = ( jet1->GetShapeProperties()->GetSoftDropPtfrac() )*( jet1->Pt() );
      else if (title=="p_{T,g,2}")
        contents[i] = ( jet2->GetShapeProperties()->GetSoftDropPtfrac() )*( jet2->Pt() );
      else if (title=="DBC_{1}")
        contents[i] = ( jet1->GetShapeProperties()->GetSoftDropDropCount() );
      else if (title=="DBC_{2}")
        contents[i] = ( jet2->GetShapeProperties()->GetSoftDropDropCount() );
      else if (title=="nsdsteps_{jet,1}")
        contents[i] = nsdsteps_jet1;
      else if (title=="nsdsteps_{jet,2}")
        contents[i] = nsdsteps_jet2;
      else if (title=="nhsplits_{jet,1}")
        contents[i] = nhsplits_jet1;
      else if (title=="nhsplits_{jet,2}")
        contents[i] = nhsplits_jet2;
      else if (title=="Z^{2}_{g,1}")    
        contents[i] = zg2_1;
      else if (title=="Z^{2}_{g,2}")
        contents[i] = zg2_2;
      else if (title=="R^{2}_{g,1}")
        contents[i] = rg2_1;
      else if (title=="R^{2}_{g,2}")
        contents[i] = rg2_2;
      else 
        AliWarning(Form("Unable to fill dimension %s!",title.Data()));
    }

    fHistMatching->Fill(contents);
  }
  
}

/**
 * Soft Drop response maker AddTask.
 */
AliAnalysisTaskSoftDropResponse* AliAnalysisTaskSoftDropResponse::AddTaskSoftDropResponse(
  const char *ntracks1,
  const char *nclusters1,
  const char *njets1,
  const char *nrho1,
  const Double_t    jetradius1,
  const char *ntracks2,
  const char *nclusters2,
  const char *njets2,
  const char *nrho2,
  const Double_t    jetradius2,
  const Double_t    jetptcut,
  const Double_t    jetareacut,
  const Double_t    jetBias,
  const Int_t       biasType,
  const AliAnalysisTaskSoftDropResponse::MatchingType matching,
  const Double_t    maxDistance1,
  const Double_t    maxDistance2,
  const char *cutType,
  const Int_t       ptHardBin,
  const Double_t    minCent,
  const Double_t    maxCent,
  const char *taskname,
  const Bool_t      biggerMatrix,
  AliAnalysisTaskSoftDropResponse* address,
  const Double_t    maxTrackPt)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskJetResponseMaker", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskJetResponseMaker", "This task requires an input event handler");
    return NULL;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("%s_%s_%s_Bias%d_BiasType%d_%s",taskname,njets1,njets2,(Int_t)floor(jetBias),biasType,cutType));

  if (minCent != -999 && maxCent != -999)
    name += Form("_Cent%d_%d", (Int_t)floor(minCent), (Int_t)floor(maxCent));

  if (ptHardBin != -999)
    name += Form("_PtHard%d", ptHardBin);

  AliAnalysisTaskSoftDropResponse* jetTask = address;
  if (jetTask) {
    new (jetTask) AliAnalysisTaskSoftDropResponse(name);
  }
  else {
    jetTask = new AliAnalysisTaskSoftDropResponse(name);
  }

  AliParticleContainer *trackCont1 = jetTask->AddParticleContainer(ntracks1);
  AliClusterContainer *clusCont1 = jetTask->AddClusterContainer(nclusters1);
  AliJetContainer *jetCont1 = jetTask->AddJetContainer(njets1, cutType, jetradius1);
  jetCont1->SetRhoName(nrho1);
  jetCont1->SetLeadingHadronType(biasType);
  jetCont1->SetPtBiasJetTrack(jetBias);
  jetCont1->SetPtBiasJetClus(jetBias);
  jetCont1->SetJetPtCut(jetptcut);
  jetCont1->SetPercAreaCut(jetareacut);
  jetCont1->SetIsParticleLevel(kFALSE);
  jetCont1->ConnectParticleContainer(trackCont1);
  jetCont1->ConnectClusterContainer(clusCont1);
  jetCont1->SetMaxTrackPt(maxTrackPt);

  AliParticleContainer *trackCont2 = jetTask->AddParticleContainer(ntracks2);
  trackCont2->SetParticlePtCut(0);
  AliClusterContainer *clusCont2 = jetTask->AddClusterContainer(nclusters2);
  AliJetContainer *jetCont2 = jetTask->AddJetContainer(njets2, cutType, jetradius2);
  jetCont2->SetRhoName(nrho2);
  jetCont2->SetLeadingHadronType(biasType);
  jetCont2->SetPtBiasJetTrack(jetBias);
  jetCont2->SetPtBiasJetClus(jetBias);
  jetCont2->SetJetPtCut(jetptcut);
  jetCont2->SetPercAreaCut(jetareacut);
  jetCont2->SetIsParticleLevel(kTRUE);
  jetCont2->ConnectParticleContainer(trackCont2);
  jetCont2->ConnectClusterContainer(clusCont2);
  jetCont2->SetMaxTrackPt(1000); // disable default 100 GeV/c track cut for particle level jets

  jetTask->SetMatching(matching, maxDistance1, maxDistance2);
  jetTask->SetVzRange(-10,10);
  jetTask->SetIsPythia(kTRUE);
  jetTask->SetPtHardBin(ptHardBin);
  jetTask->SetCentRange(minCent,maxCent);

  if (biggerMatrix) jetTask->SetHistoBins(1000,0,500);

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

void AliAnalysisTaskSoftDropResponse::Decluster(const fastjet::PseudoJet& jet) {

  fastjet::PseudoJet jet1;
  fastjet::PseudoJet jet2;

  if ( jet.has_parents(jet1, jet2) ) {

    ++fNsdsteps;    

    Float_t pt1 = jet1.pt();
    Float_t pt2 = jet2.pt();

    Float_t dr = TMath::Sqrt( jet1.plain_distance(jet2) );

    Float_t z;
    if (pt1 < pt2) z = pt1/(pt1+pt2);
    else z = pt2/(pt1+pt2);

    if (z > 0.1) {
      fZg_values.push_back(z);
      fRg_values.push_back(dr);
      fSDsteps_values.push_back(fNsdsteps);
    }

    if (pt1 > pt2) Decluster(jet1);
    else Decluster(jet2);

  }

}

//________________________________________________________________________
fastjet::ClusterSequence* AliAnalysisTaskSoftDropResponse::Recluster(const AliEmcalJet* jet) {

  std::vector<fastjet::PseudoJet> particles;
  UShort_t ntracks = jet->GetNumberOfTracks();
  for (int j = 0; j < ntracks; j++) {
    particles.push_back( fastjet::PseudoJet( jet->Track(j)->Px(), jet->Track(j)->Py(), jet->Track(j)->Pz(), jet->Track(j)->E() ) );
  }
  fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 0.4, fastjet::E_scheme);
  return ( new fastjet::ClusterSequence(particles, jet_def) );
}

void AliAnalysisTaskSoftDropResponse::FillZgRgVectors(const AliEmcalJet* jet) {
  fNsdsteps = 0;
  fZg_values.clear();
  fRg_values.clear();
  fSDsteps_values.clear();
  fastjet::ClusterSequence* cs = Recluster(jet);
  if (cs) { 
    std::vector<fastjet::PseudoJet> jetrecl = sorted_by_pt( cs->inclusive_jets() );
    if (jetrecl.size() > 0) Decluster( jetrecl[0] );
  }
  delete cs;
}

//________________________________________________________________________
void AliAnalysisTaskSoftDropResponse::UserCreateOutputObjects()
{
  // Create user objects.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  AliJetContainer *jets1 = static_cast<AliJetContainer*>(fJetCollArray.At(0));
  AliJetContainer *jets2 = static_cast<AliJetContainer*>(fJetCollArray.At(1));

  if (!jets1 || !jets2) return;

  if (jets1->GetRhoName().IsNull()) fIsJet1Rho = kFALSE;
  else fIsJet1Rho = kTRUE;
  if (jets2->GetRhoName().IsNull()) fIsJet2Rho = kFALSE;
  else fIsJet2Rho = kTRUE;

  fHistRejectionReason1 = new TH2F("fHistRejectionReason1", "fHistRejectionReason1", 32, 0, 32, 100, 0, 250);
  fHistRejectionReason1->GetXaxis()->SetTitle("Rejection reason");
  fHistRejectionReason1->GetYaxis()->SetTitle("p_{T,jet} (GeV/c)");
  fHistRejectionReason1->GetZaxis()->SetTitle("counts");
  SetRejectionReasonLabels(fHistRejectionReason1->GetXaxis());
  fOutput->Add(fHistRejectionReason1);

  fHistRejectionReason2 = new TH2F("fHistRejectionReason2", "fHistRejectionReason2", 32, 0, 32, 100, 0, 250);
  fHistRejectionReason2->GetXaxis()->SetTitle("Rejection reason");
  fHistRejectionReason2->GetYaxis()->SetTitle("p_{T,jet} (GeV/c)");
  fHistRejectionReason2->GetZaxis()->SetTitle("counts");
  SetRejectionReasonLabels(fHistRejectionReason2->GetXaxis());
  fOutput->Add(fHistRejectionReason2);

  if (fHistoType==0)
    AllocateTH2();
  else 
    AllocateTHnSparse();

  // Initialize
  const AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
  if (embeddingHelper) {
    bool res = fEmbeddingQA.Initialize();
    if (res) {
      fEmbeddingQA.AddQAPlotsToList(fOutput);
    }
  }

  PostData(1, fOutput); // Post data for ALL output slots > 0 here, to get at least an empty histogram
}


//________________________________________________________________________
AliAnalysisTaskSoftDropResponse::~AliAnalysisTaskSoftDropResponse()
{
  // Destructor
}