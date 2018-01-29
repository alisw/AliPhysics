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

  title[dim] = "CE1";
  nbins[dim] = fNbins/2;
  min[dim] = 0;
  max[dim] = 1.2;
  dim++;

  title[dim] = "CE2";
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
    title[dim] = "dR_{1}";
    nbins[dim] = 40;
    min[dim] = 0.0;
    max[dim] = 0.5;
    dim++;
    title[dim] = "dR_{2}";
    nbins[dim] = 40;
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

  fHistMatching = new THnSparseD("fHistMatching","fHistMatching",dim,nbins,min,max);

  for (Int_t i = 0; i < dim; i++)
    fHistMatching->GetAxis(i)->SetTitle(title[i]);

  fOutput->Add(fHistMatching);
}

//________________________________________________________________________
void AliAnalysisTaskSoftDropResponse::CalculateZg(AliEmcalJet* jet, const Float_t zcut, const Float_t beta, Float_t zg)
{

  zg = -1.0;

  std::vector<fastjet::PseudoJet> particles;
  UShort_t ntracks = jet->GetNumberOfTracks();
  for (int j = 0; j < ntracks; j++) {
    particles.push_back( fastjet::PseudoJet( jet->Track(j)->Px(), jet->Track(j)->Py(), jet->Track(j)->Pz(), jet->Track(j)->E() ) );
  }
  fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 0.4, fastjet::E_scheme);
  fastjet::ClusterSequence cs(particles, jet_def);
  std::vector<fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

  if (jets.size() > 0) {
    zg = AliAnalysisTaskSoftDrop::SoftDropDeclustering(jets[0], zcut, beta);
  }

}

//________________________________________________________________________
void AliAnalysisTaskSoftDropResponse::FillMatchingHistos(AliEmcalJet* jet1, AliEmcalJet* jet2, Double_t d, Double_t CE1, Double_t CE2)
{
  AliJetContainer* jets1 = GetJetContainer(0);
  AliJetContainer* jets2 = GetJetContainer(1);

  Double_t corrpt1 = jet1->Pt() - jets1->GetRhoVal() * jet1->Area();
  Double_t corrpt2 = jet2->Pt() - jets2->GetRhoVal() * jet2->Area();

  Float_t z2g1 = -1.0;
  Float_t z2g2 = -1.0;
  CalculateZg(jet1, 0.5, 1.5, z2g1);
  CalculateZg(jet2, 0.5, 1.5, z2g2);

  if (fHistoType==1) {
    Double_t contents[20]={0};

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
      else if (title=="CE1")
        contents[i] = CE1;
      else if (title=="CE2")
        contents[i] = CE2;
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
        contents[i] = z2g1;
      else if (title=="Z2_{g,2}")
        contents[i] = z2g2;
      else if (title=="dR_{1}")
        contents[i] = jet1->GetShapeProperties()->GetSoftDropdR();
      else if (title=="dR_{2}")
        contents[i] = jet2->GetShapeProperties()->GetSoftDropdR();
      else if (title=="p_{T,g,1}")
        contents[i] = ( jet1->GetShapeProperties()->GetSoftDropPtfrac() )*( jet1->Pt() );
      else if (title=="p_{T,g,2}")
        contents[i] = ( jet2->GetShapeProperties()->GetSoftDropPtfrac() )*( jet2->Pt() );
      else if (title=="DBC_{1}")
        contents[i] = ( jet1->GetShapeProperties()->GetSoftDropDropCount() );
      else if (title=="DBC_{2}")
        contents[i] = ( jet2->GetShapeProperties()->GetSoftDropDropCount() );
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


//________________________________________________________________________
AliAnalysisTaskSoftDropResponse::~AliAnalysisTaskSoftDropResponse()
{
  // Destructor
}