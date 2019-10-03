/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
//
// Jet fragmentation transverse momentum (j_T) analysis task
//
// Author: Beomkyu Kim, Beomsu Chang, Dongjo Kim, Tomas Snellman

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TVector.h>
#include <TGraphErrors.h>
#include <TGrid.h>
#include <TSystem.h>
#include <TFile.h>

#include "AliCentrality.h"



#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPicoTrack.h"
#include "AliJBaseTrack.h"
#include "AliJMCTrack.h"
#include "AliJJet.h"
#include "AliJJetTask.h"
#include "AliEmcalTrackSelectionAOD.h"


ClassImp(AliJJetTask);

//________________________________________________________________________
AliJJetTask::AliJJetTask() : 
  AliAnalysisTaskEmcalJet("AliJJetTask", kTRUE),
  fJetsCont(),
  fTrackOrMCParticle(),
  fConeSizes(),
  fJTracks("AliJBaseTrack",1000),
  fJClusters("AliJBaseTrack",1000),
  fJMCTracks("AliJMCTrack",1000),
  fJJets(),
  fTaskEntry(-1),
  fJetFinderString(),
  fTrackArrayName("nonejk"),
  fNJetFinder(0),
  debug(0),
  fIsMC(0),
  fnR(0),
  fDoFullJets(0),
  fACside(0),
  fnkt(0)

{
  // Default constructor.


  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliJJetTask::AliJJetTask(const char *name, const int nJetFinder) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fJetsCont(nJetFinder),
  fTrackOrMCParticle(nJetFinder,kJUndefined),
  fConeSizes(nJetFinder,0.0),
  fJTracks("AliJBaseTrack",1000),
  fJClusters("AliJBaseTrack",1000),
  fJMCTracks("AliJMCTrack",1000),
  fJJets(),
  fTaskEntry(-1),
  fJetFinderString(nJetFinder),
  fTrackArrayName("nonejk"),
  fNJetFinder(nJetFinder),
  debug(0),
  fIsMC(0),
  fnR(0),
  fDoFullJets(0),
  fACside(0),
  fnkt(0)
{
  SetMakeGeneralHistograms(kTRUE);
}

AliJJetTask::AliJJetTask(const AliJJetTask& ap) :
  AliAnalysisTaskEmcalJet(ap.fName, kTRUE),
  fJetsCont(ap.fJetsCont),
  //fCaloClustersCont(ap.fCaloClustersCont),
  fTrackOrMCParticle(ap.fTrackOrMCParticle),
  fConeSizes(ap.fConeSizes),
  fJTracks(ap.fJTracks),
  fJMCTracks(ap.fJMCTracks),
  fJClusters(ap.fJClusters),
  fJJets(ap.fJJets),
  fTaskEntry(ap.fTaskEntry),
  fJetFinderString(ap.fJetFinderString),
  fTrackArrayName(ap.fTrackArrayName),
  fNJetFinder(ap.fNJetFinder),
  debug(ap.debug),
  fIsMC(ap.fIsMC),
  fnR(ap.fnR),
  fDoFullJets(ap.fnR),
  fACside(ap.fACside),
  fnkt(ap.fnkt)
{

}


AliJJetTask& AliJJetTask::operator = (const AliJJetTask& ap)
{

  this->~AliJJetTask();
  new(this) AliJJetTask(ap);
  return *this;
}

//________________________________________________________________________
AliJJetTask::~AliJJetTask()
{
  // Destructor.

  //delete[] fJJets;
}




//________________________________________________________________________
void AliJJetTask::UserCreateOutputObjects()
{
  // Create user output.

  //== Init Variables
  fJJets.clear();
  fJJets.resize( fNJetFinder, TClonesArray("AliJJet",1000) );
  fJetsCont.resize( fNJetFinder ); 
  //fCaloClustersCont.resize( fNJetFinder ); 

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();


  fJetFinderString.clear();

  cout << "fNJetFinder: " << fNJetFinder << endl;
  for (int i=0; i<fNJetFinder; i++){
    fJetsCont[i]           = GetJetContainer(i);
    fJetFinderString.push_back(fJetsCont[i]->GetArrayName());
    cout << i <<"\t" << fJetFinderString[i] << endl;
  }

}

//________________________________________________________________________
Bool_t AliJJetTask::FillHistograms()
{
  // FIXME : We assume that we have only one of each. be  carefull. This must be fixed later, See for example AliEmcalJetTask.cxx FindJets()
  //== AliJTrack
  //FIXME Search by container name
  AliParticleContainer *tracksCont = GetParticleContainer("tracks");
  if(debug > 0){
    cout << "AliParticleContainer name: " << tracksCont->GetName() << endl;
    cout << "NParticles(): " << tracksCont->GetNParticles() << endl;
  }
  if( tracksCont ){
    for (int itrack = 0; itrack<tracksCont->GetNParticles(); itrack++){
      AliAODTrack *track = static_cast<AliAODTrack*>(tracksCont->GetParticle(itrack));
      if(!track) continue;
      new (fJTracks[itrack]) AliJBaseTrack(track->Px(),track->Py(), track->Pz(), track->E(), itrack,0,track->Charge());
      AliJBaseTrack * particle = static_cast<AliJBaseTrack*>(fJTracks[itrack]);
      particle->SetLabel(track->GetLabel());
      if(track->IsHybridGlobalConstrainedGlobal()){
        particle->SetFlag(1,kTRUE);
      }
    }     
  }

  //FIXME Search by container name
  if(fIsMC){
    AliMCParticleContainer * mcTracksCont = GetMCParticleContainer("mcparticles");
    if( mcTracksCont ){
      if(debug > 0){
        cout << "MCParticleContainer name: " << mcTracksCont->GetName() << endl;
        cout << "mcTracksCont->GetNParticles(): " << mcTracksCont->GetNParticles() << endl;
      }
      int tracks = 0;
      for (int itrack = 0; itrack<mcTracksCont->GetNParticles(); itrack++){
        AliAODMCParticle *track = static_cast<AliAODMCParticle*>(mcTracksCont->GetParticle(itrack));
        new (fJMCTracks[itrack]) AliJMCTrack(track->Px(),track->Py(), track->Pz(), track->E(), itrack,0,track->Charge());
        tracks++;
        AliJMCTrack * particle = static_cast<AliJMCTrack*>(fJMCTracks[itrack]);
        if(track->IsPhysicalPrimary()){
          particle->SetPrimary();
        }
        particle->SetPdgCode(track->GetPdgCode());
        particle->SetLabel(track->GetLabel());
        particle->SetMother(track->GetMother(),track->GetMother());
      }
      if(debug > 0){
        cout << "Number of accepted tracks: " << tracks << endl;
      }
    }
  }

  AliClusterContainer *clusterCont = GetClusterContainer(0);
  if( clusterCont ){
    for (int itrack = 0; itrack<clusterCont->GetNClusters(); itrack++){
      AliAODCaloCluster *track = static_cast<AliAODCaloCluster*>(clusterCont->GetCluster(itrack));
      TLorentzVector momentum;
      Double_t vertex[3] = {0,0,0}; //FIXME Get the true collision vertex
      track->GetMomentum(momentum,vertex); //FIXME Find the right way of getting the momentum
      new (fJClusters[itrack]) AliJBaseTrack(momentum.Px(),momentum.Py(), momentum.Pz(), track->E(), itrack,0,0); //No charge in AliAODCaloCluster //FIXME Particle type?
      // FIXME: AliJPhoton instead of JBaseTrack?
    }
  }

  for (int i=0; i<fNJetFinder; i++){
    if( ! fJetsCont[i] ){
      // FIXME: Give Warning!
      continue;
    }

    fJetsCont[i]->ResetCurrentID(); //Needed to reset internal iterator
    AliEmcalJet *jet = fJetsCont[i]->GetNextAcceptJet();
    int iJet =0;

    //fills fJJets[icontainer][ijet] and histograms
    while(jet){
      TClonesArray & jets = fJJets[i]; // just alias for AliJJet array
      if((fACside == 1 && jet->Eta() < 0) || (fACside == 2 && jet->Eta() > 0)){ 
        //IF ACside == 1, remove jets with negative eta (C side removed)
        //IF ACside == 2, remove jets with positive eta (A side removed)
        jet = fJetsCont[i]->GetNextAcceptJet();
        continue;
      }
      new (jets[iJet]) AliJJet(jet->Px(),jet->Py(), jet->Pz(), jet->E(), jet->GetLabel(),0,0);
      AliJJet * j = (AliJJet*) fJJets[i][iJet];
      j->SetArea( jet->Area() );

      //== TRACK or Particle
      int nTrack = jet->GetNumberOfTracks();
      for (int it=0; it<nTrack; it++){
        int iTrack = jet->TrackAt(it)%10000; //FIXME Should be 100 000?
        if( fTrackOrMCParticle[i] == kJRecoTrack ){
          j->AddConstituent(fJTracks[iTrack]); // Save as pointers
        } else {
          j->AddConstituent(fJMCTracks[iTrack]); // Save as pointers
          if(debug > 0){
            cout << "iTrack: " << iTrack << endl;
          }
        }
      }
      //Particle type is kJCluster, set particle type§
      int nCluster = jet->GetNumberOfClusters();
      for (int it=0; it<nCluster; it++){
        int iTrack = jet->ClusterAt(it);
        j->AddConstituent(fJClusters[iTrack]); // Save as pointers
      }

      Int_t nConstituents = jet->GetNumberOfConstituents();

      if (debug>0) { 
        cout
          << "    iContainer    : " << i
          << "    iJet          : " << iJet
          << "    iTrack        : " << nTrack
          << "    iCluster      : " << nCluster
          << "    nConstituents : " << nConstituents
          << endl;
      }

      //Goes to the next jet
      jet = fJetsCont[i]->GetNextAcceptJet();
      if (debug>0) {
        cout
          << "  fJJets N lists : "     << fJJets[i].GetEntries()
          << "  fJJets constituents : "<< ((AliJJet*)fJJets[i][iJet])->GetConstituents()->GetEntries()
          <<endl;
      }
      iJet++;
    }
  }
  return kTRUE;

}

//________________________________________________________________________
void AliJJetTask::ExecOnce() {

  if(debug > 0){
    cout << "AliJJetTask::ExecOnce(): " << endl;
  }
  AliAnalysisTaskEmcalJet::ExecOnce();

  for (int i=0; i<fNJetFinder; i++){
    if (fJetsCont[i] && fJetsCont[i]->GetArray() == 0) fJetsCont[i] = 0;
  }
}

//________________________________________________________________________
Bool_t AliJJetTask::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().
  fTaskEntry = fEntry; // FIXME: Comments me

  //== Clear Before Fill Histogram
  for (int i=0; i<fNJetFinder; i++) { fJJets[i]. Clear();
    fJTracks.Clear(); 
    fJMCTracks.Clear();

  }
  if (debug >0 && Entry()%1000 ==0 ) cout<<Entry()<<endl;
  return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//________________________________________________________________________
void AliJJetTask::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

