#if !defined(__CINT__) || defined(__MAKECINT__)

//Root include files 
#include <Riostream.h>
#include <TFile.h>
#include <TChain.h>
#include <TParticle.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1F.h>
#include <TVector.h>
#include <TParticle.h>
#include <TRefArray.h>
#include <TArrayS.h>

//AliRoot include files 
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODCaloCluster.h"
#include "AliAODCaloCells.h"
#include "AliPID.h"
#include "AliLog.h"

#endif

void TestAOD() {

  TFile* f = new TFile("aod.root");
  TTree* aodTree = (TTree*)f->Get("aodTree");
  
  AliAODEvent* aod = new AliAODEvent();
  aod->ReadFromTree(aodTree);

  Int_t nEvt = aodTree->GetEntries();

  for(Int_t iev = 0; iev < nEvt; iev++) {
    cout << "Event: " << iev+1 << "/" << nEvt << endl;
    aodTree->GetEvent(iev);

    TRefArray* caloClusters = new TRefArray();
    aod->GetEMCALClusters(caloClusters);

    //get reconstructed vertex position
    Double_t vertex_position[3] = { aod->GetPrimaryVertex()->GetX(),
				    aod->GetPrimaryVertex()->GetY(),
				    aod->GetPrimaryVertex()->GetZ()};

    Int_t nclus = caloClusters->GetEntries();
    for (Int_t icl = 0; icl < nclus; icl++) {

      AliAODCaloCluster* clus = (AliAODCaloCluster*)caloClusters->At(icl);
      Float_t energy = clus->E();
      TLorentzVector p;
      clus->GetMomentum(p,vertex_position);
      Int_t nMatched = clus->GetNTracksMatched();

      cout << "Cluster: " << icl+1 << "/" << nclus << " Energy: " << energy << " Phi: " << p.Phi() << " Eta: " << p.Eta() << " #Matches: " << nMatched << endl;


      
    }


  }



}
