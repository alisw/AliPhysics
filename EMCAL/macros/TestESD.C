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
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliPID.h"
#include "AliLog.h"

#endif

void TestESD() {

  TFile* f = new TFile("AliESDs.root");
  TTree* esdTree = (TTree*)f->Get("esdTree");
  
  AliESDEvent* esd = new AliESDEvent();
  esd->ReadFromTree(esdTree);

  Int_t nEvt = esdTree->GetEntries();
  Float_t pos[3] = {0.,0.,0.};

  for(Int_t iev = 0; iev < nEvt; iev++) {
    cout << "Event: " << iev+1 << "/" << nEvt << endl;
    esdTree->GetEvent(iev);

    TRefArray* caloClusters = new TRefArray();
    esd->GetEMCALClusters(caloClusters);

    //get reconstructed vertex position
    Double_t vertex_position[3];
    esd->GetVertex()->GetXYZ(vertex_position);

    //loop over clusters
    Int_t nclus = caloClusters->GetEntries();
    for (Int_t icl = 0; icl < nclus; icl++) {

      AliESDCaloCluster* clus = (AliESDCaloCluster*)caloClusters->At(icl);
      Float_t energy = clus->E();
      clus->GetPosition(pos);
      TVector3 vpos(pos[0],pos[1],pos[2]);
      TLorentzVector p;
      clus->GetMomentum(p,vertex_position);
      Double_t cphi = vpos.Phi();
      Double_t ceta = vpos.Eta();

      Int_t nMatched = clus->GetNTracksMatched();
      Int_t trackIndex = clus->GetTrackMatched();
      Int_t nLabels = clus->GetNLabels();
      Int_t labelIndex = clus->GetLabel();

      Int_t nCells = clus->GetNCells();

      //For later: ADD CHECK THAT CLUSTER IS WITHIN SM FIDUCIAL VOLUME

      cout << "Cluster: " << icl+1 << "/" << nclus << " Energy: " << energy << " Phi: " 
	   << cphi << " Eta: " << ceta << " NCells: " << nCells << " #Matches: " << nMatched 
	   << " Index: " << trackIndex << " #Labels: " << nLabels << " Index: " 
	   << labelIndex << endl;

      if(trackIndex >= 0) {
	AliESDtrack* track = esd->GetTrack(trackIndex);
	Double_t tphi = track->GetOuterParam()->Phi();
	Double_t teta = track->GetOuterParam()->Eta();
	Double_t tmom = track->GetOuterParam()->P();
	cout << "\t Track Momentum: " << tmom << " phi: " << tphi << " eta: " << teta << endl;

	Double_t deta = teta - ceta;
	Double_t dphi = tphi - cphi;
	if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
	if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
	Double_t dR = sqrt(dphi*dphi + deta*deta);

	Double_t pOverE = tmom/energy;

	if(dR < 0.02 && pOverE < 1.8 && nCells > 1) {
	  cout << "\n\t Excellent MATCH! dR = " << dR << " p/E = " << pOverE << " nCells = " << nCells << endl;
	}

      }

    }


  }



}
