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

Bool_t kPrintCaloCells    = kTRUE;
Bool_t kPrintCaloClusters = kTRUE;


void TestAOD() {

  TFile* f = new TFile("AliAOD.root");
  TTree* aodTree = (TTree*)f->Get("aodTree");
  
  AliAODEvent* aod = new AliAODEvent();
  aod->ReadFromTree(aodTree);

  Int_t nEvt = aodTree->GetEntries();

  for(Int_t iev = 0; iev < nEvt; iev++) {
    cout << "Event: " << iev+1 << "/" << nEvt << endl;
    aodTree->GetEvent(iev);

    //get reconstructed vertex position
    Double_t vertex_position[3] = { aod->GetPrimaryVertex()->GetX(),
				    aod->GetPrimaryVertex()->GetY(),
				    aod->GetPrimaryVertex()->GetZ()};

    //------------------------------------------------------
    // Clusters loop
    //------------------------------------------------------
    if(kPrintCaloClusters)
    {  
      TRefArray* caloClusters = new TRefArray();
      aod->GetEMCALClusters(caloClusters);

      Int_t nclus = caloClusters->GetEntries();
      for (Int_t icl = 0; icl < nclus; icl++) 
      {
        
        AliAODCaloCluster* clus = (AliAODCaloCluster*)caloClusters->At(icl);
        Float_t energy = clus->E();
        Float_t time   = clus->GetTOF()*1.e9;
        TLorentzVector p;
        clus->GetMomentum(p,vertex_position);
        Int_t nMatched = clus->GetNTracksMatched();
        
        cout << "Cluster: " << icl+1 << "/" << nclus << " - Energy: " << energy << "; Time "<<time
             <<"; Phi: " << p.Phi() << "; Eta: " << p.Eta() << "; #Matches: " << nMatched << endl;
        
      }
    }
    
    //------------------------------------------------------
    // Cells loop
    //------------------------------------------------------ 
    
    if(kPrintCaloCells)
    {  
      AliVCaloCells &cells= *(aod->GetEMCALCells());
      
      Int_t nTotalCells = cells.GetNumberOfCells() ;  
      //Int_t type        = cells.GetType();
      for (Int_t icell=  0; icell <  nTotalCells; icell++) {
        cout<<"Cell   : "<<icell<<"/"<<nTotalCells<<" - ID: "<<cells.GetCellNumber(icell)<<"; Amplitude: "<<cells.GetAmplitude(icell)<<"; Time: "<<cells.GetTime(icell)*1e9;
        cout << "; MC label "<<cells.GetMCLabel(icell)<<"; Embeded E fraction "<<cells.GetEFraction(icell);
        cout<<endl;	  
      }// cell loop
    }
    

  }



}
