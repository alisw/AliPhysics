#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "TFile.h"
#include <iostream>
#endif

int read_HLT_ESDs(TString file = "AliESDs.root")
{

   TH2F *clusterPosXY = new TH2F("clusterPosXY", "clusterPosXY", 1000, -500, 500, 1000, -500, 500);
   TH2F *clusterPosXZ = new TH2F("clusterPosXZ", "clusterPosXZ", 1000,  -500, 500, 1000, -500, 500);
   
   TH1F *clusterEnergies = new TH1F("clusterEnergies", "clusterEnergies", 1100, 0, 110);
   
   TH1I *clusterMult = new TH1I("clusterMult", "clusterMult", 1000, 1, 1000);
   TFile *esdFile = TFile::Open(file, "READ");
   
   TTree *esdTree = dynamic_cast<TTree*>(esdFile->Get("HLTesdTree"));
   
   if(!esdTree) return -1;
   
   std::cout << "Number of events in ESD: " << esdTree->GetEntries() << std::endl;
   
   TClonesArray *clusters = 0;
   
   esdTree->SetBranchAddress("CaloClusters", &clusters);
   
   for(int ev = 0; ev < esdTree->GetEntries(); ev++)
   {
      esdTree->GetEntry(ev);
      
      cout << "Number of clusters in event: " << clusters->GetEntries() << endl;
      
      for(int c = 0; c < clusters->GetEntries(); c++)
      {
	 AliESDCaloCluster *cl = dynamic_cast<AliESDCaloCluster*>(clusters->At(c));
      
	 cout << "\tCluster #: " << c << ",   energy: " << cl->E() << endl;
	 
	 Float_t pos[3];
	 
	 cl->GetPosition(pos);
	 
	 cout << "\tx: " << pos[0] << ", y: " << pos[1] << ", z: " << pos[2] << endl;
	 
	 clusterPosXY->Fill(pos[0], pos[1]);
	 clusterPosXZ->Fill(pos[0], pos[2]);
	 
	 clusterEnergies->Fill(cl->E());
	 clusterMult->Fill(cl->GetNCells());
	 
	 
      }
   }
   //   AliESDEvent *esdEvent = reinterpret_cast<AliESDEvent*>(esdFile->Get("
   
//   clusterPos->Draw();
   
   //clusterEnergies->Draw();

   TCanvas *c1 = new TCanvas("c1", "", 0, 0, 800, 450);
   c1->Divide(2, 1);
   c1->cd(1);
   clusterEnergies->Draw();
   c1->cd(2);
   clusterMult->Draw();
   
   
   TCanvas *c2 = new TCanvas("c2", "", 10, 10, 800, 450);
   c2->Divide(2, 1);
   c2->cd(1);
   clusterPosXY->Draw();
   c2->cd(2);
   clusterPosXZ->Draw();
	    
   return 0;
}