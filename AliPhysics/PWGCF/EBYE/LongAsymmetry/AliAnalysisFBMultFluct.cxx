#include <TChain.h>
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TFile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

//#include "AliESDEvent.h"
//#include "AliESDInputHandler.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

// Monte Carlo Event
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "AliAnalysisFBMultFluct.h"
#include "AliGenHijingEventHeader.h"

#include "AliCentrality.h"
#include "AliAODVertex.h"

#include "AliAODVZERO.h" 
#include <iostream>

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing
// Reviewed: A.Gheata (19/02/10)

using namespace std;

ClassImp(AliAnalysisFBMultFluct)

//________________________________________________________________________
/*AliAnalysisFBMultFluct::AliAnalysisFBMultFluct() 
  : AliAnalysisTaskSE(), fAOD(0), fOutputList(0) 
  , fHistVertexZ(0), fHistVertexXY(0)
  , fHistCentrality(0)
  , fHistZDCN1(0), fHistZDCN2(0)
  , fHistPt(0), fHistEta(0)//, fHistEtaAsymBin(0)
  , fHistZDCAsym(0)
  , fHistZDCAsymvsVz(0), fHistZDCAsymvsCent(0)
  , fHistLog(0)


{
  // Constructor
  // For defining multiple histograms??? RR
  // memset(fHnSigmaTPC, 0, sizeof(TH1F *) * AliPID::kSPECIES * 2);

  // Define input and output slots here
  // Input slot #0 works with a TChain
  //DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH container

  //  TH1F        *fHistEtaAsymBin[nslice]; // Eta spectrum for diff Asym Bin  
  for (Int_t i = 0; i < nslice; i++ ) fHistEtaAsymBin[i] = NULL;

   DefineOutput(1, TList::Class());
   }*/
//________________________________________________________________________
AliAnalysisFBMultFluct::AliAnalysisFBMultFluct(const char *name) 
  : AliAnalysisTaskSE(name), fAOD(0), fOutputList(0) 
  , fHistVertexZ(0), fHistVertexXY(0)
  , fHistCentrality(0)
  , fHistZDCN1(0), fHistZDCN2(0)
  , fHistTrackPt(0), fHistTrackEta(0) 
  , fHistTrackletEta(0)
//, fHistEtaAsymBin(0) histogram arrays are initalized differently
  , fHistNeventZDCAsym(0)
  , fHistZDCAsym(0)
  , fHistZDCAsymvsVz(0), fHistZDCAsymvsCent(0)
  , fHistTrackletMult(0)
  , fHistZDCAsymvsTrackletAsym(0), fHistTracklets_pos_neg(0)
  , fHistLog(0)
  , fV0SubeventAPlane(0), fV0SubeventCPlane(0)
  , fEventTree(0)


{
  // Constructor
  // For defining multiple histograms??? RR
  // memset(fHnSigmaTPC, 0, sizeof(TH1F *) * AliPID::kSPECIES * 2);

  // Define input and output slots here
  // Input slot #0 works with a TChain
  //DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH container

  //  TH1F        *fHistEtaAsymBin[nslice]; // Eta spectrum for diff Asym Bin  
  for (Int_t i = 0; i < nslice; i++ ) fHistTrackletEtaAsymBin[i] = NULL;

   DefineOutput(1, TList::Class());
   DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
void AliAnalysisFBMultFluct::UserCreateOutputObjects()
{
  // Create histograms  and add them to the respective lists

  // Called once

  fOutputList = new TList();

  fOutputList->SetOwner(kTRUE);

  //  fAnaHistList = new TList();

  Float_t AsymMin[nslice]={-10.0, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4};
  Float_t AsymMax[nslice]={-0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 10.0};
  for(Int_t i = 0;i<nslice;i++){
    asymMin[i]=AsymMin[i];
    asymMax[i]=AsymMax[i];    
  }

  fHistCentrality = new	TH1F("fHistCentrality", "Centrality of Events; Centrality Percentile (%); N",50,0.0,100.0);
  fHistVertexZ = new TH1F("fHistVertexZ","Z Coord. of primary vertex",100,-50.0,50.0);
  fHistVertexZ->GetXaxis()->SetTitle(" Vz in cm.");
  fHistVertexXY = new TH2F("fHistVertexXY","XY Coord. of primary vertex",100,-0.3,0.3,100,-0.3,0.3);
  
  fHistZDCN1 = new TH1F("fHistZDCN1","ZDC1 Neutron Energy",1000,0.0,100000); 
  fHistZDCN2 = new TH1F("fHistZDCN2","ZDC2 Neutron Energy",1000,0.0,100000); 

  fHistTrackPt = new TH1F("fHistTrackPt", "PT distribution", 15, 0.1, 3.1);
  fHistTrackPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistTrackPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistTrackPt->SetMarkerStyle(kFullCircle);

  fHistTrackEta = new TH1F("fHistTrackEta", "Eta distribution of Tracks", 40, -2.0, 2.0);
  fHistTrackEta->GetXaxis()->SetTitle("Eta ");
  fHistTrackEta->GetYaxis()->SetTitle("dN/dEta");
  fHistTrackEta->SetMarkerStyle(kFullCircle);

  fHistTrackletEta = new TH1F("fHistTrackletEta", "Eta distribution of Tracklets", 40, -2.0, 2.0);
  fHistTrackletEta->GetXaxis()->SetTitle("Eta of Tracklets ");
  fHistTrackletEta->GetYaxis()->SetTitle("dN/dEta");
  fHistTrackletEta->SetMarkerStyle(kFullCircle);


  Char_t hname[48];
  for(Int_t ibin =0;ibin<nslice;ibin++){
    //    sprintf(hname,"fHistEta_AsymBin%3.2f_%3.2f",TMath::Abs(asymMin[ibin]),TMath::Abs(asymMax[ibin]));
    sprintf(hname,"fHistEta_AsymBin%3.2f_%3.2f",asymMin[ibin],asymMax[ibin]);
    //    cout<<" hname is "<<hname<<endl;
    fHistTrackletEtaAsymBin[ibin] = new TH1F(hname,"Eta distribution of Tracklets fro ITS", 40, -2.0, 2.0);
    fHistTrackletEtaAsymBin[ibin]->GetXaxis()->SetTitle("Eta for Asym Bin");
    fHistTrackletEtaAsymBin[ibin]->GetYaxis()->SetTitle("dN/dEta");
    fHistTrackletEtaAsymBin[ibin]->SetMarkerStyle(kFullCircle);
    fOutputList->Add(fHistTrackletEtaAsymBin[ibin]);
  }

  fHistNeventZDCAsym = new TH1F("fHistNeventZDCAsym", "No. of Events in ZDC Asymmetry bin", 10, -0.5, 9.5);
  fHistNeventZDCAsym->GetXaxis()->SetTitle("ZDC Asym bin number ");
  fHistNeventZDCAsym->GetYaxis()->SetTitle("Number of Event");
  fHistNeventZDCAsym->SetMarkerColor(2);
  fHistNeventZDCAsym->SetMarkerStyle(kFullCircle);

  fHistZDCAsym = new TH1F("fHistZDCAsym", "ZDC n1 n2 Asymmetry distribution", 36, -1.8, 1.8);
  fHistZDCAsym->GetXaxis()->SetTitle("ZDC (n1-n2)/(n1+n2)  ");
  fHistZDCAsym->GetYaxis()->SetTitle("Number of Event");
  fHistZDCAsym->SetMarkerColor(2);
  fHistZDCAsym->SetMarkerStyle(kFullCircle);

  fHistZDCAsymvsVz = new TH2F("fHistZDCAsymvsVz", "ZDC n1 n2 Asymmetry vs Vertex Z", 36, -1.8, 1.8, 50,-25.0,25.0);
  fHistZDCAsymvsVz->GetXaxis()->SetTitle("ZDC (n1-n2)/(n1+n2)  ");
  fHistZDCAsymvsVz->GetYaxis()->SetTitle("Event Vertex Z in cm.");
  fHistZDCAsymvsVz->SetMarkerColor(4);
  fHistZDCAsymvsVz->SetMarkerStyle(kFullCircle);

  fHistZDCAsymvsCent = new TH2F("fHistZDCAsymvsCent", "ZDC n1 n2 Asymmetry vs EventCentrality", 36, -1.8, 1.8, 10,0.0,100.0);
  fHistZDCAsymvsCent->GetXaxis()->SetTitle("ZDC (n1-n2)/(n1+n2)  ");
  fHistZDCAsymvsCent->GetYaxis()->SetTitle("Event Centrality Percentile");
  fHistZDCAsymvsCent->SetMarkerColor(5);
  fHistZDCAsymvsCent->SetMarkerStyle(kFullCircle);

  fHistTrackletMult= new TH1F("TrackletMult","TrackletMult",100,0.0,3000);
  fHistTrackletMult->GetXaxis()->SetTitle("TrackletMult |eta|<1.0  ");
  fHistTrackletMult->GetYaxis()->SetTitle("Number of events");
  fHistTrackletMult->SetMarkerColor(6);
  fHistTrackletMult->SetMarkerStyle(kFullCircle);
  

  fHistZDCAsymvsTrackletAsym= new TProfile("ZDCAsymvsTrackletAsym","ZDCAsymvsTrackletAsym",20,-1.0,1.0);
  fHistZDCAsymvsTrackletAsym->GetXaxis()->SetTitle("ZDC (n1-n2)/(n1+n2)  ");
  fHistZDCAsymvsTrackletAsym->GetYaxis()->SetTitle("ITS (t1-t2)/(t1+t2)  ");
  fHistZDCAsymvsTrackletAsym->SetMarkerColor(6);
  fHistZDCAsymvsTrackletAsym->SetMarkerStyle(kFullCircle);
  
  fHistTracklets_pos_neg= new TH2F("Tracklet_pos_neg","Tracklet_pos_neg",100,0,2000,100,0,2000);
  fHistTracklets_pos_neg->GetXaxis()->SetTitle("ITS Tracklet Mult +eta  ");
  fHistTracklets_pos_neg->GetYaxis()->SetTitle("ITS Tracklet Mult -eta  ");
  fHistTracklets_pos_neg->SetMarkerColor(7);
  fHistTracklets_pos_neg->SetMarkerStyle(kFullCircle);
  
  fHistLog = new TH1F("fHistLog","Log Variables",100, -0.5, 99.5);

  fV0SubeventAPlane = new TH1F("fV0SubeventAPlane","Subevent A Plane from V0",8,0.0,TMath::Pi());
  fV0SubeventCPlane = new TH1F("fV0SubeventCPlane","Subevent C Plane from V0",8,0.0,TMath::Pi());
  
  fOutputList->Add(fHistVertexXY);
  fOutputList->Add(fHistVertexZ);
  fOutputList->Add(fHistCentrality);
  fOutputList->Add(fHistZDCN1);
  fOutputList->Add(fHistZDCN2);
  fOutputList->Add(fHistTrackPt);
  fOutputList->Add(fHistTrackEta);
  fOutputList->Add(fHistNeventZDCAsym);
  fOutputList->Add(fHistZDCAsym);
  fOutputList->Add(fHistZDCAsymvsVz);
  fOutputList->Add(fHistZDCAsymvsCent);
  fOutputList->Add(fHistTrackletMult);
  fOutputList->Add(fHistZDCAsymvsTrackletAsym);						   
  fOutputList->Add(fHistTracklets_pos_neg);						   
  fOutputList->Add(fHistLog);
  fOutputList->Add(fV0SubeventAPlane);
  fOutputList->Add(fV0SubeventCPlane);
  
  // Adding a Tree


  fEventTree = new TTree("fEventTree","fEventTree");
  fEventTree->Branch("fRunNumber",&fRunNumber,"fRunNumber/I");
  fEventTree->Branch("fNumberOfTracks",&fNumberOfTracks,"fNumberOfTracks/I");
  fEventTree->Branch("fNumberOfTracklets",&fNumberOfTracklets,"fNumberOfTracklets/I");
  fEventTree->Branch("fVertexX",&fVertexX,"fVertexX/F");
  fEventTree->Branch("fVertexY",&fVertexY,"fVertexY/F");
  fEventTree->Branch("fVertexZ",&fVertexZ,"fVertexZ/F");
  fEventTree->Branch("fCentPercentile",fCentPercentile,"fCentPercentile[4]/F");

  fEventTree->Branch("Qx",&Qx,"Qx/D");
  fEventTree->Branch("Qy",&Qy,"Qy/D");
  fEventTree->Branch("Q1x",&Qx,"Q1x/D");
  fEventTree->Branch("Q1y",&Qy,"Q1y/D");
  fEventTree->Branch("Q2x",&Qx,"Q2x/D");
  fEventTree->Branch("Q2y",&Qy,"Q2y/D");
  
  fEventTree->Branch("fzdcn1en",&fzdcn1en,"fzdcn1en/F");
  fEventTree->Branch("fzdcn2en",&fzdcn2en,"fzdcn2en/F");
  fEventTree->Branch("fzdcp1en",&fzdcp1en,"fzdcp1en/F");
  fEventTree->Branch("fzdcp2en",&fzdcp2en,"fzdcp2en/F");
  fEventTree->Branch("fzdcEM1en",&fzdcEM1en,"fzdcEM1en/F");
  fEventTree->Branch("fzdcEM2en",&fzdcEM2en,"fzdcEM2en/F");

  fEventTree->Branch("fTrackletEta",fTrackletEta,"fTrackletEta[fNumberOfTracklets]/F");
  fEventTree->Branch("fTrackletPhi",fTrackletPhi,"fTrackletPhi[fNumberOfTracklets]/F");

  fEventTree->Branch("fTrackPt",fTrackPt,"fTrackPt[fNumberOfTracks]/F");
  fEventTree->Branch("fTrackPhi",fTrackPhi,"fTrackPhi[fNumberOfTracks]/F");
  fEventTree->Branch("fTrackEta",fTrackEta,"fTrackEta[fNumberOfTracks]/F");
  fEventTree->Branch("fTrackCharge",fTrackCharge,"fTrackCharge[fNumberOfTracks]/I");
  fEventTree->Branch("fTrackFilterMap",fTrackFilterMap,"fTrackFilterMap[fNumberOfTracks]/I");
  fEventTree->Branch("fV0Mult",fV0Mult,"fV0Mult[64]/F");
  

  // NEW HISTO added to fOutputList here
  PostData(1, fOutputList); // Post data for ALL output slots >0 here, to get at least an empty histogram
  PostData(2, fEventTree); // Post data for ALL output slots >0 here, to get at least an empty histogram

  avfzdcn1 = 1.0;
  avfzdcn2 = 1.0;
  /*   
  TFile * fscale=NULL;
  fscale = TFile::Open("alien:///alice/cern.ch/user/r/rashmi/work0to20Cent/input/AOD.137161.input.root");
  //fscale = TFile::Open("AOD.137161.input.root");
  //if(fscale->IsOpen()==kFALSE)  
  //  fscale = TFile::Open("/alice/cern.ch/user/r/rashmi/work/AOD.137161.input.root");
  // cout<<" Is fscale open "<<fscale->IsOpen()<<endl;
  //  fscale->ls();

  TList *cscale= (TList*)fscale->Get("coutput");
  //  cscale->ls();
  if(fscale->IsOpen()==kFALSE){
    cout<<" No file for obtaining mean fZDCN1 and fZDCN2 so scalers would be 1.0"<<endl;
    avfzdcn1 = 1.0;
    avfzdcn2 = 1.0;
  }else{
    avfzdcn1 = ((TH1F*)cscale->FindObject("fHistZDCN1"))->GetMean(1);
    avfzdcn2 = ((TH1F*)cscale->FindObject("fHistZDCN2"))->GetMean(1);
    //    cout<<" Got Mean fZDC N1 and N2 as "<<avfzdcn1<<" "<<avfzdcn2<<endl;
    delete cscale; // Free memory used for TList
    delete fscale; // Free memory used for TFile 
  }

  fHistLog->SetBinContent(20,avfzdcn1);
  fHistLog->SetBinContent(21,avfzdcn2);
  */


}

//________________________________________________________________________
void AliAnalysisFBMultFluct::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  // Get INformation of the Mean of fZDCN1En andfZDCN2En for getting scaled asymmetry coeff
  
  
  // Post output data.
  /*
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }
  */
  //  printf(">>>>HERE  1\n");
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    printf("ERROR: fAOD not available\n");
    fHistLog->Fill(0);
    return;
  }

  fHistLog->Fill(1);  

  fRunNumber = fAOD->GetRunNumber();

  // For Event Centrality  
  AliCentrality* centrality = fAOD->GetCentrality();
  fCentPercentile[0] = (Float_t)centrality->GetCentralityPercentile("V0M");
  // This is zero all the time .... can be removed
  fCentPercentile[1] = (Float_t)centrality->GetCentralityPercentile("V0MEq");
  fCentPercentile[2] = (Float_t)centrality->GetCentralityPercentile("TRK");
  // This is either set upto <~44% or is 100. i.e. cent is set to 100 if 44%<cent<100
  fCentPercentile[3] = (Float_t)centrality->GetCentralityPercentile("ZEMvsZDC");
  //  if(fCentPercentile[0]<0.0 || fCentPercentile[0]>20.0) return;
  //  cout<<"fCentrality = "<<fCentPercentile[0]<<endl;
  if(fCentPercentile[0]<0.0) {
    fHistLog->Fill(3);
    return;
  }
  //  cout<<" Got Event of Centrality = "<<fCentPercentile[0]<<endl;
  
  // choose a particular trigger
  Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  if(isSelected==kFALSE) {
    //cout<<" Not selected by trigger"<<endl;
    fHistLog->Fill(13); // rejected by trigger condition
    return;
  }
  
  const AliAODVertex *vtx=fAOD->GetPrimaryVertex();
  // obtain vertex and apply vertex cuts
  //  printf(">>>>HERE  2b\n");
  if (!vtx) 
    {
      fHistLog->Fill(9); // rejected due to no vertex
      //    printf(">>>>HERE  2c\n");
      return;
    }
  //  printf(">>>>HERE  3\n");
  Double32_t fCov[6];
  vtx->GetCovarianceMatrix(fCov);
  if(vtx->GetNContributors() <= 0) { fHistLog->Fill(9); return;}
  if(fCov[5] == 0) {fHistLog->Fill(9); return;}
  
  fVertexX= (Float_t)vtx->GetX();
  fVertexY= (Float_t)vtx->GetY();
  fVertexZ= (Float_t)vtx->GetZ();
  //  cout<<" Vertex XYZ are="<<fVertexX<<" "<<fVertexY<<" "<<fVertexZ<<endl;
  
  fHistVertexZ->Fill(fVertexZ);
  fHistVertexXY->Fill(fVertexX,fVertexY);
  
  fHistLog->Fill(50);
  
  //  cout<<"Centrality Percentile is "<<fCentPercentile[0]<<" "<<fCentPercentile[1]<<" "<<fCentPercentile[2]<<" "<<fCentPercentile[3]<<endl;
  fHistCentrality->Fill(fCentPercentile[0]);

  //    printf(">>>>HERE  2\n");
  // Obtain the event plane information
  AliEventplane* Eventplane =  fAOD->GetEventplane();

  //Double_t Qx,Qy,Q1x,Q1y,Q2x,Q2y;
  // get the q-vectors by reference
  fAOD->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 8, 2, Q1x, Q1y);
  fAOD->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 9, 2, Q2x, Q2y);  
  //  cout<<"Q subevent1 x,y = "<<Q1x<<" "<<Q1y<<" angle = "<<TMath::ATan2(Q1y, Q1x)*180.0/TMath::Pi()<<endl; 
  //  cout<<"Q subevent2 x,y = "<<Q2x<<" "<<Q2y<<" angle = "<<TMath::ATan2(Q2y, Q2x)*180.0/TMath::Pi()<<endl;
  Double_t vzeroaEP = TMath::ATan2(Q1y, Q1x);           
  Double_t vzerocEP = TMath::ATan2(Q2y, Q2x);
  if(vzeroaEP<0.0) vzeroaEP = 2*TMath::Pi()+vzeroaEP;
  if(vzerocEP<0.0) vzerocEP = 2*TMath::Pi()+vzerocEP;
  
  fV0SubeventAPlane->Fill(vzeroaEP);
  fV0SubeventCPlane->Fill(vzerocEP);
  
//  printf(">>>>HERE  2plane   %p\n",Eventplane);
  if (Eventplane && Eventplane->GetQVector()) {
    TVector2 * Vevpl = Eventplane->GetQVector(); 
    Qx = Vevpl->X();
    Qy = Vevpl->Y();
    Q1x = Eventplane->GetQsub1()->X();
    Q1y = Eventplane->GetQsub1()->Y();
    Q2x = Eventplane->GetQsub2()->X();
    Q2y = Eventplane->GetQsub2()->Y();
    //    cout<<Q1x<<"+"<<Q2x<<"="<<Q1x+Q2x<<" Qx = "<<Qx<<endl;
    //    cout<<Q1y<<"+"<<Q2y<<"="<<Q1y+Q2y<<" Qy = "<<Qy<<endl;
  }
  //  else printf(">> EventPlane %p of QVector %p absent\n",Eventplane,Eventplane->GetQVector());
  //    printf(">>>>HERE  2a\n");
  //  if(TMath::Abs(Q1x+Q2x-Qx)>0.00001)
    //  if(TMath::Abs(Q1y+Q2y-Qy)>0.00001)
  // choose a particular event centrality

  // ZDC information
  fzdcn1en =  (Float_t)fAOD->GetZDCN1Energy();
  fzdcp1en =  (Float_t)fAOD->GetZDCP1Energy();
  fzdcn2en =  (Float_t)fAOD->GetZDCN2Energy();
  fzdcp2en =  (Float_t)fAOD->GetZDCP2Energy();
  fzdcEM1en =  (Float_t)fAOD->GetZDCEMEnergy(0);
  fzdcEM2en =  (Float_t)fAOD->GetZDCEMEnergy(1);
  if(fzdcn1en<=0.0||fzdcn2en<=0.0) {
    fHistLog->Fill(26);
    return;
  }
  //      printf(">>>>HERE  4\n");
  // cout<<" ZDC Energy N1,N2="<<fzdcn1en<<" "<<fzdcn2en<<endl;

  // Float_t ZDCNasym = (ZDCN1En/MeanZDCN1 - ZDCN2En/MeanZDCN2)/(ZDCN1En/MeanZDCN1 + ZDCN2En/MeanZDCN2);
  //  Double_t ZDCNasym = (fzdcn1en - fzdcn2en)/(fzdcn1en + fzdcn2en);
  // default values of average of ZDC neutron energy for both sides (A & C)

  // If average ZDC neutron energy is available, fill them here

  fHistZDCN1->Fill(fzdcn1en);
  fHistZDCN2->Fill(fzdcn2en);

  // Obtain the ZDC Asymmetry parameter and fill in the histogram or classify events accordingly.
  Double_t ZDCNasym =  GetAsymmetry(fzdcn1en,avfzdcn1,fzdcn2en,avfzdcn2);

  fHistZDCAsym->Fill(ZDCNasym);
  
  fHistZDCAsymvsVz->Fill(ZDCNasym,fVertexZ);
  fHistZDCAsymvsCent->Fill(ZDCNasym,fCentPercentile[0]);
  Int_t ibin;
  for(ibin=0;ibin<nslice;ibin++){
    if(ZDCNasym>asymMin[ibin] && ZDCNasym<=asymMax[ibin]) break;
  }
  fHistNeventZDCAsym->Fill(ibin);
  //  cout<<" ZDCNAsym is "<<ZDCNasym<<" selected ibin is "<<ibin<<endl;


  
  //Reading Tracklets in ITS from AOD-----------------------------------
  AliAODTracklets* mlt = fAOD->GetTracklets();
  if(!mlt){
    printf("Error: Could not recieve tracklets container\n");
  }
  fNumberOfTracklets = mlt->GetNumberOfTracklets();
  Int_t npos=0,nneg=0;
  for (int i=0;i<mlt->GetNumberOfTracklets();i++){
    
    fTrackletEta[i] = -log(tan(mlt->GetTheta(i)/2.0));
    fTrackletPhi[i] = mlt->GetPhi(i);

    //    printf("eta: %f  theta: %f phi: %f\n",-log(tan(mlt->GetTheta(i)/2)),mlt->GetTheta(i),mlt->GetPhi(i));
    if(fTrackletEta[i]>0.0 && fTrackletEta[i]<1.0) npos++;
    if(fTrackletEta[i]<0.0 && fTrackletEta[i]>-1.0) nneg++;
    fHistTrackletEta->Fill(fTrackletEta[i]);
    fHistTrackletEtaAsymBin[ibin]->Fill(fTrackletEta[i]);
  }
  //  if((npos+nneg)<100) 
  //cout<<" centrality = "<<fCentPercentile[0]<<" vertex  = "<<fVertexX<<" / "<<fVertexY<<" / "<<fVertexZ<<" ZDC n1/n2 = "<<fzdcn1en<<" / "<<fzdcn2en<<endl;
  fHistTracklets_pos_neg->Fill(npos,nneg);
  fHistTrackletMult->Fill(npos+nneg);
  
  Double_t TrackletAsym;
  if(npos>0 && nneg>0) {
    TrackletAsym = (1.0*npos-1.0*nneg)/(1.0*npos+1.0*nneg);
    fHistZDCAsymvsTrackletAsym->Fill(ZDCNasym,TrackletAsym);
    //    cout<<" Tracklet Asym = "<<mlt->GetNumberOfTracklets()<<" = "<<npos<<" + "<<nneg<<" asum = "<<TrackletAsym<<endl;
  }else{
    cout<<" Tracklet mult = "<<mlt->GetNumberOfTracklets()<<" = "<<npos<<" + "<<nneg<<endl;
  }
  
  // Track loop to fill a pT spectrum
  //  printf("There are %d tracks in this event\n", fAOD->GetNumberOfTracks());
  Int_t NTracks = 0; // This is the accepted track number
  //  Int_t nAllTracks = fAOD->GetNumberOfTracks();
  Int_t nAllTracks = (fAOD->GetTracks())->GetEntriesFast();
    for (Int_t iTrack = 0; iTrack < nAllTracks; iTrack++) {
    AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(iTrack);
    if (!track) {
      printf("ERROR: Could not receive track %d\n", iTrack);
      continue;
    }
// This to select Hybrid tracks in LHC10h as given in page
// https://twiki.cern.ch/twiki/bin/view/ALICE/PAGCorrelationsFAQ
//    if(track->TestFilterBit(272)==kFALSE) continue; 
//    if(track->TestFilterBit(768)==kFALSE) continue; 
    if(track->TestFilterBit(896)==kFALSE) continue; 
    fTrackPt[NTracks]  = (Float_t)track->Pt();
    fTrackPhi[NTracks] = (Float_t)track->Phi();
    fTrackEta[NTracks] = (Float_t)track->Eta();
    fTrackCharge[NTracks] = (Int_t)track->Charge();
    fTrackFilterMap[NTracks] = (Int_t)track->GetFilterMap();
    NTracks++;

    fHistTrackPt->Fill(track->Pt());
    fHistTrackEta->Fill(track->Eta());
  } //track 
  fNumberOfTracks = NTracks;
  //  cout<<" Number of AOD tracks are "<<fNumberOfTracks<<endl;
  // Getting VZERO information from AOD  
  
  AliAODVZERO * v0 = fAOD->GetVZEROData();
  for(Int_t ipm = 0;ipm<64;ipm++){
    fV0Mult[ipm]= v0->GetMultiplicity(ipm);
  }

  //  cout<<" Got all Track Info with Ntracks="<<fNumberOfTracks<<" NTracklets="<<fNumberOfTracklets<<endl;
  //  printf(">>>>HERE\n");
  fEventTree->Fill();

  PostData(1, fOutputList);
  PostData(2, fEventTree);
}      

//________________________________________________________________________
void AliAnalysisFBMultFluct::Terminate(Option_t *) 
{
}

//________________________________________________________________________
Double_t AliAnalysisFBMultFluct::GetAsymmetry( Double_t mult1, Double_t meanmult1, Double_t mult2, Double_t meanmult2)
{
  Double_t asym = ((mult1/meanmult1)-(mult2/meanmult2))/((mult1/meanmult1)+(mult2/meanmult2));
  
  return asym;
  
} 
