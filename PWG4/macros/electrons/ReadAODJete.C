void ReadAODJete(){
      gSystem->Load("libTree.so");
      gSystem->Load("libPhysics.so");
      gSystem->Load("libGeom.so");
      gSystem->Load("libVMC.so");
      gSystem->Load("libCGAL");
      gSystem->Load("libfastjet");
      gSystem->Load("libSISConePlugin");
      gSystem->Load("libSTEERBase");
      gSystem->Load("libESD");
      gSystem->Load("libAOD");
      gSystem->Load("libANALYSIS");
      gSystem->Load("libANALYSISalice");
      gSystem->Load("libJETAN");
      gSystem->Load("libFASTJETAN");
      gSystem->Load("libPWG4PartCorrBase");
      gSystem->Load("libPWG4PartCorrDep");


  TFile* fin=new TFile("aodoutput.root");

  TTree *aodTree = (TTree*)fin->Get("aodTree");
  AliAODEvent *ev = new AliAODEvent();
  ev->ReadFromTree(aodTree);
  Int_t nEvents = aodTree->GetEntries();
  cout<<nEvents<<endl;

  //header information
  AliFastJetHeaderV1* hd=( AliFastJetHeaderV1*)aodTree->GetUserInfo()->FindObject("AliJetHeader_jets");
  cout<<"Rparam = "<<hd->GetRparam()<<" BGmode = "<<hd->GetBGMode()<<endl;
  TH1F* h1=new TH1F("h1","",100,0,1000);
  TH1F* h2=new TH1F("h2","",100,0,1000);
  TH1F* h3=new TH1F("h3","",100,0,1000);
  TH1F* h4=new TH1F("h4","",100,0,1000);

  for (Int_t nEv = 0; nEv < nEvents; nEv++) {

    aodTree->GetEvent(nEv);

    //getting jets
    Int_t nJets = ev->GetNJets();
    for(Int_t iJet=0; iJet<nJets; iJet++){ 
      jet = (AliAODJet*)ev->GetJet(iJet);  
      cout<<jet->Pt()<<endl;
      cout<<"*****************************START JET PRINTOUT**************************"<<endl;
      jet->Print();
      cout<<"*****************************END   JET PRINTOUT**************************"<<endl;
    }

     //getting bkg info
    AliAODJetEventBackground* evBkg=(AliAODJetEventBackground*)ev->FindListObject("jeteventbackground");

    cout<<"*******************************evBkg = "<<evBkg<<endl;

    if(evBkg){
      cout<<"*****************************START BACKGROUND PRINTOUT********************"<<endl;
      evBkg->Print();

      Float_t bkg1=evBkg->GetBackground(0);
      Float_t bkg2=evBkg->GetBackground(1);
      Float_t bkg3=evBkg->GetBackground(2);
      Float_t bkg4=evBkg->GetBackground(3);

      cout<<"ev = "<<nEv<<" bkg FastJet ="<<bkg1<<" bkg FastJet Charg="<<bkg2<<" bkg OutOfCone ="<<bkg3<<" bkg OutOfJet ="<<bkg4<<endl;
      cout<<"*****************************END   BACKGROUND PRINTOUT********************"<<endl;
      h1->Fill(bkg1);
      h2->Fill(bkg2);
      h3->Fill(bkg3);
      h4->Fill(bkg4);
    }

  }
  h1->GetXaxis()->SetTitle("#rho");
  h2->GetXaxis()->SetTitle("#rho");
  h3->GetXaxis()->SetTitle("#rho");
  h4->GetXaxis()->SetTitle("#rho");
  h1->GetYaxis()->SetTitle("Entries");
  h2->GetYaxis()->SetTitle("Entries");
  h3->GetYaxis()->SetTitle("Entries");
  h4->GetYaxis()->SetTitle("Entries");

  h1->SetTitle("FastJet R=0.2 All");
  h2->SetTitle("FastJet R=0.4 Charged");
  h3->SetTitle("Out-of-cone");
  h4->SetTitle("Out-of-jet");


  TCanvas* c1=new TCanvas();
  c1->Divide(2,2);
  c1->cd(1);
  h1->Draw();
  c1->cd(2);
  h2->Draw();
  c1->cd(3);
  h3->Draw();
  c1->cd(4);
  h4->Draw();

}
