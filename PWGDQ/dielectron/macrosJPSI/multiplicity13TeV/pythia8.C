/// \file pythia8.C
/// \brief Simulates Pythia8 Events
///
/// Best used with the scripts in this folder
///
/// \author Julius Gronefeld <j.gronefeld@cern.ch>, GSI-Darmstadt
/// \date MAR 19, 2015
/// \see pythia8MeanPt.C
AliGenerator*  CreateGenerator();

void pythia8(Int_t taskid=1, Int_t mpiLevel = 3, Int_t nev = 1000000, Double_t energy = 13000, Int_t estimator = 1, const char* outfilename="jpsi13TeV.root")
{
  //
  //
Double_t multFactor = estimator ==1 ? 1.:4.;
  Int_t multNbins = 200;  
  Double_t binsMult[ 200 ];
  for (Int_t i=0; i<=multNbins; i++) { 
    binsMult[i] = -0.5 +  i* multFactor; 
  }
  // change binning


  
  Int_t ptNbins = 400;
  Double_t binsPt[ 400 ];
  for (Int_t i=0; i<=ptNbins; i++) { binsPt[i] = -0.5 +  0.1 * i; }
  
  Int_t yNbins = 200;
  Double_t binsY[ 200 ];
  for (Int_t i=0; i<=yNbins; i++) { binsY[i] = -1. +  0.01 * i; }
  

  TH1D* h1MCEvents = new TH1D("h1MCEvents","h1MCEvents",multNbins,binsMult);
  TH1D* histJpsiVsMult = new TH1D("histJpsiVsMult","histJpsiVsMult",multNbins,binsMult);
  TH2D* histJpsiPtVsMult = new TH2D("histJpsiPtVsMult","histJpsiPtVsMult",multNbins,binsMult,   ptNbins, binsPt);
  TH2D* histJpsiYVsMult = new TH2D("histJpsiYVsMult","histJpsiYVsMult",multNbins,binsMult,   yNbins, binsY);
  



h1MCEvents->GetXaxis()->SetTitle("multMC true in |#eta|<1.0");
  
  //  Runloader
  gSystem->Load("liblhapdf.so");    
  gSystem->Load("libEGPythia6.so"); 
  gSystem->Load("libpythia6.so");   
  gSystem->Load("libAliPythia6.so");
  //gSystem->Load("libpythia8.so");   
  gSystem->Load("libAliPythia8.so");
//  gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8175/xmldoc"));
  gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
  gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
  gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));

  AliRunLoader* rl = AliRunLoader::Open("galice.root","FASTRUN","recreate");

  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(nev);
//  rl->LoadKinematics("RECREATE");
  rl->MakeTree("E");
  gAlice->SetRunLoader(rl);

  //  Create stack
  rl->MakeStack();
  AliStack* stack      = rl->Stack();

  //  Header
  AliHeader* header = rl->GetHeader();
  //
  //  Create and Initialize Generator
  AliGenerator *gener = CreateGenerator(energy,taskid);

  if(mpiLevel > -1 && mpiLevel < 4){ 
    cout << "Setting MultipartonInteractions:processLevel to " << mpiLevel << endl;
    AliPythia8::Instance()->ReadString( Form("MultipartonInteractions:processLevel = %d", mpiLevel) );
  }
  else if(mpiLevel == -1 ){
    cout << "Switching off MPI" << mpiLevel << endl;
    AliPythia8::Instance()->ReadString( "PartonLevel:MPI = off"  ) ;
 }    


  gener->Init();
  AliPythia8::Instance()->ReadString("111:mayDecay  =  off");


  gener->SetStack(stack);

  //
  //                        Event Loop
  //
  Int_t iev;

  for (iev = 0; iev < nev; iev++) {
  //  if ( iev%100 ==0 ) cout<<"Event number: "<<iev<<" of "<<nev<<endl;
    //  Initialize event
    header->Reset(0,iev);
  //  rl->SetEventNumber(iev);
    stack->Reset();
   // rl->MakeTree("K");
    //	stack->ConnectTree();

    //  Generate event
    gener->Generate();

    //  Analysis
    Int_t npart = stack->GetNprimary();
    Int_t nTrk = stack->GetNtrack();
    
//     cout << "part: "<<npart <<endl;
//     cout << "trk: "<<nTrk <<endl;
//     continue;
    // Calculate number of charged particles
    
    Int_t multEta10 = 0; // mult in eta +-1.0
    Int_t nJpsi = 0; // mult in eta +-1.
    Int_t nPrimary = 0;
    Int_t nSecondary = 0;

    
    for (Int_t part = 0; part < npart; part++) {
      TParticle *particle = stack->Particle(part);
      if(
        stack->IsPhysicalPrimary(part)){
        
          ++nPrimary;
        
          if(
          particle->GetPDG()->Charge() != 0 && 
          TMath::Abs(particle->Eta()) < 1.0 && 
          particle->Pt() >  0. &&
          particle->Energy() > 0.
        ){
          ++multEta10;
        }
        else  if( TMath::Abs( particle->GetPdgCode() ) == 443  && 
          particle->Pt() > 0. &&
          particle->Energy() > 0. &&
          TMath::Abs(particle->Y()) < 0.9 ) {

          ++nJpsi;


        }
      }
      else{
       nSecondary++;
      }

    } 
    
    if(estimator == 1){
      h1MCEvents->Fill(multEta10);
    }
    else if( estimator ==0){
      h1MCEvents->Fill(nPrimary);
    }
     if(nJpsi){    
          if(estimator ==1){
            histJpsiVsMult->Fill(multEta10, nJpsi);
          }
          else if(estimator ==0){
            histJpsiVsMult->Fill(nPrimary, nJpsi);
          }
         // histJpsiPtVsMult->Fill(multEta10, particle->Pt(), nJpsi );
         // histJpsiYVsMult->Fill(multEta10, particle->Y(), nJpsi );
    }

  //  Finish event
    header->SetNprimary(stack->GetNprimary());
    header->SetNtrack(stack->GetNtrack());  
    //      I/O
    //	
    stack->FinishEvent();
 //   header->SetStack(stack);
   // rl->TreeE()->Fill();

//    rl->WriteKinematics("test");

  } // event loop
  TFile* fout = TFile::Open(outfilename,"RECREATE");
  h1MCEvents->Write();
  histJpsiVsMult->Write();
  histJpsiPtVsMult->Write();
  histJpsiYVsMult->Write();
  fout->Close();
  //
  //                         Termination
  //  Generator
  

//gener->FinishRun();
  //  Write file
//  rl->WriteHeader("OVERWRITE");
//  gener->Write();
//  rl->Write();


}


AliGenerator*  CreateGenerator(Double_t energy, Int_t taskid=1)
{
  AliGenPythiaPlus* gener = new AliGenPythiaPlus(AliPythia8::Instance());
  //
  //
      TDatime dt;
      UInt_t seed = dt.Get() * taskid;
  //if (gSystem->Getenv("CONFIG_SEED")) {
  //  seed = atoi(gSystem->Getenv("CONFIG_SEED"));
    cout     <<"Random seed set to "<< seed << endl;
  //} 	  


  (AliPythia8::Instance())->ReadString("Random:setSeed = on");
  (AliPythia8::Instance())->ReadString(Form("Random:seed = %ld", seed%900000000)); 
  // Set tune to Monash
  (AliPythia8::Instance())->ReadString("Tune:pp = 14");
//  gener->SetProcess(kPythia8_Monash2013); // Try this. Not working on current Aliroot
  gener->SetProcess(kPyMbDefault);
  //   Centre of mass energy 
  gener->SetEnergyCMS(energy);
  //   Initialize generator    
  //  gener->SetEventListRange(-1, 10);
  return gener;
}









