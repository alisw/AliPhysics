//
// This macro is a part of "Alice PPR fast flow analysis package"
//
// The macro does:
// 1) open existing database with Flow Pico DST ("flowPicoEvent.root")
// 2) open source data - galice.root type file with events
// 3) reconstruct psi, psiA, psiB (sub-events in eta space) and v2.
// 4) fill Flow Pico DST with reconstructed data
//
// INPUT PARAMETERS:
// inputName: file with events
// if (inputName == 0) macro tries to fetch system variable FILE_NAME
// and use its value as the name of the file. (useful for bath processing)
//
// maxN: maximum number of events to be processed, if 0 process all events
// 
// NOTE:
// file "flowPicoEvent.root" with "flowData" NTuple have to exist
// it can be created by "AliFlowCreate.C"
//
// Sylwester Radomski, GSI
// mail: S.Radomski@gsi
// 23. Oct. 2002
//
//////////////////////////////////////////////////////////////////////////////// 

//#include "TTree.h"
//#include "TNtuple.h"
//#include "TFile.h"
//#include "AliRun.h"

void AliFlowReconstruction(const char* inputName = 0, Int_t maxN = 0) {

  const char *dataName = "flowPicoEvent.root";
  
  TTree  *kine;
  TNtuple *data;

  Double_t trueV2, truePsi;
  Double_t v2, psi, psiA, psiB;
  Int_t N;
  Int_t mult;

  if (gAlice) delete gAlice;
  gAlice = 0;

  TFile *dataFile = new TFile(dataName,"UPDATE");
  data = (TNtuple *) dataFile->Get("flowData");  
  
  TFile *inputFile;
  
  if (!inputName) {
    TFile *inputFile = new TFile(gSystem->Getenv("FILE_NAME"), "READ"); 
    ::Info("AliFlowReconstruction", "Using: %s",gSystem->Getenv("FILE_NAME"));
  } else {
    TFile *inputFile = new TFile(inputName, "READ");
  }
  
  gAlice = (AliRun*) inputFile->Get("gAlice");

  N = gAlice->GetEventsPerRun();
  if (maxN != 0 && maxN < N) N = maxN;
 
  ::Info("AliFlowReconstruction", "Number of events to be processed = %d", N);

  for (Int_t i=0; i<N; i++) {
    
    ::Info("Processing event i = %d", i);

    gAlice->GetEvent(i);
    
    AliHeader *header = gAlice->GetHeader();
    AliGenGeVSimEventHeader *egHeader = (AliGenGeVSimEventHeader*)header->GenEventHeader();

    truePsi = egHeader->GetEventPlane() * 180 / TMath::Pi();
    while(truePsi > 90) truePsi -= 180;

    trueV2 = egHeader->GetEllipticFlow(211);

    mult = gAlice->GetNtrack();
    kine = gAlice->TreeK();

    psi =  FlowPsi(kine, 0);
    psiA = FlowPsi(kine, 1);
    psiB = FlowPsi(kine, -1);
    
    v2 = FlowV2(kine, psi);
    data->Fill(1, i, mult, truePsi, trueV2, psi, psiA, psiB, v2);
  }
   
  dataFile->cd();
  data->Write();

  delete data;

  delete gAlice;
  gAlice = 0;

  dataFile->Close();
  inputFile->Close();
  
  delete dataFile;
  delete inputFile;

  gSystem->Exit(0);
}

////////////////////////////////////////////////////////////////////////////////

Double_t FlowPsi(TTree *kine, Int_t etaPart) {
  //
  // this function calculates psi from list of particles
  // etaPart - division to sub-events in eta 
  // 0   all particles,
  // -1  negative eta 
  // 1   positive eta
  //

  Int_t mult = kine->GetEntries();
  Double_t *phi, *eta;
  Double_t psi;
  
  Double_t Ssin = 0, Scos = 0;

  kine->Draw("Particles->Phi():Particles->Eta()", "", "goff");
  phi = kine->GetV1();
  eta = kine->GetV2();

  for (Int_t i=0; i<mult; i++) {
    
    if ( (etaPart * eta[i]) < 0 ) continue;
  
    Ssin += sin( 2 * phi[i] );
    Scos += cos( 2 * phi[i] );
  }
  
  //psi = atan( Ssin / Scos ) / 2;
  psi = atan2 (Ssin, Scos) / 2;
  psi = psi * 180. / TMath::Pi(); 
  
  return psi;
}

////////////////////////////////////////////////////////////////////////////////

Double_t FlowV2(TTree *kine, Double_t psi) {

  Int_t mult = kine->GetEntries();
  Double_t *phi;
  Double_t V2;

  Double_t Ssin = 0, Scos = 0;
  
  kine->Draw("Particles->Phi()", "", "goff");
  phi = kine->GetV1();

  for (Int_t i=0; i<mult; i++) {
    
    Ssin += sin( 2 * (phi[i] - psi) );
    Scos += cos( 2 * (phi[i] - psi) );
  }
  
  V2 = sqrt(Ssin*Ssin + Scos*Scos) / mult;

  return V2;
}

////////////////////////////////////////////////////////////////////////////////
