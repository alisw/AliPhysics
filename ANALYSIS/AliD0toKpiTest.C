//--------------------------------------------------------------------------
// Test macro for reconstruction and analysis of D0->Kpi
//
//     Andrea Dainese, andrea.dainese@pd.infn.it
//--------------------------------------------------------------------------

void AliD0toKpiReco() {
  
  //==============  R E C O N S T R U C T I O N ==============================

  // Look for field value in galice.root
  Double_t field = 0.4;
  if (gAlice) {
    delete gAlice->GetRunLoader();
    delete gAlice; 
    gAlice=0;
  }
  if(!gSystem->AccessPathName("galice.root",kFileExists)) {
    AliRunLoader *rl = AliRunLoader::Open("galice.root");
    rl->LoadgAlice();
    field=gAlice->Field();
    Double_t bz=field->SolenoidField()/10.;
    printf("B = %3.1f T read from gAlice and set\n",bz);
    delete gAlice->GetRunLoader();
    delete gAlice; 
    gAlice=0;
  } else {
    printf(" File galice.root not found: default %3.1f T being used!\n",field);
  }

  AliD0toKpiAnalysis *analysis = new AliD0toKpiAnalysis();
  //--- set magnetic field
  analysis->SetBz(field); 
  // set simulation to take info on PDG codes from kine tree
  //analysis->SetSimulation();
  //--- set this is you want only signal candidates in output file
  //analysis->SetOnlySignal();
  //--- set this if you want to compute primary vertex D0 by D0 using 
  //    other tracks in the event (for pp, broad interaction region);
  //    it is time-consuming procedure, so it can be done after a 
  //    preselection on invariant mass
  //analysis->SetVertexOnTheFly();
  //analysis->SetMassCut(0.1); // GeV
  //--- set single-track preselections
  analysis->SetPtCut(0.5); // GeV
  analysis->Setd0Cut(50.); // micron
  //--- set cuts on D0 candidates to be written to file
  //    (see AliD0toKpiAnalysis.h for a description and for the defaults)
  analysis->SetD0Cuts(0.1,1000.,1.1,0.,0.,10000.,10000.,0.,.5);
  //analysis->SetD0Cuts();

  //--- check the settings
  analysis->PrintStatus();

  Int_t evFirst = 0;
  Int_t evLast  = 1000000;
  //analysis->SetDebug();
  //analysis->FindCandidates(evFirst,evLast,"AliD0toKpi.root");
  analysis->FindCandidatesESD(evFirst,evLast,"AliD0toKpi.root");
  delete analysis;

  return;
}
//==========================================================================
void AliD0toKpiSele() {  

  //========================  S E L E C T I O N ============================

  AliD0toKpiAnalysis *analysis = new AliD0toKpiAnalysis();
  analysis->SetSimulation();
  analysis->SetOnlySignal();
  analysis->SetD0Cuts(0.1,1000.,1.1,0.,0.,10000.,10000.,0.,.5);
  analysis->ApplySelection("AliD0toKpi.root","AliD0toKpi_sele.root");
  delete analysis;

  return;
}
//==========================================================================
 





