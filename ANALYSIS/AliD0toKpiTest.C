//--------------------------------------------------------------------------
// Test macro for reconstruction and analysis of D0->Kpi
//
//     Andrea Dainese, andrea.dainese@lnl.infn.it
//--------------------------------------------------------------------------

void AliD0toKpiReco() {
  
  gSystem->Load("libANALYSIS.so");

  //==============  R E C O N S T R U C T I O N ==============================

  Int_t evFirst = 0;
  Int_t evLast  = 1000000;

  // Get field from galice.root
  if (gAlice) {
    delete gAlice->GetRunLoader();
    delete gAlice; 
    gAlice=0;
  }  
  AliRunLoader *rl = AliRunLoader::Open("galice.root");
  if (rl == 0x0) {
    cerr<<"Can not open session"<<endl;
    return;
  }
  Int_t retval = rl->LoadgAlice();
  if (retval) {
    cerr<<"LoadgAlice returned error"<<endl;
    delete rl;
    return;
  }
  gAlice=rl->GetAliRun();
  AliMagF *fiel = (AliMagF*)gAlice->Field();
  // Set the conversion constant between curvature and Pt
  AliTracker::SetFieldMap(fiel,kTRUE);

  AliD0toKpiAnalysis *analysis = new AliD0toKpiAnalysis();
  // set simulation to take info on PDG codes from kine tree
  analysis->SetSimulation();
  rl->LoadKinematics(); 
  analysis->MakeTracksRefFile(gAlice,evFirst,evLast);
  //--- set this is you want only signal candidates in output file
  //analysis->SetOnlySignal();
  //--- set this if you want to compute primary vertex D0 by D0 using 
  //    other tracks in the event (for pp, broad interaction region);
  //    it is time-consuming procedure, so it can be done after a 
  //    preselection on invariant mass
  //analysis->SetVertexOnTheFly();
  //analysis->SetMassCut(0.1); // GeV
  //--- set single-track preselections
  analysis->SetPtCut(0.); // GeV
  analysis->Setd0Cut(0.); // micron
  //--- set cuts on D0 candidates to be written to file
  //    (see AliD0toKpiAnalysis.h for a description and for the defaults)
  //analysis->SetD0Cuts(0.1,1000.,1.1,0.,0.,10000.,10000.,0.,.5);
  analysis->SetD0Cuts();

  //--- check the settings
  analysis->PrintStatus();

  analysis->FindCandidates(evFirst,evLast,"AliD0toKpi.root");
  delete analysis;

  return;
}
//==========================================================================
void AliD0toKpiSele() {  

  gSystem->Load("libANALYSIS.so");

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
 





