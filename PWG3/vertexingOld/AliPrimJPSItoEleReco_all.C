//--------------------------------------------------------------------------
// Test macro for reconstruction and analysis of Primary J/psi into e+e-
//
//     Giuseppe Bruno, giuseppe.bruno@ba.infn.it
//      based on the for charm golden channel (D0->Kpi)
//--------------------------------------------------------------------------

void AliPrimJPSItoEleReco_all() {
  
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libAOD.so");
  gSystem->Load("libPWG3base.so");
  gSystem->Load("libPWG3.so");
  gSystem->Load("libPWG3vertexingOld.so");

  //==============  R E C O N S T R U C T I O N ==============================

  Int_t evFirst = 0;
  Int_t evLast  = 99;

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

  AliBtoJPSItoEleAnalysis *analysis = new AliBtoJPSItoEleAnalysis();
  // set simulation to take info on PDG codes from kine tree
  analysis->SetSimulation();
  rl->LoadKinematics(); 
  analysis->MakeTracksRefFile(gAlice,evFirst,evLast);
  // cout << "ci arrivo 2" << endl;
  //--- set this is you want only signal candidates in output file
   //analysis->SetOnlySignal();
   analysis->SetOnlyPrimaryJpsi();
   //analysis->SetOnlySignalAndPrimaryJpsi();
  //--- set this if you want to compute primary vertex candidate by candidate using 
  //    other tracks in the event (for pp, broad interaction region);
  //    it is time-consuming procedure, so it can be done after a 
  //    preselection on invariant mass
  //analysis->SetVertexOnTheFly();
  //analysis->SetMassCut(0.1); // GeV
  //--- set single-track preselections
  analysis->SetPtCut(0.); // GeV
  analysis->Setd0Cut(0.); // micron
  //--- set cuts on candidates to be written to file
  //    (see AliBtoJPSItoEleAnalysis.h for a description and for the defaults)
  //analysis->SetBCuts(0.1,1000.,1.1,0.,0.,10000.,10000.,0.,.5);
  analysis->SetBCuts();

  //--- check the settings
  analysis->PrintStatus();

  analysis->FindCandidates(evFirst,evLast,"AliJPSItoEle_all.root");
  delete analysis;

  return;
}
//==========================================================================
void AliBtoJPSItoEleSele() {  

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libAOD.so");
  gSystem->Load("libPWG3base.so");
  gSystem->Load("libPWG3.so");
  gSystem->Load("libPWG3vertexingOld.so");

  //========================  S E L E C T I O N ============================

  AliBtoJPSItoEleAnalysis *analysis = new AliBtoJPSItoEleAnalysis();
  analysis->SetSimulation();
  //analysis->SetOnlySignal();
  analysis->SetOnlyPrimaryJpsi();
  analysis->SetBCuts(0.1,1000.,1.1,0.,0.,10000.,10000.,0.,.5);
  analysis->ApplySelection("AliBtoJPSItoEle.root","AliBtoJPSItoEle_sele.root");
  delete analysis;

  return;
}
//==========================================================================
 





