//--------------------------------------------------------------------------
// Test macro for reconstruction of heavy-flavour vertexing candidates
//
//     Andrea Dainese, andrea.dainese@lnl.infn.it
//--------------------------------------------------------------------------

void AliAnalysisVertexingHFTest()
{


  printf("THIS MACRO IS OBSOLETE, PLEASE USE AliAnalysisTaskSEVertexingHFTest.C");
  return;

  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libPWG3base.so");
  gSystem->Load("libPWG3vertexingHF.so");

  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
  //--- switch-off candidates finding (default: all on)
  //vHF->SetD0toKpiOff();
  vHF->SetJPSItoEleOff();
  vHF->Set3ProngOff();
  vHF->Set4ProngOff();
  //--- secondary vertex with KF?
  vHF->SetSecVtxWithKF();
  //--- set cuts for single-track selection
  vHF->SetITSrefitRequired();
  vHF->SetBothSPDNotRequired();
  vHF->SetMinITSCls(5);
  vHF->SetMinPtCut(0.);
  vHF->SetMind0Cut(0.);
  //--- set cuts for candidates selection
  //vHF->SetD0toKpiCuts(); 
  //vHF->SetBtoJPSICuts(); 
  //vHF->SetDplusCuts(); 
  //--- set this if you want to reconstruct primary vertex candidate by
  //    candidate using other tracks in the event (for pp, broad 
  //    interaction region)
  //vHF->SetRecoPrimVtxSkippingTrks();
  //--- OR set this if you want to remove the candidate daughters from 
  //    the primary vertex, without recostructing it from scratch
  //vHF->SetRmTrksFromPrimVtx();

  //--- check the settings
  vHF->PrintStatus();
  //--- verbose
  //vHF->SetDebug(1);

  TFile *fout = new TFile("AliAnalysisVertexingHF.root","recreate");

  TTree **trees = 0;
  trees = new TTree*[1];
  AliAODRecoDecayHF2Prong *rd2=0;
  //AliAODRecoDecayHF3Prong *rd3=0;
  //AliAODRecoDecayHF4Prong *rd4=0;
  trees[0] = new TTree("D0toKpi","D0toKpi");
  trees[0]->Branch("D0toKpi","AliAODRecoDecayHF2Prong",&rd2);
  //trees[1] = new TTree("JPSItoEle","JPSItoEle");
  //trees[1]->Branch("JPSItoEle","AliAODRecoDecayHF2Prong",&rd2);
  //trees[2] = new TTree("Charmto3Prong","Charmto3Prong");
  //trees[2]->Branch("Charmto3Prong","AliAODRecoDecayHF3Prong",&rd3);
  //trees[3] = new TTree("Charmto4Prong","Charmto4Prong");
  //trees[3]->Branch("D0to4Prong","AliAODRecoDecayHF4Prong",&rd4);
  
  TFile *inesd = new TFile("AliESDs.root");
  AliESDEvent *esd = new AliESDEvent();
  TTree *esdTree = (TTree*)inesd->Get("esdTree");
  esd->ReadFromTree(esdTree);
  
  for(Int_t i=0; i<esdTree->GetEntries(); i++) {
    esdTree->GetEvent(i);
    vHF->FindCandidates(esd,trees);
    //if(i==evLast) i=evFirst;
  }
    
  delete esdTree;
  inesd->Close();
  delete inesd;
  
  // Write trees with candidates
  fout->cd();
  trees[0]->Write("TreeD0toKpi");
  //trees[1]->Write("TreeJPSItoEle");
  //trees[2]->Write("TreeCharm3Prong");
  //trees[3]->Write("TreeD0to4Prong");
  fout->Close();
  delete fout;

  return;
}
