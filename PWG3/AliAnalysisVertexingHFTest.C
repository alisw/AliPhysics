//--------------------------------------------------------------------------
// Test macro for reconstruction of heavy-flavour vertexing candidates
//
// mode: 0 -> input ESD, output AOD
//       1 -> input ESD, output TTree
//
//     Andrea Dainese, andrea.dainese@lnl.infn.it
//--------------------------------------------------------------------------

void AliAnalysisVertexingHFTest(Int_t mode=0)
{
  
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libPWG3base.so");

  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
  //vHF->SetUseTRef();
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

  switch(mode) {
  case 0:
    vHF->FindCandidatesESDtoAOD();
    break;
  case 1:
    TTree *trees = new TTree[1];
    AliAODRecoDecayHF2Prong *rd2=0;
    // AliAODRecoDecayHF3Prong *rd3=0;
    //AliAODRecoDecayHF4Prong *rd4=0;
    trees[0].Branch("D0toKpi","AliAODRecoDecayHF2Prong",&rd2);
    //trees[1].Branch("JPSItoEle","AliAODRecoDecayHF2Prong",&rd2);
    //trees[2].Branch("Charmto3Prong","AliAODRecoDecayHF3Prong",&rd3);
    //trees[3].Branch("D0to4Prong","AliAODRecoDecayHF4Prong",&rd4);
    
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
    TFile *fout = new TFile("AliAnalysisVertexingHF.root","recreate");
    trees[0].Write("TreeD0toKpi");
    //trees[1].Write("TreeJPSItoEle");
    //trees[2].Write("TreeCharm3Prong");
    //trees[3].Write("TreeD0to4Prong");
    fout->Close();
    delete fout;
    break;
  }

  return;
}
