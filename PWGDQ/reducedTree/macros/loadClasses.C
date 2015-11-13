void loadClasses(TString dir="../"){

  //gSystem->Load("libFlowVector.so");

  //gROOT->ProcessLine(".include "+dir+"AliQnCorrections/");
  gROOT->ProcessLine(".include "+dir);

  /*gROOT->ProcessLine(".L "+dir+"AliQnCorrections/AliQnCorrectionsAxes.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliQnCorrections/AliQnCorrectionsQnVector.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliQnCorrections/AliQnCorrectionsDataVector.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliQnCorrections/AliQnCorrectionsCuts.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliQnCorrections/AliQnCorrectionsConfiguration.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliQnCorrections/AliQnCorrectionsHistograms.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliQnCorrections/AliQnCorrectionsSteps.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliQnCorrections/AliQnCorrectionsManager.cxx+"); */
  gROOT->ProcessLine(".L "+dir+"AliReducedBaseEvent.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliReducedBaseTrack.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliReducedCaloClusterInfo.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliReducedFMDInfo.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliReducedEventInfo.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliReducedEventPlaneInfo.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliReducedPairInfo.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliReducedTrackInfo.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliReducedVarManager.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliHistogramManager.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliMixingHandler.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliReducedInfoCut.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliReducedBaseTrackCut.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliReducedEventCut.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliReducedAnalysisTaskSE.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliReducedAnalysisTest.cxx+");
  //gROOT->ProcessLine(".L "+dir+"AliReducedQnFillEvent.cxx+");
  //gROOT->ProcessLine(".L "+dir+"AliReducedAnalysisQnCorrections.cxx+");
  gROOT->ProcessLine(".L "+dir+"AliAnalysisTaskReducedTreeMaker.cxx+");
  //gROOT->ProcessLine(".L AliAnalysisTaskReducedEventProcessor.cxx+");
  
  return;
  
  AliReducedVarManager::SetDefaultVarNames();
  
  AliHistogramManager histMan("Histogram manager",AliReducedVarManager::kNVars);
  histMan.SetUseDefaultVariableNames(kTRUE);
  histMan.SetDefaultVarNames(AliReducedVarManager::fgVariableNames,AliReducedVarManager::fgVariableUnits);
  
  AliMixingHandler mh;
  mh.SetHistogramManager(&histMan);
  mh.SetPoolDepth(5);
  mh.SetMixingThreshold(1.0);
  mh.SetDownscaleEvents(1);
  mh.SetDownscaleTracks(1);
  mh.SetNParallelCuts(5);
  mh.SetEventVariables(AliReducedVarManager::kCentVZERO,AliReducedVarManager::kVtxZ,AliReducedVarManager::kVZERORP);
  
  histMan.AddHistClass("Track");
  cout << "1" << endl;
  histMan.AddHistogram("Track","pt","Track pt", kFALSE, 100, 0.0, 50.0, AliReducedVarManager::kPt);
  cout << "2" << endl;
  histMan.AddHistogram("Track","phi","Track phi", kFALSE, 100, -1.7, 1.7, AliReducedVarManager::kPhi);
  cout << "3" << endl;
  histMan.AddHistogram("Track","eta","Track eta", kFALSE, 100, -5.0, 5.0, AliReducedVarManager::kEta);
  cout << "4" << endl;
  
  TString histNames;
  for(Int_t i=0;i<mh.GetNParallelCuts();++i) {
    TString histClass = Form("Cut #%d", i);
    histNames += histClass.Data(); histNames += ";";
    histMan.AddHistClass(histClass.Data());
    histMan.AddHistogram(histClass.Data(),"minv","Invariant mass", kFALSE, 100, 0.0, 50.0, AliReducedVarManager::kMass);
    histMan.AddHistogram(histClass.Data(),"pt","Pair pt", kFALSE, 100, 0.0, 50.0, AliReducedVarManager::kPt);
    histMan.AddHistogram(histClass.Data(),"phi","Pair phi", kFALSE, 100, -1.7, 1.7, AliReducedVarManager::kPhi);
    histMan.AddHistogram(histClass.Data(),"eta","Pair eta", kFALSE, 100, -5.0, 5.0, AliReducedVarManager::kEta);
  }
  mh.SetHistClassNames(histNames.Data());
  
  AliReducedVarManager::SetUseVars(histMan.GetUsedVars());  
  
  //Float_t centLims[10] = {0.0,5.0,10.0,15.0,20.0,30.0,40.0,50.0,70.0,90.};
  Float_t centLims[2] = {0.0,100.};
  //Float_t zLims[11] = {-10.,-8.,-6.,-4.,-2.,0.,2.,4.,6.,8.,10.};
  Float_t zLims[2] = {-10.,10.};
  Float_t epLims[2] = {-1.5,1.5};
  mh.SetCentralityLimits(2,centLims);
  mh.SetEventVertexLimits(2,zLims);
  mh.SetEventPlaneLimits(2,epLims);
  
  mh.PrintMixingLists(0);
  
  TList* pos=new TList(); pos->SetOwner(kTRUE);
  TList* neg=new TList(); neg->SetOwner(kTRUE);
  TTimeStamp now;
  gRandom->SetSeed(now.GetNanoSec());
  Float_t values[AliReducedVarManager::kNVars];
  AliReducedBaseTrack* track=0x0;
  for(Int_t iev=0; iev<50; ++iev) {
    if(iev%1 == 0)
      cout << "***************************************************** Generating event #" << iev << endl;
    pos->Clear("C"); neg->Clear("C");
    Int_t ntracks = gRandom->Poisson(3.);
    values[AliReducedVarManager::kCentVZERO] = gRandom->Rndm()*100.0;
    values[AliReducedVarManager::kVtxZ] = gRandom->Gaus(0.,3.);
    values[AliReducedVarManager::kVZERORP] = 3.0*(gRandom->Rndm()-0.5);
    cout << "ntracks/cent/vtx/ep :: " << ntracks << " / " << values[AliReducedVarManager::kCentVZERO] 
         << " / " << values[AliReducedVarManager::kVtxZ] << " / " << values[AliReducedVarManager::kVZERORP] << endl;
    
    for(Int_t it=0;it<ntracks;++it) {
      track=new AliReducedBaseTrack();
      track->PxPyPz(gRandom->Landau(1.,1.), gRandom->Landau(1.,1.), gRandom->Landau(1.0,0.3));
      track->Charge((gRandom->Rndm()<0.5 ? -1 : 1));
      while(!track->GetFlags()) {
	for(UShort_t iflag=0; iflag<mh.GetNParallelCuts();++iflag)
	  if(gRandom->Rndm()<0.3) track->SetFlag(iflag);
      }
            
      AliReducedVarManager::FillBaseTrackInfo(track,values);
      histMan.FillHistClass("Track",values);
      
      cout << "track #" << it << " (pt/phi/eta/charge/flags) :: "
           << track->Pt() << "/" << track->Phi() << "/" << track->Eta() << "/" << track->Charge() << "/" << flush;
      AliReducedVarManager::PrintBits(track->GetFlags(),mh.GetNParallelCuts()); cout << endl;
      
      if(track.Charge()==-1) neg->Add(track->Clone());
      else pos->Add(track->Clone());
      delete track;
    }
    
    mh.FillEvent(pos,neg,values,1);  
    mh.PrintMixingLists(3);
  }  // end loop over events
  
  mh.RunLeftoverMixing(1);
  mh.PrintMixingLists(3);
  
  TFile* save=TFile::Open("test.root","recreate");
  histMan.WriteOutput(save);
}
