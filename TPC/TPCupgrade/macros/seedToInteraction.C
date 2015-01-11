void seedToInteraction(Bool_t bPlotOnly = kTRUE, Int_t maxPoints = -1, TString infilename = "$ALICE_ROOT/TPC/Upgrade/macros/toyMC.root"){
  //
  // dummy tracking from distorted space points
  //

  // ==================================================================================================
  // Load libraries/OCDB
  gROOT->ProcessLine(".x $ALICE_ROOT/TPC/Upgrade/macros/loadlibs.C");
  gROOT->ProcessLine(".x $ALICE_ROOT/TPC/Upgrade/macros/ConfigOCDB.C");

  // ==================================================================================================
  // New processing or only plotting
  if(!bPlotOnly){
    
    // ==================================================================================================
    // OCDB
    // TString ocdb="local://$ALICE_ROOT/OCDB/";
    // AliCDBManager::Instance()->SetDefaultStorage(ocdb);
    // AliCDBManager::Instance()->SetRun(0);   
    // TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG));  
   
    // ==================================================================================================
    // input file
    TFile * f = TFile::Open(infilename.Data(),"read");
    TTree * tree = (TTree*)f->Get("toyMCtree"); // AliToyMCEvent tree
    tree->SetMarkerStyle(25);
    tree->SetMarkerSize(0.4);
    
    tree->ls();
    tree->Print();
    
    TBranch *branchEvents = tree->GetBranch("event"); // AliToyMCEvent events
    AliToyMCEvent *events = new AliToyMCEvent();
    branchEvents->SetAddress(&events);
    
    // ==================================================================================================
    // output files
    TTreeSRedirector  * pcstream = new TTreeSRedirector("seedToInteraction.root","recreate");

    TFile *fOut     = TFile::Open("seedToInteraction_histos.root","recreate");
    TH1F  *hT0             = new TH1F("hT0","hT0",1000,0.0,0.0002);  
    TH1F  *hT0Ideal        = new TH1F("hT0Ideal","hT0Ideal",1000,0.0,0.0002);  
    TH1F  *hT0Dist         = new TH1F("hT0Dist","hT0Dist",1000,0.0,0.0002);  
    TH1F  *hT0DistCorr     = new TH1F("hT0DistCorr","hT0DistCorr",1000,0.0,0.0002);  
    TH1F  *hT0DistCorrFull = new TH1F("hT0DistCorrFull","hT0DistCorrFull",1000,0.0,0.0002);  
    TH1F  *hT0DiffIdeal        = new TH1F("hT0DiffIdeal","hT0DiffIdeal",1000,-0.0001,0.0001);  
    TH1F  *hT0DiffDist         = new TH1F("hT0DiffDist","hT0DiffDist",1000,-0.0001,0.0001);  
    TH1F  *hT0DiffDistCorr     = new TH1F("hT0DiffDistCorr","hT0DiffDistCorr",1000,-0.0001,0.0001);  
    TH1F  *hT0DiffDistCorrFull = new TH1F("hT0DiffDistCorrFull","hT0DiffDistCorrFull",1000,-0.0001,0.0001);  
    
    
    // ==================================================================================================
    //
    // Loop over all events
    //
    Int_t nEvents = tree->GetEntriesFast();
    //nEvents=1;//for testing only
    cout<<"Process "<<nEvents<<" events..."<<endl;
    
    Int_t nTracks = 0; // number of tracks per event
    Int_t nSP     = 0; // number of space points per track
    Int_t nDSP    = 0; // number of distorted space points per track
    Double_t t0   = 0.;// real time 0
    Double_t tDrift = 0.0000943;// full dift time = 94.3 us
    
    AliToyMCTrack         *mcTrack              = NULL;
    AliExternalTrackParam *recTrack             = NULL;
    AliExternalTrackParam *recTrackDist         = NULL;
    AliExternalTrackParam *recTrackDistCorr     = NULL;
    AliExternalTrackParam *recTrackDistCorrFull = NULL;
    
    for(Int_t iEvent = 0; iEvent < nEvents; iEvent++){
      tree->GetEntry(iEvent);

      // ==================================================================================================
      // Fill real time 0
      t0 = events->GetT0();      

      // ==================================================================================================
      // loop over all MC tracks
      nTracks = events->GetNumberOfTracks();
      for(Int_t iTrack = 0; iTrack < nTracks; iTrack++){
	
      // ==================================================================================================
      // retrieve the MC track
	mcTrack = events->GetTrack(iTrack); 
	nSP     = mcTrack->GetNumberOfSpacePoints();
	nDSP    = mcTrack->GetNumberOfDistSpacePoints();
	
	AliTrackPointArray pointArray(nSP);
	AliTrackPointArray pointArrayDist(nDSP);
	
	for(Int_t iSP = 0; iSP < nSP; iSP++){
	  AliTPCclusterMI *point = (AliTPCclusterMI *)(mcTrack->GetSpacePoint(iSP))->Clone();
	  point->SetZ(point->GetTimeBin());
	  pointArray.AddPoint(iSP, point->MakePoint() );     
	}
	
	for(Int_t iDSP = 0; iDSP < nDSP; iDSP++){	
	  AliTPCclusterMI *pointDist = (AliTPCclusterMI *)(mcTrack->GetDistortedSpacePoint(iDSP))->Clone();
	  pointDist->SetZ(pointDist->GetTimeBin());
	  pointArrayDist.AddPoint(iDSP, pointDist->MakePoint() );      
	}

	// for(Int_t iDSP = 0; iDSP < nDSP; iDSP++){	
	//   AliTPCclusterMI *pointDistCorr = (AliTPCclusterMI *)(mcTrack->GetDistortedSpacePoint(iDSP))->Clone();
	//   pointDistCorr->SetZ(pointDistCorr->GetTimeBin());
	//   pointArrayDistCorr.AddPoint(iDSP, pointDistCorr->MakePoint(pointDistCorr) );      
	// }

	// for(Int_t iDSP = 0; iDSP < nDSP; iDSP++){	
	//   AliTPCclusterMI *pointDistCorrFull = (AliTPCclusterMI *)(mcTrack->GetDistortedSpacePoint(iDSP))->Clone();
	//   pointDistCorrFull->SetZ(pointDistCorrFull->GetTimeBin());
	//   pointArrayDistCorrFull.AddPoint(iDSP, pointDistCorrFull->MakePoint(pointDistCorrFull) );      
	// }
	
	// ==================================================================================================
	// reconstruct
	recTrack             = SeedToVertex(&pointArray, maxPoints);
	recTrackDist         = SeedToVertex(&pointArrayDist, maxPoints);
	// recTrackDistCorr     = SeedToVertex(&pointArrayDistCorr);
	// recTrackDistCorrFull = SeedToVertex(&pointArrayDistCorrFull);

	//cout<<"ALPHA = "<< mcTrack->GetAlpha()<<" "<<recTrack->GetAlpha()<<" "<<recTrackDist->GetAlpha()<<endl;
	//cout<<mcTrack->Phi()<<" "<<recTrack->Phi()<<" "<<recTrackDist->Phi()<<endl;
  

	// // ==================================================================================================
	// // debug
	// if(recTrack->GetY()>10 || recTrack->GetY() < -1){
	//   cout<<mcTrack->Phi()<<" "<<recTrack->Phi()<<" "<<recTrackDist->Phi()<<endl;
	//   cout<<"MC    = "<<mcTrack->GetY()<<" "<<mcTrack->GetZ()<<" "<<mcTrack->GetSnp()<<" "<<mcTrack->GetTgl()<<" "<<mcTrack->GetSigned1Pt()<<endl;
	//   cout<<"Ideal = "<<recTrack->GetY()<<" "<<recTrack->GetZ()<<" "<<recTrack->GetSnp()<<" "<<recTrack->GetTgl()<<" "<<recTrack->GetSigned1Pt()<<endl;
	//   cout<<"Dist  = "<<recTrackDist->GetY()<<" "<<recTrackDist->GetZ()<<" "<<recTrackDist->GetSnp()<<" "<<recTrackDist->GetTgl()<<" "<<recTrackDist->GetSigned1Pt()<<endl;
	// }
	
	// ==================================================================================================
	// output: debug streamer
	if (pcstream) (*pcstream)<<Form("seedToInt")<<
			"trackIn.="<<mcTrack<<          //  MC track
			"trackOutIdeal.="<<recTrack<<        //  rec track from ideal space points
			"trackOut.="<<recTrackDist<<    //  rec track from dist space points
			"\n"; 

	// ==================================================================================================
	// output: histograms
	hT0->Fill(t0);
	// only for good tracks
	if(recTrack->GetY() > -50 && recTrack->GetY() < 50){
	  hT0Ideal->Fill(recTrack->GetZ());
	  hT0DiffIdeal->Fill(recTrack->GetZ()-t0-tDrift);
	}
	if(recTrackDist->GetY() > -50 && recTrackDist->GetY() < 50){
	  hT0Dist->Fill(recTrackDist->GetZ());
	  hT0DiffDist->Fill(recTrackDist->GetZ()-t0-tDrift);
	}
      }
    }

    // ==================================================================================================
    // write histograms
    fOut->cd();
    hT0->Write();
    hT0Ideal->Write();
    hT0Dist->Write();
    hT0DistCorr->Write();
    hT0DistCorrFull->Write();
    hT0DiffIdeal->Write();
    hT0DiffDist->Write();
    hT0DiffDistCorr->Write();
    hT0DiffDistCorrFull->Write();
    fOut->Close();
    
    delete pcstream;
    
  }
  
  // ==================================================================================================
  //
  // Visualize 
  //
  TFile * f = TFile::Open("seedToInteraction.root");
  TTree * tree = (TTree*)f->Get("seedToInt");
  tree->SetMarkerStyle(25);
  tree->SetMarkerSize(0.4);

  // Event with space points
  TCanvas *canvasEvent = new TCanvas("canvasEvent","canvasEvent",1200,600);
  canvasEvent->Divide(2,1);

  canvasEvent->cd(1);
  tree->SetMarkerSize(0.25);
  tree->SetMarkerColor(1);
  tree->Draw("trackIn.fSpacePoints.fY:trackIn.fSpacePoints.fX","trackIn.fSpacePoints.fX!=0","",50,0);
  tree->SetMarkerColor(2);
  tree->Draw("trackIn.fDistortedSpacePoints.fY:trackIn.fDistortedSpacePoints.fX","trackIn.fDistortedSpacePoints.fX!=0","same",50,0);
  htemp = (TH2F*)gPad->GetPrimitive("htemp");
  htemp->GetXaxis()->SetTitle("x (cm)");
  htemp->GetYaxis()->SetTitle("y (cm)");
  htemp->SetTitle("Event (x-y plane)");
  canvasEvent->cd(1)->Update();

  canvasEvent->cd(2);
  tree->SetMarkerSize(0.25);
  tree->SetMarkerColor(1);
  tree->Draw("trackIn.fSpacePoints.fY:trackIn.fSpacePoints.fZ","trackIn.fSpacePoints.fX!=0","",50,0);
  tree->SetMarkerColor(2);
  tree->Draw("trackIn.fDistortedSpacePoints.fY:trackIn.fDistortedSpacePoints.fZ","trackIn.fDistortedSpacePoints.fX!=0","same",50,0);
  htemp = (TH2F*)gPad->GetPrimitive("htemp");
  htemp->GetXaxis()->SetTitle("z (cm)");
  htemp->GetYaxis()->SetTitle("y (cm)");
  htemp->SetTitle("Event (z-y plane)");
  canvasEvent->cd(2)->Update();

  canvasEvent->SaveAs("seedToInteraction_Points.eps");

  // track (5 parameters)
  TCanvas *canvasTrack = new TCanvas("canvasTrack","canvasTrack",1200,900);
  canvasTrack->Divide(3,2);

  for(Int_t i = 0; i < 5; i++){

    canvasTrack->cd(i+1);
    tree->SetLineColor(1);
    tree->Draw(Form("trackOutIdeal.fP[%d]",i),"","");
    tree->SetLineColor(4);
    tree->Draw(Form("trackIn.AliExternalTrackParam.fP[%d]",i),"","same");
    tree->SetLineColor(2);
    tree->Draw(Form("trackOut.fP[%d]",i),"","same");
    
  }

  
  canvasTrack->cd(6);
  TLatex *l = new TLatex;
  l->SetTextColor(4);
  l->DrawLatex(0.3,0.5,"Input MC track");
  l->SetTextColor(1);
  l->DrawLatex(0.3,0.4,"rec track (ideal)");
  l->SetTextColor(2);
  l->DrawLatex(0.3,0.3,"rec track (dist)");

  canvasTrack->SaveAs("seedToInteraction_Tracks.eps");

  f->Close();


  f = TFile::Open("seedToInteraction_histos.root");
  TH1F  *hT0             = (TH1F*)f->Get("hT0");
  TH1F  *hT0Ideal        = (TH1F*)f->Get("hT0Ideal");
  TH1F  *hT0Dist         = (TH1F*)f->Get("hT0Dist");
  TH1F  *hT0DistCorr     = (TH1F*)f->Get("hT0DistCorr");
  TH1F  *hT0DistCorrFull = (TH1F*)f->Get("hT0DistCorrFull");
  TH1F  *hT0DiffIdeal        = (TH1F*)f->Get("hT0DiffIdeal");
  TH1F  *hT0DiffDist         = (TH1F*)f->Get("hT0DiffDist");
  TH1F  *hT0DiffDistCorr     = (TH1F*)f->Get("hT0DiffDistCorr");
  TH1F  *hT0DiffDistCorrFull = (TH1F*)f->Get("hT0DiffDistCorrFull");

  // track (5 parameters)
  TCanvas *canvasT0 = new TCanvas("canvasT0","canvasT0",1200,900);
  canvasT0->Divide(2,1);
 
  canvasT0->cd(1);
  hT0->SetMarkerColor(4);
  hT0->SetLineColor(4);
  hT0->DrawCopy();
  hT0Ideal->SetMarkerColor(1);
  hT0Ideal->SetLineColor(1);
  hT0Ideal->DrawCopy("same");
  hT0Dist->SetMarkerColor(2);
  hT0Dist->SetLineColor(2);
  hT0Dist->DrawCopy("same");
  hT0DistCorr->SetMarkerColor(2);
  hT0DistCorr->SetLineColor(2);
  hT0DistCorr->DrawCopy("same");
  hT0DistCorrFull->SetMarkerColor(2);
  hT0DistCorrFull->SetLineColor(2);
  hT0DistCorrFull->DrawCopy("same");

  canvasT0->cd(2);
  hT0DiffIdeal->SetMarkerColor(1);
  hT0DiffIdeal->SetLineColor(1);
  hT0DiffIdeal->DrawCopy();
  hT0DiffDist->SetMarkerColor(2);
  hT0DiffDist->SetLineColor(2);
  hT0DiffDist->DrawCopy("same");
  hT0DiffDistCorr->SetMarkerColor(2);
  hT0DiffDistCorr->SetLineColor(2);
  hT0DiffDistCorr->DrawCopy("same");
  hT0DiffDistCorrFull->SetMarkerColor(2);
  hT0DiffDistCorrFull->SetLineColor(2);
  hT0DiffDistCorrFull->DrawCopy("same");

  if(maxPoints > 0){
    canvasT0->SaveAs(Form("seedToInteraction_T0_%dSP.eps",maxPoints));
  }
  else{
    canvasT0->SaveAs(Form("seedToInteraction_T0_allSP.eps"));
  }

  f->Close();

}

AliExternalTrackParam * SeedToVertex(AliTrackPointArray *pointArray,  Int_t maxPoints = -1){
  //
  // Seeding and propagation to vertex region
  // Fixed parameters
  // dir      - direction - out=1 or in=-1
  // refX     - reference X to fit the track (needed?)
  //

  // Parameters
  Int_t dir = -1;
  Double_t refX = 0.;

  AliExternalTrackParam *track = 0x0;
  AliTrackPoint point1,point2,point3;

  Int_t ipoints[3] = {0,0,0}; // this should be a settable array later
  Int_t npoints = pointArray->GetNPoints();
  if(maxPoints > 2) npoints = maxPoints;

  // need to check if outermost points are filled
  if (dir==-1) {  // go from outermost point
    for(Int_t iPoint = 0; iPoint < npoints; iPoint++){
      pointArray->GetPoint(point1,npoints-iPoint-1);
      if(point1.GetX()!=0 && point1.GetY()!=0 && point1.GetZ()!=0){
	npoints = npoints-iPoint;
	break;
      }
    }
  }
  //
  // parameters (from AliTPCCorrection)
  //
  AliTPCROC * roc = AliTPCROC::Instance();
  const Int_t    npoints0=roc->GetNRows(0)+roc->GetNRows(36);
  const Double_t kRTPC0  =roc->GetPadRowRadii(0,0);
  const Double_t kRTPC1  =roc->GetPadRowRadii(36,roc->GetNRows(36)-1);
  const Double_t kMaxSnp = 0.85;  
  const Double_t kSigmaY=0.1;
  const Double_t kSigmaZ=0.1;
  const Double_t kMaxR=500;
  const Double_t kMaxZ=500;
  
  const Double_t kMaxZ0=220;
  const Double_t kZcut=3;
  const Double_t kMass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();

  if (dir==1) {  //make seed inner
    pointArray->GetPoint(point1,1);
    pointArray->GetPoint(point2,11);
    pointArray->GetPoint(point3,21);
    ipoints[0] = 1;
    ipoints[1] = 11;
    ipoints[2] = 21;
  }
  if (dir==-1){ //make seed outer
    pointArray->GetPoint(point1,npoints-21);
    pointArray->GetPoint(point2,npoints-11);
    pointArray->GetPoint(point3,npoints-1);
    ipoints[0] = npoints-1;
    ipoints[1] = npoints-11;
    ipoints[2] = npoints-21;
  } 
  if ((TMath::Abs(point1.GetX()-point3.GetX())+TMath::Abs(point1.GetY()-point3.GetY()))<1){
    printf("fit points not properly initialized\n");
    return 0;
  }
  track = AliTrackerBase::MakeSeed(point1, point2, point3);

  // ===========================================================================================================
  // track propagation (from AliTPCCorrection::FitDistortedTrack)
  track->ResetCovariance(10);

  for (Int_t jpoint=0; jpoint<npoints; jpoint++){

    Int_t ipoint = npoints-jpoint-1;//ipoints[jpoint]; 
    if(maxPoints==3) ipoint = ipoints[jpoint];
    //
    AliTrackPoint pIn;
    Double_t xyz[3];
    pointArray->GetPoint(pIn,ipoint);
    AliTrackPoint prot = pIn.Rotate(track->GetAlpha());   // rotate to the local frame - non distoted  point
    if (TMath::Abs(prot.GetX())<kRTPC0) continue;
    if (TMath::Abs(prot.GetX())>kRTPC1) continue;
    //
    if (!AliTrackerBase::PropagateTrackTo(track,prot.GetX(),kMass,5,kFALSE,kMaxSnp)) break;
    if (TMath::Abs(track->GetZ())>kMaxZ) break;
    if (TMath::Abs(track->GetX())>kMaxR) break;
    if (dir>0 && track->GetX()>refX) continue;
    if (dir<0 && track->GetX()<refX) continue;
    if (TMath::Abs(track->GetZ())<kZcut)continue;
    track->GetXYZ(xyz);  //track also propagated to the same reference radius
    //
    Double_t pointPos[2]={0,0};
    Double_t pointCov[3]={0,0,0};
    pointPos[0]=prot.GetY();//local y
    pointPos[1]=prot.GetZ();//local z
    pointCov[0]=prot.GetCov()[3];//simay^2
    pointCov[1]=prot.GetCov()[4];//sigmayz
    pointCov[2]=prot.GetCov()[5];//sigmaz^2
    if (!track->Update(pointPos,pointCov)) break;

  }

  AliTrackerBase::PropagateTrackTo(track,refX,kMass,5.,kTRUE,kMaxSnp);
  AliTrackerBase::PropagateTrackTo(track,refX,kMass,1.,kTRUE,kMaxSnp);
  // ===========================================================================================================


  return track;

}
