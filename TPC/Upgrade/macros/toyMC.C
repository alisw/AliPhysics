void toyMC(Bool_t bPlotOnly = kTRUE){
  //
  // dummy tracking from distorted space points
  //

  // ==================================================================================================
  // Load libraries
  gROOT->ProcessLine(".x loadlibs.C");

  // ==================================================================================================
  // New processing or only plotting
  if(!bPlotOnly){
    
    // ==================================================================================================
    // OCDB
    TString ocdb="local://$ALICE_ROOT/OCDB/";
    AliCDBManager::Instance()->SetDefaultStorage(ocdb);
    AliCDBManager::Instance()->SetRun(0);   
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG));  
    
    // ==================================================================================================
    // input file
    TFile * f = TFile::Open("/Users/physics/ALICE/TPCupgrade/ToyModel/toyMC_Martin/toymceventsten.root");
    TTree * tree = (TTree*)f->Get("ToyMC");
    tree->SetMarkerStyle(25);
    tree->SetMarkerSize(0.4);
    
    tree->ls();
    tree->Print();
    
    TBranch *branchEvents = tree->GetBranch("ToyMCEvents"); // ToyMCEvents
    ToyMCEvent *events = new ToyMCEvent();
    branchEvents->SetAddress(&events);
    
    // ==================================================================================================
    // output file
    TTreeSRedirector  * pcstream = new TTreeSRedirector("dummyTracking.root","recreate");
    
    
    // ==================================================================================================
    //
    // Loop over all events
    //
    Int_t nEvents = tree->GetEntriesFast();
    cout<<"Process "<<nEvents<<" events..."<<endl;
    
    Int_t nTracks = 0; // number of tracks per event
    Int_t nSP     = 0; // number of space points per track
    Int_t nDSP    = 0; // number of distorted space points per track
    
    ToyMCTrack            *mcTrack      = NULL;
    AliExternalTrackParam *recTrack     = NULL;
    AliExternalTrackParam *recTrackDist = NULL;
    
    for(Int_t iEvent = 0; iEvent < nEvents; iEvent++){
      tree->GetEntry(iEvent);
      cout<<iEvent<<" "<<events->GetEventNumber()<<endl;
      
      // ==================================================================================================
      // loop over all MC tracks
      nTracks = events->GetNumberOfTracks();
      for(Int_t iTrack = 0; iTrack < nTracks; iTrack++){
	
      // ==================================================================================================
      // retrieve the MC track
	mcTrack = events->GetTrack(iTrack); 
	nSP     = mcTrack->GetNumberOfSpacePoints();
	nDSP    = mcTrack->GetNumberOfDistSpacePoints();
	//cout<<iTrack<<" with "<<mcTrack->GetNumberOfSpacePoints()<<" and "<<mcTrack->GetNumberOfDistSpacePoints()<<endl;
	
	AliTrackPointArray pointArray(nSP);
	AliTrackPointArray pointArrayDist(nDSP);
	
	for(Int_t iSP = 0; iSP < nSP; iSP++){
	  AliTPCclusterMI *point = mcTrack->GetSpacePoint(iSP);
	  pointArray.AddPoint(iSP, point->MakePoint(point) );      
	}
	
	for(Int_t iDSP = 0; iDSP < nDSP; iDSP++){	
	  AliTPCclusterMI *point = mcTrack->GetDistortedSpacePoint(iDSP);
	  pointArrayDist.AddPoint(iDSP, point->MakePoint(point) );      
	}
	
	// ==================================================================================================
	// reconstruct
	recTrack     = GetDistortedTrack(&pointArray,(AliExternalTrackParam*)mcTrack);
	recTrackDist = GetDistortedTrack(&pointArrayDist,(AliExternalTrackParam*)mcTrack);
	
	// ==================================================================================================
	// output
	if (pcstream) (*pcstream)<<Form("fitDistort")<<
			"trackIn.="<<mcTrack<<          //  MC track
			"trackOutIdeal.="<<recTrack<<        //  rec track from ideal space points
			"trackOut.="<<recTrackDist<<    //  rec track from dist space points
			"\n"; 
      }
    }
    
    delete pcstream;
    
  }
  
  // ==================================================================================================
  //
  // Visualize 
  //
  TFile * f = TFile::Open("dummyTracking.root");
  TTree * tree = (TTree*)f->Get("fitDistort");
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

  canvasEvent->SaveAs("Points.eps");

  // track (5 parameters)
  TCanvas *canvasTrack = new TCanvas("canvasTrack","canvasTrack",1200,900);
  canvasTrack->Divide(3,2);

  for(Int_t i = 0; i < 5; i++){

    canvasTrack->cd(i+1);
    tree->SetLineColor(4);
    tree->Draw(Form("trackIn.AliExternalTrackParam.fP[%d]",i),"","");
    tree->SetLineColor(1);
    tree->Draw(Form("trackOutIdeal.fP[%d]",i),"","same");
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



  canvasTrack->SaveAs("Tracks.eps");


  f->Close();

}

AliExternalTrackParam * GetDistortedTrack(AliTrackPointArray *pointArray, AliExternalTrackParam *trackIn, Int_t maxPoints = -1, Int_t dir = -1, Double_t refX = 1.){
  //
  // track reconstruction (as in AliTPCCorrection for the beginning)
  // dir      - direction - out=1 or in=-1
  // refX     - reference X to fit the track
  //

  AliExternalTrackParam *track = 0x0;
  AliTrackPoint point1,point2,point3;
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
  }
  if (dir==-1){ //make seed outer
    pointArray->GetPoint(point1,npoints-21);
    pointArray->GetPoint(point2,npoints-11);
    pointArray->GetPoint(point3,npoints-1);
  } 
  if ((TMath::Abs(point1.GetX()-point3.GetX())+TMath::Abs(point1.GetY()-point3.GetY()))<1){
    printf("fit points not properly initialized\n");
    return 0;
  }
  track = AliTrackerBase::MakeSeed(point1, point2, point3);
  track->ResetCovariance(10);
  
  if (TMath::Abs(AliTrackerBase::GetBz())<0.01){
     ((Double_t*)track->GetParameter())[4] = ((Double_t*)trackIn->GetParameter())[4];    
  }
  
  for (Int_t jpoint=0; jpoint<npoints; jpoint++){

    Int_t ipoint= (dir>0) ? jpoint: npoints-1-jpoint;
    
    AliTrackPoint pIn;
    pointArray->GetPoint(pIn,ipoint);

    AliTrackPoint prot = pIn.Rotate(track->GetAlpha());   // rotate to the local frame - distorted  point

    //if(TMath::Abs(trackIn->GetTgl()) > 0.8 ) cout<<jpoint<<" "<<ipoint<<": "<<kRTPC0<<" < "<<prot.GetX()<<" < "<<kRTPC1<<" ?    "<<track->GetZ()<<" "<<track->GetX()<<" "<<track->GetSigned1Pt()<<endl;

    if (TMath::Abs(prot.GetX())<kRTPC0) continue;
    if (TMath::Abs(prot.GetX())>kRTPC1) continue;
    //
    if (!AliTrackerBase::PropagateTrackTo(track,prot.GetX(),kMass,5,kFALSE,kMaxSnp)) break;
    if (TMath::Abs(track->GetZ())>kMaxZ) break;
    if (TMath::Abs(track->GetX())>kMaxR) break;
    if (dir>0 && track->GetX()>refX) continue;
    if (dir<0 && track->GetX()<refX) continue;
    if (TMath::Abs(track->GetZ())<kZcut)continue;
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

  track->Rotate(trackIn->GetAlpha());
  AliTrackerBase::PropagateTrackTo(track,refX,kMass,5.,kTRUE,kMaxSnp);
  AliTrackerBase::PropagateTrackTo(track,refX,kMass,1.,kTRUE,kMaxSnp);

  return track;

}
