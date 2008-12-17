/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

////////////////////////////////////////////////////////////////////
//                                                                //
// Fills a set of QA histograms to check the correctness of       //
// the TRD reconstruction                                         // 
//                                                                //
////////////////////////////////////////////////////////////////////

#include "AliTRDtrackingAnalysis.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TGeoMatrix.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"

#include "AliRunLoader.h"
#include "AliTRDgeometry.h"
#include "AliRun.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliTrackReference.h"
#include "AliTracker.h"

#include "AliTRDcluster.h"
#include "AliTRDpadPlane.h"
#include "AliTRDcalibDB.h"
#include "AliTRDtracker.h"
//#include "AliTRDtracklet.h"

#include "TGeoManager.h"

//////////////////////////////////////////////////////////////////////////////////////////

AliTRDtrackingAnalysis::AliTRDtrackingAnalysis():
  TObject(),
  fPath(0),
  fRefTPC(0),
  fRefTRD(0),
  fLoader(0),
  fEsdTree(0),
  fESD(0),
  fTracker(0),
  fDeltaPt(0),
  fDeltaZ(0),
  fDeltaX(0),
  fDeltaYPos(0),
  fDeltaYNeg(0),
  fNPoints(0),
  fNGood(0),
  fRefSpace(0),
  fGeo(0),
  fClY2(0),
  fClY3(0),
  fTgPhi(0),
  fGrResTgPhi(0),
  fGrMeanTgPhi(0),
  fTrklY(0),
  fTrklZ(0),
  fClZ(0),
  fClZZ(0),
  fClYY(0),
  fClYX(0),
  fNLabels(0),
  fTestBits(0),
  fRefDx(0),
  fClZXref(0),
  fClZXcl(0),
  fClPos(0)
{

  fDeltaX = new TH1D("deltaX", ";delta X (cm)", 100, -1, 1);
  fDeltaZ = new TH1D("deltaZ", ";delta Z (cm)", 100, -2, 2);

  fDeltaYPos = new TH1D("deltaYpos", ";delta Y (mm)", 100, -1, 1);
  fDeltaYNeg = new TH1D("deltaYneg", ";delta Y (mm)", 100, -1, 1);
  
  fNPoints = new TH1D("nPoints", ";np", 40, -0.5, 39.5);
  fNGood   = new TH1D("nGood", ";np", 40, -0.5, 39.5);

  fDeltaPt = new TH1D("deltaPt", ";delta Pt/Pt (%)", 100, -10, 10);
  fRefSpace = new TH2D("refSpace", ";y;x", 120, -60, 60, 200, -4, 1);

  fTrklY = new TH1D("trklY", ";delta Y (mm)", 100, -1, 1);
  fTrklZ = new TH1D("trklZ", ";delta Z (cm)", 100, -10, 10);


  // cluster studies
  fClY2 = new TH1D("clY2", ";delta Y (mm)", 100, -10, 10);
  fClY3 = new TH1D("clY3", ";delta Y (mm)", 100, -10, 10);

  for(int i=0; i<12; i++) // bewere hidden constants in the code
    fClYTgPhi[i] = new TH1D(Form("clYtgPhi%d", i), ";delta Y (mm)", 100, -3, 3);

  fTgPhi = new TH1D("tgPhi", ";Tg(#phi)", 100, -0.3, 0.3);
  fGrResTgPhi = new TGraphErrors();
  fGrMeanTgPhi = new TGraphErrors();

  //fPullY2 = new TH1D("pullY2", ";pulls Y", 100, -5, 5);
  //fPullY3 = new TH1D("pullY3", ";pulls Y", 100, -5, 5);


  fClZ = new TH1D("clZ", ";delta Z (cm)", 200, -20, 20);
  fClZZ = new TH2D("clZZ", ";z ref;z cl", 600, -300, 300, 600, -300, 300);  

  fClYY = new TH2D("clYY", ";dY;dY", 100, -3, 3, 100, -3, 3);
  fClYX = new TH2D("clYX", ";Y;X", 250, -60, 60, 100, -4, 1);

  fNLabels = new TH1D("clLabels", ";n labels", 10, -0.5, 9.5);
  fTestBits = new TH1D("bits", ";bits", 10, -0.5, 9.5);
  
  fRefDx = new TH1D("refDX", ";delta X", 100, 0, 20);
  fClPos = new TH2D("clPos", ";z;y", 400, -400, 400, 120, -60, 60);

  fClZXref = new TH2D("clZXref", ";z;x", 36, -54, 54, 300, 280, 380);
  fClZXcl =  new TH2D("clZXcl", ";z;x", 36, -54, 54, 300, 280, 380);

  //fGeo = new AliTRDgeometry();
}

//////////////////////////////////////////////////////////////////////////////////////////

AliTRDtrackingAnalysis::AliTRDtrackingAnalysis(const AliTRDtrackingAnalysis &t):
  TObject(t),
  fPath(0),
  fRefTPC(0),
  fRefTRD(0),
  fLoader(0),
  fEsdTree(0),
  fESD(0),
  fTracker(0),
  fDeltaPt(0),
  fDeltaZ(0),
  fDeltaX(0),
  fDeltaYPos(0),
  fDeltaYNeg(0),
  fNPoints(0),
  fNGood(0),
  fRefSpace(0),
  fGeo(0),
  fClY2(0),
  fClY3(0),
  fTgPhi(0),
  fGrResTgPhi(0),
  fGrMeanTgPhi(0),
  fTrklY(0),
  fTrklZ(0),
  fClZ(0),
  fClZZ(0),
  fClYY(0),
  fClYX(0),
  fNLabels(0),
  fTestBits(0),
  fRefDx(0),
  fClZXref(0),
  fClZXcl(0),
  fClPos(0)
{

  fDeltaX = new TH1D("deltaX", ";delta X (cm)", 100, -1, 1);
  fDeltaZ = new TH1D("deltaZ", ";delta Z (cm)", 100, -2, 2);

  fDeltaYPos = new TH1D("deltaYpos", ";delta Y (mm)", 100, -1, 1);
  fDeltaYNeg = new TH1D("deltaYneg", ";delta Y (mm)", 100, -1, 1);
  
  fNPoints = new TH1D("nPoints", ";np", 40, -0.5, 39.5);
  fNGood   = new TH1D("nGood", ";np", 40, -0.5, 39.5);

  fDeltaPt = new TH1D("deltaPt", ";delta Pt/Pt (%)", 100, -10, 10);
  fRefSpace = new TH2D("refSpace", ";y;x", 120, -60, 60, 200, -4, 1);

  fTrklY = new TH1D("trklY", ";delta Y (mm)", 100, -1, 1);
  fTrklZ = new TH1D("trklZ", ";delta Z (cm)", 100, -10, 10);


  // cluster studies
  fClY2 = new TH1D("clY2", ";delta Y (mm)", 100, -10, 10);
  fClY3 = new TH1D("clY3", ";delta Y (mm)", 100, -10, 10);

  for(int i=0; i<12; i++) // bewere hidden constants in the code
    fClYTgPhi[i] = new TH1D(Form("clYtgPhi%d", i), ";delta Y (mm)", 100, -3, 3);

  fTgPhi = new TH1D("tgPhi", ";Tg(#phi)", 100, -0.3, 0.3);
  fGrResTgPhi = new TGraphErrors();
  fGrMeanTgPhi = new TGraphErrors();

  //fPullY2 = new TH1D("pullY2", ";pulls Y", 100, -5, 5);
  //fPullY3 = new TH1D("pullY3", ";pulls Y", 100, -5, 5);


  fClZ = new TH1D("clZ", ";delta Z (cm)", 200, -20, 20);
  fClZZ = new TH2D("clZZ", ";z ref;z cl", 600, -300, 300, 600, -300, 300);  

  fClYY = new TH2D("clYY", ";dY;dY", 100, -3, 3, 100, -3, 3);
  fClYX = new TH2D("clYX", ";Y;X", 250, -60, 60, 100, -4, 1);

  fNLabels = new TH1D("clLabels", ";n labels", 10, -0.5, 9.5);
  fTestBits = new TH1D("bits", ";bits", 10, -0.5, 9.5);
  
  fRefDx = new TH1D("refDX", ";delta X", 100, 0, 20);
  fClPos = new TH2D("clPos", ";z;y", 400, -400, 400, 120, -60, 60);

  fClZXref = new TH2D("clZXref", ";z;x", 36, -54, 54, 300, 280, 380);
  fClZXcl =  new TH2D("clZXcl", ";z;x", 36, -54, 54, 300, 280, 380);

  //fGeo = new AliTRDgeometry();
}

//////////////////////////////////////////////////////////////////////////////////////////

void AliTRDtrackingAnalysis::DrawResolutionPt(int startEvent, int stopEvent) 
{
  //
  // Check the pt resolution
  //

  CheckFiles();
  
  // loop over ESD events 
  int nevents = fEsdTree->GetEntries();

  for(int iEvent=startEvent; iEvent<nevents && iEvent < stopEvent; iEvent++) {

    Info("Draw", "Event = %d", iEvent);
    
    fEsdTree->GetEvent(iEvent);
    fLoader->GetEvent(iEvent);
    LoadRefs();
    
    int nTracks = fESD->GetNumberOfTracks();
    for(int iTrack=0; iTrack<nTracks; iTrack++) {
      
      //Info("Track", "Track = %d", iTrack);
      AliESDtrack *esdTrack = fESD->GetTrack(iTrack);
      if (!esdTrack->GetInnerParam()) continue;
      const AliExternalTrackParam *param = esdTrack->GetOuterParam();
      int status = esdTrack->GetStatus();

      if (!(status & AliESDtrack::kTRDout)) continue;
      if (!(status & AliESDtrack::kTRDrefit)) continue;
      if (esdTrack->GetOuterParam()->Pt() < 1.0) continue;

      int ch=0;
      while(param->GetX() > fGeo->GetTime0(ch)+2) ch++;
      fRefSpace->Fill(2.*ch+0.5, param->GetX() - fGeo->GetTime0(ch));
      //if (ch < 5) continue;

      double lastX = 0; 
      int label = abs(esdTrack->GetTRDLabel());
      int ntr = 0;
      int ngood = 0;

      for(int iPoint=GetReference(label); iPoint<fRefTRD->GetEntries(); iPoint++) {

	AliTrackReference *aRef = (AliTrackReference*)(*fRefTRD)[iPoint];
	if (aRef->GetTrack() != label) break;
	ntr++;
      
	lastX = (aRef->LocalX() < lastX)? lastX : aRef->LocalX();
	double dx = aRef->LocalX() - param->GetX();
	if (TMath::Abs(dx) > 1.) continue; 
	ngood++;
      
	double bz=fESD->GetMagneticField();
	AliExternalTrackParam out(*param);
	out.PropagateTo(aRef->LocalX(),bz);
	
	double dp = aRef->Pt() + out.GetSignedPt();
	double dy = 10. * (aRef->LocalY() - out.GetY()); // in mm

	fDeltaPt->Fill(100. * dp / aRef->Pt());
	//fDeltaPt->Fill(out.GetPt() / aRef->Pt());
	fDeltaX->Fill(dx);	

	if (esdTrack->GetSign() > 0) fDeltaYPos->Fill(dy);
	else fDeltaYNeg->Fill(dy);
      
	fDeltaZ->Fill(aRef->Z() - out.GetZ());
      }

      //if (ngood == 0) Info("X", "N = %d, X = %f, DX = %f", ntr, param->GetX(), param->GetX()-lastX);

      fNPoints->Fill(ntr);
      fNGood->Fill(ngood);
    }
  }

  new TCanvas();
  fDeltaPt->Draw();

  new TCanvas();
  fDeltaX->Draw();

  /*
  new TCanvas();
  fNPoints->Draw();
  
  new TCanvas();
  fNGood->Draw();
  */

  new TCanvas();
  fDeltaYPos->Draw();

  new TCanvas();
  fDeltaYNeg->Draw();

  new TCanvas();
  fDeltaZ->Draw();

  new TCanvas();
  fRefSpace->Draw();

}

//////////////////////////////////////////////////////////////////////////////////////////

// void  AliTRDtrackingAnalysis::DrawTrackletResolution(int startEvent, int stopEvent) {

//   LoadRecPointsFile();
  
//   TFile *file = new TFile(Form("%s/TRD.Tracklets.root", fPath), "read");
//   TTree *tree = (TTree*)file->Get("TRDtracklets");
  
//   AliTRDtracklet *tracklet = new AliTRDtracklet();
//   tree->SetBranchAddress("tracklets", &tracklet);

//   for(int ev=startEvent; ev<stopEvent; ev++) {

//     gAlice->GetEvent(ev);
//     LoadRefs();
    
//     int N = tree->GetEntries();
//     for(int i=0; i<N; i++) {
      
//       tree->GetEntry(i);
      
//       Double_t yref, zref, tgphi;
//       int stat = GetMCPosition(tracklet->GetLabel(), tracklet->GetX(), yref, zref, tgphi);
//       if (stat < 0) continue;
      
//       int plane = tracklet->GetPlane();
//       Double_t h01 = tracklet->GetTilt();
      
//       //printf("Tile = %f\tcorrection = %f um \n", h01, 1e4 * h01 * (tracklet->GetZ()-zref));
//       //double dz = zref - tracklet->GetZ() - cls->GetZ();
      
//       fTrklY->Fill(10 * (tracklet->GetY() - yref));
//       fTrklZ->Fill(tracklet->GetZ() - zref);

//       int ch=0;
//       while(tracklet->GetX() > fGeo->GetTime0(ch)+2) ch++;
//       fRefSpace->Fill(tracklet->GetY(), tracklet->GetX() - fGeo->GetTime0(ch));
//     }
//   }
  
//   new TCanvas();
//   fTrklZ->Draw();
  
//   new TCanvas();
//   fRefSpace->Draw();
  
//   gStyle->SetOptFit(1);
//   new TCanvas();
//   fTrklY->Draw();
//   fTrklY->Fit("gaus");

//   //new TCanvas();
//   //fClZ->Draw();
// }

//////////////////////////////////////////////////////////////////////////////////////////

void  AliTRDtrackingAnalysis::DrawRecPointResolution(int startEvent, int stopEvent) 
{
  //
  // Check the resolution of the reconstructed points
  //

  LoadRecPointsFile();
  TObjArray *module = new TObjArray(); 

  int nEvents = AliRunLoader::GetRunLoader()->GetNumberOfEvents();
  
  for(int ev=startEvent; ev<nEvents && ev < stopEvent; ev++) {
    
    gAlice->GetEvent(ev);
    LoadRefs();

    TTree *tree = fLoader->GetTreeR("TRD", 0);
    tree->SetBranchAddress("TRDcluster", &module);

    Info("Res", "Refs Loaded");

    Int_t nn = tree->GetEntries();
    for(int i=0; i<nn; i++) {
      
      tree->GetEntry(i);
      int m = module->GetEntries();

      for(int j=0; j<m; j++) {

	AliTRDcluster *cls = (AliTRDcluster*)module->At(j);
	if (cls->GetQ() < 10) continue;
	//fTracker->Transform(cls);
	fClPos->Fill(cls->GetZ(), cls->GetY());
		
	int layer = fGeo->GetLayer(cls->GetDetector());
	
	int nl = 0;
	for(int k=0; k<3; k++) if (cls->GetLabel(k) > -1) nl++;
	fNLabels->Fill(nl);	

	Double_t yref, zref, tgphi;
	int stat = GetMCPosition(cls->GetLabel(0), cls->GetX(), yref, zref, tgphi);
	if (stat < 0) continue;
	
	fClZXcl->Fill(cls->GetZ(), cls->GetX());
	fClZXref->Fill(zref, cls->GetX());

	AliTRDpadPlane *padPlane = fGeo->GetPadPlane(layer,0);
	Double_t h01   = TMath::Tan(-TMath::Pi() / 180.0 * padPlane->GetTiltingAngle());
	
	//double dz = zref - padPlane->GetRow0();
	double dz = zref - cls->GetZ();
	double dy = dz * h01;
	double yy = cls->GetY() - dy;
		
	if (cls->GetNPads() == 2) fClY2->Fill(10 * (yy - yref));
	if (cls->GetNPads() == 3) fClY3->Fill(10 * (yy - yref));

	int idx = GetPhiBin(tgphi);
	if (idx >= 0 && idx < 12) fClYTgPhi[idx]->Fill(10 * (yy - yref));

	fClZZ->Fill(zref, cls->GetZ());
	fClZ->Fill(dz);
	fTgPhi->Fill(tgphi);
	fClYX->Fill(cls->GetY(), cls->GetX() - fGeo->GetTime0(layer));
      }
    }    
  }

  new TCanvas();
  fClYX->Draw();

  //new TCanvas();
  //fNLabels->Draw();
  new TCanvas();
  fClZ->Draw();

  new TCanvas();
  fTgPhi->Draw();

  gStyle->SetOptFit(1);
  
  new TCanvas();
  gPad->SetLogy();
  fClY2->Draw();
  fClY2->Fit("gaus", "", "", -2*fClY2->GetRMS(), 2*fClY2->GetRMS());

  new TCanvas();
  gPad->SetLogy();
  fClY3->Draw();
  fClY3->Fit("gaus", "", "", -2*fClY3->GetRMS(), 2*fClY3->GetRMS());

  //new TCanvas();
  //fRefSpace->Draw();

  //new TCanvas();
  //fClZXcl->Draw();
  
  //new TCanvas();
  //fClZXref->Draw();

  /**/
  TCanvas *c = new TCanvas();
  c->Divide(4,3);
  
  for(int i=0; i<12; i++) {

    c->cd(i+1);
    fClYTgPhi[i]->Draw();
    if (fClYTgPhi[i]->GetSum() < 100) continue;

    double mean = fClYTgPhi[i]->GetMean();
    double rms = fClYTgPhi[i]->GetRMS();
      
    fClYTgPhi[i]->Fit("gaus", "", "", mean-2*rms, mean+2*rms);
    TF1 *f = fClYTgPhi[i]->GetFunction("gaus");
    
    int n = fGrResTgPhi->GetN();
    fGrResTgPhi->SetPoint(n, GetPhi(i), f->GetParameter(2));
    fGrResTgPhi->SetPointError(n, 0, f->GetParError(2));
    
    fGrMeanTgPhi->SetPoint(n, GetPhi(i), f->GetParameter(1));
    fGrMeanTgPhi->SetPointError(n, 0, f->GetParError(1));
  }

  //gSystem->Sleep(1000);
  gStyle->SetOptStat(0);

  c = new TCanvas();
  TH1D *dummy = new TH1D("dummy", "", 100, -0.3, 0.3);

  dummy->SetTitle(";tg(#phi);resolution (mm)");
  dummy->SetMinimum(0);
  dummy->SetMaximum(1);
  
  //c->cd();
  (dummy->Clone("dummy1"))->Draw();

  fGrResTgPhi->Draw("PL");
  fGrResTgPhi->GetHistogram()->SetTitle(";tg(#phi);resolution (mm)");
  fGrResTgPhi->SetMarkerStyle(20);

  c = new TCanvas();
  dummy->SetTitle(";tg(#phi);mean value (mm)");
  dummy->SetMinimum(-0.3);
  dummy->SetMaximum(0.3);
  
  //c->cd();
  (dummy->Clone("dummy2"))->Draw();
  
  fGrMeanTgPhi->Draw("PL");
  fGrMeanTgPhi->GetHistogram()->SetTitle(";tg(#phi);mean value (mm)");
  fGrMeanTgPhi->SetMarkerStyle(20);
  /**/
  

  //new TCanvas();
  //fClZZ->Draw("colz");

  //new TCanvas();
  //fClPos->Draw("colz");

  //new TCanvas();
  //fTestBits->Draw();

  //new TCanvas();
  //fRefDx->Draw();

}

//////////////////////////////////////////////////////////////////////////////////////////

void AliTRDtrackingAnalysis::LoadRecPointsFile() 
{
  //
  // Load the clusters from the input file
  //

  char filename[256];
  sprintf(filename, "%s/galice.root", fPath);
  
  fLoader = AliRunLoader::Open(filename);
  if (!fLoader) {
    Error("CheckFiles", "getting run loader from file %s/galice.root failed", filename);
    return;
  }
  
  fLoader->LoadgAlice();
  gAlice = fLoader->GetAliRun();
  
  if (!gAlice) {
    Error("CheckFiles", "no galice object found");
    return;
  }
  
  fLoader->LoadKinematics();
  fLoader->LoadHeader();
  fLoader->LoadTrackRefs();

  //TGeoManager::Import("/data/alice_u/radomski/condor/run_0/geometry.root");
  TGeoManager::Import(Form("%s/geometry.root", fPath));


  fLoader->CdGAFile();
  fGeo = (AliTRDgeometry*)gDirectory->Get("TRDgeometry");
  fTracker = new AliTRDtracker(gFile);

  fLoader->LoadRecPoints("TRD");

  AliTracker::SetFieldMap(gAlice->Field(), 1);
} 

//////////////////////////////////////////////////////////////////////////////////////////

void  AliTRDtrackingAnalysis::CheckFiles() 
{
  //
  // Check the presence of the input files
  //

  // MC info

  char filename[256];
  sprintf(filename, "%s/galice.root", fPath);
  
  fLoader = AliRunLoader::Open(filename);
  if (!fLoader) {
    Error("CheckFiles", "getting run loader from file %s/galice.root failed", filename);
    return;
  }
  
  fLoader->LoadgAlice();
  gAlice = fLoader->GetAliRun();
  
  if (!gAlice) {
    Error("CheckFiles", "no galice object found");
    return;
  }
  
  fLoader->LoadKinematics();
  fLoader->LoadHeader();
  fLoader->LoadTrackRefs();
  
  fLoader->CdGAFile();
  fGeo = (AliTRDgeometry*)gDirectory->Get("TRDgeometry");
  //fGeo->ReadGeoMatrices();
  
  // ESD
  
  sprintf(filename,"%s/AliESDs.root", fPath);
  TFile *esdFile = new TFile(filename, "READ");
  
  if (esdFile->IsZombie()) {
    Error("CheckFiles", "file not present: AliESDs.root");
    return;
  }
  
  fEsdTree = (TTree*)esdFile->Get("esdTree"); 
  fESD = new AliESDEvent();
  fESD->ReadFromTree(fEsdTree);
  //fEsdTree->SetBranchAddress("ESD", &fESD);
}

//////////////////////////////////////////////////////////////////////////////////////////

void  AliTRDtrackingAnalysis::LoadRefs() 
{
  //
  // Load the track references
  //

  if (fRefTPC) delete fRefTPC;
  if (fRefTRD) delete fRefTRD;
		  
  fRefTPC = new TObjArray();
  fRefTRD = new TObjArray();
    
  //fLoader->GetEvent(event);
  //AliStack* stack = gAlice->Stack();
  TTree *refTree = fLoader->TreeTR();
    
  TClonesArray *clRefs = new TClonesArray("AliTrackReference");
      
  TBranch *branch = refTree->GetBranch("TrackReferences");
  refTree->SetBranchAddress("TrackReferences",&clRefs);
    
  int nEntries = branch->GetEntries();      
  for(int iTrack = 0; iTrack < nEntries; iTrack++) {
	
    refTree->GetEvent(iTrack);
    int nPoints =  clRefs->GetEntries();
    for(int iPoint=0; iPoint<nPoints; iPoint++) {
      AliTrackReference *ref = (AliTrackReference*)clRefs->At(iPoint);
	if (ref->DetectorId() == AliTrackReference::kTPC) fRefTPC->Add(new AliTrackReference(*ref));
	if (ref->DetectorId() == AliTrackReference::kTRD) fRefTRD->Add(new AliTrackReference(*ref));	  
    }	
  }
  
  fRefTPC->Sort();
  fRefTRD->Sort();

  for(int i=0; i<fRefTRD->GetEntries(); i++) {
    AliTrackReference *ref = (AliTrackReference*)(*fRefTRD)[i]; 
    fLabels[i] = ref->GetTrack();
    
    int p=0;
    while(ref->LocalX() > fGeo->GetTime0(p)+2) p++;
    fRefSpace->Fill(ref->LocalY(), ref->LocalX()-fGeo->GetTime0(p));

    //for(int bit=0; bit<9; bit++) if (ref->TestBit(bit)) fTestBits->Fill(bit);
  }

  delete clRefs;
  Info("LoadRefs", "TPC = %d\t TRD = %d", fRefTPC->GetEntries(), fRefTRD->GetEntries());
}

//////////////////////////////////////////////////////////////////////////////////////////

Int_t AliTRDtrackingAnalysis::GetReference(Int_t label) 
{
  //
  // Sort the track references
  //
  
  int start = TMath::BinarySearch(fRefTRD->GetEntries(), fLabels, label);
  
  while (start >= 0) {
    AliTrackReference *ref = (AliTrackReference*)(*fRefTRD)[start];
    if (ref->GetTrack() != label) return start+1;
    start--;
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////

int AliTRDtrackingAnalysis::GetMCPosition(Int_t label, Double_t x, Double_t &Y, Double_t &Z, Double_t &tgphi) 
{
  //
  // Determine the MC positions from the track references
  //
  
  double lowX = 100.;
  double highX = 100.;
  int idLow = -1;
  int idHigh = -1;

  int nref= 0;
  int idx = GetReference(label);
  for(int i=idx; i<fRefTRD->GetEntries(); i++) {
    
    AliTrackReference *ref = (AliTrackReference*)(*fRefTRD)[i];
    if (ref->GetTrack() != label) break;
    nref++;

    //int p=0;
    //while(ref->LocalX() > fGeo->GetTime0(p)+2) p++;
    //if (p != layer) continue;
    
    double dX = ref->LocalX()-x;
    if ( dX > 0 ) {
      if (dX < highX) {
	idHigh = i;
	highX = dX;
      }
    } else {
      dX = TMath::Abs(dX);
      if (dX < lowX) {
	idLow = i;
	lowX = dX;
      }
    }
  }
  
  if (idLow == -1 || idHigh == -1) return -1;
  
  AliTrackReference *refI = (AliTrackReference*)(*fRefTRD)[idLow];
  AliTrackReference *refO = (AliTrackReference*)(*fRefTRD)[idHigh];

  
  double dx = refO->LocalX() - refI->LocalX();
  double dy = refO->LocalY() - refI->LocalY();
  double dz = refO->Z() - refI->Z();
  double ddx = (x - refI->LocalX())/dx;
 
  fRefDx->Fill(dx);

  Y = refI->LocalY() + ddx * dy;
  Z = refI->Z() + ddx * dz;

  tgphi = dy/dx;

  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
Int_t AliTRDtrackingAnalysis::GetPhiBin(Double_t phi) const
{
  //
  // Return the phi bin
  //
  return (int)((phi+0.3)/0.05);  
}

//////////////////////////////////////////////////////////////////////////////////////////
Double_t AliTRDtrackingAnalysis::GetPhi(Int_t bin) const
{
  //
  // Return phi for a given bin
  //
  return bin * 0.05 - 0.3 + 0.025; 
}
//////////////////////////////////////////////////////////////////////////////////////////
