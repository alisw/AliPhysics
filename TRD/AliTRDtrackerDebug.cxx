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

/* $Id: AliTRDtrackerDebug.cxx 23810 2008-02-08 09:00:27Z hristov $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Tracker debug streamer                                                   //
//                                                                           //
//  Authors:                                                                 //
//    Alex Bercuci <A.Bercuci@gsi.de>                                        //
//    Markus Fasel <M.Fasel@gsi.de>                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDtrackerDebug.h"

#include "TFile.h"
#include "TTree.h"
#include "TTreeStream.h"
#include "TLinearFitter.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"

#include "AliLog.h"
#include "AliTRDgeometry.h"
#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"
#include "AliTRDseed.h"
#include "AliTRDcluster.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDtrackerDebug)

Int_t AliTRDtrackerDebug::fgEventNumber = 0;
Int_t AliTRDtrackerDebug::fgTrackNumber = 0;
Int_t AliTRDtrackerDebug::fgCandidateNumber = 0;

//____________________________________________________
AliTRDtrackerDebug::AliTRDtrackerDebug() : AliTRDtrackerV1()
  ,fOutputStreamer(0x0)
  ,fTree(0x0)
  ,fTracklet(0x0)
  ,fTrack(0x0)
  ,fNClusters(0)
  ,fAlpha(0.)
{
        //
  // Default constructor
  //
  fOutputStreamer = new TTreeSRedirector("TRD.Debug.root");
}

//____________________________________________________
AliTRDtrackerDebug::~AliTRDtrackerDebug()
{
  // destructor
  
  delete fOutputStreamer;
}


//____________________________________________________
void AliTRDtrackerDebug::Draw(Option_t *)
{
// steer draw function
}


//____________________________________________________
Bool_t AliTRDtrackerDebug::Init()
{
// steer linking data for various debug streams	
  fTrack = new AliTRDtrackV1();
  fTree->SetBranchAddress("ncl", &fNClusters);
  fTree->SetBranchAddress("track.", &fTrack);
  return kTRUE;
}

//____________________________________________________
Bool_t AliTRDtrackerDebug::Open(const char *method)
{
  // Connect to the tracker debug file
  
  TDirectory *savedir = gDirectory; 
  TFile::Open("TRD.TrackerDebugger.root");
  fTree = (TTree*)gFile->Get(method);
  if(!fTree){
    AliInfo(Form("Can not find debug stream for the %s method.\n", method));
    savedir->cd();
    return kFALSE;
  }
  savedir->cd();
  return kTRUE;
}

//____________________________________________________
Int_t AliTRDtrackerDebug::Process()
{
// steer debug process threads
  
  for(int it = 0; it<fTree->GetEntries(); it++){
    if(!fTree->GetEntry(it)) continue;
    if(!fNClusters) continue;
    fAlpha = fTrack->GetAlpha();
    //printf("Processing track %d [%d] ...\n", it, fNClusters);
    ResidualsTrackletsTrack();

    const AliTRDseedV1 *tracklet = 0x0;
    for(int ip = 5; ip>=0; ip--){
      if(!(tracklet = fTrack->GetTracklet(ip))) continue;
      if(!tracklet->GetN()) continue;
      
      ResidualsClustersTrack(tracklet);
      ResidualsClustersTracklet(tracklet);
      ResidualsClustersParametrisation(tracklet);
    }
  }
  return kTRUE;
}


//____________________________________________________
void AliTRDtrackerDebug::ResidualsClustersTrack(const AliTRDseedV1 *tracklet)
{
// Calculate averange distances from clusters to the TRD track	
  
  Double_t x[3]; 
  AliTRDcluster *c = 0x0;
  for(int ic=0; ic<35/*AliTRDseed:knTimebins*/; ic++){
    if(!(c = tracklet->GetClusters(ic))) continue;
    Double_t xc = c->GetX(), yc = c->GetY(), zc = c->GetZ();

    // propagate track to cluster 
    PropagateToX(*fTrack, xc, 2.); 
    fTrack->GetXYZ(x);
    
    // transform to local tracking coordinates
    //Double_t xg =  x[0] * TMath::Cos(fAlpha) + x[1] * TMath::Sin(fAlpha); 
    Double_t yg = -x[0] * TMath::Sin(fAlpha) + x[1] * TMath::Cos(fAlpha);

    // apply tilt pad correction
    yc+= (zc - x[2]) * tracklet->GetTilt();
    
    Double_t dy = yc-yg;

    TTreeSRedirector &cstreamer = *fOutputStreamer;
    cstreamer << "ResidualsClustersTrack"
      << "c.="   << c
      << "dy="   << dy
      << "\n";
  }
}

//____________________________________________________
void AliTRDtrackerDebug::ResidualsClustersTracklet(const AliTRDseedV1 *tracklet) const
{
// Calculates distances from clusters to tracklets
  
  Double_t x0 = tracklet->GetX0(), 
          y0 = tracklet->GetYfit(0), 
          ys = tracklet->GetYfit(1);
          //z0 = tracklet->GetZfit(0), 
          //zs = tracklet->GetZfit(1);
  
  AliTRDcluster *c = 0x0;
  for(int ic=0; ic<35/*AliTRDseed:knTimebins*/; ic++){
    if(!(c = tracklet->GetClusters(ic))) continue;
    Double_t xc = c->GetX(), yc = c->GetY()/*, zc = c->GetZ()*/;
    Double_t dy = yc- (y0-(x0-xc)*ys);

    //To draw  use : 
    //ResidualsClustersTracklet->Draw("TMath::Abs(10.*dy):TMath::ATan(ys)*TMath::RadToDeg()>>h(20, -40, 40)", "", "prof");
    TTreeSRedirector &cstreamer = *fOutputStreamer;
    cstreamer << "ResidualsClustersTracklet"
      << "c.="   << c
      << "ys="   << ys
      << "dy="   << dy
      << "\n";
  }
}

//____________________________________________________
void AliTRDtrackerDebug::ResidualsClustersParametrisation(const AliTRDseedV1 *tracklet) const
{
// Calculates distances from clusters to Rieman fit.
  
  // store cluster positions
  Double_t x0 = tracklet->GetX0();
  AliTRDcluster *c = 0x0;
  
  Double_t x[2]; Int_t ncl, mcl, jc;
  TLinearFitter fitter(3, "hyp2");
  for(int ic=0; ic<35/*AliTRDseed:knTimebins*/; ic++){
    if(!(c = tracklet->GetClusters(ic))) continue;
    Double_t xc = c->GetX(), yc = c->GetY()/*, zc = c->GetZ()*/;
    
    jc = ic; ncl = 0; mcl=0; fitter.ClearPoints();
    while(ncl<6){
      // update index
      mcl++;
      jc = ic + ((mcl&1)?-1:1)*(mcl>>1);

      if(jc<0 || jc>=35) continue;
      if(!(c = tracklet->GetClusters(jc))) continue;

      x[0] = c->GetX()-x0;
      x[1] = x[0]*x[0];
      fitter.AddPoint(x, c->GetY(), c->GetSigmaY2());
      ncl++;
    }
    fitter.Eval();
    Double_t dy = yc - fitter.GetParameter(0) -fitter.GetParameter(1) * (xc-x0) - fitter.GetParameter(2)* (xc-x0)*(xc-x0); 
  
    TTreeSRedirector &cstreamer = *fOutputStreamer;
    cstreamer << "ResidualsClustersParametrisation"
      << "dy="   << dy
      << "\n";
  }
}


//____________________________________________________
void AliTRDtrackerDebug::ResidualsTrackletsTrack() const
{
// Calculates distances from tracklets to the TRD track.
  
  if(fTrack->GetNumberOfTracklets() < 6) return;

  // build a working copy of the tracklets attached to the track 
  // and initialize working variables fX, fY and fZ
  AliTRDseedV1 tracklet[6] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  const AliTRDseedV1 *ctracklet = 0x0;
  for(int ip = 0; ip<6; ip++){
    if(!(ctracklet = fTrack->GetTracklet(ip))) continue;
    tracklet[ip] = (*ctracklet); 
// 		Double_t x0 = tracklet[ip].GetX0();
// 		for(int ic=0; ic<AliTRDseedV1:knTimebins; ic++){
// 			if(!(c = tracklet[ip].GetClusters(ic))) continue;
// 			Double_t xc = c->GetX(), yc = c->GetY(), zc = c->GetZ();
// 			tracklet[ip].SetX(ic, xc-x0);
// 			tracklet[ip].SetY(ic, yc);
// 			tracklet[ip].SetZ(ic, zc);
// 		}
  }
  
  // Do a Rieman fit (with tilt correction) for all tracklets 
  // except the one which is tested. 
  // (Based on AliTRDseed::IsOK() return false)
  for(int ip=0; ip<6; ip++){
    // reset tracklet to be tested
    Double_t x0 = tracklet[ip].GetX0();
    new(&tracklet[ip]) AliTRDseedV1();
    tracklet[ip].SetX0(x0);

    // fit Rieman with tilt correction
    AliTRDtrackerV1::FitRiemanTilt(0x0, &tracklet[0], kTRUE);

    // make a copy of the fit result
    Double_t 
      y0   = tracklet[ip].GetYref(0),
      dydx = tracklet[ip].GetYref(1),
      z0   = tracklet[ip].GetZref(0),
      dzdx = tracklet[ip].GetZref(1);

    // restore tracklet
    tracklet[ip] = (*fTrack->GetTracklet(ip)); 
// 		for(int ic=0; ic<AliTRDseedV1:knTimebins; ic++){
// 			if(!(c = tracklet[ip].GetClusters(ic))) continue;
// 			Double_t xc = c->GetX(), yc = c->GetY(), zc = c->GetZ();
// 			tracklet[ip].SetX(ic, xc-x0);
// 			tracklet[ip].SetY(ic, yc);
// 			tracklet[ip].SetZ(ic, zc);
// 		}		
    
    // fit clusters
    AliTRDseedV1 ts(tracklet[ip]);
    ts.SetYref(0, y0); ts.SetYref(1, dydx);
    ts.SetZref(0, z0); ts.SetZref(1, dzdx);
    ts.FitMI();

    // save results for plotting
    Int_t plane   = tracklet[ip].GetPlane();
    Double_t dy   = tracklet[ip].GetYfit(0) - ts.GetYfit(0);
    Double_t tgy  = tracklet[ip].GetYfit(1);
    Double_t dtgy = (tracklet[ip].GetYfit(1) - ts.GetYfit(1))/(1. + tracklet[ip].GetYfit(1) * ts.GetYfit(1));
    Double_t dz   = tracklet[ip].GetZfit(0) - ts.GetZfit(0);
    Double_t tgz  = tracklet[ip].GetZfit(1);
    Double_t dtgz = (tracklet[ip].GetZfit(1) - ts.GetZfit(1))/(1. + tracklet[ip].GetZfit(1) * ts.GetZfit(1));
    TTreeSRedirector &cstreamer = *fOutputStreamer;
    cstreamer << "ResidualsTrackletsTrack"
      << "ref.="   << &tracklet[ip]
      << "fit.="   << &ts
      << "plane="  << plane
      << "dy="     << dy
      << "tgy="    << tgy
      << "dtgy="   << dtgy
      << "dz="     << dz
      << "tgz="    << tgz
      << "dtgz="   << dtgz
      << "\n";
  }
}

//____________________________________________________
void AliTRDtrackerDebug::AnalyseFindable(Char_t *treename){
//
// Calculates the number of findable tracklets defined as the number of tracklets
// per track candidate where the tan phi_tracklet is below 0.15 (maximum inclination
// in y-direction.
// 
// Parameters:	-the treename (this method can be used for all trees which store the
//				 tracklets
// Output:		-void
//
// A new TTree containing the number of findable tracklets and the number of clusters
// attached to the full track is stored to disk
//
  // Link the File
  TFile *debfile = TFile::Open("TRD.TrackerDebug.root");
  fTree = (TTree *)(debfile->Get(treename));
  if(!fTree){
    AliError(Form("Tree %s not found in file TRDdebug.root. Abborting!", treename));
    debfile->Close();
    return;
  }
  
  AliTRDseedV1 *tracklets[kNPlanes];
  for(Int_t iPlane = 0; iPlane < AliTRDtrackerV1::kNPlanes; iPlane++)
    tracklets[iPlane] = 0x0;
  for(Int_t iPlane = 0; iPlane < kNPlanes; iPlane++)
    fTree->SetBranchAddress(Form("S%d.", iPlane), &tracklets[iPlane]);
  fTree->SetBranchAddress("EventNumber", &fgEventNumber);
  fTree->SetBranchAddress("CandidateNumber", &fgCandidateNumber);
  
  Int_t findable = 0, nClusters = 0;
  Int_t nEntries = fTree->GetEntriesFast();
  for(Int_t iEntry = 0; iEntry < nEntries; iEntry++){
    printf("Entry %d\n", iEntry);
    fTree->GetEntry(iEntry);
    findable = 0;
    nClusters = 0;
    // Calculate Findable
    for(Int_t iPlane = 0; iPlane < kNPlanes; iPlane++){
      if (TMath::Abs(tracklets[iPlane]->GetYref(0) / tracklets[iPlane]->GetX0()) < 0.15) findable++;
      if (!tracklets[iPlane]->IsOK()) continue;
      nClusters += tracklets[iPlane]->GetN2();
    }
    
    // Fill Histogramms
    TTreeSRedirector &cstreamer = *fOutputStreamer;
    cstreamer << "AnalyseFindable"
      << "EventNumber="		<< fgEventNumber
      << "CandidateNumber="	<< fgCandidateNumber
      << "Findable="			<< findable
      << "NClusters="			<< nClusters
      << "\n";
  }
}
//____________________________________________________
void AliTRDtrackerDebug::AnalyseTiltedRiemanFit(){
//
// Creating a Data Set for the method FitTiltedRieman containing usefull variables
// Each variable can be addressed to tracks later. Data can be processed later.
//
// Parameters: -
// Output:     -
//
// TODO: Match everything with Init and Process
//
  TFile *debfile = TFile::Open("TRD.TrackerDebug.root");
  fTree = (TTree *)(debfile->Get("MakeSeeds2"));
  if(!fTree) return;
  Int_t nEntries = fTree->GetEntries();
  TLinearFitter *tiltedRiemanFitter = 0x0;
  fTree->SetBranchAddress("FitterT.", &tiltedRiemanFitter);
  fTree->SetBranchAddress("EventNumber", &fgEventNumber);
  fTree->SetBranchAddress("CandidateNumber", &fgCandidateNumber);
  for(Int_t entry = 0; entry < nEntries; entry++){
    fTree->GetEntry(entry);
    Double_t a = tiltedRiemanFitter->GetParameter(0);
    Double_t b = tiltedRiemanFitter->GetParameter(1);
    Double_t c = tiltedRiemanFitter->GetParameter(2);
    Double_t offset = tiltedRiemanFitter->GetParameter(3);
    Double_t slope  = tiltedRiemanFitter->GetParameter(4);
    Float_t radius = GetTrackRadius(a, b, c);
    Float_t curvature = GetTrackCurvature(a, b, c);
    Float_t dca = GetDCA(a, b, c);
    TTreeSRedirector &cstreamer = *fOutputStreamer;
    cstreamer << "AnalyseTiltedRiemanFit"
    << "EventNumber=" 		<< fgEventNumber
    << "CandidateNumber=" << fgCandidateNumber
    << "Radius="					<< radius
    << "Curvature="				<< curvature
    << "DCA="							<< dca
    << "Offset="					<< offset
    << "Slope="						<< slope
    << "\n";
  }
}

//____________________________________________________
Float_t AliTRDtrackerDebug::GetTrackRadius(Float_t a, Float_t b, Float_t c) const {
//
// Calculates the track radius using the parameters given by the tilted Rieman fit 
//
// Parameters: The three parameters from the Rieman fit
// Output:     The track radius
//
  Float_t radius = 0;
  if(1.0 + b*b - c*a > 0.0)
    radius = TMath::Sqrt(1.0 + b*b - c*a )/a;
  return radius;
}

//____________________________________________________
Float_t AliTRDtrackerDebug::GetTrackCurvature(Float_t a, Float_t b, Float_t c) const {
//
// Calculates the track curvature using the parameters given by the linear fitter 
//
// Parameters:	the three parameters from the tilted Rieman fitter
// Output:		the full track curvature
//
  Float_t curvature =  1.0 + b*b - c*a;
  if (curvature > 0.0) 
    curvature  =  a / TMath::Sqrt(curvature);
  return curvature;
}

//____________________________________________________
Float_t AliTRDtrackerDebug::GetDCA(Float_t a, Float_t b, Float_t c) const {
//
// Calculates the Distance to Clostest Approach for the Vertex using the paramters
// given by the tilted Rieman fit 
//
// Parameters: the three parameters from the tilted Rieman fitter
// Output:     the Distance to Closest Approach
//
  Float_t dca  =  0.0;
  if (1.0 + b*b - c*a > 0.0) {
    dca = -c / (TMath::Sqrt(1.0 + b*b - c*a) + TMath::Sqrt(1.0 + b*b));
  }
  return dca;
}

//____________________________________________________
void AliTRDtrackerDebug::AnalyseMinMax()
{
//
  TFile *debfile = TFile::Open("TRD.TrackerDebug.root");
  if(!debfile){
    AliError("File TRD.TrackerDebug.root not found!");
    return; 
  }
  fTree = (TTree *)(debfile->Get("MakeSeeds0"));
  if(!fTree){
    AliError("Tree MakeSeeds0 not found in File TRD.TrackerDebug.root.");
    return;
  }
  AliTRDseedV1 *cseed[4] = {0x0, 0x0, 0x0, 0x0};
  AliTRDcluster *c[4] = {0x0, 0x0, 0x0, 0x0};
  for(Int_t il = 0; il < 4; il++){
    fTree->SetBranchAddress(Form("Seed%d.", il),	&cseed[il]);
    fTree->SetBranchAddress(Form("c%d.",il), &c[il]);
  }
  fTree->SetBranchAddress("CandidateNumber",	&fgCandidateNumber);
  fTree->SetBranchAddress("EventNumber",	&fgEventNumber);
  Int_t entries = fTree->GetEntries();
  for(Int_t ientry = 0; ientry < entries; ientry++){
    fTree->GetEntry(ientry);
    Float_t minmax[2] = { -100.0,  100.0 };
    for (Int_t iLayer = 0; iLayer < 4; iLayer++) {
      Float_t max = c[iLayer]->GetZ() + cseed[iLayer]->GetPadLength() * 0.5 + 1.0 - cseed[iLayer]->GetZref(0);
      if (max < minmax[1]) minmax[1] = max;
      Float_t min = c[iLayer]->GetZ()-cseed[iLayer]->GetPadLength() * 0.5 - 1.0 - cseed[iLayer]->GetZref(0);
      if (min > minmax[0]) minmax[0] = min;
    }
    TTreeSRedirector &cstreamer = *fOutputStreamer;
    cstreamer << "AnalyseMinMaxLayer"
    << "EventNumber="				<< fgEventNumber
    << "CandidateNumber="		<< fgCandidateNumber
    << "Min="								<< minmax[0]
    << "Max="								<< minmax[1]
    << "\n";
  }
}

//____________________________________________________
TCanvas* AliTRDtrackerDebug::PlotSeedingConfiguration(const Char_t *direction, Int_t event, Int_t candidate){
//
// Plots the four seeding clusters, the helix fit and the reference Points for
// a given combination consisting of eventnr. and candidatenr.
//
// Parameters: 	-direction (y or z)
//				-Event Nr
//            	-Candidate that has to be plotted
//
  const Float_t kxmin = 280;
  const Float_t kxmax = 380;
  const Float_t kxdelta = (kxmax - kxmin)/1000;
  
  if((strcmp(direction, "y") != 0) && (strcmp(direction, "z") != 0)){
    AliError(Form("Direction %s does not exist. Abborting!", direction));
    return 0x0;
  }

  TFile *debfile = TFile::Open("TRD.TrackerDebug.root");
  if(!debfile){
    AliError("File TRD.TrackerDebug.root not found!");
    return 0x0; 
  }
  fTree = (TTree *)(debfile->Get("MakeSeeds0"));
  if(!fTree){
    AliError("Tree MakeSeeds0 not found in File TRD.TrackerDebug.root.");
    return 0x0;
  }
  
  TGraph *seedcl = new TGraph(4);
  TGraph *seedRef = new TGraph(4);
  TGraph *riemanFit = new TGraph(1000);
  seedcl->SetMarkerStyle(20);
  seedcl->SetMarkerColor(kRed);
  seedRef->SetMarkerStyle(2);

  AliTRDcluster *c[4] = {0x0, 0x0, 0x0, 0x0};
  AliRieman *rim = 0x0;
  Bool_t found = kFALSE;
  for(Int_t il = 0; il < 4; il++) fTree->SetBranchAddress(Form("c%d.",il), &c[il]);
  fTree->SetBranchAddress("EventNumber", &fgEventNumber);
  fTree->SetBranchAddress("CandidateNumber", &fgCandidateNumber);
  fTree->SetBranchAddress("RiemanFitter.", &rim);
  Int_t entries = fTree->GetEntries();
  for(Int_t entry = 0; entry < entries; entry++){
    fTree->GetEntry(entry);
    if(fgEventNumber < event) continue;
    if(fgEventNumber > event) break;
    // EventNumber matching: Do the same for the candidate number
    if(fgCandidateNumber < candidate) continue;
    if(fgCandidateNumber > candidate) break;
    found = kTRUE;
    Int_t nPoints = 0;
    for(Int_t il = 0; il < 4; il++){
      Float_t cluster = 0.0, reference = 0.0;
      if(!strcmp(direction, "y")){
        cluster = c[il]->GetY();
        reference = rim->GetYat(c[il]->GetX());
      }
      else{
        cluster = c[il]->GetZ();
        reference = rim->GetZat(c[il]->GetX());
      }
      seedcl->SetPoint(nPoints, cluster, c[il]->GetX());
      seedRef->SetPoint(nPoints, reference , c[il]->GetX());
      nPoints++;
    }
    // evaluate the fitter Function numerically
    nPoints = 0;
    for(Int_t ipt = 0; ipt < 1000; ipt++){
      Float_t x = kxmin + ipt * kxdelta;
      Float_t point = 0.0;
      if(!strcmp(direction, "y"))
        point = rim->GetYat(x);
      else
        point = rim->GetZat(x);
      riemanFit->SetPoint(nPoints++, point, x);
    }
    // We reached the End: break
    break;
  }
  if(found){
    seedcl->SetTitle(Form("Event %d, Candidate %d\n", fgEventNumber, fgCandidateNumber));
    seedRef->SetTitle(Form("Event %d, Candidate %d\n", fgEventNumber, fgCandidateNumber));
    riemanFit->SetTitle(Form("Event %d, Candidate %d\n", fgEventNumber, fgCandidateNumber));
    TCanvas *c1 = new TCanvas();
    seedcl->Draw("ap");
    seedRef->Draw("psame");
    riemanFit->Draw("lpsame");
    return c1;
  }
  else{
    AliError(Form("Combination consisting of event %d and candidate %d not found", event, candidate));
    delete seedcl;
    delete seedRef;
    delete riemanFit;
    return 0x0;
  }
}

//____________________________________________________
TCanvas* AliTRDtrackerDebug::PlotFullTrackFit(Int_t event, Int_t candidate, Int_t iteration, const Char_t *direction){
//
// Plots the tracklets (clusters and reference in y direction) and the fitted function for several iterations
// in the function ImproveSeedQuality (default is before ImproveSeedQuality)
// 
// Parameters: -Event Number
//             -Candidate Number
//             -Iteration Number in ImproveSeedQuality (default: -1 = before ImproveSeedQuality)
//			   -direction (default: y)
// Output:     -TCanvas (containing the Picture);
//
  const Float_t kxmin = 280;
  const Float_t kxmax = 380;
  const Float_t kxdelta = (kxmax - kxmin)/1000;
  
  if(strcmp(direction, "y") && strcmp(direction, "z")){
    AliError(Form("Direction %s does not exist. Abborting!", direction));
    return 0x0;
  }

  TFile *debfile = TFile::Open("TRD.TrackerDebug.root");
  if(!debfile){
    AliError("File TRD.TrackerDebug.root not found.");
    return 0x0;
  }
  TString *treename = 0x0;
  if(iteration > -1)
    treename = new TString("ImproveSeedQuality");
  else
    treename = new TString("MakeSeeds1");
  fTree = (TTree *)(debfile->Get(treename->Data()));
  if(!fTree){
    AliError(Form("Tree %s not found in File TRD.TrackerDebug.root.", treename->Data()));
    delete treename;
    return 0x0;
  }
  delete treename;

  TGraph *fitfun = new TGraph(1000);
  // Prepare containers
  Float_t x0[AliTRDtrackerV1::kNPlanes],
      refP[AliTRDtrackerV1::kNPlanes],
      clx[AliTRDtrackerV1::kNPlanes * AliTRDtrackerV1::kNTimeBins],
      clp[AliTRDtrackerV1::kNPlanes * AliTRDtrackerV1::kNTimeBins];
  Int_t nLayers = 0, ncls = 0;
  
  TLinearFitter *fitter = 0x0;
  AliTRDseedV1 *tracklet[6] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  for(Int_t iLayer = 0; iLayer < 6; iLayer++)
    fTree->SetBranchAddress(Form("S%d.", iLayer), &tracklet[iLayer]);
  fTree->SetBranchAddress("FitterT.", &fitter);
  fTree->SetBranchAddress("EventNumber", &fgEventNumber);
  fTree->SetBranchAddress("CandidateNumber", &fgCandidateNumber);
  
  Int_t nEntries = fTree->GetEntriesFast();
  Bool_t found = kFALSE;
  for(Int_t entry = 0; entry < nEntries; entry++){
    fTree->GetEntry(entry);
    if(fgEventNumber < event) continue;
    if(fgEventNumber > event) break;
    // EventNumber matching: Do the same for the candidate number
    if(fgCandidateNumber < candidate) continue;
    if(fgCandidateNumber > candidate) break;
    found = kTRUE;
    
    for(Int_t iLayer = 0; iLayer < 6; iLayer++){
      if(!tracklet[iLayer]->IsOK()) continue;
      x0[nLayers] = tracklet[iLayer]->GetX0();
      if(!strcmp(direction, "y"))
        refP[nLayers] = tracklet[iLayer]->GetYref(0);
      else
        refP[nLayers] = tracklet[iLayer]->GetZref(0);
      nLayers++;
      for(Int_t itb = 0; itb < 30; itb++){
        if(!tracklet[iLayer]->IsUsable(itb)) continue;
        AliTRDcluster *cl = tracklet[iLayer]->GetClusters(itb);
        if(!strcmp(direction, "y"))
          clp[ncls] = cl->GetY();
        else
          clp[ncls] = cl->GetZ();
        clx[ncls] = cl->GetX();
        ncls++;
      }
    }
    // Add function derived by the tilted Rieman fit (Defined by the curvature)
    Int_t nPoints = 0;
    if(!strcmp(direction, "y")){
      Double_t a = fitter->GetParameter(0);
      Double_t b = fitter->GetParameter(1);
      Double_t c = fitter->GetParameter(2);
      Double_t curvature =  1.0 + b*b - c*a;
      if (curvature > 0.0) {
        curvature  =  a / TMath::Sqrt(curvature);
      }
      // Numerical evaluation of the function:
      for(Int_t ipt = 0; ipt < 1000; ipt++){
        Float_t x = kxmin + ipt * kxdelta;
        Double_t res = (x * a + b);								// = (x - x0)/y0
        res *= res;
        res  = 1.0 - c * a + b * b - res;					// = (R^2 - (x - x0)^2)/y0^2
        Double_t y = 0.;
        if (res >= 0) {
          res = TMath::Sqrt(res);
          y    = (1.0 - res) / a;
        }
        fitfun->SetPoint(nPoints++, y, x);
      }
    }
    else{
      Double_t offset	= fitter->GetParameter(3);
      Double_t slope	= fitter->GetParameter(4);	 
      // calculate the reference x (defined as medium between layer 2 and layer 3)
      // same procedure as in the tracker code
      Float_t medx = 0, xref = 0;
      Int_t startIndex = 5, nDistances = 0;
      for(Int_t iLayer = 5; iLayer > 0; iLayer--){
        if(tracklet[iLayer]->IsOK() && tracklet[iLayer - 1]->IsOK()){
          medx += tracklet[iLayer]->GetX0() - tracklet[iLayer - 1]->GetX0();
          startIndex = iLayer - 1;
          nDistances++;
        }
      }
      if(nDistances){
        medx /= nDistances;
      }
      else{
        Float_t xpos[2];	memset(xpos, 0, sizeof(Float_t) * 2);
        Int_t ien = 0, idiff = 0;
        for(Int_t iLayer = 5; iLayer > 0; iLayer--){
          if(tracklet[iLayer]->IsOK()){
            xpos[ien++] = tracklet[iLayer]->GetX0();
            startIndex = iLayer;
          }
          if(ien)
            idiff++;
          if(ien >=2)
            break;
        }
        medx = (xpos[0] - xpos[1])/idiff;
      }
      xref = tracklet[startIndex]->GetX0() + medx * (2.5 - startIndex) - 0.5 * (AliTRDgeometry::AmThick() + AliTRDgeometry::DrThick());

      for(Int_t ipt = 0; ipt < 1000; ipt++){
        Float_t x = kxmin + ipt * kxdelta;
        Float_t z = offset + slope * (x - xref);
        fitfun->SetPoint(nPoints++, z, x);
      }
    }
    break;
  }
  if(found){
    TGraph *trGraph		= new TGraph(ncls);
    TGraph *refPoints	= new TGraph(nLayers);
    trGraph->SetMarkerStyle(20);
    trGraph->SetMarkerColor(kRed);
    refPoints->SetMarkerStyle(21);
    refPoints->SetMarkerColor(kBlue);
    // fill the graphs
    for(Int_t iLayer = 0; iLayer < nLayers; iLayer++)
      refPoints->SetPoint(iLayer, refP[iLayer], x0[iLayer]);
    for(Int_t icls = 0; icls < ncls; icls++)
      trGraph->SetPoint(icls, clp[icls], clx[icls]);
    TCanvas *c1 = new TCanvas();
    trGraph->SetTitle(Form("Event %d, Candidate %d\n", fgEventNumber, fgCandidateNumber));
    refPoints->SetTitle(Form("Event %d, Candidate %d\n", fgEventNumber, fgCandidateNumber));
    fitfun->SetTitle(Form("Event %d, Candidate %d\n", fgEventNumber, fgCandidateNumber));
    trGraph->Draw("ap");
    refPoints->Draw("psame");
    fitfun->Draw("lpsame");
    return c1;
  }
  else{
    AliError(Form("Combination consisting of event %d and candidate %d not found", event, candidate));
    delete fitfun;
    return 0x0;
  }
}

