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

#include "AliLog.h"
#include "AliTRDgeometry.h"
#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"
#include "AliTRDseed.h"
#include "AliTRDcluster.h"

ClassImp(AliTRDtrackerDebug)

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
void AliTRDtrackerDebug::Draw(const Option_t *)
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
	AliTRDcluster *c = 0x0;
	AliTRDseedV1 tracklet[6] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
	const AliTRDseedV1 *ctracklet = 0x0;
	for(int ip = 0; ip<6; ip++){
		if(!(ctracklet = fTrack->GetTracklet(ip))) continue;
		tracklet[ip] = (*ctracklet); 
		Double_t x0 = tracklet[ip].GetX0();
		for(int ic=0; ic<35/*AliTRDseed:knTimebins*/; ic++){
			if(!(c = tracklet[ip].GetClusters(ic))) continue;
			Double_t xc = c->GetX(), yc = c->GetY(), zc = c->GetZ();
			tracklet[ip].SetX(ic, xc-x0);
			tracklet[ip].SetY(ic, yc);
			tracklet[ip].SetZ(ic, zc);
		}
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
		AliTRDseedV1::FitRiemanTilt(&tracklet[0], kTRUE);

		// make a copy of the fit result
		Double_t 
			y0   = tracklet[ip].GetYref(0),
			dydx = tracklet[ip].GetYref(1),
			z0   = tracklet[ip].GetZref(0),
			dzdx = tracklet[ip].GetZref(1);

		// restore tracklet
		tracklet[ip] = (*fTrack->GetTracklet(ip)); 
		for(int ic=0; ic<35/*AliTRDseed:knTimebins*/; ic++){
			if(!(c = tracklet[ip].GetClusters(ic))) continue;
			Double_t xc = c->GetX(), yc = c->GetY(), zc = c->GetZ();
			tracklet[ip].SetX(ic, xc-x0);
			tracklet[ip].SetY(ic, yc);
			tracklet[ip].SetZ(ic, zc);
		}		
		
		// fit clusters
		AliTRDseedV1 ts(tracklet[ip]);
		ts.SetYref(0, y0); ts.SetYref(1, dydx);
		ts.SetZref(0, z0); ts.SetZref(1, dzdx);
		ts.Update();

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

