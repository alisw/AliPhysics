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

///////////////////////////////////////////////////////////////////////////////
//
//
//  TRD calibration class for building reference data for PID
//  - 2D reference histograms (responsible A.Bercuci) 
//  - 3D reference histograms (not yet implemented) (responsible A.Bercuci)
//  - Neural Network (responsible A.Wilk)
//
//   Origin
//   Alex Bercuci  (A.Bercuci@gsi.de)
//
///////////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH2I.h>
#include <TH3D.h>
#include <TParticle.h>
#include <TParticle.h>
#include <TPrincipal.h>
#include <TVector3.h>
#include <TLinearFitter.h>
#include <TVectorT.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TMarker.h>

#include "AliLog.h"
#include "AliPID.h"
#include "AliESD.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliESDtrack.h"

#include "AliTRDCalPIDRefMaker.h"
#include "AliTRDCalPID.h"
#include "AliTRDcalibDB.h"
#include "AliTRDgeometry.h"
#include "AliTRDtrack.h"

#include <vector>

ClassImp(AliTRDCalPIDRefMaker)

TLinearFitter *AliTRDCalPIDRefMaker::fFitter2D2 = 0x0;
TLinearFitter *AliTRDCalPIDRefMaker::fFitter2D1 = 0x0;

//__________________________________________________________________
AliTRDCalPIDRefMaker::AliTRDCalPIDRefMaker()
  :TObject()
{
  //
  // AliTRDCalPIDRefMaker default constructor
  //

	// histogram settings
	for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
		fH2dEdx[ispec] = 0x0;
		fPrinc[ispec] = 0x0;
	}
}

//__________________________________________________________________
AliTRDCalPIDRefMaker::AliTRDCalPIDRefMaker(const AliTRDCalPIDRefMaker &ref)
  :TObject()
{
  //
  // AliTRDCalPIDRefMaker copy constructor
  // 

	for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
		if(ref.fH2dEdx[ispec]){
		       fH2dEdx[ispec] = new  TH2D((TH2D&)*(ref.fH2dEdx[ispec]));
		} else fH2dEdx[ispec] = 0x0;
		fPrinc[ispec] = 0x0;
	}
}

//__________________________________________________________________
AliTRDCalPIDRefMaker::~AliTRDCalPIDRefMaker()
{
  //
  // AliTRDCalPIDQRef destructor
  //

	for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
		if(fH2dEdx[ispec]) delete fH2dEdx[ispec];
		if(fPrinc[ispec]) delete fPrinc[ispec]; 
	}	
	if(fFitter2D1){ delete fFitter2D1; fFitter2D1 = 0x0;}
	if(fFitter2D2){ delete fFitter2D2; fFitter2D2 = 0x0;}

}

//__________________________________________________________________
AliTRDCalPIDRefMaker& AliTRDCalPIDRefMaker::operator=(const AliTRDCalPIDRefMaker &ref)
{
  //
  // Assignment operator
  //

	for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
		if(ref.fH2dEdx[ispec]) (ref.fH2dEdx[ispec])->Copy(*fH2dEdx[ispec]);
		fPrinc[ispec] = 0x0;
	}
	return *this;
}


//__________________________________________________________________
Bool_t AliTRDCalPIDRefMaker::BuildLQReferences(const Char_t *File, const Char_t *dir)
{
	// Build, Fill and write to file the histograms used for PID.
	// The simulations are looked in the
	// directories with the general form Form("p%3.1f", momentum)
	// starting from dir (default .). Here momentum belongs to the list
	// of known momentum to PID (see trackMomentum).
	// The output histograms are
	// written to the file "File" in cwd (default
	// TRDPIDHistograms.root). In order to build a DB entry
	// consider running $ALICE_ROOT/Cal/AliTRDpidCDB.C
	// 
	// Author:
	// Alex Bercuci (A.Bercuci@gsi.de)

  Int_t partCode[AliPID::kSPECIES] =
    {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton};
		
	// check and retrive number of directories in the production
	Int_t nBatches;
	if(!(nBatches = CheckProdDirTree(dir))) return kFALSE;

			
	// Number of Time bins
	Int_t nTimeBins;
	AliTRDcalibDB *calibration = AliTRDcalibDB::Instance();
	if(!calibration){
	  AliError("No AliTRDcalibDB available.");
	  return kFALSE;
	} else {
		nTimeBins = calibration->GetNumberOfTimeBins();
// 		if (calibration->GetRun() > -1) nTimeBins = calibration->GetNumberOfTimeBins();
// 		else {
// 			AliError("No run number set.");
// 	    return kFALSE;
// 	  }
	}

	// Build PID reference/working objects
	fH1TB[0] = new TH1F("h1TB_0", "", nTimeBins, -.5, nTimeBins-.5);
	fH1TB[0]->SetLineColor(4);
	fH1TB[1] = new TH1F("h1TB_1", "", nTimeBins, -.5, nTimeBins-.5);
	fH1TB[1]->SetLineColor(2);
	for(Int_t ispec=0; ispec<AliPID::kSPECIES; ispec++) if(!fPrinc[ispec]) fPrinc[ispec] = new TPrincipal(2, "ND");

	// Print statistics header
	Int_t nPart[AliPID::kSPECIES], nTotPart;
	printf("P[GeV/c] ");
	for(Int_t ispec=0; ispec<AliPID::kSPECIES; ispec++) printf(" %s[%%] ", AliTRDCalPID::GetPartSymb(ispec));
	printf("\n-----------------------------------------------\n");
	
	Float_t trackMomentum[AliTRDCalPID::kNMom];
        //Float_t trackSegLength[AliTRDCalPID::kNLength];
	for(int i=0; i<AliTRDCalPID::kNMom; i++) trackMomentum[i] = AliTRDCalPID::GetMomentum(i);
	AliRunLoader *fRunLoader = 0x0;
	TFile *esdFile = 0x0;
	TTree *esdTree = 0x0;
	AliESD *esd = 0x0;
	AliESDtrack *esdTrack = 0x0;
	
	//
	// Momentum loop
	for (Int_t imom = 0; imom < AliTRDCalPID::kNMom; imom++) {
		Reset();
		
		// init statistics for momentum
		nTotPart = 0;
		for(Int_t ispec=0; ispec<AliPID::kSPECIES; ispec++) nPart[ispec] = 0;

		// loop over production directories
		for(Int_t ibatch = 0; ibatch<nBatches; ibatch++){
			// open run loader and load gAlice, kinematics and header
			fRunLoader = AliRunLoader::Open(Form("%s/%3.1fGeV/%03d/galice.root", dir, trackMomentum[imom], ibatch));
			if (!fRunLoader) {
				AliError(Form("Getting run loader for momentum %3.1f GeV/c batch %03d failed.", trackMomentum[imom], ibatch));
				return kFALSE;
			}
			TString s; s.Form("%s/%3.1fGeV/%03d/", dir, trackMomentum[imom], ibatch);
			fRunLoader->SetDirName(s);
			fRunLoader->LoadgAlice();
			gAlice = fRunLoader->GetAliRun();
			if (!gAlice) {
				AliError(Form("galice object not found for momentum %3.1f GeV/c batch %03d.", trackMomentum[imom], ibatch));
				return kFALSE;
			}
			fRunLoader->LoadKinematics();
			fRunLoader->LoadHeader();
	
			// open the ESD file
			esdFile = TFile::Open(Form("%s/%3.1fGeV/%03d/AliESDs.root", dir, trackMomentum[imom], ibatch));
			if (!esdFile || esdFile->IsZombie()) {
				AliError(Form("Opening ESD file failed for momentum  %3.1f GeV/c batch %03d.", trackMomentum[imom], ibatch));
				return kFALSE;
			}
			esdTree = (TTree*)esdFile->Get("esdTree");
			if (!esdTree) {
				AliError(Form("ESD tree not found for momentum %3.1f GeV/c batch %03d.", trackMomentum[imom], ibatch));
				return kFALSE;
			}
			esd = new AliESD;
			esdTree->SetBranchAddress("ESD", &esd);
	
			
			// Event loop
			for (Int_t iEvent = 0; iEvent < fRunLoader->GetNumberOfEvents(); iEvent++) {
				
				// Load MC info
				fRunLoader->GetEvent(iEvent);
				AliStack* stack = AliRunLoader::GetRunLoader()->Stack();
				TArrayF vertex(3);
				fRunLoader->GetHeader()->GenEventHeader()->PrimaryVertex(vertex);
							
				// Load event summary data
				esdTree->GetEvent(iEvent);
				if (!esd) {
					AliWarning(Form("ESD object not found for event %d. (@ momentum %3.1f GeV/c, batch %03d)", iEvent, trackMomentum[imom], ibatch));
					continue;
				}
	
				// Track loop
				for(Int_t iTrack=0; iTrack<esd->GetNumberOfTracks(); iTrack++){
					esdTrack = esd->GetTrack(iTrack);
	
					//if(!AliTRDpidESD::CheckTrack(esdTrack)) continue;

					if((esdTrack->GetStatus() & AliESDtrack::kITSrefit) == 0) continue;
					if(esdTrack->GetConstrainedChi2() > 1E9) continue;
					if ((esdTrack->GetStatus() & AliESDtrack::kESDpid) == 0) continue;
					if (esdTrack->GetTRDsignal() == 0.) continue;
	
					// read MC info
					Int_t label = esdTrack->GetLabel();
					if(label<0) continue;
					if (label > stack->GetNtrack()) continue;     // background
					TParticle* particle = stack->Particle(label);
					if(!particle){
						AliWarning(Form("Retriving particle with index %d from AliStack failed. [@ momentum %3.1f batch %03d event %d track %d]", label, trackMomentum[imom], ibatch, iEvent, iTrack));
						continue;
					}
					if(particle->Pt() < 1.E-3) continue;
					//      if (TMath::Abs(particle->Eta()) > 0.3) continue;
					TVector3 dVertex(particle->Vx() - vertex[0],
										particle->Vy() - vertex[1],
										particle->Vz() - vertex[2]);
					if (dVertex.Mag() > 1.E-4){
						//AliInfo(Form("Particle with index %d generated too far from vertex. Skip from analysis. Details follows. [@ event %d track %d]", label, iEvent, iTrack));
						//particle->Print();
						continue;
					}
					Int_t iGen = -1;
					for (Int_t ispec=0; ispec<AliPID::kSPECIES; ispec++)
						if(TMath::Abs(particle->GetPdgCode()) == partCode[ispec]){
							iGen = ispec;
							break;
						}
					if(iGen<0) continue;
	
					nPart[iGen]++; nTotPart++;
					
					Float_t mom;
                                        //Float_t length;
					Double_t dedx[AliTRDtrack::kNslice], dEdx;
					Int_t timebin;
					for (Int_t iLayer=0; iLayer<AliTRDgeometry::kNlayer; iLayer++){
						// read data for track segment
						for(int iSlice=0; iSlice<AliTRDtrack::kNslice; iSlice++)
							dedx[iSlice] = esdTrack->GetTRDslice(iLayer, iSlice);
						dEdx    = esdTrack->GetTRDslice(iLayer, -1);
						timebin = esdTrack->GetTRDTimBin(iLayer);
			
						// check data
						if ((dEdx <=  0.) || (timebin <= -1.)) continue;
			
						// retrive kinematic info for this track segment
						//if(!AliTRDpidESD::RecalculateTrackSegmentKine(esdTrack, iLayer, mom, length)) continue;
						mom = esdTrack->GetOuterParam()->GetP();
						
						// find segment length and momentum bin
						Int_t jmom = 1, refMom = -1;
						while(jmom<AliTRDCalPID::kNMom-1 && mom>trackMomentum[jmom]) jmom++;
						if(TMath::Abs(trackMomentum[jmom-1] - mom) < trackMomentum[jmom-1] * .2) refMom = jmom-1;
						else if(TMath::Abs(trackMomentum[jmom] - mom) < trackMomentum[jmom] * .2) refMom = jmom;
						if(refMom<0){
							AliInfo(Form("Momentum at plane %d entrance not in momentum window. [@ momentum %3.1f batch %03d event %d track %d]", iLayer, trackMomentum[imom], ibatch, iEvent, iTrack));
							continue;
						}
						/*while(jleng<AliTRDCalPID::kNLength-1 && length>trackSegLength[jleng]) jleng++;*/
						
						// this track segment has fulfilled all requierments
						//nPlanePID++;

						if(dedx[0] > 0. && dedx[1] > 0.){
							dedx[0] = log(dedx[0]); dedx[1] = log(dedx[1]);
							fPrinc[iGen]->AddRow(dedx);
						}
						fH1TB[(iGen>0)?1:0]->Fill(timebin);
					} // end plane loop
				} // end track loop
			} // end events loop
			
			delete esd; esd = 0x0;
			esdFile->Close();
			delete esdFile; esdFile = 0x0;
	
			fRunLoader->UnloadHeader();
			fRunLoader->UnloadKinematics();
			delete fRunLoader; fRunLoader = 0x0;
		} // end directory loop
		
		// use data to prepare references
		Prepare2D();
		// save references
		SaveReferences(imom, File);

			
		// print momentum statistics
		printf("  %3.1f  ", trackMomentum[imom]);
		for(Int_t ispec=0; ispec<AliPID::kSPECIES; ispec++) printf(" %5.2f ", 100.*nPart[ispec]/nTotPart);
		printf("\n");
	} // end momentum loop
	
	TFile *fSave = 0x0;
	TListIter it((TList*)gROOT->GetListOfFiles());
	while((fSave=(TFile*)it.Next()))
		if(strcmp(File, fSave->GetName())==0) break;

	fSave->cd();
	fSave->Close();
	delete fSave;

	return kTRUE;
}

//__________________________________________________________________
Bool_t AliTRDCalPIDRefMaker::BuildNNReferences(const Char_t* /*File*/, const Char_t* /*dir*/)
{
	return kTRUE;
}

//__________________________________________________________________
Double_t AliTRDCalPIDRefMaker::Estimate2D2(TH2 *h, Float_t &x, Float_t &y)
{
  //
  // Linear interpolation of data point with a parabolic expresion using
  // the logarithm of the bin content in a close neighbourhood. It is
  // assumed that the bin content of h is in number of events !
  //
  // Observation:
  // This function have to be used when the statistics of each bin is
  // sufficient and all bins are populated. For cases of sparse data
  // please refere to Estimate2D1().
  //
  // Author : Alex Bercuci (A.Bercuci@gsi.de)
  //

	if(!h){
		AliErrorGeneral("AliTRDCalPIDRefMaker::Estimate2D2()", "No histogram defined.");
		return 0.;
	}

	TAxis *ax = h->GetXaxis(), *ay = h->GetYaxis();
	Int_t binx   = ax->FindBin(x);
	Int_t biny   = ay->FindBin(y);
	Int_t nbinsx = ax->GetNbins();
	Int_t nbinsy = ay->GetNbins();
	Double_t p[2];
	Double_t entries;

	// late construction of fitter
	if(!fFitter2D2) fFitter2D2 = new TLinearFitter(6, "1++x++y++x*x++y*y++x*y");
		
	fFitter2D2->ClearPoints();
	Int_t npoints=0;
	Int_t binx0, binx1, biny0, biny1;
	for(int bin=0; bin<5; bin++){
		binx0 = TMath::Max(1, binx-bin);
		binx1 = TMath::Min(nbinsx, binx+bin);
		for(int ibin=binx0; ibin<=binx1; ibin++){
			biny0 = TMath::Max(1, biny-bin);
			biny1 = TMath::Min(nbinsy, biny+bin);
			for(int jbin=biny0; jbin<=biny1; jbin++){
				if(ibin != binx0 && ibin != binx1 && jbin != biny0 && jbin != biny1) continue;
				if((entries = h->GetBinContent(ibin, jbin)) == 0.) continue;
				p[0] = ax->GetBinCenter(ibin);
				p[1] = ay->GetBinCenter(jbin);
				fFitter2D2->AddPoint(p, log(entries), 1./sqrt(entries));
				npoints++;
			}
		}
		if(npoints>=25) break;
	}
	if(fFitter2D2->Eval() == 1){
		printf("<I2> x = %9.4f y = %9.4f\n", x, y);
		printf("\tbinx %d biny %d\n", binx, biny);
		printf("\tpoints %d\n", npoints);

		return 0.;
	}
	TVectorD vec(6);
	fFitter2D2->GetParameters(vec);
	Double_t result = vec[0] + x*vec[1] + y*vec[2] + x*x*vec[3] + y*y*vec[4] + x*y*vec[5];
	return exp(result);

}

//__________________________________________________________________
Double_t AliTRDCalPIDRefMaker::Estimate2D1(TH2 *h, Float_t &x, Float_t &y
                                      , Float_t &dCT, Float_t &rmin
                                      , Float_t &rmax)
{
  //
  // Linear interpolation of data point with a plane using
  // the logarithm of the bin content in the area defined by the
  // d(cos(phi)) and dr=(rmin, rmax). It is assumed that the bin content
  // of h is number of events !
  //

	if(!h){
		AliErrorGeneral("AliTRDCalPIDRefMaker::Estimate2D1()", "No histogram defined.");
		return 0.;
	}

	TAxis *ax = h->GetXaxis(), *ay = h->GetYaxis();
// 	Int_t binx   = ax->FindBin(x);
// 	Int_t biny   = ay->FindBin(y);
	Int_t nbinsx = ax->GetNbins();
	Int_t nbinsy = ay->GetNbins();
	Double_t p[2];
	Double_t entries;
	Double_t rxy = sqrt(x*x + y*y), rpxy;

	// late construction of fitter	
	if(!fFitter2D1) fFitter2D1 = new TLinearFitter(3, "1++x++y");
	
	fFitter2D1->ClearPoints();
	Int_t npoints=0;
	for(int ibin=1; ibin<=nbinsx; ibin++){
		for(int jbin=1; jbin<=nbinsy; jbin++){
			if((entries = h->GetBinContent(ibin, jbin)) == 0.) continue;
			p[0] = ax->GetBinCenter(ibin);
			p[1] = ay->GetBinCenter(jbin);
			rpxy = sqrt(p[0]*p[0] + p[1]*p[1]);
			if((x*p[0] + y*p[1])/rxy/rpxy < dCT) continue;
			if(rpxy<rmin || rpxy > rmax) continue;
			
			fFitter2D1->AddPoint(p, log(entries), 1./sqrt(entries));
			npoints++;
		}
	}
	if(npoints<15) return 0.;
	if(fFitter2D1->Eval() == 1){
		printf("<O2> x = %9.4f y = %9.4f\n", x, y);
		printf("\tpoints %d\n", npoints);
		return 0.;
	}
	TVectorD vec(3);
	fFitter2D1->GetParameters(vec);
	Double_t result = vec[0] + x*vec[1] + y*vec[2];
	return exp(result);
}

//__________________________________________________________________
// Double_t	AliTRDCalPIDRefMaker::Estimate3D2(TH3 *h, Float_t &x, Float_t &y, Float_t &z)
// {
// 	// Author Alex Bercuci (A.Bercuci@gsi.de)
// 	return 0.;
// }



/////////////  Private functions ///////////////////////////////////

//__________________________________________________________________
Int_t AliTRDCalPIDRefMaker::CheckProdDirTree(const Char_t *dir)
{
	// Scan directory tree for momenta. Returns the smallest number of
	// batches found in all directories or 0 if one momentum is missing.

	const char *pwd = gSystem->pwd();

	if(!gSystem->ChangeDirectory(dir)){
		AliError(Form("Couldn't access production root directory %s.", dir));
		return 0;
	}

	Int_t iDir;
	Int_t nDir = Int_t(1.E6);
	for(int imom=0; imom<AliTRDCalPID::kNMom; imom++){
		if(!gSystem->ChangeDirectory(Form("%3.1fGeV", AliTRDCalPID::GetMomentum(imom)))){
			AliError(Form("Couldn't find data for momentum %3.1f GeV/c.", AliTRDCalPID::GetMomentum(imom)));
			return 0;	
		}
		
		iDir = 0;
		while(gSystem->ChangeDirectory(Form("%03d", iDir))){
			iDir++;
			gSystem->ChangeDirectory("..");
		}
		if(iDir < nDir) nDir = iDir;
		gSystem->ChangeDirectory(dir);
	}

	gSystem->ChangeDirectory(pwd);

	return nDir;
}


//__________________________________________________________________
void  AliTRDCalPIDRefMaker::Prepare2D()
{
  //
  // Prepares the 2-dimensional reference histograms
  //
	
	// histogram settings
	Float_t xmin, xmax, ymin, ymax;
	Int_t nbinsx, nbinsy, nBinsSector, nSectors;
	const Int_t color[] = {4, 3, 2, 7, 6};
	for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
		// check PCA data
		if(!fPrinc[ispec]){
			AliError(Form("No data defined for %s.", AliTRDCalPID::GetPartName(ispec)));
			return;
		}
		// build reference histograms
		nBinsSector = 10;
		nSectors = 20;
		nbinsx = nbinsy = nBinsSector * nSectors;
		xmin = ymin = 0.;
		xmax = 8000.; ymax = 6000.;
		if(!fH2dEdx[ispec]){
			fH2dEdx[ispec] = new  TH2D(Form("h2%s", AliTRDCalPID::GetPartSymb(ispec)), "", nbinsx, xmin, xmax, nbinsy, ymin, ymax);
			fH2dEdx[ispec]->SetLineColor(color[ispec]);
		}
	}

	// build transformed rotated histograms
	nbinsx = nbinsy = 100;
	xmin = ymin = -6.;
	xmax = 10.; ymax = 6.;
	TH2I *hProj = 0x0;
	if(!(hProj = (TH2I*)gDirectory->Get("hProj"))) hProj = new  TH2I("hProj", "", nbinsx, xmin, xmax, nbinsy, ymin, ymax);
	else hProj->Reset();

	printf("Doing interpolation and invertion ... "); fflush(stdout);
	Bool_t kDebugPlot = kTRUE;
	//Bool_t kDebugPrint = kFALSE;

	TCanvas *c = 0x0;
	TEllipse *ellipse = 0x0;
	TMarker *mark = 0x0;
	if(kDebugPlot){
		c=new TCanvas("c2", "Interpolation 2D", 10, 10, 500, 500);
		ellipse = new TEllipse();
		ellipse->SetFillStyle(0); ellipse->SetLineColor(2);
		mark = new TMarker();
		mark->SetMarkerColor(2); mark->SetMarkerSize(2); mark->SetMarkerStyle(2);
	}
	
	// define observable variables
	TAxis *ax, *ay;
	Double_t xy[2], lxy[2], rxy[2];
	Double_t estimate, position;
	const TVectorD *eValues;
	
	// define radial sectors
	const Int_t   nPhi      = 36;
	const Float_t dPhi      = TMath::TwoPi()/nPhi;
	//const Float_t dPhiRange = .1;
	Int_t nPoints[nPhi], nFitPoints, binStart, binStop;
	TLinearFitter refsFitter[nPhi], refsLongFitter(6, "1++x++y++x*x++y*y++x*y");
	Float_t fFitterRange[nPhi];
	Bool_t kFitterStatus[nPhi];
	for(int iphi=0; iphi<nPhi; iphi++){
		refsFitter[iphi].SetDim(3);
		refsFitter[iphi].SetFormula("1++x++y");//++x*x++y*y++x*y");
		fFitterRange[iphi] = .8;
		kFitterStatus[iphi] = kFALSE;
	}
	std::vector<UShort_t> storeX[nPhi], storeY[nPhi];

	// define working variables
	Float_t x0, y0, rx, ry;
	//Float_t rc, rmin, rmax, dr, dCT;
	Double_t Phi, r;
	Int_t iPhi;
	Double_t entries;
	for(int ispec=0; ispec<5; ispec++){
		hProj->Reset();
		for(int iphi=0; iphi<nPhi; iphi++){
			nPoints[iphi] = 0;
			refsFitter[iphi].ClearPoints();
			fFitterRange[iphi] = .8;
			kFitterStatus[iphi] = kFALSE;
			storeX[iphi].clear();
			storeY[iphi].clear();
		}
		
		// calculate covariance ellipse
		fPrinc[ispec]->MakePrincipals();
		eValues  = fPrinc[ispec]->GetEigenValues();
		x0  = 0.;
		y0  = 0.;
		rx  = 3.5*sqrt((*eValues)[0]);
		ry  = 3.5*sqrt((*eValues)[1]);

		// rotate to principal axis
		Int_t irow = 0;
		const Double_t *xx;
		printf("filling for spec %d ...\n", ispec);
		while((xx = fPrinc[ispec]->GetRow(irow++))){
			fPrinc[ispec]->X2P(xx, rxy);
			hProj->Fill(rxy[0], rxy[1]);
		}
		
		
//  	// debug plot
// 		if(kDebugPlot){
// 			hProj->Draw();
// 			ellipse->DrawEllipse(x0, y0, rx, ry, 0., 360., 0.);
// 			mark->DrawMarker(x0, y0);
// 			gPad->Modified(); gPad->Update();
// 		}
				
		// define radial sectors
		ax=hProj->GetXaxis();
		ay=hProj->GetYaxis();
		for(int ibin=1; ibin<=ax->GetNbins(); ibin++){
			rxy[0] = ax->GetBinCenter(ibin);
			for(int jbin=1; jbin<=ay->GetNbins(); jbin++){
				rxy[1] = ay->GetBinCenter(jbin);

				if((entries = hProj->GetBinContent(ibin, jbin)) == 0) continue;

				position = rxy[0]*rxy[0]/rx/rx + rxy[1]*rxy[1]/ry/ry;
				if(position < 1.) continue;
				
				r = sqrt(rxy[0]*rxy[0] + rxy[1]*rxy[1]);
				Phi   = ((rxy[1] > 0.) ? 1. : -1.) * TMath::ACos(rxy[0]/r); // [-pi, pi]
				iPhi = nPhi/2 + Int_t(Phi/dPhi) - ((Phi/dPhi > 0.) ? 0 : 1);
				
				refsFitter[iPhi].AddPoint(rxy, log(entries), 1./sqrt(entries));
				nPoints[iPhi]++;
			}
		}
		
		// do interpolation
		ax=fH2dEdx[ispec]->GetXaxis();
		ay=fH2dEdx[ispec]->GetYaxis();
		for(int ibin=1; ibin<=ax->GetNbins(); ibin++){
			xy[0]  = ax->GetBinCenter(ibin);
			lxy[0] = log(xy[0]);
			for(int jbin=1; jbin<=ay->GetNbins(); jbin++){
				xy[1]  = ay->GetBinCenter(jbin);
				lxy[1] = log(xy[1]);

				// rotate to PCA
				fPrinc[ispec]->X2P(lxy, rxy);

				// calculate border of covariance ellipse
				position = rxy[0]*rxy[0]/rx/rx + rxy[1]*rxy[1]/ry/ry;

				// interpolation inside the covariance ellipse
				if(position < 1.){
				  Float_t xTemp = rxy[0];
                                  Float_t yTemp = rxy[1];
					estimate = Estimate2D2((TH2 *) hProj, xTemp, yTemp); 
					rxy[0] = xTemp;
                                        rxy[1] = yTemp;
					//					estimate = Estimate2D2((TH2 *) hProj, (Float_t)rxy[0], (Float_t)rxy[1]);
					fH2dEdx[ispec]->SetBinContent(ibin, jbin, estimate/xy[0]/xy[1]);
				} else { // interpolation outside the covariance ellipse
					r = sqrt(rxy[0]*rxy[0] + rxy[1]*rxy[1]);
					Phi   = ((rxy[1] > 0.) ? 1. : -1.) * TMath::ACos(rxy[0]/r); // [-pi, pi]
					iPhi = nPhi/2 + Int_t(Phi/dPhi) - ((Phi/dPhi > 0.) ? 0 : 1);
	
					storeX[iPhi].push_back(ibin);
					storeY[iPhi].push_back(jbin);
				}
			}
		}
		
		// Fill outliers
		// Radial fit on transformed rotated
		TVectorD vec(3);
		Int_t xbin, ybin;
		for(int iphi=0; iphi<nPhi; iphi++){
			Phi = iphi * dPhi - TMath::Pi();
			if(TMath::Abs(TMath::Abs(Phi)-TMath::Pi()) < 100.*TMath::DegToRad()) continue;
				
			
			refsFitter[iphi].Eval();
			refsFitter[iphi].GetParameters(vec);
			for(UInt_t ipoint=0; ipoint<storeX[iphi].size(); ipoint++){
				xbin = storeX[iphi].at(ipoint);
				ybin = storeY[iphi].at(ipoint);
				xy[0] = ax->GetBinCenter(xbin); lxy[0] = log(xy[0]);
				xy[1] = ay->GetBinCenter(ybin); lxy[1] = log(xy[1]);

				// rotate to PCA
				fPrinc[ispec]->X2P(lxy, rxy);
				estimate = exp(vec[0] + rxy[0]*vec[1] + rxy[1]*vec[2]);//+ rxy[0]*rxy[0]*vec[3]+ rxy[1]*rxy[1]*vec[4]+ rxy[0]*rxy[1]*vec[5]);
				fH2dEdx[ispec]->SetBinContent(xbin, ybin, estimate/xy[0]/xy[1]);
			}
		}
		
		// Longitudinal fit on ref histo in y direction
		for(int isector=1; isector<nSectors; isector++){
			// define sectors
			binStart = ((isector>0)?isector-1:isector)*nBinsSector + 1;
			binStop  = (isector+1)*nBinsSector;
			
			nPoints[0] = 0;
			refsLongFitter.ClearPoints();
			storeX[0].clear();
			storeY[0].clear();
	
			for(int ibin=1; ibin<=ax->GetNbins(); ibin++){
				xy[0] = ax->GetBinCenter(ibin);
				for(int jbin=binStart; jbin<=binStop; jbin++){
					xy[1] = ay->GetBinCenter(jbin);
					if((entries = fH2dEdx[ispec]->GetBinContent(ibin, jbin)) > 0.){
						refsLongFitter.AddPoint(xy, log(entries), 1.e-7);
						nPoints[0]++;
					} else {
						storeX[0].push_back(ibin);
						storeY[0].push_back(jbin);
					}
				}
				if(nPoints[0] >= ((isector==0)?60:420)) break;
			}
	
			
			if(refsLongFitter.Eval() == 1){
				printf("Error in Y sector %d\n", isector);
				continue;
			}
			refsLongFitter.GetParameters(vec);
			for(UInt_t ipoint=0; ipoint<storeX[0].size(); ipoint++){
				xbin = storeX[0].at(ipoint);
				ybin = storeY[0].at(ipoint);
				xy[0] = ax->GetBinCenter(xbin);
				xy[1] = ay->GetBinCenter(ybin);
				
				estimate = vec[0] + xy[0]*vec[1] + xy[1]*vec[2] + xy[0]*xy[0]*vec[3]+ xy[1]*xy[1]*vec[4]+ xy[0]*xy[1]*vec[5];
				fH2dEdx[ispec]->SetBinContent(xbin, ybin, exp(estimate));
			}
		}

		// Longitudinal fit on ref histo in x direction
		for(int isector=1; isector<nSectors; isector++){
			binStart = ((isector>0)?isector-1:isector)*nBinsSector + 1;
			binStop  = (isector+1)*nBinsSector;
			
			nPoints[0] = 0;
			nFitPoints = 0;
			refsLongFitter.ClearPoints();
			storeX[0].clear();
			storeY[0].clear();
	
			for(int jbin=1; jbin<=ay->GetNbins(); jbin++){
				xy[1] = ay->GetBinCenter(jbin);
				for(int ibin=binStart; ibin<=binStop; ibin++){
					xy[0] = ax->GetBinCenter(ibin);
					if((entries = fH2dEdx[ispec]->GetBinContent(ibin, jbin)) > 0.){
						refsLongFitter.AddPoint(xy, log(entries), 1.e-7);
						nPoints[0]++;
					} else {
						storeX[0].push_back(ibin);
						storeY[0].push_back(jbin);
						nFitPoints++;
					}
				}
				if(nPoints[0] >= ((isector==0)?60:420)) break;
			}
			if(nFitPoints == 0) continue;
	
			
			if(refsLongFitter.Eval() == 1){
				printf("Error in X sector %d\n", isector);
				continue;
			}
			refsLongFitter.GetParameters(vec);
			for(UInt_t ipoint=0; ipoint<storeX[0].size(); ipoint++){
				xbin = storeX[0].at(ipoint);
				ybin = storeY[0].at(ipoint);
				xy[0] = ax->GetBinCenter(xbin);
				xy[1] = ay->GetBinCenter(ybin);
				
				estimate = vec[0] + xy[0]*vec[1] + xy[1]*vec[2] + xy[0]*xy[0]*vec[3]+ xy[1]*xy[1]*vec[4]+ xy[0]*xy[1]*vec[5];
				fH2dEdx[ispec]->SetBinContent(xbin, ybin, exp(estimate));
			}
		}

		
		fH2dEdx[ispec]->Scale(1./fH2dEdx[ispec]->Integral());
		
		// debug plot
		if(kDebugPlot){
			fH2dEdx[ispec]->Draw("cont1"); gPad->SetLogz();
			gPad->Modified(); gPad->Update();
		}
	}
	AliInfo("Finish interpolation.");
}

//__________________________________________________________________
void	AliTRDCalPIDRefMaker::Reset()
{
  //
  // Reset reference histograms
  //

	for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
		if(fH2dEdx[ispec]) fH2dEdx[ispec]->Reset();
		fPrinc[ispec]->Clear();
	}	
	if(fH1TB[0]) fH1TB[0]->Reset(); 
	if(fH1TB[1]) fH1TB[1]->Reset();
}

//__________________________________________________________________
void  AliTRDCalPIDRefMaker::SaveReferences(const Int_t mom, const char *fn)
{
  //
  // Save the reference histograms
  //

	TFile *fSave = 0x0;
	TListIter it((TList*)gROOT->GetListOfFiles());
	Bool_t kFOUND = kFALSE;
	TDirectory *pwd = gDirectory;
	while((fSave=(TFile*)it.Next()))
		if(strcmp(fn, fSave->GetName())==0){
			kFOUND = kTRUE;
			break;
		}
	if(!kFOUND) fSave = TFile::Open(fn, "RECREATE");
	fSave->cd();

	Float_t fmom = AliTRDCalPID::GetMomentum(mom);
	
	// save dE/dx references
	TH2 *h2 = 0x0;
	for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
		h2 = (TH2D*)fH2dEdx[ispec]->Clone(Form("h2dEdx%s%d", AliTRDCalPID::GetPartSymb(ispec), mom));
		h2->SetTitle(Form("2D dEdx for particle %s @ %d", AliTRDCalPID::GetPartName(ispec), mom));
		h2->GetXaxis()->SetTitle("dE/dx_{TRD}^{amplif} [au]");
		h2->GetYaxis()->SetTitle("dE/dx_{TRD}^{drift} [au]");
		h2->GetZaxis()->SetTitle("Entries");
		h2->Write();
	}

	// save maximum time bin references 
	TH1 *h1 = 0x0;
	h1 = (TH1F*)fH1TB[0]->Clone(Form("h1MaxTBEL%02d", mom));
	h1->SetTitle(Form("Maximum Time Bin distribution for electrons @ %4.1f GeV", fmom));
	h1->GetXaxis()->SetTitle("time [100 ns]");
	h1->GetYaxis()->SetTitle("Probability");
	h1->Write();

	h1 = (TH1F*)fH1TB[1]->Clone(Form("h1MaxTBPI%02d", mom));
	h1->SetTitle(Form("Maximum Time Bin distribution for pions @ %4.1f GeV", fmom));
	h1->GetXaxis()->SetTitle("time [100 ns]");
	h1->GetYaxis()->SetTitle("Probability");
	h1->Write();


	pwd->cd();
}

