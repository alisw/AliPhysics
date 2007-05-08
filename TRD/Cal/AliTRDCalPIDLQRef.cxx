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
//                                                                           //
//  TRD calibration class for 2-dim PID reference histograms                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH2I.h>
#include <TH3D.h>
#include <TPrincipal.h>
#include <TLinearFitter.h>
#include <TVectorT.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TMarker.h>

#include "AliLog.h"
#include "AliPID.h"
#include "AliESD.h"

#include "AliTRDCalPIDLQ.h"
#include "AliTRDCalPIDLQRef.h"

ClassImp(AliTRDCalPIDLQRef)

//__________________________________________________________________
AliTRDCalPIDLQRef::AliTRDCalPIDLQRef() 
  :TObject()
  ,fFitter2D2(0x0)
  ,fFitter2D1(0x0)
{
  //
  // AliTRDCalPIDLQRef default constructor
  //

	// histogram settings
	for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
		h2dEdx[ispec] = 0x0;
		fPrinc[ispec] = 0x0;
	}
	// build parabolic 2D fitter
	fFitter2D2 = new TLinearFitter(6, "1++x++y++x*x++y*y++x*y");
	fFitter2D1 = new TLinearFitter(3, "1++x++y");

}

//__________________________________________________________________
AliTRDCalPIDLQRef::AliTRDCalPIDLQRef(const AliTRDCalPIDLQRef &ref) 
  :TObject()
  ,fFitter2D2(0x0)
  ,fFitter2D1(0x0)
{
  //
  // AliTRDCalPIDLQRef copy constructor
  // 

	for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
		if(ref.h2dEdx[ispec]){
			h2dEdx[ispec] = new  TH2D((TH2D&)*(ref.h2dEdx[ispec]));
		} else h2dEdx[ispec] = 0x0;
		fPrinc[ispec] = 0x0;
	}
	
	// use the argument constructor because copy constructor is not
	// working as for the ROOT version 5.15/07
	fFitter2D2 = new TLinearFitter(6, "1++x++y++x*x++y*y++x*y");
	fFitter2D1 = new TLinearFitter(3, "1++x++y");

}

//__________________________________________________________________
AliTRDCalPIDLQRef::~AliTRDCalPIDLQRef()
{
  //
  // AliTRDCalPIDQRef destructor
  //

	for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
		if(h2dEdx[ispec]) delete h2dEdx[ispec];
		if(fPrinc[ispec]) delete fPrinc[ispec]; 
	}	
	delete fFitter2D1;
	delete fFitter2D2;

}

//__________________________________________________________________
AliTRDCalPIDLQRef& AliTRDCalPIDLQRef::operator=(const AliTRDCalPIDLQRef &ref)
{
  //
  // Assignment operator
  //

	for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
		if(ref.h2dEdx[ispec]) (ref.h2dEdx[ispec])->Copy(*h2dEdx[ispec]);
		fPrinc[ispec] = 0x0;
	}
	return *this;

}

//__________________________________________________________________
Double_t AliTRDCalPIDLQRef::Estimate2D2(TH2 *h, Float_t &x, Float_t &y)
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
		AliError("No histogram defined.");
		return 0.;
	}

	TAxis *ax = h->GetXaxis(), *ay = h->GetYaxis();
	Int_t binx   = ax->FindBin(x);
	Int_t biny   = ay->FindBin(y);
	Int_t nbinsx = ax->GetNbins();
	Int_t nbinsy = ay->GetNbins();
	Double_t p[2];
	Double_t entries;
		
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
Double_t AliTRDCalPIDLQRef::Estimate2D1(TH2 *h, Float_t &x, Float_t &y
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
		AliError("No histogram defined.");
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
// Double_t	AliTRDCalPIDLQRef::Estimate3D2(TH3 *h, Float_t &x, Float_t &y, Float_t &z)
// {
// 	// Author Alex Bercuci (A.Bercuci@gsi.de)
// 	return 0.;
// }

//__________________________________________________________________
void  AliTRDCalPIDLQRef::Prepare2DReferences()
{
  //
  // Prepares the 2-dimensional reference histograms
  //
	
	// histogram settings
	Float_t xmin, xmax, ymin, ymax;
	Int_t nbinsx, nbinsy;
	const Int_t color[] = {4, 3, 2, 7, 6};
	for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
		// build reference histograms
		nbinsx = nbinsy = 500;
		xmin = ymin = 0.;
		xmax = ymax = 500.;
		if(!h2dEdx[ispec]){
			h2dEdx[ispec] = new  TH2D(Form("h2%s", AliTRDCalPIDLQ::fpartSymb[ispec]), "", nbinsx, xmin, xmax, nbinsy, ymin, ymax);
			h2dEdx[ispec]->SetLineColor(color[ispec]);
		}
		// build PCA
		if(fPrinc[ispec]) fPrinc[ispec] = new TPrincipal(2, "ND");
	}

	// build transformed rotated histograms
	nbinsx = nbinsy = 100;
	xmin = ymin = -6.;
	xmax = 10.; ymax = 6.;
	TH2I *hProj = 0x0;
	if(!(hProj = (TH2I*)gDirectory->Get("hProj"))) hProj = new  TH2I("hProj", "", nbinsx, xmin, xmax, nbinsy, ymin, ymax);
	else hProj->Reset();
	
	// build transformed smoothed histogram
	nbinsx = nbinsy = 200;
	xmin = -1.5; xmax = 7.5;
	ymin = -3.; ymax = 7.;
	TH2D *hSmooth = 0x0;
	if(!(hSmooth = (TH2D*)gDirectory->Get("hSmooth"))) hSmooth = new TH2D("hSmooth", "", nbinsx, xmin, xmax, nbinsy, ymin, ymax);
	else hSmooth->Reset();
	

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
	
	TAxis *ax, *ay;
	Double_t xy[2], lxy[2];
	Double_t estimate, position;
	const TVectorD *eValues;
	Float_t x0, y0, rx, ry, rc, rmin, rmax, dr, dCT;
	for(int ispec=0; ispec<5; ispec++){
		hProj->Reset(); hSmooth->Reset();
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
		while((xx = fPrinc[ispec]->GetRow(irow++))){
			fPrinc[ispec]->X2P(xx, lxy);
			hProj->Fill(lxy[0], lxy[1]);
		}

 		// debug plot
		if(kDebugPlot){
			hProj->Draw();
			ellipse->DrawEllipse(x0, y0, rx, ry, 0., 360., 0.);
			mark->DrawMarker(x0, y0);
			gPad->Modified(); gPad->Update();
		}
		
				
		// do interpolation
		ax=hSmooth->GetXaxis();
		ay=hSmooth->GetYaxis();
		for(int ibin=1; ibin<=ax->GetNbins(); ibin++){
			xy[0] = ax->GetBinCenter(ibin); 
			for(int jbin=1; jbin<=ay->GetNbins(); jbin++){
				xy[1] = ay->GetBinCenter(jbin);

				// rotate to PCA
				fPrinc[ispec]->X2P(xy, lxy);

				// calculate border of covariance ellipse
				position = lxy[0]*lxy[0]/rx/rx+lxy[1]*lxy[1]/ry/ry;

				// interpolation inside the covariance ellipse
				if(position < 1.) estimate = Estimate2D2(hProj, (Float_t&)lxy[0], (Float_t&)lxy[1]);
				else{ // interpolation outside the covariance ellipse
					dCT  = .9977;
					dr = .5;
					rc = sqrt(lxy[0]*lxy[0]+lxy[1]*lxy[1]);
					rmin = rc - dr*sqrt(rc);
					rmax = rc + dr/rc;
					while((estimate = Estimate2D1(hProj, (Float_t&)lxy[0], (Float_t&)lxy[1], dCT, rmin, rmax)) < 1.E-8){
						rmin -= dr;
						dCT -= .001; 
					}
				}
				hSmooth->SetBinContent(ibin, jbin, estimate);
			}
		}
		
		// debug plot
		if(kDebugPlot){
			hSmooth->Draw("lego2 fb"); gPad->SetLogz();
			gPad->Modified(); gPad->Update();
		}
				
		// doing invertion
		ax=h2dEdx[ispec]->GetXaxis();
		ay=h2dEdx[ispec]->GetYaxis();
		for(int ibin=1; ibin<=ax->GetNbins(); ibin++){
			xy[0] = ax->GetBinCenter(ibin); 
			for(int jbin=1; jbin<=ay->GetNbins(); jbin++){
				xy[1] = ay->GetBinCenter(jbin);
				estimate = (ispec ==0 && (ibin==1 || jbin==1)) ? 1.E-10 : hSmooth->GetBinContent(hSmooth->FindBin(log(xy[0]), log(xy[1])))/xy[0]/xy[1];
				h2dEdx[ispec]->SetBinContent(ibin, jbin, (estimate>0.) ? estimate : 1.E-10);
			}
		}
		h2dEdx[ispec]->Scale(1./h2dEdx[ispec]->Integral());
	}
	printf("Done \n");
}

//__________________________________________________________________
void	AliTRDCalPIDLQRef::Reset()
{
  //
  // Resets the reference histograms
  //

	for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
		h2dEdx[ispec]->Reset();
		fPrinc[ispec]->Clear();
	}	
}

//__________________________________________________________________
void  AliTRDCalPIDLQRef::SaveReferences(const Int_t mom, const char *fn)
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
	if(!kFOUND) fSave = new TFile(fn, "RECREATE");
	fSave->cd();

	TH2 *h;
	for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
		h = (TH2D*)h2dEdx[ispec]->Clone(Form("h2dEdx%s%d", AliTRDCalPIDLQ::fpartSymb[ispec], mom));
		h->SetTitle(Form("2D dEdx for particle %s @ %d", AliTRDCalPIDLQ::fpartName[ispec], mom));
		h->GetXaxis()->SetTitle("dE/dx_{TRD}^{amplif} [au]");
		h->GetYaxis()->SetTitle("dE/dx_{TRD}^{drift} [au]");
		h->GetZaxis()->SetTitle("Entries");
		h->Write();
	}
		
	pwd->cd();
}
