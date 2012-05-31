// #include "TSystem.h"
// #include "TFile.h"
// #include "TTree.h"
// #include "TProfile.h"
// #include "TF1.h"
// #include "TF2.h"
// #include "TF3.h"
// #include "TLegend.h"
// #include "TRandom.h"
// #include "TString.h"
// #include "TGraph.h"
// #include "TMarker.h"
// #include "TStopwatch.h"
// #include "../src/TKDPDF.h"
// #include "../src/TKDInterpolator.h"

void interpolTest(const Int_t ndim = 1);
Double_t interpolation(const int nstat, const int ndim, const int first);
void build(const int ndim, const int npoints);

TF1 *f = 0x0;

//__________________________________________________________
void interpolTest(const Int_t ndim)
{
// Test interpolator on a variable dimension (ndim<=5) uncorrelated
// Gauss model. The macro displays the chi2 between interpolation and
// model for various statistics of experimental data.
//
// Check the worker function interpolation() to fine tune the
// algorithm.
// 
// Attention:
// The macro works in compiled mode !

	// define model
	switch(ndim){
	case 1:
		f=new TF1("f1", "gaus(0);x;f(x)", -5., 5.);
		f->SetParameter(0, 1.);
		f->SetParameter(1, 0.);
		f->SetParameter(2, 1.);
		f->SetParameter(0, 1./f->Integral(-5., 5.));
		break;
	case 2:
		f=new TF2("f2", "[0]*exp(-0.5*((x-[1])/[2])**2)*exp(-0.5*((y-[3])/[4])**2)", -5., 5., -5., 5.);
		f->SetParameter(0, 1.);
		f->SetParameter(1, 0.);
		f->SetParameter(2, 1.);
		f->SetParameter(3, 0.);
		f->SetParameter(4, 1.);
		f->SetParameter(0, 1./f->Integral(-5., 5., -5., 5.));
		break;
	case 3:
	case 4:
	case 5:
		f=new TF3("f3", "[0]*exp(-0.5*((x-[1])/[2])**2)*exp(-0.5*((y-[3])/[4])**2)*exp(-0.5*((z-[5])/[6])**2)", -5., 5., -5., 5., -5., 5.);
		f->SetParameter(0, 1.);
		f->SetParameter(1, 0.);
		f->SetParameter(2, 1.);
		f->SetParameter(3, 0.);
		f->SetParameter(4, 1.);
		f->SetParameter(5, 0.);
		f->SetParameter(6, 1.);
		f->SetParameter(0, 1./f->Integral(-5., 5., -5., 5., -5., 5.));
		if(ndim>3) printf("Warning : chi2 results are unreliable !\n");
		break;
	default:
		printf("%dD not supported in this macro.\n", ndim);
		return;
	}

	Double_t chi2;
	const Int_t nchi2 = 18;
	const Int_t nfirst = 1;
 	//TProfile *hp = new TProfile("hp", "", 100, 1.E4, 1.E6);
 	//hp->SetMarkerStyle(24);
	TGraph *gChi2 = new TGraph(nchi2);
 	gChi2->SetMarkerStyle(24);
	Int_t ip = 0, nstat, first;
	for(int ifirst=0; ifirst<nfirst; ifirst++){
		first = 0;//Int_t(gRandom->Uniform(1.E4));
		for(int i=2; i<10; i++){
			nstat = 10000*i;
			chi2 = interpolation(nstat, ndim, first);
			gChi2->SetPoint(ip++, TMath::Log10(float(nstat)), chi2);
			//hp->Fill(float(nstat), chi2);
		}
		for(int i=1; i<10; i++){
			nstat = 100000*i;
			chi2 = interpolation(nstat, ndim, first);
			gChi2->SetPoint(ip++, TMath::Log10(float(nstat)), chi2);
			//hp->Fill(float(nstat), chi2);
		}
		gChi2->Draw("apl");
	}
	//hp->Draw("pl");
}

void test()
{
	TClonesArray ar("TKDNodeInfo", 10);
}

//__________________________________________________________
Double_t interpolation(const int nstat, const int ndim, const int first)
{
// Steer interpolation of a nD (n<=5) uncorrelated Gauss distribution.
// The user is suppose to give the experimental statistics (nstat) the
// dimension of the experimental point (ndim). The entry from which the
// data base is being read is steered by "first".

 	const int bs = 400;
 	const int npoints = 100;
	// switch on chi2 calculation between interpolation and model.
	// Default calculates distance.
	const Bool_t kCHI2 = kFALSE; 
	const Bool_t kINT  = kFALSE; // switch on integral interpolator

	// open data base
	TString fn("5D_Gauss.root");
	if(!gSystem->FindFile(".", fn)) build(5, 1000000);
	TFile::Open("5D_Gauss.root");
	TTree *t = (TTree*)gFile->Get("db");

	// build interpolator
	TString var = "x0"; for(int idim=1; idim<ndim; idim++) var += Form(":x%d", idim);
	TKDPDF pdf(t, var.Data(), "", bs, nstat, first);
	printf("\n\nFINISH BUILDING PDF\n\n");
	
	Double_t fx;
	Int_t inode = 0;
	TKDNodeInfo *node = 0x0;
	TKDInterpolator in(ndim, pdf.GetNTNodes());
/*	while(node = pdf.GetNodeInfo(inode)){
		if(node->fVal[0] > 0.){
			fx = node->fVal[0];
			node->fVal[0] = TMath::Log(fx);
			node->fVal[1] = node->fVal[1]/fx;
		}
		in.SetNode(inode, *node);
		inode++;
	}*/
	//pdf.SetInterpolationMethod(kINT);
	pdf.SetAlpha(2.);
	//in.SetStore();


	Double_t x[2], r, e, chi2;
	TH2 *h2 = new TH2F("h2", "", 50, -4., 4., 50, -4., 4.);
	TAxis *ax = h2->GetXaxis(), *ay = h2->GetYaxis();
	for(int ix=1; ix<=ax->GetNbins(); ix++){
		x[0] = ax->GetBinCenter(ix);
		for(int iy=1; iy<=ay->GetNbins(); iy++){
			x[1] = ay->GetBinCenter(iy);

			chi2 = pdf.Eval(x, r, e);
			printf("x[%2d] x[%2d] r[%f] e[%f] chi2[%f]\n", ix, iy, r, e, chi2);
			h2->SetBinContent(ix, iy, r/*TMath::Exp(r)*/);
		}
	}
	h2->Draw("lego2");
	return;
	
	// define testing grid
	Double_t x[5], r, e;
	Double_t dx[5] = {4., 4., 4., 4., 4.};
	Double_t chi2, theor, result, err;
	for(int idim=0; idim<ndim; idim++){
		dx[idim] = (ndim<<2)/float(npoints);
		x[idim]  = -2.+dx[idim]/2.;
	}

	
	// setup the integral interpolator and do memory test
	Float_t c[5], v, ve;
	for(int ip=0; ip<in.GetNTNodes(); ip++){
		in.GetCOGPoint(ip, c, v, ve);
		for(int idim=0; idim<ndim; idim++) x[idim] = (Double_t)c[idim];
		in.Eval(x, result, err);
	}
	printf("\nInterpolating %d data points in %dD ...\n", nstat, ndim);
	printf("\tbucket size %d\n", bs);
	printf("\tstarting interpolation ...\n");
	return;

	
	// evaluate relative error
	for(int idim=0; idim<ndim; idim++) x[idim] = 0.;
	in.Eval(x, result, err);
	Double_t threshold = err/result;
	
	// starting interpolation over the testing grid
	chi2 = 0.; Int_t nsample = 0;
	x[4]  = -2.+dx[4]/2.;
	do{
		x[3]  = -2.+dx[3]/2.;
		do{
			x[2]  = -2.+dx[2]/2.;
			do{
				x[1]  = -2.+dx[1]/2.;
				do{
					x[0]  = -2.+dx[0]/2.;
					do{
						in.Eval(x, result, err);
						if(err/TMath::Abs(result) < 1.1*threshold){
							theor = f->Eval(x[0], x[1], x[2]);
							Double_t tmp = result - theor; tmp *= tmp;
							chi2 += tmp/(kCHI2 ? err*err : theor);
							nsample++;
						}
					} while((x[0] += dx[0]) < 2.);
				} while((x[1] += dx[1]) < 2.);
			} while((x[2] += dx[2]) < 2.);
		} while((x[3] += dx[3]) < 2.);
	} while((x[4] += dx[4]) < 2.);
	chi2 /= float(nsample);
	printf("\tinterpolation quality (chi2) %f\n", chi2);
	
	return kCHI2 ? chi2 : TMath::Sqrt(chi2);
}


//___________________________________________________________
void build(const int ndim, const int npoints)
{
	printf("Building DB ...\n");
	Float_t data[ndim];
	TFile::Open(Form("%dD_Gauss.root", ndim), "RECREATE");
	TTree *t = new TTree("db", Form("%dD data base for kD statistics", ndim));
	for(int idim=0; idim<ndim; idim++) t->Branch(Form("x%d", idim), &data[idim], Form("x%d/F", idim));

	for (Int_t ip=0; ip<npoints; ip++){
		for (Int_t id=0; id<ndim; id++) data[id]= gRandom->Gaus();
		t->Fill();
	}

	t->Write();
	gFile->Close();
	delete gFile;
}
