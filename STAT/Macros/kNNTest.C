#include <malloc.h>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2I.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TString.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TStopwatch.h"
#include "../src/TKDPDF.h"
#include "../src/TKDTree.h"

void kNNTest(const int np = 10000, const int ndim = 2);
void kNNDraw(const Float_t *p, const int kNN=20);
void build(const Int_t ndim = 2, const Int_t nstat = 1000000);
Float_t Mem();


//______________________________________________________________
void kNNTest(const int np, const int ndim)
{
// Test speed and quality of nearest neighbors search.
// The results are compared with a simple loop search through the data.
// The macro should run in compiled mode.

	const int ntimes = 1000; // times to repeate kNN search
	const Int_t kNN = 1; // number of neighbors

	Float_t **x = new Float_t*[ndim];
	for(int idim =0 ; idim<ndim; idim++){
		x[idim] = new Float_t[np];
		for(int ip=0; ip<np; ip++) x[idim][ip] = gRandom->Gaus();
	}
	TKDTreeIF nnFinder(np, ndim, 1, x);
	
	// generate test sample
	Float_t **p = new Float_t*[ntimes];
	for(int ip =0 ; ip<ntimes; ip++){
		p[ip] = new Float_t[ndim];
		for(int idim=0; idim<ndim; idim++) p[ip][idim] = gRandom->Gaus();
	}


	Int_t *index;
	Float_t *d;
	TStopwatch timer;
	timer.Start(kTRUE);
	for(int itime=0; itime<ntimes; itime++) nnFinder.FindNearestNeighbors(p[itime], kNN, index, d);
	timer.Stop();
	printf("kDTree NN calculation.\n");
	printf("cpu = %f [ms] real = %f [ms]\n", 1.E3*timer.CpuTime()/ntimes, 1.E3*timer.RealTime()/ntimes);
	printf("Points indexes in increasing order of their distance to the target point.\n");
	for(int i=0; i<kNN; i++){
		printf("%5d [%7.4f] ", index[i], d[i]);
		if(i%5==4) printf("\n");
	}
	printf("\n");
	
	// draw kNN
	TLegend *leg = new TLegend(.7, .7, .9, .9);
	leg->SetBorderSize(1);
	leg->SetHeader("NN finders");
	TH2 *h2 = new TH2I("h2", "", 100, -2., 2., 100, -2., 2.);
	h2->Draw(); h2->SetStats(kFALSE);
	TMarker *mp = new TMarker(p[ntimes-1][0], p[ntimes-1][1], 3);
	mp->SetMarkerColor(4);
	mp->Draw(); leg->AddEntry(mp, "Target", "p");
	TGraph *gKD = new TGraph(kNN);
	gKD->SetMarkerStyle(24);
	gKD->SetMarkerColor(2);
	gKD->SetMarkerSize(.5);
	for(int ip=0; ip<kNN; ip++) gKD->SetPoint(ip, x[0][index[ip]], x[1][index[ip]]);
	gKD->Draw("p"); leg->AddEntry(gKD, "kDTree", "p");
	

		
	// STAND ALONE
	Float_t ftmp, gtmp, dist;
	Int_t itmp, jtmp;
	Int_t   fkNN[kNN];
	Float_t fkNNdist[kNN];
	for(int i=0; i<kNN; i++) fkNNdist[i] = 9999.;

	// calculate
	timer.Start(kTRUE);
	for(int idx=0; idx<np; idx++){
		// calculate distance in the L1 metric
		dist = 0.;
		for(int idim=0; idim<ndim; idim++) dist += TMath::Abs(p[ntimes-1][idim] - x[idim][idx]);
		if(dist >= fkNNdist[kNN-1]) continue;

		//insert neighbor
		int iNN=0;
		while(dist >= fkNNdist[iNN]) iNN++;
		itmp = fkNN[iNN]; ftmp = fkNNdist[iNN];
		fkNN[iNN]     = idx;
		fkNNdist[iNN] = dist;
		for(int ireplace=iNN+1; ireplace<kNN; ireplace++){
			jtmp = fkNN[ireplace]; gtmp = fkNNdist[ireplace];
			fkNN[ireplace] = itmp; fkNNdist[ireplace] = ftmp;
			itmp = jtmp; ftmp = gtmp;
			if(ftmp == 9999.) break;
		}
	}
	timer.Stop();
	printf("\nStand Alone NN calculation.\n");
	printf("cpu = %f [ms] real = %f [ms]\n", 1.E3*timer.CpuTime(), 1.E3*timer.RealTime());
	printf("Points indexes in increasing order of their distance to the target point.\n");
	for(int i=0; i<kNN; i++){
		printf("%5d [%7.4f] ", fkNN[i], fkNNdist[i]);
		if(i%5==4) printf("\n");
	}
	printf("\n");	

	TGraph *gSA = new TGraph(kNN);
	gSA->SetMarkerStyle(24);
	for(int ip=0; ip<kNN; ip++) gSA->SetPoint(ip, x[0][fkNN[ip]], x[1][fkNN[ip]]);

	gSA->Draw("p"); leg->AddEntry(gSA, "Stand Alone", "p");
	leg->Draw();
	
	for(int ip=0; ip<ntimes; ip++) delete [] p[ip];
	delete [] p;
	for(int idim=0; idim<ndim; idim++) delete [] x[idim];
	delete [] x;
}


//______________________________________________________________
Float_t p[]={1.4, -.6}; //}{1.7, -.4};
void kNNDraw(const Float_t *p, const int kNN)
{
// Test memory consumption and draw "kNN" nearest neighbours of point "p".
// The distance (in the L1 metric) is encoded in the color code.

	const Int_t npoints = 10000;
	Float_t *data0 =  new Float_t[npoints*2];
  Float_t *data[2];
  data[0] = &data0[0];
  data[1] = &data0[npoints];
  for (Int_t i=0;i<npoints;i++) {
    data[1][i]= gRandom->Gaus();
    data[0][i]= gRandom->Gaus();
  }

	TKDPDF pdf(npoints, 2, 100, data);
	pdf.DrawNode(pdf.GetNodeIndex(p));
	
	TMarker *mp = new TMarker(p[0], p[1], 3);
	mp->SetMarkerColor(4);
	mp->Draw();
		
	Int_t *index, color;
	Float_t *d, d0, pknn[2];
	Float_t start = Mem();
	pdf.FindNearestNeighbors(p, kNN, index, d);
	Float_t end = Mem();
	printf("FindNearestNeighbors memory usage %fKB\n", end-start);
	d0 = d[kNN-1];
	TMarker *ma = new TMarker[kNN];
	for(int ip=0; ip<kNN; ip++){
		pdf.GetDataPoint(index[ip], pknn);
		color = 101 - Int_t(50. * d[ip]/d0);
		ma[ip].SetMarkerStyle(4);
		ma[ip].SetMarkerColor(color);
		ma[ip].DrawMarker(pknn[0], pknn[1]);
	}
}


//______________________________________________________________
Float_t Mem()
{
  // size in KB
  static struct mallinfo memdebug; 
  memdebug = mallinfo();
  return memdebug.uordblks/1024.;
}

//______________________________________________________________
void build(const Int_t ndim, const Int_t nstat)
{
// Build "nstat" data points in "ndim" dimensions taken from an
// uncorrelated 2D Gauss distribution.
	

	printf("build data ... \n");
	Double_t pntTrain[ndim];
	TFile *f = TFile::Open(Form("%dD_Gauss.root", ndim), "RECREATE");
	TTree *t = new TTree("db", "gauss database");
	for(int idim=0; idim<ndim; idim++) t->Branch(Form("x%d", idim), &pntTrain[idim], Form("x%d/D", idim));
	
	for(int istat=0; istat<nstat; istat++){
		for(int idim=0; idim<ndim; idim++) pntTrain[idim] = gRandom->Gaus();
		t->Fill();
	}
	f->cd();
	t->Write();
	f->Close();
	delete f;
}

