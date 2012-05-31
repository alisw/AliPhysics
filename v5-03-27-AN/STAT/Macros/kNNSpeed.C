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
#include "../src/TKDInterpolator.h"
#include "../src/TKDTree.h"

const int ntimes = 500000; // times to repeate kNN search
TStopwatch timer;

Float_t kdTreeNN(Float_t **p, const int ndim = 5, const int np = 10000);
Float_t standAloneNN(Float_t **p, const int ndim = 5, const int np = 10000);

//______________________________________________________________
TTree *db = 0x0;
void kNNSpeed(const Int_t method = 0)
{
// Estimate kNN algorithm speed:
// method = 0 : run kDTree kNN search
// method = 1 : run kNN loop search

	const Int_t ndim  = 5;
	const Int_t nstat = 9;


	// generate test sample
	Float_t **p = new Float_t*[ntimes];
	for(int ip =0 ; ip<ntimes; ip++){
		p[ip] = new Float_t[ndim];
		for(int idim=0; idim<ndim; idim++) p[ip][idim] = gRandom->Uniform(-.5, .5);
	}

	Float_t time;
	Int_t stat;
	TLegend *leg = new TLegend(.7, .7, .9, .9);
	TH1 *h1 = new TH1I("h1", "", 100, 9.5, 20.5);
	h1->SetMaximum(1.E2);
	h1->SetMinimum(1.E-1);
	h1->GetXaxis()->SetTitle("Statistics");
	h1->GetYaxis()->SetTitle("CPU [#mus]");
	h1->SetStats(0);
	h1->Draw();
	TGraph **g=new TGraph*[ndim];
	for(int idim=2; idim<3/*ndim*/; idim++){
		g[idim] = new TGraph(nstat);
		g[idim]->SetMarkerStyle(20+idim);
		g[idim]->SetMarkerColor(idim+1);
		leg->AddEntry(g[idim], Form("%dD", idim+1), "pl");
		stat = 16384;//1024;
		for(int istat = 0; istat<nstat; istat++){
			time = kdTreeNN(p, idim+1, stat);
			//time = standAloneNN(p, idim+1, stat);
			g[idim]->SetPoint(istat, TMath::Log2(float(stat)), time);
			stat *= 2;
		}
		g[idim]->Draw("pl");
	}
	leg->Draw();
	
	for(int ip=0; ip<ntimes; ip++) delete [] p[ip];
	delete [] p;
}

//______________________________________________________________
Float_t kdTreeNN(Float_t **p, const int ndim, const int np)
{
	Float_t **x = new Float_t*[ndim];
	for(int idim =0 ; idim<ndim; idim++){
		x[idim] = new Float_t[np];
		for(int ip=0; ip<np; ip++) x[idim][ip] = gRandom->Gaus();
	}
	TKDTreeIF nnFinder(np, ndim, 1, x);
	nnFinder.MakeBoundaries();
	
	Int_t *index, fNN = 1;
	Float_t *d;
	timer.Start(kTRUE);
	for(int itime=0; itime<ntimes; itime++){
		nnFinder.FindNearestNeighbors(p[itime], fNN, index, d);
	}
	timer.Stop();
	Float_t time = 1.E6*timer.CpuTime()/float(ntimes);

	printf("np[%7d] nd[%d] time = %5.2f[mus]\n", np, ndim, time);
	
	// clean
	for(int idim=0; idim<ndim; idim++) delete [] x[idim];
	delete [] x;
	

	return time;
}

//______________________________________________________________
Float_t standAloneNN(Float_t **p, const int ndim, const int np)
{	
	// STAND ALONE
	const Int_t kNN = 1;
	Float_t ftmp, gtmp, dist;
	Int_t itmp, jtmp;
	Int_t   fkNN[kNN];
	Float_t fkNNdist[kNN];
	for(int i=0; i<kNN; i++) fkNNdist[i] = 9999.;

	Float_t **x = new Float_t*[ndim];
	for(int idim =0 ; idim<ndim; idim++){
		x[idim] = new Float_t[np];
		for(int ip=0; ip<np; ip++) x[idim][ip] = gRandom->Gaus();
	}
	Int_t npoints = np;
	
	// calculate
	Int_t ntimes = 100;
	timer.Start(kTRUE);
	for(int it=0; it<ntimes; it++){
		for(int idx=0; idx<npoints; idx++){
			// calculate distance in the L1 metric
			dist = 0.;
			for(int idim=0; idim<ndim; idim++) dist += TMath::Abs(p[it*10][idim] - x[idim][idx]);
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
	}
	timer.Stop();
	
	Float_t time = timer.CpuTime()/float(ntimes);
	printf("np[%5d] nd[%d] time[%f] %6.2f\n", np, ndim, time, timer.CpuTime());

	return time;
}

