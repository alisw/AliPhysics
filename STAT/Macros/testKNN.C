const Float_t p[]={1.4, -.6}; //}{1.7, -.4};
void testKNN(const Float_t *p, const int kNN=20)
{
// Draw "kNN" nearest neighbours of point "p". The distance (in the L1
// metric) is encoded in the color code.
// To build the data refere to function build().

	TFile::Open("2D_Gauss.root");
	TKDInterpolator in(db, "x0:x1", "x0>-1.5&&x0<2.&&x1>-2.&&x1<2.", 300);
	in.DrawNode(in.FindNode(p)-in.GetNNodes());
	
	TMarker *mp = new TMarker(p[0], p[1], 3);
	mp->SetMarkerColor(4);
	mp->Draw();
		
	Int_t *index, color;
	Float_t d, d0, pknn[2];
	in.FindNearestNeighbors(p, kNN, index, d0);
	TMarker *ma = new TMarker[kNN];
	for(int ip=0; ip<kNN; ip++){
		in.GetDataPoint(index[ip], pknn);
		d = TMath::Abs(p[0]-pknn[0]) + TMath::Abs(p[1]-pknn[1]);
		ma[ip].SetMarkerStyle(4);
		color = 101 - Int_t(50. * d/d0);
		ma[ip].SetMarkerColor(color);
		ma[ip].DrawMarker(pknn[0], pknn[1]);
	}
}

void build(const Int_t ndim = 2, const Int_t nstat = 100000)
{
// Build "nstat" data points in "ndim" dimensions taken from an
// uncorrelated 2D Gauss distribution.
	

	printf("build data ... \n");
	Double_t pntTrain[ndim];
	f = TFile::Open(Form("%dD_Gauss.root", ndim), "RECREATE");
	t = new TTree("db", "gauss database");
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