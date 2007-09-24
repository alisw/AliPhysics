Float_t p[]={1.4, -.6}; //}{1.7, -.4};
//______________________________________________________________
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

void performanceKNN(const int np = 10000)
{
	const int ntimes = 1000; // times to repeate kNN search
	const int ndim = 5;      // no of dimensions for data

	if(!gSystem->FindFile(".", Form("%dD_Gauss.root", ndim))) build(ndim);
	TFile::Open(Form("%dD_Gauss.root", ndim));
	TString var="x0";
	for(int i=1; i<ndim; i++) var+=Form(":x%d", i);
	TKDInterpolator in(db, var.Data(), "", 30, np);
	
	Int_t *index, fNN = 20;
	Float_t *d, p[ndim];
	for(int idim=0; idim<ndim; idim++) p[idim] = gRandom->Gaus();
	TStopwatch time;
	time.Start();
	for(int itime=0; itime<ntimes; itime++) in.FindNearestNeighbors(p, fNN, index, d);
	time.Stop();
	printf("cpu = %f [mus] real = %f [mus]\n", 1.E6*time.CpuTime()/ntimes, 1.E6*time.RealTime()/ntimes);
	for(int i=0; i<fNN; i++){
		printf("%5d ", index[i]);
		if(i%5==4) printf("\n");
	}
	
	// draw kNN
	TH2 *h2 = new TH2I("h2", "", 100, -2., 2., 100, -2., 2.);
	h2->Draw();
	TGraph *gKD = new TGraph(fNN);
	gKD->SetMarkerStyle(24);
	gKD->SetMarkerColor(2);
	gKD->SetMarkerSize(.5);
	Float_t pknn[ndim];
	for(int ip=0; ip<fNN; ip++){
		in.GetDataPoint(index[ip], pknn);
		gKD->SetPoint(ip, pknn[0], pknn[1]);
	}
	gKD->Draw("p");
	TMarker *mp = new TMarker(p[0], p[1], 3);
	mp->SetMarkerColor(4);
	mp->Draw();

	
	// STAND ALONE
	const Int_t kNN = fNN;
	Float_t ftmp, gtmp, dist;
	Int_t itmp, jtmp;
	Int_t   fkNN[kNN];
	Float_t fkNNdist[kNN];
	for(int i=0; i<kNN; i++) fkNNdist[i] = 9999.;

	// read data in memory
	Int_t npoints = db->Draw("x0", "", "goff", np);
	Float_t **x = new Float_t*[ndim];
	Double_t *v;
	for(int idim=0; idim<ndim; idim++){
		x[idim] = new Float_t[npoints];
		db->Draw(Form("x%d", idim), "", "goff", np); v = db->GetV1();
		for(int i=0; i<npoints; i++) x[idim][i] = v[i];
	}
	// calculate
	printf("stand alone calculation for %d neighbors in %d points\n", kNN, npoints);
	time.Start(kTRUE);
	for(int idx=0; idx<npoints; idx++){
		// calculate distance in the L1 metric
		dist = 0.;
		for(int idim=0; idim<ndim; idim++) dist += TMath::Abs(p[idim] - x[idim][idx]);
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
	time.Stop();
	printf("cpu = %f [s] real = %f [s]\n", time.CpuTime(), time.RealTime());
	
	for(int i=0; i<kNN; i++){
		printf("%5d ", fkNN[i]);
		if(i%5==4) printf("\n");
	}
	
	TGraph *gSA = new TGraph(kNN);
	gSA->SetMarkerStyle(24);
	for(int ip=0; ip<kNN; ip++) gSA->SetPoint(ip, x[0][fkNN[ip]], x[1][fkNN[ip]]);

	gSA->Draw("p");
}

//______________________________________________________________
void build(const Int_t ndim = 2, const Int_t nstat = 1000000)
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