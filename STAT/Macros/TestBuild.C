Int_t npoints = 301;
Int_t bsize = 10;

Float_t *data0 =  new Float_t[npoints*2];
Float_t *data[2];
data[0] = &data0[0];
data[1] = &data0[npoints];
Float_t dataf[] = {.54, .54};
TKDTreeIF* TestBuild(Int_t k = 5)
{
	gStyle->SetOptStat(0);
	
	for (Int_t i=0;i<npoints;i++) {
		data[1][i]= gRandom->Gaus(.5, .1);
		data[0][i]= gRandom->Gaus(.5, .1);
		//data[1][i]= gRandom->Rndm();
		//data[0][i]= gRandom->Rndm();
	}
	//TKDTreeIF *tree = new TKDTreeIF(npoints, 2, bsize, data);
	//TKDSpline *spline = new TKDSpline(npoints, 2, bsize, data);
	//spline->DrawNodes(0, 1, depth);
	TKDInterpolator *s = new TKDInterpolator(npoints, 2, bsize, data);
	s->DrawNodes(-1);
	TMarker *m = new TMarker(dataf[0], dataf[1], 20);
	m->Draw();
	TGraph *g=new TGraph(npoints);
	g->SetMarkerStyle(7);
	for(int ip=0; ip<npoints; ip++) g->SetPoint(ip, data[0][ip], data[1][ip]);
	g->Draw("p");
	
	Int_t *index = 0x0;
	Float_t dist;
	s->FindNearestNeighbors(dataf, k, index, dist);
	
	TGraph *gNN=new TGraph(k);
	gNN->SetMarkerStyle(24);
	gNN->SetMarkerColor(2);
	for(int i=0; i<k; i++){
		//printf("%d x %d y %d\n", i, index[i], index[i]);
		gNN->SetPoint(i, data[0][index[i]], data[1][index[i]]);
	}
	gNN->Draw("p");

	return s;	
}
