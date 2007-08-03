void testInterpolator()
{
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	
	Int_t npoints = 301;
	Int_t bsize = 10;
	
	Float_t *data0 =  new Float_t[npoints*2];
	Float_t *data[2];
	data[0] = &data0[0];
	data[1] = &data0[npoints];
	for (Int_t i=0;i<npoints;i++) {
		data[1][i]= gRandom->Gaus(.5, .1);
		data[0][i]= gRandom->Gaus(.5, .1);
	}
	TKDInterpolator *in = new TKDInterpolator(npoints, 2, bsize, data);

	TH2 *hS = new TH2F("hS", "", 10, .3, .7, 10, .3, .7);
	TAxis *ax = hS->GetXaxis(), *ay = hS->GetYaxis();
	Float_t p[2], eval;
	for(int ix=1; ix<=ax->GetNbins(); ix++){
		p[0] = ax->GetBinCenter(ix);
		for(int iy=1; iy<=ay->GetNbins(); iy++){
			p[1] = ay->GetBinCenter(iy);
			eval = in->Eval(p);
			//printf("x %f y %f eval %f [%d]\n", p[0], p[1], eval, TMath::IsNaN(eval));
			if(!TMath::IsNaN(eval)) hS->SetBinContent(ix, iy, eval);
		}
	}
	hS->Draw("lego2");
}