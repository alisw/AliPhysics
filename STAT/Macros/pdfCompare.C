void pdfCompare(const Int_t method = 0, const int nstat = 40000)
{
// Compare interpolation methods.
// method = 0 : Select COG method
// method = 1 : Select INT method

	const int bucket = nstat/100; //400;
	const int nevals = 100;
	Float_t *data[1];
	Float_t data1[nstat];
	for(int istat=0; istat<nstat; istat++) data1[istat] = gRandom->Gaus();
	data[0] = &data1[0];

	TLegend *leg = new TLegend(0.7157258,0.8280255,0.9637097,0.9596603);
	leg->SetBorderSize(1);
	TF1 *f=new TF1("f1", "gaus(0);x;f(x)", -5., 5.);
	f->SetLineWidth(1);
	f->SetLineColor(2);
	f->SetParameter(0, 1.);
	f->SetParameter(1, 0.);
	f->SetParameter(2, 1.);
	f->SetParameter(0, 1./f->Integral(-10., 10.));
	f->Draw(); 	leg->AddEntry(f, "model", "l");
	
	
	Float_t cog, val, val_err;
	Double_t x, chi2, CHI2, result, err;
	Int_t npoints;
	TKDPDF pdf(nstat, 1, bucket, data);
	pdf.SetInterpolationMethod(method);
	pdf.SetStore();
	pdf.SetAlpha(1.);
	Float_t *bounds;// = in.GetBoundary(0);
	Double_t dx = 10./nevals;

	TGraph *gRE = new TGraph(npoints);
	gRE->SetMarkerStyle(4);
	gRE->SetMarkerColor(4);
	gRE->SetLineColor(4);
	
	// do a preevaluation for INT method on COG points
	Float_t *pcog;
	for(int icog=0; icog<pdf.GetNTNodes(); icog++){
		pdf.GetCOGPoint(icog, pcog, val, val_err);
		x = pcog[0];
		chi2 = pdf.Eval(&x, result, err);
	}
	//pdf.GetStatus();

	TGraphErrors *gINT = new TGraphErrors(npoints);
	gINT->SetMarkerStyle(7);
	x = 5.-dx/2.;
 	CHI2 = 0.; npoints = 0;
	for(int ip=0; ip<nevals; ip++){
		chi2 = pdf.Eval(&x, result, err);
		gINT->SetPoint(ip, x, result);
		gINT->SetPointError(ip, 0., err);
		gRE->SetPoint(ip, x, .01*err/TMath::Abs(result));
		if(TMath::Abs(x)<2.){
			CHI2 += (result - f->Eval(x))*(result - f->Eval(x))/err/err;
			npoints++;
		}
		x -= dx;
	}
	gINT->Draw("p"); leg->AddEntry(gINT, Form("interpolation [%s]", method ? "INT" : "COG"), "pl");
	
	gRE->Draw("pl"); leg->AddEntry(gRE, "1% relative errors", "pl");
	CHI2 /= npoints;
	printf("chi2(%d) = %f\n", npoints, CHI2);

	leg->Draw();
}



