void interpolCompare(const Int_t method = 0, const int nstat = 40000)
{
// Compare interpolation methods.
// method = 0 : Select COG method
// method = 1 : Select INT method

	const Float_t alpha = 20.;
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
	Double_t x, chi2, chi2pdf, chi2int, result, err;
	Int_t npPDF, npINT;
	TKDPDF pdf(nstat, 1, bucket, data);
	pdf.SetInterpolationMethod(method);
	pdf.SetAlpha(alpha);
	
	TKDInterpolator in(1, pdf.GetNTNodes());
	in.SetInterpolationMethod(method);
	in.SetAlpha(alpha);
	
	// change to the reprezentation which is closer to the local
	// interpolating polynom.
	Double_t fx;
	Int_t inode = 0;
	TKDNodeInfo *node = 0x0;
	while(node = pdf.GetNodeInfo(inode)){
		if(node->fVal[0] > 0.){
			in.SetNode(inode, *node);
			node  = in.GetNodeInfo(inode);
			fx = node->fVal[0];
			node->fVal[0] = TMath::Log(fx);
			node->fVal[1] = node->fVal[1]/fx;
		}
		inode++;
	}

	// Do preprocessing for INT
	Float_t *pcog;
	if(method) for(int icog=0; icog<pdf.GetNTNodes(); icog++){
		pdf.GetCOGPoint(icog, pcog, val, val_err);
		x = pcog[0];
		pdf.Eval(&x, result, err);
		in.Eval(&x, result, err);
	}


	
	// prepare drawing containers
	TGraphErrors *gPDF = new TGraphErrors(nevals);
	gPDF->SetMarkerStyle(7);
	TGraph *gREP = new TGraph(nevals);
	gREP->SetMarkerStyle(7);
	gREP->SetMarkerColor(4);
	gREP->SetLineColor(4);
	
	TGraphErrors *gINT = new TGraphErrors(nevals);
	gINT->SetMarkerStyle(24);
	gINT->SetMarkerSize(.6);
	gINT->SetMarkerColor(3);
	gINT->SetLineColor(3);
	TGraph *gREI = new TGraph(nevals);
	gREI->SetMarkerStyle(24);
	gREI->SetMarkerColor(4);
	gREI->SetLineColor(4);
	
	// Do the interpolation with PDF and Interpolator
	Double_t dx = 10./nevals;
	x = 5.-dx/2.;
 	chi2pdf = 0.; npPDF = 0;
 	chi2int = 0.; npINT = 0;
	for(int ip=0; ip<nevals; ip++){
		chi2 = pdf.Eval(&x, result, err);
		gPDF->SetPoint(ip, x, result);
		gPDF->SetPointError(ip, 0., err);
		if(TMath::Abs(x)<3.){
			chi2pdf += (result - f->Eval(x))*(result - f->Eval(x))/f->Eval(x)/f->Eval(x);
			npPDF++;
		}
		gREP->SetPoint(ip, x, .1*err/TMath::Abs(result));
		//printf("%2d x[%f] r[%f +- %f] chi2[%e]\n", ip, x, result, err, chi2);

		chi2 = in.Eval(&x, result, err);
		result = TMath::Exp(result);
		err = err*result;
		gINT->SetPoint(ip, x, result);
		gINT->SetPointError(ip, 0., err);
		if(TMath::Abs(x)<3.){
			chi2int += (result - f->Eval(x))*(result - f->Eval(x))/f->Eval(x)/f->Eval(x);
			npINT++;
		}
		gREI->SetPoint(ip, x, .1*err/TMath::Abs(result));
		//printf("%2d x[%f] r[%f +- %f] chi2[%e]\n\n", ip, x, result, err, chi2);

		x -= dx;
	}
	gINT->Draw("p"); leg->AddEntry(gINT, Form("Interpolator [%s]", method ? "INT" : "COG"), "pl");
	gPDF->Draw("p"); leg->AddEntry(gPDF, Form("PDF [%s]", method ? "INT" : "COG"), "pl");	
 	gREI->Draw("pl"); leg->AddEntry(gREI, "1% #epsilon Interpolator", "pl");
 	gREP->Draw("pl"); leg->AddEntry(gREP, "1% #epsilon PDF", "pl");

	chi2pdf /= npPDF;
	printf("PDF : chi2(%d) = %f\n", npPDF, chi2pdf);
	chi2int /= npINT;
	printf("INT : chi2(%d) = %f\n", npINT, chi2int);

	leg->Draw();
}



