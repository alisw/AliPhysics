const Int_t ndim = 2;
Double_t testInterpolator(const Int_t nstat = 100000)
{
// Macro for testing the TKDInterpolator.
// 
// The function which it is interpolated is an uncorrelated landau
// distribution in "ndim" dimensions. The function returns the chi2 of
// the interpolation.
// 
// Parameters
// 	nstat - number of points to be used for training
// 	kBuild - on/off generate data
// 	kTransform - on/off outliers compresion
// 	kPCA - on/off pricipal component analysis 

	const Bool_t kBuild = 1;
	const Bool_t kTransform = 1;
	const Bool_t kPCA = 1;
	gStyle->SetPalette(1);
	
	Double_t pntTrain[ndim], pntTest[ndim], pntRotate[ndim];
	Double_t pdf;
	TFile *f = 0x0, *fEval = 0x0;
	TTree *t = 0x0, *tEval = 0x0;


	// build data
	if(kBuild){
		printf("build data ... \n");
		f = TFile::Open(Form("%dD_LL.root", ndim), "RECREATE");
		t = new TTree("db", "Log-Log database");
		for(int idim=0; idim<ndim; idim++) t->Branch(Form("x%d", idim), &pntTrain[idim], Form("x%d/D", idim));
		
		fEval = TFile::Open(Form("%dD_Eval.root", ndim), "RECREATE");
		tEval = new TTree("db", "Evaluation database");
		for(int idim=0; idim<ndim; idim++) tEval->Branch(Form("x%d", idim), &pntTest[idim], Form("x%d/D", idim));
		
		for(int istat=0; istat<nstat; istat++){
			for(int idim=0; idim<ndim; idim++) pntTrain[idim] = gRandom->Landau(5.);
			if(!(istat%3)){ // one third of the statistics is for testing
				memcpy(pntTest, pntTrain, ndim*sizeof(Double_t));
				tEval->Fill();
				continue;
			}
			if(kTransform)
				for(int idim=0; idim<ndim; idim++)
					if(pntTrain[idim] > 0.) pntTrain[idim] = TMath::Log(pntTrain[idim]);
					else pntTrain[idim] = 0.;
			t->Fill();
		}
		f->cd();
		t->Write();
		f->Flush();
		
		fEval->cd();
		tEval->Write();
		fEval->Flush();
	} else {// link data
		printf("link data ... \n");
		f = TFile::Open(Form("%dD_LL.root", ndim));
		t = (TTree*)f->Get("db");
		for(int idim=0; idim<ndim; idim++) t->SetBranchAddress(Form("x%d", idim), &pntTrain[idim]);
		
		fEval = TFile::Open(Form("%dD_Eval.root", ndim));
		tEval = (TTree*)fEval->Get("db");
		for(int idim=0; idim<ndim; idim++) tEval->SetBranchAddress(Form("x%d", idim), &pntTest[idim]);
	}

	
	// do principal component analysis (PCA)
	TPrincipal princ(ndim, "N");
	if(kPCA && kBuild){
		printf("do principal component analysis (PCA) ... \n");
		f->cd();
		TTree *tt = new TTree("db1", "PCA database");
		for(int idim=0; idim<ndim; idim++) tt->Branch(Form("x%d", idim), &pntRotate[idim], Form("x%d/D", idim));
		for(int ientry=0; ientry<t->GetEntries(); ientry++){
			t->GetEntry(ientry);
			princ.AddRow(pntTrain);
		}
		princ.MakePrincipals();
		for(int ientry=0; ientry<t->GetEntries(); ientry++){
			t->GetEntry(ientry);
			princ.X2P(pntTrain, pntRotate);
			tt->Fill();
		}
		tt->Write();
		f->Flush();
		for(int idim=0; idim<ndim; idim++) tt->SetBranchAddress(Form("x%d", idim), &pntTrain[idim]);
		t = tt;
	}
	gROOT->cd();
	
	// do interpolation
	printf("do interpolation ... \n");
	Double_t pdf, pdf_estimate, chi2;
	TString vl = "x0";
	for(int idim=1; idim<ndim; idim++) vl+=Form(":x%d", idim);
	TKDInterpolator in(t, vl.Data(), "", 200.);
	chi2 = 0.;
/*	for(int ip=0; ip<tEval->GetEntries(); ip++){
		tEval->GetEntry(ip);
		printf("\nEval %d\n", ip);*/
	TH1 *h1 = new TH2F("h1", "", 50, 0., 100., 50, 0., 100.);
	TH1 *h2 = new TH2F("h2", "", 50, 0., 100., 50, 0., 100.);
	TAxis *ax = h2->GetXaxis(), *ay = h2->GetYaxis();
	for(int ix=2; ix<ax->GetNbins(); ix++){
		pntTest[0] = ax->GetBinCenter(ix);
	for(int iy=2; iy<ay->GetNbins(); iy++){
		pntTest[1] = ay->GetBinCenter(iy);
		memcpy(pntTrain, pntTest, ndim*sizeof(Double_t));

		if(kTransform)
			for(int idim=0; idim<ndim; idim++)
				if(pntTrain[idim] > 0.) pntTrain[idim] = TMath::Log(pntTrain[idim]);
				else pntTrain[idim] = 0.;
		
		if(kPCA){
			princ.X2P(pntTrain, pntRotate);
			memcpy(pntTrain, pntRotate, ndim*sizeof(Double_t));
		}

		pdf_estimate = in.Eval(pntTrain, 30);
		// calculate chi2
		if(kTransform)
			for(int idim=0; idim<ndim; idim++)
				if(pntTest[idim] > 0.) pdf_estimate /= pntTest[idim];
				else continue; 
		
		h1->SetBinContent(ix, iy, pdf_estimate);
		
		pdf = 1.; for(int idim=0; idim<ndim; idim++) pdf *= TMath::Landau(pntTest[idim], 5.);
		h2->SetBinContent(ix, iy, pdf);
		pdf_estimate -= pdf;
		chi2 += pdf_estimate*pdf_estimate/pdf;
	}}
	f->Close(); delete f;
	fEval->Close(); delete fEval;
	
	// results presentation
	printf("chi2 = %f\n", chi2);
	TCanvas *c = 0x0;
	if(!(c = (TCanvas*)gROOT->FindObject("c"))){
		c = new TCanvas("c", "", 10, 10, 900, 500);
		c->Divide(2, 1);
	}
	c->cd(1);
	h1->Draw("lego2"); h1->GetZaxis()->SetRangeUser(1.e-9, 5.e-2); gPad->SetLogz(); gPad->Modified(); gPad->Update();
	
	c->cd(2);
	h2->Draw("lego2"); h2->GetZaxis()->SetRangeUser(1.e-9, 5.e-2); gPad->SetLogz(); gPad->Modified(); gPad->Update();
	return chi2;
}

