TStopwatch timer;
const Int_t nt = 100;
//___________________________________________________________
void pdfIO(const Int_t method = 0)
{
// Test interpolator IO response for several data dimensions.
// method = 0 : use COG interpolator
// method = 1 : use INT interpolator

	Float_t tw, tr;
	TGraph *gw = new TGraph(5);
	gw->SetMarkerStyle(24);gw->SetMarkerColor(2);
	TGraph *gr = new TGraph(5);
	gr->SetMarkerStyle(25);gr->SetMarkerColor(4);
	for(int idim = 1; idim<6; idim++){
		tw = interpolWrite(method, idim);
		gw->SetPoint(idim-1, idim, tw);
		tr = interpolRead(idim);
		gr->SetPoint(idim-1, idim, tr);
	}
	gw->Draw("apl");
	gr->Draw("pl");
}

//___________________________________________________________
Float_t interpolWrite(const Int_t method, const Int_t ndim)
{
	if(!gSystem->FindFile(".", "5D_Gauss.root")) build(5, 1000000);
	TFile::Open("5D_Gauss.root");
	printf("\nWriting %dD interpolator ...\n", ndim);
	TTree *db = (TTree*)gFile->Get("db");
	TString var = "x0"; for(int idim=1; idim<ndim; idim++) var += Form(":x%d", idim);
	printf("\tBuilding interpolator ...\n");
	TKDPDF pdf(db, var.Data(), "", 400, 100000);
	pdf.SetInterpolationMethod(method);
	pdf.SetStore();

	
	printf("\tSetting up the interpolator ...\n");
	Float_t *c, val, v_err;
	Double_t *p = new Double_t[ndim], res, r_err;
	timer.Start(kTRUE);
	for(int inode=0; inode<pdf.GetNTNodes(); inode++){
		pdf.GetCOGPoint(inode, c, val, v_err);
		for(int idim=0; idim<ndim; idim++) p[idim] = (Double_t)c[idim];
// 		printf("Evaluate for node [%d] ", inode);
// 		for(int idim=0; idim<ndim; idim++) printf("%5.3f ", p[idim]); printf("\n");
		for(int it=0; it<nt; it++) pdf.Eval(p, res, r_err, kTRUE);
//		printf("R = %6.4f +- %6.4f\n", res, r_err);
	}
	timer.Stop();
	
	Float_t time = 1.E3*timer.CpuTime()/float(nt);
	printf("\tFinish interpolation in %6.2f [ms]\n", time);
	
	printf("\tSaving interpolator ...\n");
	TFile *fi = TFile::Open(Form("%dD_Interpolator.root", ndim), "RECREATE");
	pdf.Write(Form("%dDgauss", ndim));
	fi->Close();
	delete fi;
	
	delete [] p;

	return time;
}

//___________________________________________________________
Float_t interpolRead(const Int_t ndim)
{
	printf("Reading %dD interpolator ...\n", ndim);
	TFile::Open(Form("%dD_Interpolator.root", ndim));
	TKDPDF *pdf = (TKDPDF*)gFile->Get(Form("%dDgauss", ndim));
	gFile->Close();
	
	printf("\tDoing interpolation ...\n");
	Float_t *c, v, ve;
	Double_t *p = new Double_t[ndim], r, re;
	timer.Start(kTRUE);	
	for(int ip=0; ip<pdf->GetNTNodes(); ip++){
		pdf->GetCOGPoint(ip, c, v, ve);
		for(int idim=0; idim<ndim; idim++) p[idim] = (Double_t)c[idim];
		for(int it=0; it<nt; it++) pdf->Eval(p, r, re);
	}
	timer.Stop();
	Float_t time = 1.E3*timer.CpuTime()/float(nt);
	printf("\tFinish interpolation in %6.2f [ms]\n", time);
	delete [] p;
	return time;
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


