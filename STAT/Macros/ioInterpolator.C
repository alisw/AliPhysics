const int ndim = 2;

//___________________________________________________________
void writeInterpolator()
{
	if(!gSystem->FindFile(".", "6Ddb.root")) build(6, 1000000);
	TFile::Open("6Ddb.root");
	TKDInterpolator in(db, "x0:x1", "", 400, 100000);
	in.SetIntInterpolation();
	in.SetSetStore();

	Float_t c[ndim], val, v_err;
	Double_t p[ndim], res, r_err;
	TGraph2DErrors *g= new TGraph2DErrors(in.GetNTerminalNodes());g->SetMarkerStyle(7);
	for(int inode=0; inode<in.GetNTerminalNodes(); inode++){
		in.GetCOGPoint(inode, c, val, v_err);
		for(int idim=0; idim<ndim; idim++){
			p[idim] = (Double_t)c[idim];
			//printf("%f ", p[idim]);
		}
		//printf("\n");
		in.Eval(p, res, r_err);
		g->SetPoint(inode, p[0], p[1], res);
		g->SetPointError(inode, 0., 0., r_err);
	}
	g->Draw("ap");

	TFile *fi = TFile::Open(Form("%dD_interpolator.root", ndim), "RECREATE");
	in.Write(Form("%dDgauss", ndim));
	fi->Close();
	delete fi;
}

//___________________________________________________________
void readInterpolator()
{
	TFile::Open(Form("%dD_interpolator.root", ndim));
	TKDInterpolator *in = (TKDInterpolator*)gFile->Get(Form("%dDgauss", ndim));
	in->Dump();
	//in->GetStatus();
}

//___________________________________________________________
void build(const int ndim, const int npoints)
{
	printf("building db ...\n");
	Float_t data[ndim];
	TFile *f = new TFile(Form("%dDdb.root", ndim), "RECREATE");
	TTree *t = new TTree("db", Form("%dD data base for kD statistics", ndim));
	for(int idim=0; idim<ndim; idim++) t->Branch(Form("x%d", idim), &data[idim], Form("x%d/F", idim));

	for (Int_t ip=0; ip<npoints; ip++){
		for (Int_t id=0; id<ndim; id++) data[id]= gRandom->Gaus();
		t->Fill();
	}

	t->Write();
	f->Close();
	delete f;
}
