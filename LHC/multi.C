Double_t fs1(Double_t *x, Double_t* dum);
Double_t fs2(Double_t *x, Double_t* dum);
Double_t fs3(Double_t *x, Double_t* dum);

static Int_t n;        // array size
static Double_t* as1;  // dsigma/db
static Double_t* as2;  // <N>(b)
static Double_t* ab;   // b
// |  PbPb           |    geom.   |      hard
// |  0   -  5   fm |   0 - 10 % |   0    -  42.8 %
// |  5   -  8.6 fm |  10 - 30 % |  42.8  -  83.5 %
// |  8.6 - 11.2 fm |  30 - 50 % |  83.5  -  96.8 %
// | 11.2 - 13.2 fm |  50 - 70 % |  96.8  -  99.6 %
// | 13.2 - 15.0 fm |  70 - 90 % |  99.6  -  99.9 %
// | 15.0 -  oo  fm |  90 -100 % |  99.9  - 100.0 %
 
void multi()
{
    const Float_t drift = 90;
    
//  open file and get graph objects
    TFile*  file = new TFile("DsigmaDb.root", "read");
    TGraph* gs1 =  gAlice = (TGraph*)(file->Get("Graph;1"));
    TGraph* gs2 =  gAlice = (TGraph*)(file->Get("Graph;2"));
// initialize arrays
    n   = gs1->GetN();
    ab  = gs1->GetX();
    as1 = gs1->GetY();
    as2 = gs2->GetY();
//
// function
//
// dsigma_mb/db 
    TF1* tf1 = new TF1("tf1", fs1, 0., ab[n-1], 0);
// dsigma_hard/db
    TF1* tf2 = new TF1("tf2", fs3, 0., ab[n-1], 0);

    Float_t db = 20./100.;
    for (Int_t i = 0; i< 100; i++) 
    {
	Float_t bb = Float_t(i)*db;
	Float_t s  = tf1->Integral(0,bb);
	Float_t sh = tf2->Integral(0,bb);
	printf("\n %d %f %f %f %f" , i, bb, s, s/7.769*100., sh/4.327*100.);
    }
    
    
// histogram for random distribution
    TH1F* hs1 = new TH1F("hs1", "dSigma/db", n, 0,ab[n-1]);
    hs1->Eval(tf1);
//    hs1->Draw();
    TH1F* hmult  = new TH1F("hmult",  "Multiplicity", 100, 0, 10000);
    TH1F* hmulte = new TH1F("hmulte", "Extra Multiplicity", 100, 0, 10000);
//     
    Float_t t1, t2;
    Double_t b, mult, dm, multe;
    Int_t i;
    Double_t *par;
    Double_t scale = 8000./as2[0];
    Double_t dT    = 125./2.;
    
    for (Int_t i=0; i < 1000000; i++)
    {
	// impact parameters
	b = hs1->GetRandom();
	// multiplicity
        mult = scale*(Float_t) fs2(&b,par);
	mult += gRandom->Gaus(0, TMath::Sqrt(mult));
	
	hmult->Fill(Float_t(mult));
	t1 = 0;
	t2 = 0;
	multe = 0;
	Double_t arg;
/*	
	while (1) {
	    arg=0.;
	    while(arg==0.) arg = gRandom->Rndm();
	    t1 -= dT*TMath::Log(arg); // (musec)   
	    if (t1 > drift) break;
	    b = hs1->GetRandom();
	    mult = scale*(Float_t) fs2(&b,par);
	    mult += gRandom->Gaus(0, TMath::Sqrt(mult));
//	    multe+=mult*(drift-t1)/drift;
	    multe+=mult;
	}

	while (1) {
	    arg=0.;
	    while(arg==0.) arg = gRandom->Rndm();
	    t2 -= dT*TMath::Log(arg); // (musec)   
	    if (t2 > drift) break;
	    b = hs1->GetRandom();
	    mult = scale*(Float_t) fs2(&b,par);
	    mult += gRandom->Gaus(0, TMath::Sqrt(mult));
//	    multe+=mult*(drift-t2)/drift;
	    multe+=mult;
	}
*/
	hmulte->Fill(Float_t(multe));
    }
    Double_t scale = 1./hmult->Integral();
    hmult->Scale(scale/100.);
    printf("\n scale %f \n ", scale);
    
    
    TCanvas *c1 = new TCanvas("c1","Canvas 1",400,10,600,700);
//    c1->Divide(1,2);
//    c1->cd(1);
    hmult->Draw();
//    c1->cd(2);
//   hmulte->Draw();
    
    
    TCanvas *c2 = new TCanvas("c2","Canvas 2",400,10,600,700);
    c2->Divide(1,2);
    c2->cd(1);
    tf1->Draw();
    tf1->GetHistogram()->SetXTitle("b (fm)");
    tf1->GetHistogram()->SetYTitle("d#sigma /db (barn/fm)");    
    c2->cd(2);
    tf2->Draw();
    tf2->GetHistogram()->SetXTitle("b (fm)");
    tf2->GetHistogram()->SetYTitle("d#sigma /db (a.u.)");    
    c2->Update();
    out = new TFile("mult.root", "recreate");
    hmult->Write();
    out->Close();
    
}


Double_t fs1(Double_t* x, Double_t* dum)
{
    Double_t xx = x[0];
    Int_t i;
    Double_t y;
    
    for (i=0; i<n; i++) {
	
	if (xx <ab[i]) {
	    y = as1[i-1]+(xx-ab[i-1])/(ab[i]-ab[i-1])*(as1[i]-as1[i-1]);
	    break;
	}
    }
    return y;
}


Double_t fs2(Double_t *x, Double_t* dum)
{
    Double_t xx = x[0];
    Int_t i;
    for (i=0; i<n; i++) {
	if (xx <ab[i]) {
	    Double_t y = as2[i-1]+(xx-ab[i-1])/(ab[i]-ab[i-1])*(as2[i]-as2[i-1]);
	    break;
	}
    }
    return y;
}

Double_t fs3(Double_t *x, Double_t* dum)
{
    Double_t xx = x[0];
    Int_t i;
    Double_t y1, y2;
    
    for (i=0; i<n; i++) {
	if (xx <ab[i]) {
	    y1 = as1[i-1]+(xx-ab[i-1])/(ab[i]-ab[i-1])*(as1[i]-as1[i-1]);
	    break;
	}
    }

    for (i=0; i<n; i++) {
	if (xx <ab[i]) {
	    Double_t y2 = as2[i-1]+(xx-ab[i-1])/(ab[i]-ab[i-1])*(as2[i]-as2[i-1]);
	    break;
	}
    }


    return y1*y2;
}


