{
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Wed Jun 16 15:35:08 1999 by ROOT version 2.21/07)
//   from TTree ntuple/Reconst ntuple
//   found on file: reconst.root
//////////////////////////////////////////////////////////

//Reset ROOT and connect tree file
gROOT->Reset();

TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("reconst.root");
if (!f) {
     f = new TFile("reconst.root");
}

TTree *ntuple = (TTree*)gDirectory->Get("ntuple");

TPostScript ps("mass.ps",112);
//ps->Open("mass.ps");
//Declaration of leaves types
Int_t           ntrack=500;

Int_t           ievr;
Int_t           ntrackr;
Int_t           istatr[ntrack];
Int_t           isignr[ntrack];
Float_t         pxr[ntrack];
Float_t         pyr[ntrack];
Float_t         pzr[ntrack];
Float_t         zvr[ntrack];
Float_t         chi2r[ntrack];

//Set branch addresses
ntuple->SetBranchAddress("ievr",&ievr);
ntuple->SetBranchAddress("ntrackr",&ntrackr);
ntuple->SetBranchAddress("istatr",&istatr);
ntuple->SetBranchAddress("isignr",&isignr);
ntuple->SetBranchAddress("pxr",&pxr);
ntuple->SetBranchAddress("pyr",&pyr);
ntuple->SetBranchAddress("pzr",&pzr);
ntuple->SetBranchAddress("zvr",&zvr);
ntuple->SetBranchAddress("chi2r",&chi2r);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// ntuple->SetBranchStatus("*",0);  // disable all branches
// ntuple->SetBranchStatus("branchname",1);  // activate branchname

TH1F *h1 = new TH1F("h1","mass",240,0.,12.);
TH1F *h2 = new TH1F("h2","mass good muon",240,0.,12.);
// Upsilon  mass cut
TH1F *h5 = new TH1F("h5"," Mass Upsilon",70,8.,11.);
// J/psi  mass cut
//TH1F *h5 = new TH1F("h5"," Mass J/psi",120,1.,5.);
TH1F *h3 = new TH1F("h3","Pt",100,0.,20.);
TH1F *h4 = new TH1F("h4","Pt bon muon",100,0.,20.);

gStyle->SetOptFit();
//gStyle->SetOptStat(0000);
gStyle->SetOptStat(10001110);

Float_t mass = 1.;
Float_t amassmu = 0.10566;
Float_t Chi2Cut = 100.;
Float_t PtCut = 0.;
//Float_t MassMin = 8.96;
//Float_t MassMax = 9.96;
// Upsilon
Float_t MassMin = 8.5;
Float_t MassMax = 10.5;
// J/psi
//Float_t MassMin = 2.5;
//Float_t MassMax = 3.8;

//Float_t Chi2Cut = 999.;

Int_t nentries = ntuple->GetEntries();

Int_t nbytes = 0;
Int_t nbtrack[50]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

Int_t compt = 0;
Int_t Nres = 0;

//    Good muons
for (Int_t iev=0; iev<nentries; iev++) {
    nbytes += ntuple->GetEvent(iev);
//    ntuple->SetBranchStatus("main",1);
    ntuple->SetBranchStatus("*",1);

    //printf("\n");
    //  printf("ntrackr=%d\n",ntrackr);
    
    if(ntrackr < 50)   nbtrack[ntrackr]=nbtrack[ntrackr]+1 ;
    
     
    for (Int_t i=0; i<(ntrackr-1); i++) {

      Float_t px1 =  pxr[i] ;
      Float_t py1 =  pyr[i] ;
      Float_t pz1 =  pzr[i] ;
      Float_t Pt1 =  TMath::Sqrt(px1*px1+py1*py1) ;
      
//      if (chi2r[i] < Chi2Cut && istatr[i] == 1 && Pt1 > PtCut) {
      h4->Fill(Pt1);
      
      for (Int_t j=i+1; j<ntrackr; j++) {
	Float_t px2 =  pxr[j] ;
	Float_t py2 =  pyr[j] ;
	Float_t pz2 =  pzr[j] ;
	Float_t Pt2 =  TMath::Sqrt(px2*px2+py2*py2) ;
	if (chi2r[j]<Chi2Cut && istatr[j] == 1 && Pt2 > PtCut ) {
	    if (isignr[i]!=isignr[j]) {
		Float_t p12 = px1*px1+py1*py1+pz1*pz1;
		Float_t p22 = px2*px2+py2*py2+pz2*pz2;
		Float_t e1 = TMath::Sqrt(amassmu*amassmu+p12);
		Float_t e2 = TMath::Sqrt(amassmu*amassmu+p22);
		Float_t amass = TMath::Sqrt((e1+e2)*(e1+e2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2));
//		printf("amass=%f\n",amass);
		h2->Fill(amass);
		if (amass<11. && amass>8.) compt++;
	    }
	}
      }
    }




    //   All tracks
    for (Int_t i=0; i<(ntrackr-1); i++) {
      Float_t px1 =  pxr[i] ;
      Float_t py1 =  pyr[i] ;
      Float_t pz1 =  pzr[i] ;
      Float_t Pt1 =  TMath::Sqrt(px1*px1+py1*py1) ;
      if (chi2r[i] < Chi2Cut && Pt1 > PtCut) {
	h3->Fill(Pt1);
	for (Int_t j=i+1; j<ntrackr; j++) {
	  Float_t px2 =  pxr[j] ;
	  Float_t py2 =  pyr[j] ;
	  Float_t pz2 =  pzr[j] ;
	  Float_t Pt2 =  TMath::Sqrt(px2*px2+py2*py2) ;
	  if (chi2r[j]<Chi2Cut && Pt2 > PtCut) {
	    if (isignr[i]!=isignr[j]) {
	      Float_t p12 = px1*px1+py1*py1+pz1*pz1;
	      Float_t p22 = px2*px2+py2*py2+pz2*pz2;
	      Float_t e1 = TMath::Sqrt(amassmu*amassmu+p12);
	      Float_t e2 = TMath::Sqrt(amassmu*amassmu+p22);
	      Float_t amass = TMath::Sqrt((e1+e2)*(e1+e2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2));
	      //			printf("amass=%f\n",amass);
	      h1->Fill(amass);
	      if (amass > MassMin && amass < MassMax) {
		//		if (ntrackr==2) h5->Fill(amass);
		h5->Fill(amass);
//		if (amass < 9.715 && amass > 9.205) 
		if (amass < 9.766 && amass > 9.154) 
		    {
			printf("amass=%f\n",amass);
			Nres++;
		    }
	      }
	    }
	  }
	}
      }
    }
}

//    printf("compt= %d\n",compt);

for (Int_t i = 0 ; i<50 ; i++) 
{
    
  printf(" nb of events= %d with %d tracks \n",nbtrack[i], i) ;
  
}
g1= new TF1("g1","gaus",9.30,9.75) ;
//g1= new TF1("g1","gaus",8.7,10.2) ;
//g1= new TF1("g1","gaus",2.95,3.25) ;
h5->Fit("g1","R0")      ;
h5->SetXTitle("Mass (GeV/c^2!)");
h5->Draw();
printf("Nombre de resonances : %d\n",Nres);

ps->Close();
}








