void MUONtestsim(Int_t ncham, Int_t evNumber1=0,Int_t evNumber2=0) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////

// Dynamically link some shared libs                    


    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("loadlibs.C");
	loadlibs();
    }

// Connect the Root Galice file containing Geometry, Kine and Hits

    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");

    
    if (!file) {
	printf("\n Creating galice.root \n");
	file = new TFile("galice.root");
    } else {
	printf("\n galice.root found in file list");
    }
    file->ls();
    
// Get AliRun object from file or create it if not on file
    if (!gAlice) {
	gAlice = (AliRun*)(file->Get("gAlice"));
	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) {
	    printf("\n create new gAlice object");
	    gAlice = new AliRun("gAlice","Alice test program");
	}
    }
    
//  Create some histograms


   TH1F *ccharge1 = new TH1F("ccharge1","Cluster Charge distribution"
			   ,100,0.,5000.);
   TH1F *pcharge1 = new TH1F("pcharge1","Pad      Charge distribution"
			   ,100,0.,200.);
   TH1F *xresid1  = new TH1F("xresid1","x-Residuals"
			   ,100,-10.,10);
   TH1F *yresid1  = new TH1F("yresid1","y-Residuals"
			   ,100,-0.2,0.2);
   TH1F *npads1   = new TH1F("npads1" ,"Pads per Hit"
			   ,20,-0.5,19.5);
   TH2F *xresys1  = new TH2F("xresys1","x-Residuals systematics"
			   ,50,-2.0,2.0,100,-1.0,1.0);
   TH2F *yresys1  = new TH2F("yresys1","y-Residuals systematics"
			   ,50,-0.4,0.4,100,-0.1,0.1);

   TH1F *ccharge2 = new TH1F("ccharge2","Cluster Charge distribution"
			   ,100,0.,5000.);
   TH1F *pcharge2 = new TH1F("pcharge2","Pad      Charge distribution"
			   ,100,0.,200.);
   TH1F *xresid2  = new TH1F("xresid2","x-Residuals"
			   ,100,-1,1);
   TH1F *yresid2  = new TH1F("yresid2","y-Residuals"
			   ,100,-10, 10.);
   TH1F *npads2   = new TH1F("npads2" ,"Pads per Hit"
			   ,20,-0.5,19.5);
   TH2F *xresys2  = new TH2F("xresys2","x-Residuals systematics"
			   ,50,-2.0,2.0,100,-1.0,1.0);
   TH2F *yresys2  = new TH2F("yresys2","y-Residuals systematics"
			   ,50,-0.4,0.4,100,-0.1,0.1);


   AliMUONChamber*  iChamber;
//
//   Loop over events 
//
   Int_t Nh=0;
   Int_t Nh1=0;
   for (Int_t nev=0; nev<= evNumber2; nev++) {
       cout << "nev         " << nev <<endl;
       Int_t nparticles = gAlice->GetEvent(nev);
       cout << "nparticles  " << nparticles <<endl;
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;

       AliMUON *MUON  = (AliMUON*) gAlice->GetDetector("MUON");
	   printf("\n track %d %d \n ", nev, MUON);

       TTree *TH = gAlice->TreeH();
       Int_t ntracks = TH->GetEntries();
       Int_t Nc=0;
       Int_t   npad[2];
       Float_t Q[2], xmean[2],ymean[2],xres[2],yres[2], xonpad[2], yonpad[2];
//
//   Loop over events
//
       for (Int_t track=0; track<ntracks;track++) {
	   
	   gAlice->ResetHits();
	   Int_t nbytes += TH->GetEvent(track);
	   if (MUON)  {
	       for(AliMUONHit* mHit=(AliMUONHit*)MUON->FirstHit(-1); 
		   mHit;
		   mHit=(AliMUONHit*)MUON->NextHit()) 
	       {
		   Int_t   nch   = mHit->fChamber;  // chamber number
		   Float_t x     = mHit->X();        // x-pos of hit
		   Float_t y     = mHit->Y();        // y-pos
		   Float_t z     = mHit->Z();
		   
//
//
		   iChamber = & MUON->Chamber(nch-1);
		   response=iChamber->ResponseModel();
//
//
		   if (nch == ncham) {
		       for (Int_t i=0; i<2; i++) {
			   xmean[i]=0;
			   ymean[i]=0;
			   Q[i]=0;
			   npad[i]=0;
		       }
		       
		       for (AliMUONPadHit* mPad=
				(AliMUONPadHit*)MUON->FirstPad(mHit,MUON->PadHits());
			    mPad;
			    mPad=(AliMUONPadHit*)MUON->NextPad(MUON->PadHits()))
		       {

			   Int_t nseg     = mPad->fCathode;   // segmentation 
			   Int_t nhit     = mPad->fHitNumber; // hit number
			   Int_t qtot     = mPad->fQ;         // charge
			   Int_t ipx      = mPad->fPadX;      // pad number on X
			   Int_t ipy      = mPad->fPadY;      // pad number on Y
			   Int_t iqpad    = mPad->fQpad;      // charge per pad
			   Int_t secpad   = mPad->fRSec;      // r-pos of pad
//
//   
			   segmentation=iChamber->SegmentationModel(nseg);
			   Int_t ind=nseg-1;
			   Float_t xg,yg,zg;
			   segmentation->GetPadC(ipx, ipy, xg, yg, zg);
			   Int_t sec = segmentation->Sector(ipx, ipy);
			   Float_t dpx=segmentation->Dpx(sec);
			   Float_t dpy=segmentation->Dpy(sec);
			   
			   if (nseg == 1) {
			       pcharge1->Fill(iqpad,(float) 1);
			       ccharge1->Fill(qtot ,(float) 1);
			   }  else {
			       pcharge2->Fill(iqpad,(float) 1);
			       ccharge2->Fill(qtot ,(float) 1);
			   }
// Calculate centre of gravity
//
			   if (iqpad == 0 && ind==0) {
			       printf("\n attention iqpad 0");
			   }
			   
			   if (iqpad >  0) {
			       Float_t xc,yc,zc;
			       npad[ind]++;
			       segmentation->GetPadC(ipx,ipy,xc,yc,zc);
			       xmean[ind]+=Float_t(iqpad*xc);
			       ymean[ind]+=Float_t(iqpad*yc);
			       Q[ind]+=iqpad;
			   }
			   
		       } //pad hit loop
		       for (Int_t i=0; i<2; i++) {
			   segmentation = iChamber->SegmentationModel(i+1);
			   if (Q[i] >0) {
			       xmean[i] = xmean[i]/Q[i];
			       xres[i]  = xmean[i]-x;
			       ymean[i] = ymean[i]/Q[i];
			       yres[i]  = ymean[i]-y;
// Systematic Error
//
			       Int_t icx, icy;
			       Float_t zonpad;
			       
			       segmentation->
				   GetPadI(x,y,z,icx,icy);
			       segmentation->
				   GetPadC(icx,icy,xonpad[i],yonpad[i],zonpad);
			       xonpad[i]-=x;
			       yonpad[i]-=y;
			   } // charge not 0
		       } // plane loop
		       xresid1->Fill(xres[0],(float) 1);
		       yresid1->Fill(yres[0],(float) 1);
		       npads1->Fill(npad[0],(float) 1);
		       if (npad[0] >=2 && Q[0] > 6 ) {
			   xresys1->Fill(xonpad[0],xres[0],(float) 1);
			   yresys1->Fill(yonpad[0],yres[0],(float) 1);
		       }
		       
		       xresid2->Fill(xres[1],(float) 1);
		       yresid2->Fill(yres[1],(float) 1);
		       npads2->Fill(npad[1],(float) 1);
		       if (npad[1] >=2 && Q[1] > 6) {
			   xresys2->Fill(xonpad[1],xres[1],(float) 1);
			   yresys2->Fill(yonpad[1],yres[1],(float) 1);
		       }
		   } // chamber 1
	       } // hit loop
	   } // if MUON
       } // track loop
   } // event loop 
   
//Create a canvas, set the view range, show histograms
   TCanvas *c1 = new TCanvas("c1","Charge and Residuals (1)",400,10,600,700);
   pad11 = new TPad("pad11"," ",0.01,0.51,0.49,0.99);
   pad12 = new TPad("pad12"," ",0.51,0.51,0.99,0.99);
   pad13 = new TPad("pad13"," ",0.01,0.01,0.49,0.49);
   pad14 = new TPad("pad14"," ",0.51,0.01,0.99,0.49);
   pad11->SetFillColor(11);
   pad12->SetFillColor(11);
   pad13->SetFillColor(11);
   pad14->SetFillColor(11);
   pad11->Draw();
   pad12->Draw();
   pad13->Draw();
   pad14->Draw();
   
   pad11->cd();
   ccharge1->SetFillColor(42);
   ccharge1->SetXTitle("ADC units");
   ccharge1->Draw();

   pad12->cd();
   pcharge1->SetFillColor(42);
   pcharge1->SetXTitle("ADC units");
   pcharge1->Draw();

   pad13->cd();
   xresid1->SetFillColor(42);
   xresid1->Draw();

   pad14->cd();
   yresid1->SetFillColor(42);
   yresid1->Draw();

   TCanvas *c2 = new TCanvas("c2","Cluster Size (1)",400,10,600,700);
   pad21 = new TPad("pad21"," ",0.01,0.51,0.49,0.99);
   pad22 = new TPad("pad22"," ",0.51,0.51,0.99,0.99);
   pad23 = new TPad("pad23"," ",0.01,0.01,0.49,0.49);
   pad24 = new TPad("pad24"," ",0.51,0.01,0.99,0.49);
   pad21->SetFillColor(11);
   pad22->SetFillColor(11);
   pad23->SetFillColor(11);
   pad24->SetFillColor(11);
   pad21->Draw();
   pad22->Draw();
   pad23->Draw();
   pad24->Draw();
   
   pad21->cd();
   npads1->SetFillColor(42);
   npads1->SetXTitle("Cluster Size");
   npads1->Draw();

   pad23->cd();
   xresys1->SetXTitle("x on pad");
   xresys1->SetYTitle("x-xcog");
   xresys1->Draw();

   pad24->cd();
   yresys1->SetXTitle("y on pad");
   yresys1->SetYTitle("y-ycog");
   yresys1->Draw();

   TCanvas *c3 = new TCanvas("c3","Charge and Residuals (2)",400,10,600,700);
   pad31 = new TPad("pad31"," ",0.01,0.51,0.49,0.99);
   pad32 = new TPad("pad32"," ",0.51,0.51,0.99,0.99);
   pad33 = new TPad("pad33"," ",0.01,0.01,0.49,0.49);
   pad34 = new TPad("pad34"," ",0.51,0.01,0.99,0.49);
   pad31->SetFillColor(11);
   pad32->SetFillColor(11);
   pad33->SetFillColor(11);
   pad34->SetFillColor(11);
   pad31->Draw();
   pad32->Draw();
   pad33->Draw();
   pad34->Draw();
   
   pad31->cd();
   ccharge2->SetFillColor(42);
   ccharge2->SetXTitle("ADC units");
   ccharge2->Draw();

   pad32->cd();
   pcharge2->SetFillColor(42);
   pcharge2->SetXTitle("ADC units");
   pcharge2->Draw();

   pad33->cd();
   xresid2->SetFillColor(42);
   xresid2->Draw();

   pad34->cd();
   yresid2->SetFillColor(42);
   yresid2->Draw();

   TCanvas *c4 = new TCanvas("c4","Cluster Size (2)",400,10,600,700);
   pad41 = new TPad("pad41"," ",0.01,0.51,0.49,0.99);
   pad42 = new TPad("pad42"," ",0.51,0.51,0.99,0.99);
   pad43 = new TPad("pad43"," ",0.01,0.01,0.49,0.49);
   pad44 = new TPad("pad44"," ",0.51,0.01,0.99,0.49);
   pad41->SetFillColor(11);
   pad42->SetFillColor(11);
   pad43->SetFillColor(11);
   pad44->SetFillColor(11);
   pad41->Draw();
   pad42->Draw();
   pad43->Draw();
   pad44->Draw();
   
   pad41->cd();
   npads2->SetFillColor(42);
   npads2->SetXTitle("Cluster Size");
   npads2->Draw();

   pad43->cd();
   xresys2->SetXTitle("x on pad");
   xresys2->SetYTitle("x-xcog");
   xresys2->Draw();

   pad44->cd();
   yresys2->SetXTitle("y on pad");
   yresys2->SetYTitle("y-ycog");
   yresys2->Draw();
   
}
