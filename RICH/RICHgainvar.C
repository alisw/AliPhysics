void RICHgainvar (Int_t evNumber1=0,Int_t evNumber2=0) 
{

  gClassTable->GetID("AliRun");


  // Dynamically link some shared libs
 
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }else {
    delete gAlice;
    gAlice = 0;
  }
  
  gAlice=0;
  
  // Connect the Root Galice file containing Geometry, Kine and Hits
  
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
  if (!file) file = new TFile("galice.root","UPDATE");
  
  // Get AliRun object from file or create it if not on file
  
  if (!gAlice) {
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
  }
  else {
    delete gAlice;
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
  }
  
  AliRICH *RICH  = (AliRICH*)gAlice->GetDetector("RICH");
  
  TH1F *gainX = new TH1F("gainX","Gain along the wires",40,0,40);
  TH1F *gainIn = new TH1F("gainIn","Gain near the center of photocathode",300,0,300);
  TH1F *gainOut = new TH1F("gainOut","Gain near the edge of photocathode",300,0,300);
  TH1F *gainPhi1 = new TH1F("gainPhi1","Gain variation along wires, 0-30 degrees",300,0,300);
  TH1F *gainPhi2 = new TH1F("gainPhi2","Gain variation along wires, 30-60 degrees",300,0,300);
  TH1F *gainPhi3 = new TH1F("gainPhi3","Gain variation along wires, 60-90 degrees",300,0,300);
  TH1F *gainPhi4 = new TH1F("gainPhi4","Gain variation along wires, 90-120 degrees",300,0,300);
  TH1F *gainPhi5 = new TH1F("gainPhi5","Gain variation along wires, 120-150 degrees",300,0,300);
  TH1F *gainPhi6 = new TH1F("gainPhi6","Gain variation along wires, 150-180 degrees",300,0,300);
  TH1F *gainPhi7 = new TH1F("gainPhi7","Gain variation along wires, 180-210 degrees",300,0,300);
  TH1F *gainPhi8 = new TH1F("gainPhi8","Gain variation along wires, 210-240 degrees",300,0,300);
  TH1F *gainPhi9 = new TH1F("gainPhi9","Gain variation along wires, 240-270 degrees",300,0,300);
  TH1F *gainPhi10 = new TH1F("gainPhi10","Gain variation along wires, 270-300 degrees",300,0,300);
  TH1F *gainPhi11 = new TH1F("gainPhi11","Gain variation along wires, 300-330 degrees",300,0,300);
  TH1F *gainPhi12 = new TH1F("gainPhi12","Gain variation along wires, 330-360 degrees",300,0,300);
  TH1F *gainPhi = new TH1F("gainPhi","Gain variation along phi",36,0,360);

  //   Start loop over events 

  for (int nev=0; nev<= evNumber2; nev++) {
    Int_t nparticles = gAlice->GetEvent(nev);
    
    
    printf ("\n**********************************\nProcessing Event: %d\n",nev);
    printf ("Particles       : %d\n\n",nparticles);
    if (nev < evNumber1) continue;
    if (nparticles <= 0) return;
       
    // Get pointers to RICH detector and Hits containers
       
    
    TTree *TH = gAlice->TreeH(); 
    Int_t ntracks = TH->GetEntries();
    Int_t nent=(Int_t)gAlice->TreeR()->GetEntries();
    gAlice->TreeR()->GetEvent(nent-1);
    TClonesArray *Rawclusters = RICH->RawClustAddress(2);    //  Raw clusters branch
    //printf ("Rawclusters:%p",Rawclusters);
    Int_t nrawclusters = Rawclusters->GetEntriesFast();
    //printf (" nrawclusters:%d\n",nrawclusters);
    
    
    // Start loop on tracks in the hits containers
    Int_t Nc=0;
    for (Int_t track=0; track<ntracks;track++) {
      printf ("\nProcessing Track: %d\n",track);
      gAlice->ResetHits();
      Int_t nbytes += TH->GetEvent(track);
      if (RICH)
	TClonesArray *Hits = RICH->Hits();            // Hits branch
      //see hits distribution
      
      Int_t nhits = Hits->GetEntriesFast();
      printf("Hits            : %d\n",nhits);
      for (Int_t hit=0;hit<nhits;hit++) {
	mHit = (AliRICHHit*) Hits->UncheckedAt(hit);
	Float_t mx  = mHit->X();                    // x-pos of hit
	Float_t my  = mHit->Z();                    // y-pos
      }


      if (nrawclusters) {
	printf("Raw Clusters    : %d\n",nrawclusters);
	for (Int_t hit=0;hit<nrawclusters;hit++) {
	  rcHit = (AliRICHRawCluster*) Rawclusters->UncheckedAt(hit);
	  //Int_t nchamber = rcHit->fChamber;     // chamber number
	  //Int_t nhit = cHit->fHitNumber;        // hit number
	  Int_t qtot = rcHit->fQ;                 // charge
	  Int_t fx  =  rcHit->fX;                 // x-position
	  Int_t fy  =  rcHit->fY;                 // y-position
	  Int_t mult = rcHit->fMultiplicity;      // How many pads form the cluster
	  		    
	  Float_t radfid = TMath::Sqrt((TMath::Power(fx-mx,2)+(TMath::Power(fy-my,2))));
	  //printf("Radius:%f\n", radfid);
	  
	  
	  if(qtot<200 && radfid>=9 && radfid<=15)
	    {
	      if (fy>=0)
		gainX->Fill(fy-40,(float) qtot);
	      else
		gainX->Fill(fy+40,(float) qtot);
	      
	      if(fy<=0 && fy>=-4)
		gainOut->Fill(qtot,(float) 1);
	      if(fy<=-11 && fy>=-15)
		gainIn->Fill(qtot,(float) 1);
	      
	      Float_t phi = TMath::ATan2(fx-mx, fy-my);
	      Float_t phi_deg = phi/TMath::Pi()*180;
	      
	      if ((fx-mx)<0)
		{
		  phi_deg = 360+phi_deg;
		  //printf("X:%f, Y:%f, Phi:%f\n",-fy+my,fx-mx,phi_deg);
		}
	      else
		{
		  //printf("X:%f, Y:%f, Phi:%f\n",-fy+my,fx-mx,phi_deg);
		}
	      
	      
	      if (phi_deg>=0 && phi_deg<30)
		gainPhi1->Fill(qtot,(float) 1);
	      if (phi_deg>=30 && phi_deg<60)
		gainPhi2->Fill(qtot,(float) 1);
	      if (phi_deg>=60 && phi_deg<90)
		gainPhi3->Fill(qtot,(float) 1);
	      if (phi_deg>=90 && phi_deg<120)
		gainPhi4->Fill(qtot,(float) 1);
	      if (phi_deg>=120 && phi_deg<150)
		gainPhi5->Fill(qtot,(float) 1);
	      if (phi_deg>=150 && phi_deg<180)
		gainPhi6->Fill(qtot,(float) 1);
	      if (phi_deg>=180 && phi_deg<210)
		gainPhi7->Fill(qtot,(float) 1);
	      if (phi_deg>=210 && phi_deg<240)
		gainPhi8->Fill(qtot,(float) 1);
	      if (phi_deg>=240 && phi_deg<270)
		gainPhi9->Fill(qtot,(float) 1);
	      if (phi_deg>=270 && phi_deg<300)
		gainPhi10->Fill(qtot,(float) 1);
	      if (phi_deg>=300 && phi_deg<330)
		gainPhi11->Fill(qtot,(float) 1);
	      if (phi_deg>=330 && phi_deg<360)
		gainPhi12->Fill(qtot,(float) 1);
	      
	    }
	  
	  //printf("Y:%d, Q:%d %d in %d\n",fy, qtot,hit,nrawclusters);
	}
      }
    }
  }

  
  TCanvas *c16 = new TCanvas("c16","Gain Variation over Phi",50,50,1200,700);
  c16->Divide(4,3);
  //c16->SetFillColor(42);
  
  c16->cd(1);
  c16_1->SetLogy();
  gainPhi1->SetFillColor(5);
  gainPhi1->SetXTitle("ADC counts");
  gainPhi1->Fit("expo","M");
  expo->SetLineColor(2);
  expo->SetLineWidth(3);
  Float_t inv_slope1 = -1/expo->GetParameter(1);
  gainPhi->Fill(15, (float) inv_slope1);
  gainPhi1->Draw();
  
  c16->cd(2);
  c16_2->SetLogy();
  gainPhi2->SetFillColor(5);
  gainPhi2->SetXTitle("ADC counts");
  gainPhi2->Fit("expo","M");
  expo->SetLineColor(2);
  expo->SetLineWidth(3);
  Float_t inv_slope2 = -1/expo->GetParameter(1);
  gainPhi->Fill(45, (float) inv_slope2);
  gainPhi2->Draw();
  
  c16->cd(3);
  c16_3->SetLogy();
  gainPhi3->SetFillColor(5);
  gainPhi3->SetXTitle("ADC counts");
  gainPhi3->Fit("expo","M");
  expo->SetLineColor(2);
  expo->SetLineWidth(3);
  Float_t inv_slope3 = -1/expo->GetParameter(1);
  gainPhi->Fill(75, (float) inv_slope3);
  gainPhi3->Draw();
  
  c16->cd(4);
  c16_4->SetLogy();
  gainPhi4->SetFillColor(5);
  gainPhi4->SetXTitle("ADC counts");
  gainPhi4->Fit("expo","M");
  expo->SetLineColor(2);
  expo->SetLineWidth(3);
  Float_t inv_slope4 = -1/expo->GetParameter(1);
  gainPhi->Fill(105, (float) inv_slope4);
  gainPhi4->Draw();
  
  c16->cd(5);
  c16_5->SetLogy();
  gainPhi5->SetFillColor(5);
  gainPhi5->SetXTitle("ADC counts");
  gainPhi5->Fit("expo","M");
  expo->SetLineColor(2);
  expo->SetLineWidth(3);
  Float_t inv_slope5 = -1/expo->GetParameter(1);
  gainPhi->Fill(135, (float) inv_slope5);
  gainPhi5->Draw();
  
  c16->cd(6);
  c16_6->SetLogy();
  gainPhi6->SetFillColor(5);
  gainPhi6->SetXTitle("ADC counts");
  gainPhi6->Fit("expo","M");
  expo->SetLineColor(2);
  expo->SetLineWidth(3);
  Float_t inv_slope6 = -1/expo->GetParameter(1);
  gainPhi->Fill(165, (float) inv_slope6);
  gainPhi6->Draw();
  
  c16->cd(7);
  c16_7->SetLogy();
  gainPhi7->SetFillColor(5);
  gainPhi7->SetXTitle("ADC counts");
  gainPhi7->Fit("expo","M");
  expo->SetLineColor(2);
  expo->SetLineWidth(3);
  Float_t inv_slope7 = -1/expo->GetParameter(1);
  gainPhi->Fill(195, (float) inv_slope7);
  gainPhi7->Draw();
  
  c16->cd(8);
  c16_8->SetLogy();
  gainPhi8->SetFillColor(5);
  gainPhi8->SetXTitle("ADC counts");
  gainPhi8->Fit("expo","M");
  expo->SetLineColor(2);
  expo->SetLineWidth(3);
  Float_t inv_slope8 = -1/expo->GetParameter(1);
  gainPhi->Fill(225, (float) inv_slope8);
  gainPhi8->Draw();
  
  c16->cd(9);
  c16_9->SetLogy();
  gainPhi9->SetFillColor(5);
  gainPhi9->SetXTitle("ADC counts");
  gainPhi9->Fit("expo","M");
  expo->SetLineColor(2);
  expo->SetLineWidth(3);
  Float_t inv_slope9 = -1/expo->GetParameter(1);
  gainPhi->Fill(255, (float) inv_slope9);
  c16->Update();
  gainPhi9->Draw();
  
  
  c16->cd(10);
  c16_10->SetLogy();
  gainPhi10->SetFillColor(5);
  gainPhi10->SetXTitle("ADC counts");
  gainPhi10->Fit("expo","M");
  expo->SetLineColor(2);
  expo->SetLineWidth(3);
  Float_t inv_slope10 = -1/expo->GetParameter(1);
  gainPhi->Fill(285, (float) inv_slope10);
  gainPhi10->Draw();
  
  c16->cd(11);
  c16_11->SetLogy();
  gainPhi11->SetFillColor(5);
  gainPhi11->SetXTitle("ADC counts");
  gainPhi11->Fit("expo","M");
  expo->SetLineColor(2);
  expo->SetLineWidth(3);
  Float_t inv_slope11 = -1/expo->GetParameter(1);
  gainPhi->Fill(315, (float) inv_slope11);
  gainPhi11->Draw();
  
  c16->cd(12);
  c16_12->SetLogy();
  gainPhi12->SetFillColor(5);
  gainPhi12->SetXTitle("ADC counts");
  gainPhi12->Fit("expo","M");
  expo->SetLineColor(2);
  expo->SetLineWidth(3);
  Float_t inv_slope12 = -1/expo->GetParameter(1);
  gainPhi->Fill(335, (float) inv_slope12);
  gainPhi12->Draw();
  
  TCanvas *c17 = new TCanvas("c17","Gain Variation",50,50,600,700);
  c17->Divide(2,2);
  
  //c11->SetFillColor(42);
  
  c17->cd(1);
  c17_1->SetLogy();
  gainIn->SetFillColor(5);
  gainIn->SetXTitle("ADC counts");
  gainIn->Fit("expo","M");
  expo->SetLineColor(2);
  expo->SetLineWidth(3);
  gainIn->Draw();
  
  c17->cd(2);
  c17_2->SetLogy();
  gainOut->SetFillColor(5);
  gainOut->SetXTitle("ADC counts");
  gainOut->Fit("expo","M");
  expo->SetLineColor(2);
  expo->SetLineWidth(3);
  gainOut->Draw();
  
  c17->cd(3);
  gainX->SetFillColor(5);
  gainX->SetYTitle("Total ADC counts");
  gainX->SetXTitle("(pads)");
  gainX->Draw();
  
  c17->cd(4);
  gainPhi->SetFillColor(5);
  gainPhi->SetYTitle("Gain (ADC)");
  gainPhi->SetXTitle("phi (degrees)");
  gainPhi->Draw();

  printf("\nEnd of macro\n");
  printf("**********************************\n");
  
}
