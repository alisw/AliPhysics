void abso_dedx()
{
    TGeant3 *geant3 = (TGeant3*)gMC;
    Float_t   tkin[100], valPb[100], valW[100], val1[100], val2[100];
    Float_t   val3[100], val4[100], val5[100], pcut[5];
    Int_t ipart=5;
    char chmeca[4];
    strcpy(chmeca,"LOSS");
    Int_t ixst, i;
    Int_t kdim=100;
    
    for (i=0; i< kdim; i++) {
	tkin[i]=Float_t(i+1)*2.;
    }
//  Carbon
    geant3->Gftmat( 4, ipart, chmeca, kdim, tkin, val1, pcut, ixst);
//  Concrete
    geant3->Gftmat(25, ipart, chmeca, kdim, tkin, val2, pcut, ixst);
//  Ch2
    geant3->Gftmat(28, ipart, chmeca, kdim, tkin, val3, pcut, ixst);
//  Lead
    geant3->Gftmat(16, ipart, chmeca, kdim, tkin, val4, pcut, ixst);
//  W
    geant3->Gftmat(13, ipart, chmeca, kdim, tkin, val5, pcut, ixst);

    for (i=0; i< kdim; i++) {
	valPb[i]=(225.*val1[i]+153.*val2[i]+15.*val3[i]+20*val4[i])/1000.;
    }

    for (i=0; i< kdim; i++) {
	valW[i]=(225.*val1[i]+153.*val2[i]+35*val5[i])/1000.;
    }


    TGraph*  dedx1 = new TGraph(kdim, tkin, valPb);
    TGraph*  dedx2 = new TGraph(kdim, tkin, valW);

    TCanvas *c1=new TCanvas("c1","dedx",400,10,600,700);
    dedx1->SetFillColor(42);
    dedx1->SetMarkerColor(4);
    dedx1->SetMarkerStyle(21);
    dedx1->Draw("AC");
    dedx1->GetHistogram()->SetXTitle("Kinetic Energy (GeV)");
    dedx1->GetHistogram()->SetYTitle("Mean Energy Loss (GeV)"); 

    TCanvas *c2=new TCanvas("c2","dedx",400,10,600,700);
    dedx2->SetFillColor(42);
    dedx2->SetMarkerColor(4);
    dedx2->SetMarkerStyle(21);
    dedx2->Draw("AC");
    dedx2->GetHistogram()->SetXTitle("Kinetic Energy (GeV)");
    dedx2->GetHistogram()->SetYTitle("Mean Energy Loss (GeV)"); 

 }
