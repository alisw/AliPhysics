void dedx(Int_t imate)
{
    TGeant3 *geant3 = (TGeant3*)gMC;
    Float_t   tkin[100], value[100], pcut[5];
    Int_t ipart=5;
    char chmeca[4];
    strcpy(chmeca,"LOSS");
    Int_t ixst, i;
    Int_t kdim=100;
    
    for (i=0; i< kdim; i++) {
	tkin[i]=Float_t(i)*1.+1;
    }
    geant3->Gftmat(imate, ipart, chmeca, kdim, tkin, value, pcut, ixst);
    for (i=0; i< kdim; i++) {
	printf("\n Energy %f dE/dx %f", tkin[i], value[i]);
    }
    TGraph*  dedx = new TGraph(kdim, tkin, value);
    TCanvas *c1=new TCanvas("c1","dedx",400,10,600,700);
    dedx->SetFillColor(42);
    dedx->SetMarkerColor(4);
    dedx->SetMarkerStyle(21);
    dedx->Draw("AC");
    dedx->GetHistogram()->SetXTitle("Kinetic Energy (GeV)");
    dedx->GetHistogram()->SetYTitle("dE/dx   "); 
 }










