void PlotPressure()
{
    FILE* file;
    file=fopen("gasPressure.dat","r");
    Float_t z1[20];
    Float_t g1[20], g2[20], g3[20], g4[20], g5[20];
    Float_t z2[21];
    Float_t h1[21], h2[21], h3[21], h4[21], h5[21];

    char c[45];
    Float_t z;
    
    for (Int_t i = 0; i < 20; i++)
    {
	fscanf(file, "%f %f %f %f %f %f", &z, 
	       &g1[i], &g2[i], &g3[i], &g4[i], &g5[i]);
//	printf("%d %f %f %f %f %f %f \n", i, z1[i], 
//		       g1[i][0], g1[i][1], g1[i][2], g1[i][3], g1[i][4]);
	if (i > 0) {
	    z1[i] = z1[i-1] + z;
	} else {
	    z1[i] = 20.;
	}
    }

    
    for (Int_t i = 0; i < 21; i++)
    {
	fscanf(file, "%f %f %f %f %f %f", &z, 
	       &h1[i], &h2[i], &h3[i], &h4[i], &h5[i]);
	if (i > 0) {
	    z2[i] = z2[i-1] + z;
	} else {
	    z2[i] = 20.;
	}
    }

//
// 
    TCanvas *c1 = new TCanvas("c1","Gas Pressure Beam 1", 200, 10, 700, 500);
    gPad->SetLogy();
    
    TGraph* gr1  = new TGraph(20, z1, g1);
    gr1->SetMaximum(1e17);
    gr1->SetMinimum(1e11);    
    gr1->SetLineColor(1);
    gr1->SetTitle("Ring 1: Beginning of Run");
    
    TGraph* gr2  = new TGraph(20, z1, g3);
    gr2->SetLineColor(2);

    TGraph* gr3  = new TGraph(20, z1, g5);
    gr3->SetLineColor(4);

    gr1->Draw("AL");
    gr2->Draw("L");
//    gr3->Draw("L");
    text();

//
// 
    TCanvas *c2 = new TCanvas("c2","Gas Pressure Beam 1", 200, 10, 700, 500);
    gPad->SetLogy();
    
    TGraph* gr4  = new TGraph(20, z1, g2);
    gr4->SetMaximum(1e17);
    gr4->SetMinimum(1e11);    
    gr4->SetLineColor(1);
    gr4->SetTitle("Ring 1");
    
    TGraph* gr5  = new TGraph(20, z1, g4);
    gr5->SetLineColor(2);

    TGraph* gr6  = new TGraph(20, z1, g5);
    gr6->SetLineColor(4);

    gr4->Draw("AL");
    gr5->Draw("L");
    gr6->Draw("L");
    text();

//
// 
    TCanvas *c3 = new TCanvas("c3","Gas Pressure Beam 2", 200, 10, 700, 500);
    gPad->SetLogy();
    
    TGraph* hr1  = new TGraph(21, z2, h1);
    hr1->SetMaximum(1e17);
    hr1->SetMinimum(1e11);    
    hr1->SetLineColor(1);
    hr1->SetTitle("Ring 2: Beginning of Run");
    
    TGraph* hr2  = new TGraph(21, z2, h3);
    hr2->SetLineColor(2);

    TGraph* hr3  = new TGraph(21, z2, h5);
    hr3->SetLineColor(4);

    hr1->Draw("AL");
    hr2->Draw("L");
//    hr3->Draw("L");
    text();

//
// 
    TCanvas *c4 = new TCanvas("c4","Gas Pressure Beam 2", 200, 10, 700, 500);
    gPad->SetLogy();
    
    TGraph* hr4  = new TGraph(21, z2, h2);
    hr4->SetMaximum(1e17);
    hr4->SetMinimum(1e11);    
    hr4->SetLineColor(1);
    hr4->SetTitle("Ring 2");
    
    TGraph* hr5  = new TGraph(21, z2, h4);
    hr5->SetLineColor(2);

    TGraph* hr6  = new TGraph(21, z2, h5);
    hr6->SetLineColor(4);

    hr4->Draw("AL");
    hr5->Draw("L");
    hr6->Draw("L");
    text();
//
// One example of interaction rate calculation
//
// 3rd year, full intensity
//
    Float_t r[20];
    const Float_t crossSection = 0.094e-28;   //  m^2
    const Float_t pFlux        = 1.e11/25.e-9; // 1/s
    
    for (Int_t i= 0; i < 20; i++) 
    {
	r[i] =  g5[i] * crossSection * pFlux; // 1/m/s 
    }

    TCanvas *c5 = new TCanvas("c5","Interaction Rate Beam 1, 3rd year", 200, 10, 700, 500);
    gPad->SetLogy();
    
    TGraph* rr1  = new TGraph(20, z1, r);
    rr1->SetMaximum(1e6);
    rr1->SetMinimum(1e1);    
    rr1->SetLineColor(4);
    rr1->SetTitle("Ring 1");
    rr1->Draw("AL");
}

void text()
{
    
    TPave *pave = new TPave(194.619,11.0495,293.103,12.6485,4,"br");
    pave->SetFillColor(18);
    pave->Draw();
    TLine *line = new TLine(199.702,11.8109,232.106,11.8109);
    line->SetLineColor(2);
    line->Draw();
    line = new TLine(198.431,12.2868,232.106,12.2868);
    line->Draw();
    line = new TLine(199.066,11.3731,232.742,11.3731);
    line->SetLineColor(4);
    line->Draw();
    line = new TLine(215.586,13.3338,215.586,13.3147);
    line->Draw();
    tex = new TLatex(239.096,11.6777,"2nd year");
    tex->SetTextSize(0.05);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(236.554,12.1536,"1st year");
    tex->SetTextSize(0.0507614);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(239.096,11.2589,"3rd year");
    tex->SetTextSize(0.05);
    tex->SetLineWidth(2);
    tex->Draw();
    c1->Modified();
    c1->cd();
}
