{
    gSystem->Load("libPWG2LRC.so");
 	TH1D *gaus = new TH1D("gaus","source hist",1000,-5,5);
    gaus->FillRandom("gaus",10000);
    TH2D *source = new TH2D("NN","NN Test",100,0,100,100,0,100);
	Double_t x, y;
    for (Int_t i=0;i<10000;i++) 
	{   
	    x = 10 * (gaus->GetRandom() + 5);
	    y = 0.678 * (x + gaus->GetRandom() * 15);
	    //Create source 2D histogram with correlation coefficient 0.678
		source->Fill(x,y);
    }
    //Create TNN class encapsulated NN correlation algorithms
    //2D histogram pass into TNN constructor
	AliLRCNN final1("name", source);

  TCanvas *c1 = new TCanvas("c1","c1",800,1000);
	c1->Divide(3);
	c1->cd(1);
	source->Draw();
	c1->cd(2);
	final1.Draw_abs();
	c1->cd(3);
	final1.Draw_rel();
	c1->cd();

	AliLRCNN final2;
        final2.MakeHistogramm("name", source);

   TCanvas *c2 = new TCanvas("c2","c2",800,1000);
	c2->Divide(3);
	c2->cd(1);
	source->Draw();
	c2->cd(2);
	final2.Draw_abs();
	c2->cd(3);
	final2.Draw_rel();
	c2->cd();

}