void drawMCPesd()
{
  /*
 Float_t xc[28] = {-8.95, -3.05, 3.05, 8.95, -14.85,
                    -8.95, -3.05, 3.05, 8.95,  14.85,
                    -14.85, -8.95, 8.95, 14.85, -14.85, 
                    -8.95, 8.95, 14.85, -14.85,-8.95,
                    -3.05, 3.05, 8.95, 14.85, -8.95,
                    -3.05, 3.05, 8.95 };
  
  Float_t yc[28] = {14.85, 14.85,14.85,14.85,8.95,
                    8.95,8.95,8.95,8.95,8.95, 
                    3.05, 3.05,3.05,3.05,-3.05,
                    -3.05,-3.05,-3.05,-8.95,-8.95,
                    -8.95,-8.95,-8.95,-8.95,-14.85, 
                    -14.85,-14.85,-14.85};
  */
  Double_t crad = 82.;		//define concave c-side radius here

  Double_t dP = 3.31735114408; 	// Work in Progress side length

  Float_t pstartC[3] = {20, 20 ,5};
  //uniform angle between detector faces==
  Double_t btta = 2*TMath::ATan(dP/crad);
 
  //get noncompensated translation data
  Double_t grdin[6] = {-3, -2, -1, 1, 2, 3};  
  Double_t gridpoints[6];
  for(Int_t i = 0; i < 6; i++){
    gridpoints[i] = crad*TMath::Sin((1 - 1/(2*TMath::Abs(grdin[i])))*grdin[i]*btta);
  } 

  std::vector<Double_t> xi,yi;
    
  for(Int_t j = 5; j >= 0; j--){
      for(Int_t i = 0; i < 6; i++){
          if(!(((j == 5 || j == 0) && (i == 5 || i == 0)) || 
               ((j == 3 || j == 2) && (i == 3 || i == 2))))
          {
              xi.push_back(gridpoints[i]);
              yi.push_back(gridpoints[j]);
          }
      }
      
  }
  
  Double_t zi[28];
  for(Int_t i = 0; i < 28; i++) {
    zi[i] = TMath::Sqrt(TMath::Power(crad, 2) - TMath::Power(xi[i], 2) - TMath::Power(yi[i], 2));
  }

  //get rotation data
  Double_t ac[28], bc[28], gc[28];
  for(Int_t i = 0; i < 28; i++) {
    ac[i] = TMath::ATan(yi[i]/xi[i]) - TMath::Pi()/2 + 2*TMath::Pi();
    if(xi[i] < 0){
      bc[i] = TMath::ACos(zi[i]/crad);
    }
    else {
      bc[i] = -1 * TMath::ACos(zi[i]/crad);
    }
  }
  Double_t xc[28], yc[28], zc[28];
    
  //compensation based on node position within individual detector geometries
  //determine compensated radius
  Double_t rcomp = crad + pstartC[2] / 2.0; //
  for(Int_t i = 0; i < 28; i++) {
    //Get compensated translation data
    xc[i] = rcomp*TMath::Cos(ac[i] + TMath::Pi()/2)*TMath::Sin(-1*bc[i]);
    yc[i] = rcomp*TMath::Sin(ac[i] + TMath::Pi()/2)*TMath::Sin(-1*bc[i]);
    cout<<i<<" "<<xc[i]<<" "<<yc[i]<<endl;
  }

  
 Float_t xa[24] = {-11.8, -5.9, 0, 5.9, 11.8, 
		    -11.8, -5.9, 0, 5.9, 11.8,
		    -12.8, -6.9, 6.9, 12.8, 
		    -11.8, -5.9, 0, 5.9, 11.8,
		    -11.8, -5.9, 0, 5.9, 11.8 };
  
  Float_t ya[24] = { 11.9, 11.9, 12.9, 11.9, 11.9,
		     6.0,   6.0,  7.0, 6.0,  6.0,
		    -0., -0., 0., 0.,
		    -6.0, -6.0, -7.0, -6.0, -6.0,
		     -11.9, -11.9, -12.9,  -11.9, -11.9 };  
 
  Float_t xpmt[4] = {-1.32, -1.32, 1.32, 1.32};
  Float_t ypmt[4] = {-1.32, 1.32, 1.32, -1.32};
  Float_t xma[208], yma[208];
  Int_t na=0;
  TString num[208];
  Int_t time[208];
  TCanvas *c1 = new TCanvas("c1");
  c1->Range(-20,-20,20,20);
  TFile * f1 = new TFile("FigESDfit.ruben.root");
  f1->ls();
  TH1F *hTime[208];
  for (int ipmt=0; ipmt<208; ipmt++){
    hTime[ipmt]=(TH1F*)f1->Get(Form("hTime%i",ipmt));
    time[ipmt]=hTime[ipmt]->GetEntries();
    num[ipmt] = Form("%i",time[ipmt]);
  }
  // TH2F *mcpA = (TH2F*)f1->Get("hXYa");
  // TH2F *mcpC = (TH2F*)f1->Get("hXYc");
  c1->cd();
  // mcpA->Draw("COL");
  TCanvas *c2 = new TCanvas("c2");
  c2->Range(-20,-20,20,20);
  c2->cd();
  // mcpC->Draw("COL");
  
  TText t(0,0,"a");
  t.SetTextFont(62);
  t.SetTextSize(0.025);
  t.SetTextAlign(12);
  c1->cd();
 
 for (int imcp=0; imcp<24; imcp++) {
      for (int ix=0; ix<4; ix++) {
	xma[na]=xa[imcp] + xpmt[ix];
	yma[na]=ya[imcp] + ypmt[ix];
	//	cout<<" A "<<na<<" "<<xma[na]<<", "<<yma[na]<<endl;;
	//	num[na] = Form("%i",na);  
	t.DrawText(xma[na],yma[na],num[na].Data());
	na++;
      }
    } // A side
 
   c2->cd();
  for (int imcp=24; imcp<52; imcp++) {
    for (int ix=0; ix<4; ix++) {
      xma[na]=xc[imcp-24] + xpmt[ix];
      yma[na]=yc[imcp-24] + ypmt[ix];
      //   cout<<"C "<<na<<" "<<xma[na]<<", "<<yma[na]<<endl;
      //  num[na] = Form("%i",na);  
      t.DrawText(xma[na],yma[na],num[na]);
     na++;
    }
  } //C side
   
 
  
}
