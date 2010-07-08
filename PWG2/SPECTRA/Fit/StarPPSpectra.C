TGraphErrors ** StarPPSpectra() {

  TGraphErrors ** gStar = new TGraphErrors*[6];
  const Int_t kNMaxPoints = 19;
  Double_t pt[kNMaxPoints];
  Double_t errpt[kNMaxPoints] = {0};
  Double_t pim [kNMaxPoints];
  Double_t pime[kNMaxPoints];
  Double_t pip [kNMaxPoints];
  Double_t pipe[kNMaxPoints];
  Double_t km [kNMaxPoints];
  Double_t kme[kNMaxPoints];
  Double_t kp [kNMaxPoints];
  Double_t kpe[kNMaxPoints];
  Double_t pm [kNMaxPoints];
  Double_t pme[kNMaxPoints];
  Double_t pp [kNMaxPoints];
  Double_t ppe[kNMaxPoints];

  Int_t ipoint = 0;
  pt[ipoint++] = 0.225;
  pt[ipoint++] = 0.275;
  pt[ipoint++] = 0.325;
  pt[ipoint++] = 0.375;
  pt[ipoint++] = 0.425;
  pt[ipoint++] = 0.475;
  pt[ipoint++] = 0.525;
  pt[ipoint++] = 0.575;
  pt[ipoint++] = 0.625;
  pt[ipoint++] = 0.675;
  pt[ipoint++] = 0.725;
  pt[ipoint++] = 0.775;
  pt[ipoint++] = 0.825;
  pt[ipoint++] = 0.875;
  pt[ipoint++] = 0.925;
  pt[ipoint++] = 0.975;
  pt[ipoint++] = 1.025;
  pt[ipoint++] = 1.075;
  pt[ipoint++] = 1.125;


  // pi-
  ipoint = 0;
  pim[ipoint++]  = 2.02;
  pim[ipoint++]  = 1.52;
  pim[ipoint++]  = 1.13;
  pim[ipoint++]  = 8.44*1e-1;
  pim[ipoint++]  = 6.35*1e-1;
  pim[ipoint++]  = 4.69*1e-1;
  pim[ipoint++]  = 3.54*1e-1;
  pim[ipoint++]  = 2.67*1e-1;
  pim[ipoint++]  = 2.02*1e-1;
  pim[ipoint++]  = 1.53*1e-1;
  pim[ipoint++]  = 1.16*1e-1;
  pim[ipoint++]  = 8.65*1e-2;

  ipoint = 0;
  pime[ipoint++] = 0.06;
  pime[ipoint++] = 0.03;
  pime[ipoint++] = 0.02;
  pime[ipoint++] = 0.09*1e-1;
  pime[ipoint++] = 0.07*1e-1;
  pime[ipoint++] = 0.05*1e-1;
  pime[ipoint++] = 0.04*1e-1;
  pime[ipoint++] = 0.03*1e-1;
  pime[ipoint++] = 0.04*1e-1;
  pime[ipoint++] = 0.03*1e-1;
  pime[ipoint++] = 0.03*1e-1;
  pime[ipoint++] = 0.30*1e-2;

  gStar[0] = new TGraphErrors(ipoint,pt,pim,errpt,pime);

  // pi+
  ipoint = 0;
  pip[ipoint++]  = 2.07;
  pip[ipoint++]  = 1.54;
  pip[ipoint++]  = 1.14;
  pip[ipoint++]  = 8.57*1e-1;
  pip[ipoint++]  = 6.38*1e-1;
  pip[ipoint++]  = 4.76*1e-1;
  pip[ipoint++]  = 3.59*1e-1;
  pip[ipoint++]  = 2.73*1e-1;
  pip[ipoint++]  = 2.07*1e-1;
  pip[ipoint++]  = 1.55*1e-1;
  pip[ipoint++]  = 1.16*1e-1;
  pip[ipoint++]  = 8.95*1e-2;

  ipoint = 0;
  pipe[ipoint++] = 0.06;
  pipe[ipoint++] = 0.03;
  pipe[ipoint++] = 0.02;
  pipe[ipoint++] = 0.09*1e-1;
  pipe[ipoint++] = 0.07*1e-1;
  pipe[ipoint++] = 0.05*1e-1;
  pipe[ipoint++] = 0.04*1e-1;
  pipe[ipoint++] = 0.03*1e-1;
  pipe[ipoint++] = 0.04*1e-1;
  pipe[ipoint++] = 0.03*1e-1;
  pipe[ipoint++] = 0.03*1e-1;
  pipe[ipoint++] = 0.31*1e-2;

  gStar[1] = new TGraphErrors(ipoint,pt,pip,errpt,pipe);

  // k-
  ipoint = 0;
  km[ipoint++]  = 1.43*1e-1;
  km[ipoint++]  = 1.26*1e-1;
  km[ipoint++]  = 1.08*1e-1;
  km[ipoint++]  = 8.77*1e-2;
  km[ipoint++]  = 7.34*1e-2;
  km[ipoint++]  = 6.17*1e-2;
  km[ipoint++]  = 4.87*1e-2;
  km[ipoint++]  = 4.11*1e-2;
  km[ipoint++]  = 3.70*1e-2;
  km[ipoint++]  = 2.86*1e-2;
  km[ipoint++]  = 2.42*1e-2;

  ipoint = 0;
  kme[ipoint++] = 0.11*1e-1;
  kme[ipoint++] = 0.05*1e-1;
  kme[ipoint++] = 0.02*1e-1;
  kme[ipoint++] = 0.24*1e-2;
  kme[ipoint++] = 0.26*1e-2;
  kme[ipoint++] = 0.63*1e-2;
  kme[ipoint++] = 0.50*1e-2;
  kme[ipoint++] = 0.42*1e-2;
  kme[ipoint++] = 0.28*1e-2;
  kme[ipoint++] = 0.22*1e-2;
  kme[ipoint++] = 0.27*1e-2;

  gStar[2] = new TGraphErrors(ipoint,pt,km,errpt,kme);


  // k+
  ipoint = 0;
  kp[ipoint++]  = 1.52*1e-1;
  kp[ipoint++]  = 1.30*1e-1;
  kp[ipoint++]  = 1.08*1e-1;
  kp[ipoint++]  = 9.16*1e-2;
  kp[ipoint++]  = 7.47*1e-2;
  kp[ipoint++]  = 6.26*1e-2;
  kp[ipoint++]  = 5.26*1e-2;
  kp[ipoint++]  = 4.41*1e-2;
  kp[ipoint++]  = 3.81*1e-2;
  kp[ipoint++]  = 3.06*1e-2;
  kp[ipoint++]  = 2.49*1e-2;

  ipoint = 0;
  kpe[ipoint++] = 0.11*1e-1;
  kpe[ipoint++] = 0.05*1e-1;
  kpe[ipoint++] = 0.02*1e-1;
  kpe[ipoint++] = 0.25*1e-2;
  kpe[ipoint++] = 0.27*1e-2;
  kpe[ipoint++] = 0.64*1e-2;
  kpe[ipoint++] = 0.54*1e-2;
  kpe[ipoint++] = 0.45*1e-2;
  kpe[ipoint++] = 0.28*1e-2;
  kpe[ipoint++] = 0.24*1e-2;
  kpe[ipoint++] = 0.28*1e-2;

  gStar[3] = new TGraphErrors(ipoint,pt,kp,errpt,kpe);


  // pbar
  ipoint = 0;
  pm[ipoint++]  = 0;
  pm[ipoint++]  = 0;
  pm[ipoint++]  = 0;
  pm[ipoint++]  = 5.54*1e-2;
  pm[ipoint++]  = 4.81*1e-2;
  pm[ipoint++]  = 4.25*1e-2;
  pm[ipoint++]  = 3.77*1e-2;
  pm[ipoint++]  = 3.38*1e-2;
  pm[ipoint++]  = 2.78*1e-2;
  pm[ipoint++]  = 2.49*1e-2;
  pm[ipoint++]  = 2.03*1e-2;
  pm[ipoint++]  = 1.75*1e-2;
  pm[ipoint++]  = 1.52*1e-2;
  pm[ipoint++]  = 1.27*1e-2;
  pm[ipoint++]  = 1.05*1e-2;
  pm[ipoint++]  = 8.95*1e-3;
  pm[ipoint++]  = 7.36*1e-3;
  pm[ipoint++]  = 6.59*1e-3;
  pm[ipoint++]  = 5.21*1e-3;

  ipoint = 0;
  pme[ipoint++] = 0;
  pme[ipoint++] = 0;
  pme[ipoint++] = 0;
  pme[ipoint++] = 0.13*1e-2;
  pme[ipoint++] = 0.11*1e-2;
  pme[ipoint++] = 0.10*1e-2;
  pme[ipoint++] = 0.09*1e-2;
  pme[ipoint++] = 0.08*1e-2;
  pme[ipoint++] = 0.07*1e-2;
  pme[ipoint++] = 0.06*1e-2;
  pme[ipoint++] = 0.06*1e-2;
  pme[ipoint++] = 0.06*1e-2;
  pme[ipoint++] = 0.05*1e-2;
  pme[ipoint++] = 0.04*1e-2;
  pme[ipoint++] = 0.04*1e-2;
  pme[ipoint++] = 0.34*1e-3;
  pme[ipoint++] = 0.34*1e-3;
  pme[ipoint++] = 0.32*1e-3;
  pme[ipoint++] = 0.29*1e-3;

  gStar[4] =  new TGraphErrors(ipoint,pt,pm,errpt,pme);


  // p
  ipoint = 0;
  pp[ipoint++]  = 0;
  pp[ipoint++]  = 0;
  pp[ipoint++]  = 0;
  pp[ipoint++]  = 0;
  pp[ipoint++]  = 0;
  pp[ipoint++]  = 5.07*1e-2;
  pp[ipoint++]  = 4.69*1e-2;
  pp[ipoint++]  = 4.08*1e-2;
  pp[ipoint++]  = 3.42*1e-2;
  pp[ipoint++]  = 2.87*1e-2;
  pp[ipoint++]  = 2.41*1e-2;
  pp[ipoint++]  = 2.13*1e-2;
  pp[ipoint++]  = 1.82*1e-2;
  pp[ipoint++]  = 1.54*1e-2;
  pp[ipoint++]  = 1.31*1e-2;
  pp[ipoint++]  = 1.11*1e-2;
  pp[ipoint++]  = 9.78*1e-3;
  pp[ipoint++]  = 8.56*1e-3;
  pp[ipoint++]  = 7.38*1e-3;

  ipoint = 0;
  ppe[ipoint++] = 0;
  ppe[ipoint++] = 0;
  ppe[ipoint++] = 0;
  ppe[ipoint++] = 0;
  ppe[ipoint++] = 0;
  ppe[ipoint++] = 0.25*1e-2;
  ppe[ipoint++] = 0.19*1e-2;
  ppe[ipoint++] = 0.14*1e-2;
  ppe[ipoint++] = 0.10*1e-2;
  ppe[ipoint++] = 0.07*1e-2;
  ppe[ipoint++] = 0.07*1e-2;
  ppe[ipoint++] = 0.06*1e-2;
  ppe[ipoint++] = 0.05*1e-2;
  ppe[ipoint++] = 0.05*1e-2;
  ppe[ipoint++] = 0.04*1e-2;
  ppe[ipoint++] = 0.04*1e-2;
  ppe[ipoint++] = 0.40*1e-3;
  ppe[ipoint++] = 0.37*1e-3;
  ppe[ipoint++] = 0.38*1e-3;



  gStar[5] = new TGraphErrors(ipoint,pt,pp,errpt,ppe);

  gStar[0]->SetMarkerStyle(kOpenTriangleUp);
  gStar[1]->SetMarkerStyle(kFullTriangleUp);
  gStar[2]->SetMarkerStyle(kOpenCircle);
  gStar[3]->SetMarkerStyle(kFullCircle);
  gStar[4]->SetMarkerStyle(kOpenSquare);
  gStar[5]->SetMarkerStyle(kFullSquare);

  gStar[0]->SetMarkerColor(kBlack);
  gStar[1]->SetMarkerColor(kBlack);
  gStar[2]->SetMarkerColor(kRed);
  gStar[3]->SetMarkerColor(kRed);
  gStar[4]->SetMarkerColor(kBlue);
  gStar[5]->SetMarkerColor(kBlue);

  gStar[0]->SetLineColor (kBlack);
  gStar[1]->SetLineColor (kBlack);
  gStar[2]->SetLineColor (kRed);
  gStar[3]->SetLineColor (kRed);
  gStar[4]->SetLineColor (kBlue);
  gStar[5]->SetLineColor (kBlue);

  // Scale all for 2*pi and pt in order to go to dNdpt
  for(Int_t istar = 0; istar < 6; istar++){
    Int_t npoint = gStar[istar]->GetN();
    for(Int_t ipoint = 0; ipoint < npoint; ipoint++){
      Double_t ptbin = gStar[istar]->GetX()[ipoint];
      gStar[istar]->SetPoint(ipoint,
			     gStar[istar]->GetX()[ipoint],
			     gStar[istar]->GetY()[ipoint]*2*TMath::Pi()*ptbin);
      gStar[istar]->SetPointError(ipoint,
				  gStar[istar]->GetEX()[ipoint],
				  gStar[istar]->GetEY()[ipoint]*2*TMath::Pi()*ptbin);
    }
    
  }
  


//   gStar[0]->Draw("AP");
//   gStar[1]->Draw("P");
//   gStar[2]->Draw("P");
//   gStar[3]->Draw("P");
//   gStar[4]->Draw("P");
//   gStar[5]->Draw("P");


  return gStar;


  
  // table from PRC79  79 (2009) 34909
//         pi-                   pi+               k-                 K+               pbar              p
//  0.225  2.02+-0.06            2.07+-0.06        (1.43+-0.11)*10-1 (1.52+-0.11)*10-1
//  0.275  1.52+-0.03            1.54+-0.03        (1.26+-0.05)*10-1 (1.30+-0.05)*10-1
//  0.325  1.13+-0.02            1.14+-0.02        (1.08+-0.02)*10-1 (1.08+-0.02)*10-1
//  0.375  (8.44+-0.09)*10-1     (8.57+-0.09)*10-1 (8.77+-0.24)*10-2 (9.16+-0.25)*10-2 (5.54+-0.13)*10-2
//  0.425  (6.35+-0.07)*10-1     (6.38+-0.07)*10-1 (7.34+-0.26)*10-2 (7.47+-0.27)*10-2 (4.81+-0.11)*10-2
//  0.475  (4.69+-0.05)*10-1     (4.76+-0.05)*10-1 (6.17+-0.63)*10-2 (6.26+-0.64)*10-2 (4.25+-0.10)*10-2 (5.07+-0.25)*10-2
//  0.525  (3.54+-0.04)*10-1     (3.59+-0.04)*10-1 (4.87+-0.50)*10-2 (5.26+-0.54)*10-2 (3.77+-0.09)*10-2 (4.69+-0.19)*10-2
//  0.575  (2.67+-0.03)*10-1     (2.73+-0.03)*10-1 (4.11+-0.42)*10-2 (4.41+-0.45)*10-2 (3.38+-0.08)*10-2 (4.08+-0.14)*10-2
//  0.625  (2.02+-0.04)*10-1     (2.07+-0.04)*10-1 (3.70+-0.28)*10-2 (3.81+-0.28)*10-2 (2.78+-0.07)*10-2 (3.42+-0.10)*10-2
//  0.675  (1.53+-0.03)*10-1     (1.55+-0.03)*10-1 (2.86+-0.22)*10-2 (3.06+-0.24)*10-2 (2.49+-0.06)*10-2 (2.87+-0.07)*10-2
//  0.725  (1.16+-0.03)*10-1     (1.16+-0.03)*10-1 (2.42+-0.27)*10-2 (2.49+-0.28)*10-2 (2.03+-0.06)*10-2 (2.41+-0.07)*10-2
//  0.775  (8.65+-0.30)*10-2     (8.95+-0.31)*10-2                                     (1.75+-0.06)*10-2 (2.13+-0.06)*10-2
//  0.825                                                                              (1.52+-0.05)*10-2 (1.82+-0.05)*10-2
//  0.875                                                                              (1.27+-0.04)*10-2 (1.54+-0.05)*10-2
//  0.925                                                                              (1.05+-0.04)*10-2 (1.31+-0.04)*10-2
//  0.975                                                                              (8.95+-0.34)*10-3 (1.11+-0.04)*10-2
//  1.025                                                                              (7.36+-0.34)*10-3 (9.78+-0.40)*10-3
//  1.075                                                                              (6.59+-0.32)*10-3 (8.56+-0.37)*10-3
//  1.125                                                                              (5.21+-0.29)*10-3 (7.38+-0.38)*10-3


  }
