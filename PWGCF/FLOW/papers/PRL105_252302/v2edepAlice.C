// S. Voloshin   03-oct-2010
// modified from plotV2edepBevalacRhicSqrtS.C - from Art   
{
  float m_p=0.940;
  
  //  gROOT->Reset();
 
  //  gROOT->SetStyle("Bold");  
  //gROOT->SetStyle("Plain");  
  //  gStyle->SetTitleTextColor(kBlack);
  //  gStyle->SetLabelColor(kBlack,"X");
  //  gStyle->SetLabelColor(kBlack,"Y");

  int cCERES=kOrange+7;
  int cPHENIX=kYellow+3;
  int cFOPI=kGreen+2;
  int cE895=kMagenta+2;
  int cAlice=2;
  
  TCanvas *canvas = new TCanvas("v2edep","v2edep",10,10,750,600);
  canvas->cd();

  hist = new TH1F("v2 vs sqrt(s)","v2 vs sqrt(s)", 1, 1., 10000.);
  hist->SetLineColor(0);
  
  TAxis *axis = hist->GetXaxis();
  axis->SetTitle("#sqrt{s_{NN}} (GeV)");
  //axis->CenterTitle(kTRUE);
  axis->SetTitleSize(0.05);
  axis->SetLabelSize(0.04);
  axis->SetTitleOffset(1.5);
  axis = hist->GetYaxis();
  axis->SetTitle("v_{2}(y^{*}=0)");
  axis->SetTitle("v_{2}");
  axis->SetTitleOffset(1.6);
  //axis->CenterTitle(kTRUE);
  axis->SetTitleSize(0.05);
  axis->SetLabelSize(0.04);
 
  hist->SetStats(0); 
  //  hist->SetTitle("Elliptic Flow");
  hist->SetTitle("");
  hist->SetMaximum(0.08); 
  hist->SetMinimum(-0.085); 
  //hist->SetLabelOffset(-0.01,"X");
  gPad->SetTopMargin(.1);
  gPad->SetBottomMargin(.2);
  gPad->SetLeftMargin(.2);
  gPad->SetRightMargin(.1);
  hist->Draw(); 
  canvas->SetLogx(); 
  hist->Draw(); 

  TLine *tline=new TLine(1.,0.,10000.,0.);
  tline->SetLineWidth(1);
  tline->SetLineStyle(7); //wide dash
  tline->Draw("same");
  
  
  
  // CERES ----------------------------------------------------------------------------
  float xxr[3]={8.7,12.3,17.}; //sqrt(s)
  //float xxr[3]={9.57,18.7}; 
  //float yyr[3]={0.028,0.031}; 
  //float yyr[3]={0.028,0.031}; // taken from QM01 and INPC01
  //float yyr[3]={0.028,0.03,0.031}; // QM01 and INPC01, 80GeV interpolated
  //float yyr[3]={0.025,0.035,0.04}; // values integrated over phase space
  float yyr[3]={0.026,0.035,0.04}; // changed to those values after request by H. Appleshaeuser
  //CERES data from yugoslavian conference
  //float xxr[3]={9.57,13.5,18.7}; 
  //float yyr[3]={0.028,0.033,0.036}; 
   
  float exr[3]={0.,0.,0.}; 
  float exr[3]={0.,0.,0.}; 
  float eyr[3]={0.,0.,0.}; 
  //float eyr[3]={0.0025,0.0025,0.0025}; 
  float eyr[3]={0.005,0.005,0.005};  // changed to those values after request by H. Appleshaeuser

  //float eyr[3]={0.004,0.005,0.004}; 

  ceres = new TGraphAsymmErrors(3,xxr,yyr,exr,exr,eyr,eyr); 
  ceres->SetTitle(""); 
  ceres->SetMarkerColor(cCERES); 
  ceres->SetMarkerStyle(kOpenCircle); 
  ceres->SetMarkerSize(1.2);
  ceres->SetLineWidth(2); 
  ceres->SetLineColor(cCERES); 
  ceres->Draw("p"); 
  
  // NA49 --------------------------------------------------------------------------------
  float xxd[4]={8.7,17.,8.7,17.}; //sqrt(s)
  //  float yyd[3]={0.029867,0.0324}; //mod method with Full Acceptance (40)
  //float yyd[4]={0.023676,0.0324,0.023676,0.0324}; //mod method with Cos Cut (40)
  float yyd[4]={0.023676,0.0309,0.023676,0.0309}; //mod method with Cos Cut (40) '96 and '00 data (158)
  //float yyd[4]={0.0452,0.0379,0.0452,0.0379}; //mod method with Cos Cut (40) centrality 5;
  //float yyd[4]={0.0411,0.0406,0.0411,0.0406}; //mod method with Cos Cut (40) centrality 6;
  //float yyd[4]={0.0265,0.0341,0.0265,0.0341}; //mod method with Cos Cut (40) centrality 4;
  //float yyd[4]={0.0206,0.0276,0.0206,0.0276}; //mod method with Cos Cut (40) centrality 3;
  //float yyd[4]={-0.00495,0.0189,-0.00495,0.0189}; //mod method with Cos Cut (40) centrality 2;
  //hist->SetTitle("std method: centrality 2");
  float exd[4]={0.,0.,0.}; 
  //  float eyd[3]={0.001734,0.000716}; //mod method with Full Acceptance (40)
  //float eyd[4]={0.00276,0.000716,0.005,0.005}; //mod method with Cos Cut (40)
  float eyd[4]={0.00276,0.00054,0.005,0.005}; //mod method with Cos Cut (40) '96 and '00 data (158)
  data = new TGraphAsymmErrors(2,xxd,yyd,exd,exd,eyd,eyd); 
  data->SetTitle(""); 
  data->SetMarkerColor(kRed); 
  data->SetMarkerStyle(kFullSquare); 
  data->SetMarkerSize(1.2); 
  data->SetLineWidth(2); 
  data->SetLineColor(kRed); 
  //data->Draw("p"); 
  
  // NA49 
  float xxd[4]={8.7,17.}; //sqrt(s)
  float yyd[4]={0.023676,0.0324}; //mod method with Cos Cut (40)
  float exd[4]={0.,0.,0.}; 
  float eyd[4]={0.01,0.005}; //mod method with Cos Cut (40)
  dataSys = new TGraphAsymmErrors(2,xxd,yyd,exd,exd,eyd,eyd); 
  dataSys->SetTitle(""); 
  //dataSys->SetMarkerColor(kRed); 
  //dataSys->SetMarkerStyle(kFullSquare); 
  dataSys->SetMarkerSize(1.2); 
  dataSys->SetLineWidth(2); 
  dataSys->SetLineColor(kRed); 
  //dataSys->Draw("[]"); 
  
  
  // NA49 cumulant
  float xxd[2]={7.91,15.45};  //sqrt(s)
  float yyd[2]={0.02184,0.0291746}; //order 2
  //float yyd[2]={0.034883,0.027746}; //order 4
  float exd[3]={0.,0.,0.}; 
  //float eyd[3]={0.0105334,0.0026345}; //order 4
  float eyd[3]={0.001173,0.000526268}; //order 2
  cumul = new TGraphAsymmErrors(2,xxd,yyd,exd,exd,eyd,eyd); 
  cumul->SetTitle(""); 
  cumul->SetMarkerColor(kRed); 
  cumul->SetMarkerStyle(kFullSquare); 
  cumul->SetMarkerSize(1.2); 
  cumul->SetLineWidth(2); 
  cumul->SetLineColor(kBlue); 
  cumul->Draw("p"); 
  
  // STAR ------------------------------------------------------------------------------
  float xxd[2]={130.,200.}; //sqrt(s)
  float yyd[2]={0.0426,0.0478}; 
  float exd[2]={0.}; 
  float exd[2]={0.}; 
  float eyd[2]={0.}; 
  float eyd[2]={0.0026,0.0026}; 
  star = new TGraphAsymmErrors(2,xxd,yyd,exd,exd,eyd,eyd); 
  star->SetTitle(""); 
  star->SetMarkerColor(kRed); 
  star->SetMarkerStyle(30); 
  star->SetMarkerSize(1.6); 
  star->SetLineWidth(2); 
  star->SetLineColor(kRed); 
  star->Draw("p"); 
  
  
  // Phenix --------------------------------------------------------------------------
  float xxd[1]={220}; //sqrt(s) 
  float yyd[1]={0.054}; 
  float exd[1]={0.}; 
  float exd[1]={0.}; 
  float eyd[1]={0.}; 
  float eyd[1]={0.0041}; 
  phenix = new TGraphAsymmErrors(1,xxd,yyd,exd,exd,eyd,eyd); 
  phenix->SetTitle(""); 
  phenix->SetMarkerColor(cPHENIX); 
  phenix->SetMarkerStyle(kOpenSquare); 
  phenix->SetMarkerSize(1.2); 
  phenix->SetLineWidth(2); 
  phenix->SetLineColor(cPHENIX); 
  phenix->Draw("p"); 
  

  // Phobos -----------------------------------------------------------------------
  float xxd[2]={117.,180.}; //sqrt(s)
  float yyd[2]={0.048,0.051}; 
  float exd[2]={0.}; 
  float exd[2]={0.}; 
  float eyd[2]={0.}; 
  float eyd[2]={0.005,0.005}; 
  phobos = new TGraphAsymmErrors(2,xxd,yyd,exd,exd,eyd,eyd); 
  phobos->SetTitle(""); 
  phobos->SetMarkerColor(kBlue); 
  phobos->SetMarkerStyle(28); 
  phobos->SetMarkerSize(1.2); 
  phobos->SetLineWidth(2); 
  phobos->SetLineColor(kBlue); 
  phobos->Draw("p"); 


  //E877 -------------------------------------------------------------------------
  //float xxhh[4]={4.75,17.,130.,200.}; 
  //float yyhh[4]={8.2,9.3,10.,9.4}; 
  float xxhh[1]={4.75}; // sqrt(s)
  float yyhh[1]={0.019}; 
  float exhh[1]={0.}; 
  float exhh[1]={0.}; 
  float eyhh[1]={0.}; 
  float eyhh[1]={0.002}; 
  e877 = new TGraphAsymmErrors(1,xxhh,yyhh,exhh,exhh,eyhh,eyhh); 
  e877->SetTitle(""); 
  e877->SetMarkerColor(kBlue); 
  e877->SetMarkerStyle(34); 
  e877->SetMarkerSize(1.2); 
  e877->SetLineWidth(2); 
  e877->SetLineColor(kBlue); 
  e877->Draw("p"); 
  
  //E895 ------------------------------------------------------------------------
  //float xxhq[4]={4.75,17.,130.,200.}; 
  //float yyhq[4]={9.1.,10.,9.8,8.4}; 
  float xxhq[4]={2.68,3.32,3.83,4.24}; // sqrt(s)
  float yyhq[4]={-0.05,-0.005,0.01,0.015}; 
  float exhq[4]={0.,0.,0.}; 
  float eyhq[4]={0.004,0.003,0.003,0.004}; 
  e895 = new TGraphAsymmErrors(4,xxhq,yyhq,exhq,exhq,eyhq,eyhq); 
  e895->SetTitle(""); 
  e895->SetMarkerColor(cE895); 
  e895->SetMarkerStyle(kFullTriangleUp); 
  e895->SetMarkerSize(1.2); 
  e895->SetLineWidth(2); 
  e895->SetLineColor(cE895); 
  e895->Draw("p"); 
  //Plastic Ball ----------------------------------------------------
  //float xxhq[4]={4.75,17.,130.,200.}; 
  //float yyhq[4]={9.1.,10.,9.8,8.4}; 
  float xxPB[6]={.15,.2,.4,.6,.8,1.2}; 
  float yyPB[6]={-.035,-.08,-.09,-.06,-.04,-.03}; 
  float exPB[6]={0.,0.,0.}; 
  float eyPB[6]={0.,0.,0.}; 
  plasticBall = new TGraphAsymmErrors(6,xxPB,yyPB,exPB,exPB,eyPB,eyPB); 
  plasticBall->SetTitle(""); 
  plasticBall->SetMarkerColor(kBlue); 
  plasticBall->SetMarkerStyle(kFullTriangleDown); 
  plasticBall->SetMarkerSize(1.2); 
  plasticBall->SetLineWidth(1); 
  plasticBall->SetLineColor(kBlue); 
  //plasticBall->Draw("p");


  //Fopi------------------------------------------------------------
  float xxFopi[10]={0.09,0.12,0.15,0.25,0.4,0.6,0.8,1.0,1.2,1.49}; 
  float yyFopi[10]={0.07456,0.02847,-0.00774,-0.05784,-0.08200,-0.07087,-0.06845,-0.06327,-0.05523,-0.04344}; 
  float exFopi[10]={0.,0.,0.}; 
  float eyFopi[10]={0.00746,0.00285,0.00077,0.00578,0.00656,0.00850,0.00684,0.00823,0.00828,0.00956}; 
  for (int i=0;i<10;i++) xxFopi[i] = sqrt((xxFopi[i]+m_p)*m_p*2. +2.*m_p*m_p);
    
  fopi = new TGraphAsymmErrors(10,xxFopi,yyFopi,exFopi,exFopi,eyFopi,eyFopi); 
  fopi->SetTitle(""); 
  fopi->SetMarkerColor(cFOPI); 
  fopi->SetMarkerStyle(kFullTriangleDown); 
  fopi->SetMarkerSize(1.2); 
  fopi->SetLineWidth(2); 
  fopi->SetLineColor(cFOPI);
  fopi->Draw("p"); 


  //EOS --------------------------------------------------------------------
  //float xxhq[4]={4.75,17.,130.,200.}; 
  //float yyhq[4]={9.1.,10.,9.8,8.4}; 
  float xxhq[4]={2.35}; 
  float yyhq[4]={-0.065}; 
  float exhq[4]={0.}; 
  float eyhq[4]={0.007}; 
  eos = new TGraphAsymmErrors(1,xxhq,yyhq,exhq,exhq,eyhq,eyhq); 
  eos->SetTitle(""); 
  eos->SetMarkerColor(kRed); 
  eos->SetMarkerStyle(3); 
  eos->SetMarkerSize(1.2); 
  eos->SetLineWidth(2); 
  eos->SetLineColor(kRed); 
  eos->Draw("p"); 



// ALICE measured v2{4} + statistical error
  float v24ALICE = 0.073;
  float v24eALICE = 0.001;
  float corrpt = 0.88;
  float xxhh[1]={2760.}; // sqrt(s)
  float yyhh[1]={v24ALICE*corrpt}; 
  float exhh[1]={0.};  


// ALICE systematic error
// hijing 15.3 %
// therminator 6.1%
  float v24eALICE2 = 0.004;
  float eyhh2[1]={v24eALICE2};

  //draw with systematic error
  //alice = new TGraphAsymmErrors(1,xxhh,yyhh,exhh,exhh,eyhh,eyhh); 
  alice = new TGraphAsymmErrors(1,xxhh,yyhh,exhh,exhh,eyhh2,eyhh2); 
  alice->SetTitle(""); 
  alice->SetMarkerColor(cAlice); 
  alice->SetMarkerStyle(20); 
  alice->SetMarkerSize(1.2); 
  alice->SetLineWidth(2); 
  alice->SetLineColor(cAlice); 
  alice->Draw("p"); 
  
  //------------------------------------------------------------------  
  TLegend *legend = new TLegend(.68,.24,.86,.66);
  // TLegend *legend = new TLegend(.45,.23,.62,.52);
  legend->Clear();
  legend->SetBorderSize(1);
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.028);
  //legend->AddEntry(e895,"E895","P");
//   legend->AddEntry(e877,"E877","P");
//   legend->AddEntry(data,"NA49 pions","P");
//   legend->AddEntry(ceres,"CERES charged particles","P");
//   legend->AddEntry(star,"STAR charged particles","P");
  //legend->SetHeader("A+A:");
  //legend->AddEntry(data,"NA49 std/mod","P");
  legend->AddEntry(alice,"ALICE","p");
  legend->AddEntry(star,"STAR","P");
  legend->AddEntry(phobos,"PHOBOS","P");
  legend->AddEntry(phenix,"PHENIX","P");
  legend->AddEntry(cumul,"NA49","P");
  legend->AddEntry(ceres,"CERES","P");
  legend->AddEntry(e877,"E877","P");
  legend->AddEntry(eos,"EOS","P");
  legend->AddEntry(e895,"E895","P");
  legend->AddEntry(fopi,"FOPI","P");
  legend->SetTextFont(22); // 22 = Times New Roman (bold)

  legend->Draw();
  //--------------------------- insert logo -----------
  TImage *ps = TImage::Open("macros/ALICElogo.png");
  float xlogo=0.75;
  float ylogo=0.76;
  float wlogo=0.08;
  TPad *aliceLogo = new TPad("aliceLogo", "aliceLogo",xlogo,ylogo,xlogo+wlogo,ylogo+wlogo);
  //  aliceLogo->Draw();
  aliceLogo->cd();
  //  ps->Draw("same"); 
  //-------------------------- print -----------
  canvas->Print("v2edep.png");
  canvas->Print("v2edep.pdf");
  canvas->Print("v2edep.eps");
  
  
} 
 
