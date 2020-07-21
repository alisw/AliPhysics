TFile* data1;

TFile* syst_std_001[4]; //for the syst uncertainties
TFile* syst_BDflat_001[4]; //to propagate the FD uncertainty under BD vs cent flat assumption
TFile* syst_BDextr_001[4]; //to propagate the FD uncertainty under BD vs cent extreme assumption
TFile* syst_std_0110[4]; 
TFile* syst_BDflat_0110[4];
TFile* syst_BDextr_0110[4]; 
TFile* syst_std_1030[4]; 
TFile* syst_BDflat_1030[4]; 
TFile* syst_BDextr_1030[4];
TFile* syst_std_30100[4];
TFile* syst_BDflat_30100[4]; //to propagate the FD uncertainty under BD vs cent flat assumption
TFile* syst_BDextr_30100[4]; //to propagate the FD uncertainty under BD vs cent extreme assumption

TFile* totalsyst_NSy_001[4]; //to get baseline var fit unc and v2 fit unc
TFile* totalsyst_NSw_001[4]; //to get baseline var fit unc and v2 fit unc
TFile* totalsyst_NSy_0110[4]; 
TFile* totalsyst_NSw_0110[4]; 
TFile* totalsyst_NSy_1030[4]; 
TFile* totalsyst_NSw_1030[4]; 
TFile* totalsyst_NSy_30100[4]; 
TFile* totalsyst_NSw_30100[4]; 

TFile* dPhisyst_001[4][4]; //to get yield unc, bkg shape unc
TFile* dPhisyst_0110[4][4]; 
TFile* dPhisyst_1030[4][4]; 
TFile* dPhisyst_30100[4][4]; 

void Ratio_Centralities_NS_Final(){
 
 data1 = TFile::Open("ComparePPbVsCentFitResults_4x2Bins.root");

 LoadFilesForSystPropagation();

 TCanvas *c1 = (TCanvas*)data1->Get("cPPvsPPbFitResultsFinalPaperStyle");
 c1->ls();
 
 TCanvas * cfinal = new TCanvas("cout","coutput",1000,1000); 

 TPad * ptmp = (TPad*)c1->FindObject("cPPvsPPbFitResultsFinalPaperStyle_1");
 ptmp->SetName("pad1");
 TPad * ptmp = (TPad*)c1->FindObject("cPPvsPPbFitResultsFinalPaperStyle_2");
 ptmp->SetName("pad2");
 TPad * ptmp = (TPad*)c1->FindObject("cPPvsPPbFitResultsFinalPaperStyle_3");
 ptmp->SetName("pad3");
 TPad * ptmp = (TPad*)c1->FindObject("cPPvsPPbFitResultsFinalPaperStyle_4");
 ptmp->SetName("pad4");
 TPad * ptmp = (TPad*)c1->FindObject("cPPvsPPbFitResultsFinalPaperStyle_5");
 ptmp->SetName("pad5");
 TPad * ptmp = (TPad*)c1->FindObject("cPPvsPPbFitResultsFinalPaperStyle_6");
 ptmp->SetName("pad6");
 TPad * ptmp = (TPad*)c1->FindObject("cPPvsPPbFitResultsFinalPaperStyle_7");
 ptmp->SetName("pad7");
 TPad * ptmp = (TPad*)c1->FindObject("cPPvsPPbFitResultsFinalPaperStyle_8");
 ptmp->SetName("pad8");

 TPad * p1 = (TPad*)c1->FindObject("pad1");
 TPad * p2 = (TPad*)c1->FindObject("pad2");
 TPad * p3 = (TPad*)c1->FindObject("pad3");
 TPad * p4 = (TPad*)c1->FindObject("pad4");

 TPad * p5 = (TPad*)(c1->FindObject("pad1")->Clone("pad1b"));
 TPad * p6 = (TPad*)(c1->FindObject("pad2")->Clone("pad2b"));
 TPad * p7 = (TPad*)(c1->FindObject("pad3")->Clone("pad3b"));
 TPad * p8 = (TPad*)(c1->FindObject("pad4")->Clone("pad4b"));

 TPad * p12 = (TPad*)c1->FindObject("pad5");
 TPad * p22 = (TPad*)c1->FindObject("pad6");
 TPad * p32 = (TPad*)c1->FindObject("pad7");
 TPad * p42 = (TPad*)c1->FindObject("pad8");

 TPad * p52 = (TPad*)(c1->FindObject("pad5")->Clone("pad5b"));
 TPad * p62 = (TPad*)(c1->FindObject("pad6")->Clone("pad6b"));
 TPad * p72 = (TPad*)(c1->FindObject("pad7")->Clone("pad7b"));
 TPad * p82 = (TPad*)(c1->FindObject("pad8")->Clone("pad8b"));
 
//=========Pad1==========================
 printf("*** PROCESSING PAD 1 ***\n");
  cfinal->cd();
 p1->SetFillStyle(4000);
 p1->SetFrameFillStyle(4000);
 p1->SetPad(0.,0.7408+0.0092,0.33+0.0048,1.0);  // 0.02 Margin = 0.005133 pad dimension
 p1->SetMargin(0.30,0.,0.03,0.02);
 TList *lc=p1->GetListOfPrimitives();

 TLatex *pt1=(TLatex*)lc->At(18);
 pt1->SetY(pt1->GetY()+0.02);
 pt1->SetX(pt1->GetX()+0.08);
 pt1->SetTextSize(22.0);

 TLatex *pt1=(TLatex*)lc->At(19);
 pt1->SetX(pt1->GetX()+0.03);
 pt1->SetTextSize(22.0);

 TLatex *pt1=(TLatex*)lc->At(20);
 pt1->SetX(pt1->GetX()+0.06);
 pt1->SetY(pt1->GetY()-0.02);
 pt1->SetTextSize(20.0);

 TH1D *h=(TH1D*)p1->FindObject("hDraw0");
 h->GetYaxis()->SetTitleSize(25);
 h->GetYaxis()->SetTitleOffset(3);
 h->GetYaxis()->SetLabelSize(20);
 h->GetYaxis()->SetRangeUser(0,3.7);

  ((TGraphAsymmErrors*)lc->At(2))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(2))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(3))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(3))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(4))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(4))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(5))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(5))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(6))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(6))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(7))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(7))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(8))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(8))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(9))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(9))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(10))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(10))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(11))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(11))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(12))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(12))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(13))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(13))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(14))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(14))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(15))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(15))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(16))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(16))->GetMarkerSize()-0.4);          
  ((TGraphAsymmErrors*)lc->At(17))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(17))->GetMarkerSize()-0.4);

 p1->Draw();

//=========Pad2==========================
 printf("*** PROCESSING PAD 2 ***\n");
 cfinal->cd();
 p2->SetFillStyle(4000);
 p2->SetFrameFillStyle(4000);
 p2->SetPad(0.33,0.7408+0.0092,0.55+0.0048,1.0);
 p2->SetMargin(0.02,0.,0.03,0.02);
 TList *lc=p2->GetListOfPrimitives();

 TH1D *h=(TH1D*)p2->FindObject("hDraw1");
 h->GetYaxis()->SetRangeUser(0,3.7);

 TLegend *lg=(TLegend*)lc->At(18);
 lg->SetBorderSize(0.);
 lg->SetTextSize(16);
 lg->SetY1(lg->GetY1()-1.51);
 lg->SetY2(lg->GetY2()-2.11);
 lg->DeleteEntry();
 lg->DeleteEntry();
 lg->DeleteEntry();
 lg->DeleteEntry();
 lg->AddEntry((TGraphAsymmErrors*)lc->At(10),"0#minus0.1\% V0M","lp");
 lg->AddEntry((TGraphAsymmErrors*)lc->At(12),"0.1#minus10\% V0M","lp");
 lg->AddEntry((TGraphAsymmErrors*)lc->At(14),"10#minus30\% V0M","lp");
 lg->AddEntry((TGraphAsymmErrors*)lc->At(16),"30#minus100\% V0M","lp");

 TLatex *pt1=(TLatex*)lc->At(19);
 pt1->SetX(pt1->GetX()+0.03);
 pt1->SetTextSize(22.0);

 TLatex *pt1=(TLatex*)lc->At(20);
 pt1->SetX(pt1->GetX()+0.06);
 pt1->SetY(pt1->GetY()-0.02);
 pt1->SetTextSize(20.0);

  ((TGraphAsymmErrors*)lc->At(2))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(2))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(3))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(3))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(4))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(4))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(5))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(5))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(6))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(6))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(7))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(7))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(8))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(8))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(9))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(9))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(10))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(10))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(11))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(11))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(12))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(12))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(13))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(13))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(14))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(14))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(15))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(15))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(16))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(16))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(17))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(17))->GetMarkerSize()-0.4);

 p2->Draw();

//=========Pad3==========================
 printf("*** PROCESSING PAD 3 ***\n");
 p3->SetFillStyle(4000);
 p3->SetFrameFillStyle(4000);
 p3->SetPad(0.55,0.7408+0.0092,0.76+0.0048,1.0);
 p3->SetMargin(0.02,0.,0.03,0.02);
 TList *lc=p3->GetListOfPrimitives();

 TH1D *h=(TH1D*)p3->FindObject("hDraw2");
 h->GetYaxis()->SetRangeUser(0,3.7);

 TLatex *pt1=(TLatex*)lc->At(18);
 pt1->SetX(pt1->GetX()-0.05);
 pt1->SetTextSize(22.0);

 TLatex *pt1=(TLatex*)lc->At(19);
 pt1->SetX(pt1->GetX()+0.06);
 pt1->SetY(pt1->GetY()-0.02);
 pt1->SetTextSize(20.0);

  ((TGraphAsymmErrors*)lc->At(2))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(2))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(3))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(3))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(4))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(4))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(5))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(5))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(6))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(6))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(7))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(7))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(8))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(8))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(9))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(9))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(10))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(10))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(11))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(11))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(12))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(12))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(13))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(13))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(14))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(14))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(15))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(15))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(16))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(16))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(17))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(17))->GetMarkerSize()-0.4);

 p3->Draw();

//=========Pad4==========================
 printf("*** PROCESSING PAD 4 ***\n");
 cfinal->cd();
 p4->SetFillStyle(4000);
 p4->SetFrameFillStyle(4000);
 p4->SetPad(0.76,0.7408+0.0092,0.99,1.0);
 p4->SetMargin(0.02,0.,0.03,0.02);
 TList *lc=p4->GetListOfPrimitives();

 TH1D *h=(TH1D*)p4->FindObject("hDraw3");
 h->GetYaxis()->SetRangeUser(0,3.7);

 TLatex *pt1=(TLatex*)lc->At(18);
 pt1->SetX(pt1->GetX()-0.05);
 pt1->SetTextSize(22.0);

  ((TGraphAsymmErrors*)lc->At(2))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(2))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(3))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(3))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(4))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(4))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(5))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(5))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(6))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(6))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(7))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(7))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(8))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(8))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(9))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(9))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(10))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(10))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(11))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(11))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(12))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(12))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(13))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(13))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(14))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(14))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(15))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(15))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(16))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(16))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(17))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(17))->GetMarkerSize()-0.4);

 p4->Draw();

//=========Pad5==========================
 printf("*** PROCESSING PAD 5 ***\n");
 cfinal->cd();
 p5->SetFillStyle(4000);
 p5->SetFrameFillStyle(4000);
 p5->SetPad(0.,0.5122+0.0092,0.33+0.0048,0.7408+0.017);  // 0.02 Margin = 0.005133 pad dimension
 p5->SetMargin(0.30,0.,0.03,0.0);
 TList *lc=p5->GetListOfPrimitives();

  TGraphAsymmErrors* grRatio1 = (TGraphAsymmErrors*)lc->At(12);
  TGraphAsymmErrors* grRatiosyst1 = (TGraphAsymmErrors*)lc->At(4);
  TGraphAsymmErrors* grRatio2 = (TGraphAsymmErrors*)lc->At(14);
  TGraphAsymmErrors* grRatiosyst2 = (TGraphAsymmErrors*)lc->At(6);
  TGraphAsymmErrors* grRatio3 = (TGraphAsymmErrors*)lc->At(16);
  TGraphAsymmErrors* grRatiosyst3 = (TGraphAsymmErrors*)lc->At(8);

 //pick legend from next panel
 TList *lcLEG=p6->GetListOfPrimitives();
 TLegend *lg=(TLegend*)lcLEG->At(18);
 lg->SetBorderSize(0.);
 lg->SetTextSize(16);
 lg->SetX1NDC(0.11);
 lg->SetX2NDC(0.3);
 lg->SetY1NDC(0.535);
 lg->SetY2NDC(0.59);
 //lg->SetY1(lg->GetY1()-1.92);
 //lg->SetY2(lg->GetY2()-2.52);
 //lg->SetX1(lg->GetX1()-0.58);
 lg->DeleteEntry();
 lg->DeleteEntry();
 lg->DeleteEntry();
 lg->DeleteEntry();
 lg->AddEntry(grRatio1,"0.1#minus10\ / 0#minus0.1\% V0M","lp");
 lg->AddEntry(grRatio2,"10#minus30\ / 0#minus0.1\% V0M","lp");
 lg->AddEntry(grRatio3,"30#minus100\ / 0#minus0.1\% V0M","lp");
 lg->Draw();

  TH1D *h=(TH1D*)lc->At(1);
  h->GetYaxis()->SetTitleSize(25);
  h->GetYaxis()->SetTitle("Yield ratio to 0#minus0.1%");
  h->GetYaxis()->SetTitleOffset(3);
  h->GetYaxis()->SetLabelSize(20);
  h->GetYaxis()->SetRangeUser(-2.2,1.6);

 for(int i=20;i>=2;i--) {
   TObject *obj = lc->At(i);
   if(obj && (i!=4 && i!=6 && i!=8 && i!=12 && i!=14 && i!=16)) lc->Remove(obj); //remove everything except second, thrid and fourth data points
 }

  TGraphAsymmErrors *gr001 = (TGraphAsymmErrors*)(p1->GetListOfPrimitives()->At(10));
  TGraphAsymmErrors *gr0110 = (TGraphAsymmErrors*)(p1->GetListOfPrimitives()->At(12));
  TGraphAsymmErrors *gr1030 = (TGraphAsymmErrors*)(p1->GetListOfPrimitives()->At(14));
  TGraphAsymmErrors *gr30100 = (TGraphAsymmErrors*)(p1->GetListOfPrimitives()->At(16));
  TGraphAsymmErrors *gr001syst = (TGraphAsymmErrors*)(p1->GetListOfPrimitives()->At(2));
  TGraphAsymmErrors *gr0110syst = (TGraphAsymmErrors*)(p1->GetListOfPrimitives()->At(4));
  TGraphAsymmErrors *gr1030syst = (TGraphAsymmErrors*)(p1->GetListOfPrimitives()->At(6));
  TGraphAsymmErrors *gr30100syst = (TGraphAsymmErrors*)(p1->GetListOfPrimitives()->At(8));

  //1st arg: number of panel, 0 to 7 (NSy*4 ranges, NSw*4 ranges)
  //2nd and 3rd arguments: index of kMin,kMax variations (corresponding to feed-down lower and upwards [+ a contribution from MCclosure also])
  EvaluateRatioPointsAndUnc(0,10,11,gr001,gr0110,gr1030,gr30100,gr001syst,gr0110syst,gr1030syst,gr30100syst,grRatio1,grRatiosyst1,grRatio2,grRatiosyst2,grRatio3,grRatiosyst3);

////DEBUG PROVA POINTS DOPO
  Double_t x, y, e, esl, esh;
    for(int i=0; i<4; i++) {
    grRatio1->GetPoint(i,x,y);
    e = grRatio1->GetErrorYlow(i);
    esl = grRatiosyst1->GetErrorYlow(i);
    esh = grRatiosyst1->GetErrorYhigh(i);
    printf("DOPO il punto ratio1 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);
    grRatio2->GetPoint(i,x,y);
    e = grRatio2->GetErrorYlow(i);
    esl = grRatiosyst2->GetErrorYlow(i);
    esh = grRatiosyst2->GetErrorYhigh(i);
    printf("DOPO il punto ratio2 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);    
    grRatio3->GetPoint(i,x,y);
    e = grRatio3->GetErrorYlow(i);
    esl = grRatiosyst3->GetErrorYlow(i);
    esh = grRatiosyst3->GetErrorYhigh(i);
    printf("DOPO il punto ratio3 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);
  }

  p5->Update();

  grRatio1->SetMarkerSize(grRatio1->GetMarkerSize()-0.4);
  grRatio2->SetMarkerSize(grRatio2->GetMarkerSize()-0.4);
  grRatio3->SetMarkerSize(grRatio3->GetMarkerSize()-0.4);  

  TLine *line = new TLine(0.1,1,27.2,1);
  line->SetLineColor(kGray+1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->Draw("same");

  p5->Draw();

//=========Pad6==========================
 printf("*** PROCESSING PAD 6 ***\n");
 cfinal->cd();
 p6->SetFillStyle(4000);
 p6->SetFrameFillStyle(4000);
 p6->SetPad(0.33,0.5122+0.0092,0.55+0.0048,0.7408+0.017);  // 0.02 Margin = 0.005133 pad dimension
 p6->SetMargin(0.02,0.,0.03,0.0);
 TList *lc=p6->GetListOfPrimitives();

  TGraphAsymmErrors* grRatio1 = (TGraphAsymmErrors*)lc->At(12);
  TGraphAsymmErrors* grRatiosyst1 = (TGraphAsymmErrors*)lc->At(4);
  TGraphAsymmErrors* grRatio2 = (TGraphAsymmErrors*)lc->At(14);
  TGraphAsymmErrors* grRatiosyst2 = (TGraphAsymmErrors*)lc->At(6);
  TGraphAsymmErrors* grRatio3 = (TGraphAsymmErrors*)lc->At(16);
  TGraphAsymmErrors* grRatiosyst3 = (TGraphAsymmErrors*)lc->At(8);

  TH1D *h=(TH1D*)lc->At(1);
  h->GetYaxis()->SetTitleSize(25);
  h->GetYaxis()->SetTitle("");
  h->GetYaxis()->SetTitleOffset(3);
  h->GetYaxis()->SetLabelSize(0);
  h->GetYaxis()->SetRangeUser(-2.2,1.6);

 for(int i=20;i>=2;i--) {
   TObject *obj = lc->At(i);
   if(obj && (i!=4 && i!=6 && i!=8 && i!=12 && i!=14 && i!=16)) lc->Remove(obj); //remove everything except second, thrid and fourth data points
 }

  TGraphAsymmErrors *gr001 = (TGraphAsymmErrors*)(p2->GetListOfPrimitives()->At(10));
  TGraphAsymmErrors *gr0110 = (TGraphAsymmErrors*)(p2->GetListOfPrimitives()->At(12));
  TGraphAsymmErrors *gr1030 = (TGraphAsymmErrors*)(p2->GetListOfPrimitives()->At(14));
  TGraphAsymmErrors *gr30100 = (TGraphAsymmErrors*)(p2->GetListOfPrimitives()->At(16));
  TGraphAsymmErrors *gr001syst = (TGraphAsymmErrors*)(p2->GetListOfPrimitives()->At(2));
  TGraphAsymmErrors *gr0110syst = (TGraphAsymmErrors*)(p2->GetListOfPrimitives()->At(4));
  TGraphAsymmErrors *gr1030syst = (TGraphAsymmErrors*)(p2->GetListOfPrimitives()->At(6));
  TGraphAsymmErrors *gr30100syst = (TGraphAsymmErrors*)(p2->GetListOfPrimitives()->At(8));

  //1st arg: number of panel, 0 to 7 (NSy*4 ranges, NSw*4 ranges)
  //2nd and 3rd arguments: index of kMin,kMax variations (corresponding to feed-down lower and upwards [+ a contribution from MCclosure also])
  EvaluateRatioPointsAndUnc(1,10,11,gr001,gr0110,gr1030,gr30100,gr001syst,gr0110syst,gr1030syst,gr30100syst,grRatio1,grRatiosyst1,grRatio2,grRatiosyst2,grRatio3,grRatiosyst3);

////DEBUG PROVA POINTS DOPO
  Double_t x, y, e, esl, esh;
    for(int i=0; i<4; i++) {
    grRatio1->GetPoint(i,x,y);
    e = grRatio1->GetErrorYlow(i);
    esl = grRatiosyst1->GetErrorYlow(i);
    esh = grRatiosyst1->GetErrorYhigh(i);
    printf("DOPO il punto ratio1 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);
    grRatio2->GetPoint(i,x,y);
    e = grRatio2->GetErrorYlow(i);
    esl = grRatiosyst2->GetErrorYlow(i);
    esh = grRatiosyst2->GetErrorYhigh(i);
    printf("DOPO il punto ratio2 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);    
    grRatio3->GetPoint(i,x,y);
    e = grRatio3->GetErrorYlow(i);
    esl = grRatiosyst3->GetErrorYlow(i);
    esh = grRatiosyst3->GetErrorYhigh(i);
    printf("DOPO il punto ratio3 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);
  }

  p5->Update();

  grRatio1->SetMarkerSize(grRatio1->GetMarkerSize()-0.4);
  grRatio2->SetMarkerSize(grRatio2->GetMarkerSize()-0.4);
  grRatio3->SetMarkerSize(grRatio3->GetMarkerSize()-0.4);  

  TLine *line = new TLine(0.1,1,27.2,1);
  line->SetLineColor(kGray+1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->Draw("same");

 p6->Draw();

//=========Pad7==========================
 printf("*** PROCESSING PAD 7 ***\n");
 cfinal->cd();
 p7->SetFillStyle(4000);
 p7->SetFrameFillStyle(4000);
 p7->SetPad(0.55,0.5122+0.0092,0.76+0.0048,0.7408+0.017);  // 0.02 Margin = 0.005133 pad dimension
 p7->SetMargin(0.02,0.,0.03,0.0);
 TList *lc=p7->GetListOfPrimitives();

  TGraphAsymmErrors* grRatio1 = (TGraphAsymmErrors*)lc->At(12);
  TGraphAsymmErrors* grRatiosyst1 = (TGraphAsymmErrors*)lc->At(4);
  TGraphAsymmErrors* grRatio2 = (TGraphAsymmErrors*)lc->At(14);
  TGraphAsymmErrors* grRatiosyst2 = (TGraphAsymmErrors*)lc->At(6);
  TGraphAsymmErrors* grRatio3 = (TGraphAsymmErrors*)lc->At(16);
  TGraphAsymmErrors* grRatiosyst3 = (TGraphAsymmErrors*)lc->At(8);

  TH1D *h=(TH1D*)lc->At(1);
  h->GetYaxis()->SetTitleSize(25);
  h->GetYaxis()->SetTitle("");
  h->GetYaxis()->SetTitleOffset(3);
  h->GetYaxis()->SetLabelSize(0);
  h->GetYaxis()->SetRangeUser(-2.2,1.6);

 for(int i=20;i>=2;i--) {
   TObject *obj = lc->At(i);
   if(obj && (i!=4 && i!=6 && i!=8 && i!=12 && i!=14 && i!=16)) lc->Remove(obj); //remove everything except second, thrid and fourth data points
 }

  TGraphAsymmErrors *gr001 = (TGraphAsymmErrors*)(p3->GetListOfPrimitives()->At(10));
  TGraphAsymmErrors *gr0110 = (TGraphAsymmErrors*)(p3->GetListOfPrimitives()->At(12));
  TGraphAsymmErrors *gr1030 = (TGraphAsymmErrors*)(p3->GetListOfPrimitives()->At(14));
  TGraphAsymmErrors *gr30100 = (TGraphAsymmErrors*)(p3->GetListOfPrimitives()->At(16));
  TGraphAsymmErrors *gr001syst = (TGraphAsymmErrors*)(p3->GetListOfPrimitives()->At(2));
  TGraphAsymmErrors *gr0110syst = (TGraphAsymmErrors*)(p3->GetListOfPrimitives()->At(4));
  TGraphAsymmErrors *gr1030syst = (TGraphAsymmErrors*)(p3->GetListOfPrimitives()->At(6));
  TGraphAsymmErrors *gr30100syst = (TGraphAsymmErrors*)(p3->GetListOfPrimitives()->At(8));

  //1st arg: number of panel, 0 to 7 (NSy*4 ranges, NSw*4 ranges)
  //2nd and 3rd arguments: index of kMin,kMax variations (corresponding to feed-down lower and upwards [+ a contribution from MCclosure also])
  EvaluateRatioPointsAndUnc(2,10,11,gr001,gr0110,gr1030,gr30100,gr001syst,gr0110syst,gr1030syst,gr30100syst,grRatio1,grRatiosyst1,grRatio2,grRatiosyst2,grRatio3,grRatiosyst3);

////DEBUG PROVA POINTS DOPO
  Double_t x, y, e, esl, esh;
    for(int i=0; i<4; i++) {
    grRatio1->GetPoint(i,x,y);
    e = grRatio1->GetErrorYlow(i);
    esl = grRatiosyst1->GetErrorYlow(i);
    esh = grRatiosyst1->GetErrorYhigh(i);
    printf("DOPO il punto ratio1 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);
    grRatio2->GetPoint(i,x,y);
    e = grRatio2->GetErrorYlow(i);
    esl = grRatiosyst2->GetErrorYlow(i);
    esh = grRatiosyst2->GetErrorYhigh(i);
    printf("DOPO il punto ratio2 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);    
    grRatio3->GetPoint(i,x,y);
    e = grRatio3->GetErrorYlow(i);
    esl = grRatiosyst3->GetErrorYlow(i);
    esh = grRatiosyst3->GetErrorYhigh(i);
    printf("DOPO il punto ratio3 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);
  }

  p5->Update();

  grRatio1->SetMarkerSize(grRatio1->GetMarkerSize()-0.4);
  grRatio2->SetMarkerSize(grRatio2->GetMarkerSize()-0.4);
  grRatio3->SetMarkerSize(grRatio3->GetMarkerSize()-0.4);  

  TLine *line = new TLine(0.1,1,27.2,1);
  line->SetLineColor(kGray+1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->Draw("same");

 p7->Draw();

//=========Pad8==========================
 printf("*** PROCESSING PAD 8 ***\n");
 cfinal->cd();
 p8->SetFillStyle(4000);
 p8->SetFrameFillStyle(4000);
 p8->SetPad(0.76,0.5122+0.0092,0.99,0.7408+0.017);  // 0.02 Margin = 0.005133 pad dimension
 p8->SetMargin(0.02,0.,0.03,0.0);
 TList *lc=p8->GetListOfPrimitives();

  TGraphAsymmErrors* grRatio1 = (TGraphAsymmErrors*)lc->At(12);
  TGraphAsymmErrors* grRatiosyst1 = (TGraphAsymmErrors*)lc->At(4);
  TGraphAsymmErrors* grRatio2 = (TGraphAsymmErrors*)lc->At(14);
  TGraphAsymmErrors* grRatiosyst2 = (TGraphAsymmErrors*)lc->At(6);
  TGraphAsymmErrors* grRatio3 = (TGraphAsymmErrors*)lc->At(16);
  TGraphAsymmErrors* grRatiosyst3 = (TGraphAsymmErrors*)lc->At(8);

  TH1D *h=(TH1D*)lc->At(1);
  h->GetYaxis()->SetTitleSize(25);
  h->GetYaxis()->SetTitle("");
  h->GetYaxis()->SetTitleOffset(3);
  h->GetYaxis()->SetLabelSize(0);
  h->GetYaxis()->SetRangeUser(-2.2,1.6);

 for(int i=20;i>=2;i--) {
   TObject *obj = lc->At(i);
   if(obj && (i!=4 && i!=6 && i!=8 && i!=12 && i!=14 && i!=16)) lc->Remove(obj); //remove everything except second, thrid and fourth data points
 }

  TGraphAsymmErrors *gr001 = (TGraphAsymmErrors*)(p4->GetListOfPrimitives()->At(10));
  TGraphAsymmErrors *gr0110 = (TGraphAsymmErrors*)(p4->GetListOfPrimitives()->At(12));
  TGraphAsymmErrors *gr1030 = (TGraphAsymmErrors*)(p4->GetListOfPrimitives()->At(14));
  TGraphAsymmErrors *gr30100 = (TGraphAsymmErrors*)(p4->GetListOfPrimitives()->At(16));
  TGraphAsymmErrors *gr001syst = (TGraphAsymmErrors*)(p4->GetListOfPrimitives()->At(2));
  TGraphAsymmErrors *gr0110syst = (TGraphAsymmErrors*)(p4->GetListOfPrimitives()->At(4));
  TGraphAsymmErrors *gr1030syst = (TGraphAsymmErrors*)(p4->GetListOfPrimitives()->At(6));
  TGraphAsymmErrors *gr30100syst = (TGraphAsymmErrors*)(p4->GetListOfPrimitives()->At(8));

  //1st arg: number of panel, 0 to 7 (NSy*4 ranges, NSw*4 ranges)
  //2nd and 3rd arguments: index of kMin,kMax variations (corresponding to feed-down lower and upwards [+ a contribution from MCclosure also])
  EvaluateRatioPointsAndUnc(3,10,11,gr001,gr0110,gr1030,gr30100,gr001syst,gr0110syst,gr1030syst,gr30100syst,grRatio1,grRatiosyst1,grRatio2,grRatiosyst2,grRatio3,grRatiosyst3);

////DEBUG PROVA POINTS DOPO
  Double_t x, y, e, esl, esh;
    for(int i=0; i<4; i++) {
    grRatio1->GetPoint(i,x,y);
    e = grRatio1->GetErrorYlow(i);
    esl = grRatiosyst1->GetErrorYlow(i);
    esh = grRatiosyst1->GetErrorYhigh(i);
    printf("DOPO il punto ratio1 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);
    grRatio2->GetPoint(i,x,y);
    e = grRatio2->GetErrorYlow(i);
    esl = grRatiosyst2->GetErrorYlow(i);
    esh = grRatiosyst2->GetErrorYhigh(i);
    printf("DOPO il punto ratio2 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);    
    grRatio3->GetPoint(i,x,y);
    e = grRatio3->GetErrorYlow(i);
    esl = grRatiosyst3->GetErrorYlow(i);
    esh = grRatiosyst3->GetErrorYhigh(i);
    printf("DOPO il punto ratio3 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);
  }

  p5->Update();

  grRatio1->SetMarkerSize(grRatio1->GetMarkerSize()-0.4);
  grRatio2->SetMarkerSize(grRatio2->GetMarkerSize()-0.4);
  grRatio3->SetMarkerSize(grRatio3->GetMarkerSize()-0.4);  

  TLine *line = new TLine(0.1,1,27.2,1);
  line->SetLineColor(kGray+1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->Draw("same");

 p8->Draw();

//=========Pad9==========================
 printf("*** PROCESSING PAD 9 ***\n");
 cfinal->cd();
 p12->SetFillStyle(4000);
 p12->SetFrameFillStyle(4000);
 p12->SetPad(0.,0.2836+0.0092,0.33+0.0048,0.5122+0.016);  // 0.02 Margin = 0.005133 pad dimension
 p12->SetMargin(0.30,0.,0.03,0.0);
 TList *lc=p12->GetListOfPrimitives();
p12->ls();

 TH1D *h=(TH1D*)p12->FindObject("hDraw10");
 h->GetYaxis()->SetTitleSize(25);
 h->GetYaxis()->SetTitleOffset(3);
 h->GetYaxis()->SetTitle("Associated width");
 h->GetYaxis()->SetLabelSize(20);
 h->GetXaxis()->SetLabelSize(0);
 h->GetXaxis()->SetTitle("");

  ((TGraphAsymmErrors*)lc->At(2))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(2))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(3))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(3))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(4))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(4))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(5))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(5))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(6))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(6))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(7))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(7))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(8))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(8))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(9))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(9))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(10))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(10))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(11))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(11))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(12))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(12))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(13))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(13))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(14))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(14))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(15))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(15))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(16))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(16))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(17))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(17))->GetMarkerSize()-0.4);  

 p12->Draw();

//=========Pad10==========================
 printf("*** PROCESSING PAD 10 ***\n");
 cfinal->cd();
 p22->SetFillStyle(4000);
 p22->SetFrameFillStyle(4000);
 p22->SetPad(0.33,0.2836+0.0092,0.55+0.0048,0.5122+0.016);  // 0.02 Margin = 0.005133 pad dimension
 p22->SetMargin(0.02,0.,0.03,0.0);
 TList *lc=p22->GetListOfPrimitives();

 TH1D *h=(TH1D*)p22->FindObject("hDraw11");
 h->GetXaxis()->SetLabelSize(0);
 h->GetXaxis()->SetTitle("");
 h->GetYaxis()->SetLabelSize(0); 

  ((TGraphAsymmErrors*)lc->At(2))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(2))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(3))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(3))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(4))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(4))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(5))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(5))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(6))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(6))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(7))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(7))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(8))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(8))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(9))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(9))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(10))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(10))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(11))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(11))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(12))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(12))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(13))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(13))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(14))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(14))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(15))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(15))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(16))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(16))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(17))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(17))->GetMarkerSize()-0.4);  

 p22->Draw();

//=========Pad11==========================
 printf("*** PROCESSING PAD 11 ***\n");
 cfinal->cd();
 p32->SetFillStyle(4000);
 p32->SetFrameFillStyle(4000);
 p32->SetPad(0.55,0.2836+0.0092,0.76+0.0048,0.5122+0.016);  // 0.02 Margin = 0.005133 pad dimension
 p32->SetMargin(0.02,0.,0.03,0.0);
 TList *lc=p32->GetListOfPrimitives();

 TH1D *h=(TH1D*)p32->FindObject("hDraw12");
 h->GetXaxis()->SetLabelSize(0);
 h->GetXaxis()->SetTitle("");
 h->GetYaxis()->SetLabelSize(0); 

  ((TGraphAsymmErrors*)lc->At(2))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(2))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(3))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(3))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(4))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(4))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(5))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(5))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(6))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(6))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(7))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(7))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(8))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(8))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(9))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(9))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(10))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(10))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(11))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(11))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(12))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(12))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(13))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(13))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(14))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(14))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(15))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(15))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(16))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(16))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(17))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(17))->GetMarkerSize()-0.4);  

 p32->Draw();

//=========Pad12==========================
 printf("*** PROCESSING PAD 12 ***\n");
 cfinal->cd();
 p42->SetFillStyle(4000);
 p42->SetFrameFillStyle(4000);
 p42->SetPad(0.76,0.2836+0.0092,0.99,0.5122+0.016);  // 0.02 Margin = 0.005133 pad dimension
 p42->SetMargin(0.02,0.,0.03,0.0);
 TList *lc=p42->GetListOfPrimitives();

 TH1D *h=(TH1D*)p42->FindObject("hDraw13");
 h->GetXaxis()->SetLabelSize(0);
 h->GetXaxis()->SetTitle("");
 h->GetYaxis()->SetLabelSize(0); 

  ((TGraphAsymmErrors*)lc->At(2))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(2))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(3))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(3))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(4))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(4))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(5))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(5))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(6))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(6))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(7))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(7))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(8))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(8))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(9))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(9))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(10))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(10))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(11))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(11))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(12))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(12))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(13))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(13))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(14))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(14))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(15))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(15))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(16))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(16))->GetMarkerSize()-0.4);
  ((TGraphAsymmErrors*)lc->At(17))->SetMarkerSize(((TGraphAsymmErrors*)lc->At(17))->GetMarkerSize()-0.4);  

 p42->Draw();

//=========Pad13==========================
 printf("*** PROCESSING PAD 13 ***\n");
 cfinal->cd();
 p52->SetFillStyle(4000);
 p52->SetFrameFillStyle(4000);
 p52->SetPad(0.,0.0,0.33+0.0048,0.2836+0.016);  // 0.02 Margin = 0.005133 pad dimension
 p52->SetMargin(0.30,0.,0.18,0.0);
 TList *lc=p52->GetListOfPrimitives();  

  TGraphAsymmErrors* grRatio1 = (TGraphAsymmErrors*)lc->At(12);
  TGraphAsymmErrors* grRatiosyst1 = (TGraphAsymmErrors*)lc->At(4);
  TGraphAsymmErrors* grRatio2 = (TGraphAsymmErrors*)lc->At(14);
  TGraphAsymmErrors* grRatiosyst2 = (TGraphAsymmErrors*)lc->At(6);
  TGraphAsymmErrors* grRatio3 = (TGraphAsymmErrors*)lc->At(16);
  TGraphAsymmErrors* grRatiosyst3 = (TGraphAsymmErrors*)lc->At(8);

  TH1D *h=(TH1D*)lc->At(1);
  h->GetXaxis()->SetTitleSize(22);  
  h->GetXaxis()->SetTitleOffset(3.5);
  h->GetXaxis()->SetLabelSize(20);  
  h->GetYaxis()->SetTitleSize(25);
  h->GetYaxis()->SetTitle("Width ratio to 0#minus0.1%");
  h->GetYaxis()->SetTitleOffset(3);
  h->GetYaxis()->SetLabelSize(20);
  h->GetYaxis()->SetRangeUser(-2.2,2.2);

 for(int i=20;i>=2;i--) {
   TObject *obj = lc->At(i);
   if(obj && (i!=4 && i!=6 && i!=8 && i!=12 && i!=14 && i!=16)) lc->Remove(obj); //remove everything except second, thrid and fourth data points
 }

  TGraphAsymmErrors *gr001 = (TGraphAsymmErrors*)(p12->GetListOfPrimitives()->At(10));
  TGraphAsymmErrors *gr0110 = (TGraphAsymmErrors*)(p12->GetListOfPrimitives()->At(12));
  TGraphAsymmErrors *gr1030 = (TGraphAsymmErrors*)(p12->GetListOfPrimitives()->At(14));
  TGraphAsymmErrors *gr30100 = (TGraphAsymmErrors*)(p12->GetListOfPrimitives()->At(16));
  TGraphAsymmErrors *gr001syst = (TGraphAsymmErrors*)(p12->GetListOfPrimitives()->At(2));
  TGraphAsymmErrors *gr0110syst = (TGraphAsymmErrors*)(p12->GetListOfPrimitives()->At(4));
  TGraphAsymmErrors *gr1030syst = (TGraphAsymmErrors*)(p12->GetListOfPrimitives()->At(6));
  TGraphAsymmErrors *gr30100syst = (TGraphAsymmErrors*)(p12->GetListOfPrimitives()->At(8));

  //1st arg: number of panel, 0 to 7 (NSy*4 ranges, NSw*4 ranges)
  //2nd and 3rd arguments: index of kMin,kMax variations (corresponding to feed-down lower and upwards [+ a contribution from MCclosure also])
  EvaluateRatioPointsAndUnc(4,10,11,gr001,gr0110,gr1030,gr30100,gr001syst,gr0110syst,gr1030syst,gr30100syst,grRatio1,grRatiosyst1,grRatio2,grRatiosyst2,grRatio3,grRatiosyst3);

////DEBUG PROVA POINTS DOPO
  Double_t x, y, e, esl, esh;
    for(int i=0; i<4; i++) {
    grRatio1->GetPoint(i,x,y);
    e = grRatio1->GetErrorYlow(i);
    esl = grRatiosyst1->GetErrorYlow(i);
    esh = grRatiosyst1->GetErrorYhigh(i);
    printf("DOPO il punto ratio1 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);
    grRatio2->GetPoint(i,x,y);
    e = grRatio2->GetErrorYlow(i);
    esl = grRatiosyst2->GetErrorYlow(i);
    esh = grRatiosyst2->GetErrorYhigh(i);
    printf("DOPO il punto ratio2 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);    
    grRatio3->GetPoint(i,x,y);
    e = grRatio3->GetErrorYlow(i);
    esl = grRatiosyst3->GetErrorYlow(i);
    esh = grRatiosyst3->GetErrorYhigh(i);
    printf("DOPO il punto ratio3 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);
  }

  p5->Update();

  grRatio1->SetMarkerSize(grRatio1->GetMarkerSize()-0.4);
  grRatio2->SetMarkerSize(grRatio2->GetMarkerSize()-0.4);
  grRatio3->SetMarkerSize(grRatio3->GetMarkerSize()-0.4);  

  TLine *line = new TLine(0.1,1,27.2,1);
  line->SetLineColor(kGray+1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->Draw("same");


 p52->Draw();

//=========Pad14==========================
 printf("*** PROCESSING PAD 14 ***\n");
 cfinal->cd();
 p62->SetFillStyle(4000);
 p62->SetFrameFillStyle(4000);
 p62->SetPad(0.33,0.0,0.55+0.0048,0.2836+0.016);  // 0.02 Margin = 0.005133 pad dimension
 p62->SetMargin(0.02,0.,0.18,0.0);
 TList *lc=p62->GetListOfPrimitives();  

  TGraphAsymmErrors* grRatio1 = (TGraphAsymmErrors*)lc->At(12);
  TGraphAsymmErrors* grRatiosyst1 = (TGraphAsymmErrors*)lc->At(4);
  TGraphAsymmErrors* grRatio2 = (TGraphAsymmErrors*)lc->At(14);
  TGraphAsymmErrors* grRatiosyst2 = (TGraphAsymmErrors*)lc->At(6);
  TGraphAsymmErrors* grRatio3 = (TGraphAsymmErrors*)lc->At(16);
  TGraphAsymmErrors* grRatiosyst3 = (TGraphAsymmErrors*)lc->At(8);

  TH1D *h=(TH1D*)lc->At(1);
  h->GetXaxis()->SetTitleSize(22);  
  h->GetXaxis()->SetTitleOffset(3.5);
  h->GetXaxis()->SetLabelSize(20);  
  h->GetYaxis()->SetTitleSize(25);
  h->GetYaxis()->SetTitle("");
  h->GetYaxis()->SetTitleOffset(3);
  h->GetYaxis()->SetLabelSize(0);
  h->GetYaxis()->SetRangeUser(-2.2,2.2);

 for(int i=20;i>=2;i--) {
   TObject *obj = lc->At(i);
   if(obj && (i!=4 && i!=6 && i!=8 && i!=12 && i!=14 && i!=16)) lc->Remove(obj); //remove everything except second, thrid and fourth data points
 }

  TGraphAsymmErrors *gr001 = (TGraphAsymmErrors*)(p22->GetListOfPrimitives()->At(10));
  TGraphAsymmErrors *gr0110 = (TGraphAsymmErrors*)(p22->GetListOfPrimitives()->At(12));
  TGraphAsymmErrors *gr1030 = (TGraphAsymmErrors*)(p22->GetListOfPrimitives()->At(14));
  TGraphAsymmErrors *gr30100 = (TGraphAsymmErrors*)(p22->GetListOfPrimitives()->At(16));
  TGraphAsymmErrors *gr001syst = (TGraphAsymmErrors*)(p22->GetListOfPrimitives()->At(2));
  TGraphAsymmErrors *gr0110syst = (TGraphAsymmErrors*)(p22->GetListOfPrimitives()->At(4));
  TGraphAsymmErrors *gr1030syst = (TGraphAsymmErrors*)(p22->GetListOfPrimitives()->At(6));
  TGraphAsymmErrors *gr30100syst = (TGraphAsymmErrors*)(p22->GetListOfPrimitives()->At(8));

  //1st arg: number of panel, 0 to 7 (NSy*4 ranges, NSw*4 ranges)
  //2nd and 3rd arguments: index of kMin,kMax variations (corresponding to feed-down lower and upwards [+ a contribution from MCclosure also])
  EvaluateRatioPointsAndUnc(5,10,11,gr001,gr0110,gr1030,gr30100,gr001syst,gr0110syst,gr1030syst,gr30100syst,grRatio1,grRatiosyst1,grRatio2,grRatiosyst2,grRatio3,grRatiosyst3);

////DEBUG PROVA POINTS DOPO
  Double_t x, y, e, esl, esh;
    for(int i=0; i<4; i++) {
    grRatio1->GetPoint(i,x,y);
    e = grRatio1->GetErrorYlow(i);
    esl = grRatiosyst1->GetErrorYlow(i);
    esh = grRatiosyst1->GetErrorYhigh(i);
    printf("DOPO il punto ratio1 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);
    grRatio2->GetPoint(i,x,y);
    e = grRatio2->GetErrorYlow(i);
    esl = grRatiosyst2->GetErrorYlow(i);
    esh = grRatiosyst2->GetErrorYhigh(i);
    printf("DOPO il punto ratio2 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);    
    grRatio3->GetPoint(i,x,y);
    e = grRatio3->GetErrorYlow(i);
    esl = grRatiosyst3->GetErrorYlow(i);
    esh = grRatiosyst3->GetErrorYhigh(i);
    printf("DOPO il punto ratio3 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);
  }

  p5->Update();

  grRatio1->SetMarkerSize(grRatio1->GetMarkerSize()-0.4);
  grRatio2->SetMarkerSize(grRatio2->GetMarkerSize()-0.4);
  grRatio3->SetMarkerSize(grRatio3->GetMarkerSize()-0.4);  

  TLine *line = new TLine(0.1,1,27.2,1);
  line->SetLineColor(kGray+1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->Draw("same");


 p62->Draw();

//=========Pad15==========================
 printf("*** PROCESSING PAD 15 ***\n");
 cfinal->cd();
 p72->SetFillStyle(4000);
 p72->SetFrameFillStyle(4000);
 p72->SetPad(0.55,0.0,0.76+0.0048,0.2836+0.016);  // 0.02 Margin = 0.005133 pad dimension
 p72->SetMargin(0.02,0.,0.18,0.0);
 TList *lc=p72->GetListOfPrimitives();  

  TGraphAsymmErrors* grRatio1 = (TGraphAsymmErrors*)lc->At(12);
  TGraphAsymmErrors* grRatiosyst1 = (TGraphAsymmErrors*)lc->At(4);
  TGraphAsymmErrors* grRatio2 = (TGraphAsymmErrors*)lc->At(14);
  TGraphAsymmErrors* grRatiosyst2 = (TGraphAsymmErrors*)lc->At(6);
  TGraphAsymmErrors* grRatio3 = (TGraphAsymmErrors*)lc->At(16);
  TGraphAsymmErrors* grRatiosyst3 = (TGraphAsymmErrors*)lc->At(8);

  TH1D *h=(TH1D*)lc->At(1);
  h->GetXaxis()->SetTitleSize(22);  
  h->GetXaxis()->SetTitleOffset(3.5);
  h->GetXaxis()->SetLabelSize(20);
  h->GetYaxis()->SetTitleSize(25);
  h->GetYaxis()->SetTitle("");
  h->GetYaxis()->SetTitleOffset(3);
  h->GetYaxis()->SetLabelSize(0);
  h->GetYaxis()->SetRangeUser(-2.2,2.2);

 for(int i=20;i>=2;i--) {
   TObject *obj = lc->At(i);
   if(obj && (i!=4 && i!=6 && i!=8 && i!=12 && i!=14 && i!=16)) lc->Remove(obj); //remove everything except second, thrid and fourth data points
 }

  TGraphAsymmErrors *gr001 = (TGraphAsymmErrors*)(p32->GetListOfPrimitives()->At(10));
  TGraphAsymmErrors *gr0110 = (TGraphAsymmErrors*)(p32->GetListOfPrimitives()->At(12));
  TGraphAsymmErrors *gr1030 = (TGraphAsymmErrors*)(p32->GetListOfPrimitives()->At(14));
  TGraphAsymmErrors *gr30100 = (TGraphAsymmErrors*)(p32->GetListOfPrimitives()->At(16));
  TGraphAsymmErrors *gr001syst = (TGraphAsymmErrors*)(p32->GetListOfPrimitives()->At(2));
  TGraphAsymmErrors *gr0110syst = (TGraphAsymmErrors*)(p32->GetListOfPrimitives()->At(4));
  TGraphAsymmErrors *gr1030syst = (TGraphAsymmErrors*)(p32->GetListOfPrimitives()->At(6));
  TGraphAsymmErrors *gr30100syst = (TGraphAsymmErrors*)(p32->GetListOfPrimitives()->At(8));

  //1st arg: number of panel, 0 to 7 (NSy*4 ranges, NSw*4 ranges)
  //2nd and 3rd arguments: index of kMin,kMax variations (corresponding to feed-down lower and upwards [+ a contribution from MCclosure also])
  EvaluateRatioPointsAndUnc(6,10,11,gr001,gr0110,gr1030,gr30100,gr001syst,gr0110syst,gr1030syst,gr30100syst,grRatio1,grRatiosyst1,grRatio2,grRatiosyst2,grRatio3,grRatiosyst3);

////DEBUG PROVA POINTS DOPO
  Double_t x, y, e, esl, esh;
    for(int i=0; i<4; i++) {
    grRatio1->GetPoint(i,x,y);
    e = grRatio1->GetErrorYlow(i);
    esl = grRatiosyst1->GetErrorYlow(i);
    esh = grRatiosyst1->GetErrorYhigh(i);
    printf("DOPO il punto ratio1 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);
    grRatio2->GetPoint(i,x,y);
    e = grRatio2->GetErrorYlow(i);
    esl = grRatiosyst2->GetErrorYlow(i);
    esh = grRatiosyst2->GetErrorYhigh(i);
    printf("DOPO il punto ratio2 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);    
    grRatio3->GetPoint(i,x,y);
    e = grRatio3->GetErrorYlow(i);
    esl = grRatiosyst3->GetErrorYlow(i);
    esh = grRatiosyst3->GetErrorYhigh(i);
    printf("DOPO il punto ratio3 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);
  }

  p5->Update();

  grRatio1->SetMarkerSize(grRatio1->GetMarkerSize()-0.4);
  grRatio2->SetMarkerSize(grRatio2->GetMarkerSize()-0.4);
  grRatio3->SetMarkerSize(grRatio3->GetMarkerSize()-0.4);  

  TLine *line = new TLine(0.1,1,27.2,1);
  line->SetLineColor(kGray+1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->Draw("same");

 p72->Draw();

//=========Pad16==========================

 printf("*** PROCESSING PAD 16 ***\n");
 cfinal->cd();
 p82->SetFillStyle(4000);
 p82->SetFrameFillStyle(4000);
 p82->SetPad(0.76,0.0,0.99,0.2836+0.016);  // 0.02 Margin = 0.005133 pad dimension
 p82->SetMargin(0.02,0.,0.18,0.0);
 TList *lc=p82->GetListOfPrimitives();  

  TGraphAsymmErrors* grRatio1 = (TGraphAsymmErrors*)lc->At(12);
  TGraphAsymmErrors* grRatiosyst1 = (TGraphAsymmErrors*)lc->At(4);
  TGraphAsymmErrors* grRatio2 = (TGraphAsymmErrors*)lc->At(14);
  TGraphAsymmErrors* grRatiosyst2 = (TGraphAsymmErrors*)lc->At(6);
  TGraphAsymmErrors* grRatio3 = (TGraphAsymmErrors*)lc->At(16);
  TGraphAsymmErrors* grRatiosyst3 = (TGraphAsymmErrors*)lc->At(8);

  TH1D *h=(TH1D*)lc->At(1);
  h->GetXaxis()->SetTitleSize(22);  
  h->GetXaxis()->SetTitleOffset(3.5);  
  h->GetXaxis()->SetLabelSize(20);
  h->GetYaxis()->SetTitleSize(25);
  h->GetYaxis()->SetTitle("");
  h->GetYaxis()->SetTitleOffset(3);
  h->GetYaxis()->SetLabelSize(0);
  h->GetYaxis()->SetRangeUser(-2.2,2.2);

 for(int i=20;i>=2;i--) {
   TObject *obj = lc->At(i);
   if(obj && (i!=4 && i!=6 && i!=8 && i!=12 && i!=14 && i!=16)) lc->Remove(obj); //remove everything except second, thrid and fourth data points
 }

  TGraphAsymmErrors *gr001 = (TGraphAsymmErrors*)(p42->GetListOfPrimitives()->At(10));
  TGraphAsymmErrors *gr0110 = (TGraphAsymmErrors*)(p42->GetListOfPrimitives()->At(12));
  TGraphAsymmErrors *gr1030 = (TGraphAsymmErrors*)(p42->GetListOfPrimitives()->At(14));
  TGraphAsymmErrors *gr30100 = (TGraphAsymmErrors*)(p42->GetListOfPrimitives()->At(16));
  TGraphAsymmErrors *gr001syst = (TGraphAsymmErrors*)(p42->GetListOfPrimitives()->At(2));
  TGraphAsymmErrors *gr0110syst = (TGraphAsymmErrors*)(p42->GetListOfPrimitives()->At(4));
  TGraphAsymmErrors *gr1030syst = (TGraphAsymmErrors*)(p42->GetListOfPrimitives()->At(6));
  TGraphAsymmErrors *gr30100syst = (TGraphAsymmErrors*)(p42->GetListOfPrimitives()->At(8));

  //1st arg: number of panel, 0 to 7 (NSy*4 ranges, NSw*4 ranges)
  //2nd and 3rd arguments: index of kMin,kMax variations (corresponding to feed-down lower and upwards [+ a contribution from MCclosure also])
  EvaluateRatioPointsAndUnc(7,10,11,gr001,gr0110,gr1030,gr30100,gr001syst,gr0110syst,gr1030syst,gr30100syst,grRatio1,grRatiosyst1,grRatio2,grRatiosyst2,grRatio3,grRatiosyst3);

////DEBUG PROVA POINTS DOPO
  Double_t x, y, e, esl, esh;
    for(int i=0; i<4; i++) {
    grRatio1->GetPoint(i,x,y);
    e = grRatio1->GetErrorYlow(i);
    esl = grRatiosyst1->GetErrorYlow(i);
    esh = grRatiosyst1->GetErrorYhigh(i);
    printf("DOPO il punto ratio1 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);
    grRatio2->GetPoint(i,x,y);
    e = grRatio2->GetErrorYlow(i);
    esl = grRatiosyst2->GetErrorYlow(i);
    esh = grRatiosyst2->GetErrorYhigh(i);
    printf("DOPO il punto ratio2 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);    
    grRatio3->GetPoint(i,x,y);
    e = grRatio3->GetErrorYlow(i);
    esl = grRatiosyst3->GetErrorYlow(i);
    esh = grRatiosyst3->GetErrorYhigh(i);
    printf("DOPO il punto ratio3 ha (pTbin=%.1f): %f, %f - err stat %f - err syst %f, %f\n",i,x,y,e,esl,esh);
  }

  p5->Update();

  grRatio1->SetMarkerSize(grRatio1->GetMarkerSize()-0.4);
  grRatio2->SetMarkerSize(grRatio2->GetMarkerSize()-0.4);
  grRatio3->SetMarkerSize(grRatio3->GetMarkerSize()-0.4);  

  TLine *line = new TLine(0.1,1,27.2,1);
  line->SetLineColor(kGray+1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->Draw("same");

  p82->Draw();

 //Final update and save of file

 cfinal->Modified();
 cfinal->Update();

 cfinal->SaveAs("Centrality_Ratios_NS_Final.eps");
 cfinal->SaveAs("Centrality_Ratios_NS_Final.png");
 cfinal->SaveAs("Centrality_Ratios_NS_Final.pdf");
 cfinal->SaveAs("Centrality_Ratios_NS_Final.root");

}


void EvaluateRatioPointsAndUnc(Int_t ratiopanel, Int_t kMinInd, Int_t kMaxInd, TGraphAsymmErrors *gr001, TGraphAsymmErrors *gr0110, TGraphAsymmErrors *gr1030, TGraphAsymmErrors *gr30100, TGraphAsymmErrors *gr001syst, TGraphAsymmErrors *gr0110syst, TGraphAsymmErrors *gr1030syst, TGraphAsymmErrors *gr30100syst, TGraphAsymmErrors *grRatio1, TGraphAsymmErrors *grRatiosyst1, TGraphAsymmErrors *grRatio2, TGraphAsymmErrors *grRatiosyst2, TGraphAsymmErrors *grRatio3, TGraphAsymmErrors *grRatiosyst3) {

Double_t x001[4], y001[4], x0110[4], y0110[4], x1030[4], y1030[4], x30100[4], y30100[4], ex[4], ey001stat[4], ey0110stat[4], ey1030stat[4], ey30100stat[4], ey001systlow[4], ey0110systlow[4], ey1030systlow[4], ey30100systlow[4], ey001systhigh[4], ey0110systhigh[4], ey1030systhigh[4], ey30100systhigh[4];
  for(int iD=0; iD<4; iD++) {
    gr001->GetPoint(iD,x001[iD],y001[iD]);
    gr0110->GetPoint(iD,x0110[iD],y0110[iD]);
    gr1030->GetPoint(iD,x1030[iD],y1030[iD]);
    gr30100->GetPoint(iD,x30100[iD],y30100[iD]);
    ex[iD] = gr0110->GetErrorX(iD);
    ey001stat[iD] = gr001->GetErrorY(iD);
    ey0110stat[iD] = gr0110->GetErrorY(iD);
    ey1030stat[iD] = gr1030->GetErrorY(iD);
    ey30100stat[iD] = gr30100->GetErrorY(iD);
    ey001systlow[iD] = gr001syst->GetErrorYlow(iD);
    ey001systhigh[iD] = gr001syst->GetErrorYhigh(iD);
    ey0110systlow[iD] = gr0110syst->GetErrorYlow(iD);
    ey0110systhigh[iD] = gr0110syst->GetErrorYhigh(iD);
    ey1030systlow[iD] = gr1030syst->GetErrorYlow(iD);
    ey1030systhigh[iD] = gr1030syst->GetErrorYhigh(iD);
    ey30100systlow[iD] = gr30100syst->GetErrorYlow(iD);
    ey30100systhigh[iD] = gr30100syst->GetErrorYhigh(iD);    
    printf("001   point %d (x %f) = %f, stat %f, systl %f, systh %f\n",iD,x001[iD],y001[iD],ey001stat[iD],ey001systlow[iD],ey001systhigh[iD]);
    printf("0110  point %d (x %f) = %f, stat %f, systl %f, systh %f\n",iD,x0110[iD],y0110[iD],ey0110stat[iD],ey0110systlow[iD],ey0110systhigh[iD]);
    printf("1030 point %d (x %f) = %f, stat %f, systl %f, systh %f\n",iD,x1030[iD],y1030[iD],ey1030stat[iD],ey1030systlow[iD],ey1030systhigh[iD]);
    printf("30100 point %d (x %f) = %f, stat %f, systl %f, systh %f\n",iD,x30100[iD],y30100[iD],ey30100stat[iD],ey30100systlow[iD],ey30100systhigh[iD]);
  }
  
  Double_t ratio1, ratio2, ratio3; //trivially, the three ratios of central points 0110/001, 1030/001, 30100/001
  Double_t ratiostat1, ratiostat2, ratiostat3; //stat error to ratios 0110/001, 1030/001, 30100/001
  Double_t ratiosystUp1, ratiosystUp2, ratiosystUp3; //syst error to ratios 0110/001, 1030/001, 30100/001
  Double_t ratiosystDown1, ratiosystDown2, ratiosystDown3; //syst error to ratios 0110/001, 1030/001, 30100/001

  //group 1: the symmetric, dPhi-related uncertainties (excluded for widths)
  Double_t yield_001, yield_0110, yield_1030, yield_30100;
  Double_t bkgshape_001, bkgshape_0110, bkgshape_1030, bkgshape_30100;
  Double_t fitvar_001, fitvar_0110, fitvar_1030, fitvar_30100;
  Double_t sum_GROUP1_001, sum_GROUP1_0110, sum_GROUP1_1030, sum_GROUP1_30100; //sum in quadrature of the above three, for each cent
  Double_t rel_GROUP1_Ratio1, rel_GROUP1_Ratio2, rel_GROUP1_Ratio3; //the propagation af the above sum to the three ratios!

  //group 2: the v2 uncertainty - envelope of no-v2-hyp ratio (standard) and v2-hyp ratio
  Double_t v2relmin_001, v2relmin_0110, v2relmin_1030, v2relmin_30100;
  Double_t v2relmax_001, v2relmax_0110, v2relmax_1030, v2relmax_30100;
  Double_t v2absmin_001, v2absmin_0110, v2absmin_1030, v2absmin_30100;
  Double_t v2absmax_001, v2absmax_0110, v2absmax_1030, v2absmax_30100;
  Double_t v2absmin_Ratio1, v2absmin_Ratio2, v2absmin_Ratio3; //ratios of the v2 hypotheses
  Double_t v2absmax_Ratio1, v2absmax_Ratio2, v2absmax_Ratio3;
  Double_t v2_ENVLOW_abs_Ratio1, v2_ENVLOW_abs_Ratio2, v2_ENVLOW_abs_Ratio3; //downward envelope (absolute)
  Double_t v2_ENVHIG_abs_Ratio1, v2_ENVHIG_abs_Ratio2, v2_ENVHIG_abs_Ratio3; //upward envelope (absolute)
  Double_t v2_ENVLOW_rel_Ratio1, v2_ENVLOW_rel_Ratio2, v2_ENVLOW_rel_Ratio3; //downward envelope (relative, for further propagation)
  Double_t v2_ENVHIG_rel_Ratio1, v2_ENVHIG_rel_Ratio2, v2_ENVHIG_rel_Ratio3; //upward envelope (relative, for further propagation)    
 
  //group 3: the feed-down uncertainty (includes also MC-closure, same treatment) 
  // ->envelope of minvar/minvar and maxvar/maxvar ratio with flat and extreme B/D vs cent hypotheses
  Double_t FDmin_BDflat_abs_001, FDmin_BDflat_abs_0110, FDmin_BDflat_abs_1030, FDmin_BDflat_abs_30100; //is the absolute value of kMin case in BDflat scenario
  Double_t FDmax_BDflat_abs_001, FDmax_BDflat_abs_0110, FDmax_BDflat_abs_1030, FDmax_BDflat_abs_30100; //is the absolute value of kMax case in BDflat scenario
  Double_t FDmin_BDextr_abs_001, FDmin_BDextr_abs_0110, FDmin_BDextr_abs_1030, FDmin_BDextr_abs_30100; //is the absolute value of kMin case in BDextr scenario
  Double_t FDmax_BDextr_abs_001, FDmax_BDextr_abs_0110, FDmax_BDextr_abs_1030, FDmax_BDextr_abs_30100; //is the absolute value of kMax case in BDextr scenario
  Double_t FDmin_BDflat_abs_Ratio1, FDmin_BDflat_abs_Ratio2, FDmin_BDflat_abs_Ratio3; //ratios of centralities
  Double_t FDmax_BDflat_abs_Ratio1, FDmax_BDflat_abs_Ratio2, FDmax_BDflat_abs_Ratio3;
  Double_t FDmin_BDextr_abs_Ratio1, FDmin_BDextr_abs_Ratio2, FDmin_BDextr_abs_Ratio3; //ratios of centralities
  Double_t FDmax_BDextr_abs_Ratio1, FDmax_BDextr_abs_Ratio2, FDmax_BDextr_abs_Ratio3;
  Double_t FD_ENVLOW_abs_Ratio1, FD_ENVLOW_abs_Ratio2, FD_ENVLOW_abs_Ratio3; //downward envelope (absolute)
  Double_t FD_ENVHIG_abs_Ratio1, FD_ENVHIG_abs_Ratio2, FD_ENVHIG_abs_Ratio3; //upward envelope (absolute)
  Double_t FD_ENVLOW_rel_Ratio1, FD_ENVLOW_rel_Ratio2, FD_ENVLOW_rel_Ratio3; //downward envelope (relative, for further propagation)
  Double_t FD_ENVHIG_rel_Ratio1, FD_ENVHIG_rel_Ratio2, FD_ENVHIG_rel_Ratio3; //upward envelope (relative, for further propagation)      

  //support ingredients (overwritten at each loop)
  AliHFDhadronCorrSystUnc *a_001;
  AliHFDhadronCorrSystUnc *a_0110;
  AliHFDhadronCorrSystUnc *a_1030;
  AliHFDhadronCorrSystUnc *a_30100;
  TCanvas *c_001;
  TCanvas *c_0110;
  TCanvas *c_1030;
  TCanvas *c_30100;
  TGraphAsymmErrors *gr_fitvar_001;
  TGraphAsymmErrors *gr_fitvar_0110;
  TGraphAsymmErrors *gr_fitvar_1030;
  TGraphAsymmErrors *gr_fitvar_30100;
  TGraphAsymmErrors *gr_v2var_001;
  TGraphAsymmErrors *gr_v2var_0110;
  TGraphAsymmErrors *gr_v2var_1030;
  TGraphAsymmErrors *gr_v2var_30100;
  TH1D *hist_FDmin_BDflat_001;
  TH1D *hist_FDmax_BDflat_001;
  TH1D *hist_FDmin_BDextr_001;
  TH1D *hist_FDmax_BDextr_001;
  TH1D *hist_FDmin_BDflat_0110;
  TH1D *hist_FDmax_BDflat_0110;
  TH1D *hist_FDmin_BDextr_0110;
  TH1D *hist_FDmax_BDextr_0110;
  TH1D *hist_FDmin_BDflat_1030;
  TH1D *hist_FDmax_BDflat_1030;
  TH1D *hist_FDmin_BDextr_1030;
  TH1D *hist_FDmax_BDextr_1030;
  TH1D *hist_FDmin_BDflat_30100;
  TH1D *hist_FDmax_BDflat_30100;
  TH1D *hist_FDmin_BDextr_30100;
  TH1D *hist_FDmax_BDextr_30100;

  for(int iD=0; iD<4; iD++) {
    printf("************************************\n");
    printf("****** Panel %d, pT(D) bin %d ******\n", ratiopanel,iD);
    printf("************************************\n");

    //central ratio and stat error
    ratio1 = y0110[iD]/y001[iD];
    ratio2 = y1030[iD]/y001[iD];
    ratio3 = y30100[iD]/y001[iD];
    ratiostat1 = ratio1*TMath::Sqrt(ey0110stat[iD]*ey0110stat[iD]/(y0110[iD]*y0110[iD]) + ey001stat[iD]*ey001stat[iD]/(y001[iD]*y001[iD])); //absolute error, which is the one to insert in the tgraph bars
    ratiostat2 = ratio2*TMath::Sqrt(ey1030stat[iD]*ey1030stat[iD]/(y1030[iD]*y1030[iD]) + ey001stat[iD]*ey001stat[iD]/(y001[iD]*y001[iD])); //absolute error, which is the one to insert in the tgraph bars
    ratiostat3 = ratio3*TMath::Sqrt(ey30100stat[iD]*ey30100stat[iD]/(y30100[iD]*y30100[iD]) + ey001stat[iD]*ey001stat[iD]/(y001[iD]*y001[iD])); //absolute error, which is the one to insert in the tgraph bars

    //syst group 1: extraction
    Int_t iAss = ratiopanel;
    if(iAss > 3) iAss -= 4;
    a_001 = (AliHFDhadronCorrSystUnc*)dPhisyst_001[iAss][iD]->Get("AverageSystematicUncertainty");
    a_0110 = (AliHFDhadronCorrSystUnc*)dPhisyst_0110[iAss][iD]->Get("AverageSystematicUncertainty");
    a_1030 = (AliHFDhadronCorrSystUnc*)dPhisyst_1030[iAss][iD]->Get("AverageSystematicUncertainty");
    a_30100 = (AliHFDhadronCorrSystUnc*)dPhisyst_30100[iAss][iD]->Get("AverageSystematicUncertainty");
    yield_001 = a_001->GetHistoYieldUnc()->GetBinContent(1); //histos are flat in dPhi, all bins are equal
    yield_0110 = a_0110->GetHistoYieldUnc()->GetBinContent(1); 
    yield_1030 = a_1030->GetHistoYieldUnc()->GetBinContent(1); 
    yield_30100 = a_30100->GetHistoYieldUnc()->GetBinContent(1); 
    bkgshape_001 = a_001->GetHistoBackSubUncMax()->GetBinContent(1); //histos are flat in dPhi, all bins are equal
    bkgshape_0110 = a_0110->GetHistoBackSubUncMax()->GetBinContent(1);
    bkgshape_1030 = a_1030->GetHistoBackSubUncMax()->GetBinContent(1);
    bkgshape_30100 = a_30100->GetHistoBackSubUncMax()->GetBinContent(1);

    if(ratiopanel < 4) { //NS yield
      c_001 = (TCanvas*)totalsyst_NSy_001[iAss]->Get("TotalSystematicSourcesNSYield");
      c_0110 = (TCanvas*)totalsyst_NSy_0110[iAss]->Get("TotalSystematicSourcesNSYield");
      c_1030 = (TCanvas*)totalsyst_NSy_1030[iAss]->Get("TotalSystematicSourcesNSYield");
      c_30100 = (TCanvas*)totalsyst_NSy_30100[iAss]->Get("TotalSystematicSourcesNSYield");
      gr_fitvar_001 = (TGraphAsymmErrors*)c_001->FindObject("grSystBasNSYield");
      gr_fitvar_0110 = (TGraphAsymmErrors*)c_0110->FindObject("grSystBasNSYield");
      gr_fitvar_1030 = (TGraphAsymmErrors*)c_1030->FindObject("grSystBasNSYield");
      gr_fitvar_30100 = (TGraphAsymmErrors*)c_30100->FindObject("grSystBasNSYield");
    } else { //NS width
      c_001 = (TCanvas*)totalsyst_NSw_001[iAss]->Get("TotalSystematicSourcesNSSigma");
      c_0110 = (TCanvas*)totalsyst_NSw_0110[iAss]->Get("TotalSystematicSourcesNSSigma");
      c_1030 = (TCanvas*)totalsyst_NSw_1030[iAss]->Get("TotalSystematicSourcesNSSigma");
      c_30100 = (TCanvas*)totalsyst_NSw_30100[iAss]->Get("TotalSystematicSourcesNSSigma");      
      gr_fitvar_001 = (TGraphAsymmErrors*)c_001->FindObject("grSystBasNSSigma");
      gr_fitvar_0110 = (TGraphAsymmErrors*)c_0110->FindObject("grSystBasNSSigma");
      gr_fitvar_1030 = (TGraphAsymmErrors*)c_1030->FindObject("grSystBasNSSigma");
      gr_fitvar_30100 = (TGraphAsymmErrors*)c_30100->FindObject("grSystBasNSSigma");
    }
    fitvar_001 = gr_fitvar_001->GetEYhigh()[iD]; //picks the upward (it's symm) baseline fit variations from point iD, which is the relative uncertainty
    fitvar_0110 = gr_fitvar_0110->GetEYhigh()[iD]; 
    fitvar_1030 = gr_fitvar_1030->GetEYhigh()[iD]; 
    fitvar_30100 = gr_fitvar_30100->GetEYhigh()[iD]; 

    //syst group 1: sum in quadrature of contributions
    if(ratiopanel > 3) { //For NS widths PUT THE first two to zero!!!
      yield_001 = 0; yield_0110 = 0; yield_1030 = 0; yield_30100 = 0;
      bkgshape_001 = 0; bkgshape_0110 = 0; bkgshape_1030 = 0; bkgshape_30100 = 0;
    }

    sum_GROUP1_001 = TMath::Sqrt(yield_001*yield_001 + bkgshape_001*bkgshape_001 + fitvar_001*fitvar_001);
    sum_GROUP1_0110 = TMath::Sqrt(yield_0110*yield_0110 + bkgshape_0110*bkgshape_0110 + fitvar_0110*fitvar_0110);
    sum_GROUP1_1030 = TMath::Sqrt(yield_1030*yield_1030 + bkgshape_1030*bkgshape_1030 + fitvar_1030*fitvar_1030);
    sum_GROUP1_30100 = TMath::Sqrt(yield_30100*yield_30100 + bkgshape_30100*bkgshape_30100 + fitvar_30100*fitvar_30100);

    rel_GROUP1_Ratio1 = TMath::Sqrt(sum_GROUP1_0110*sum_GROUP1_0110 + sum_GROUP1_001*sum_GROUP1_001);
    rel_GROUP1_Ratio2 = TMath::Sqrt(sum_GROUP1_1030*sum_GROUP1_1030 + sum_GROUP1_001*sum_GROUP1_001);
    rel_GROUP1_Ratio3 = TMath::Sqrt(sum_GROUP1_30100*sum_GROUP1_30100 + sum_GROUP1_001*sum_GROUP1_001);

    printf("001   - Yield %.2f, Bkgshape %.2f, FitVar %.2f --> Sum %.2f\n",yield_001,bkgshape_001,fitvar_001,sum_GROUP1_001);
    printf("0110  - Yield %.2f, Bkgshape %.2f, FitVar %.2f --> Sum %.2f\n",yield_0110,bkgshape_0110,fitvar_0110,sum_GROUP1_0110);
    printf("1030  - Yield %.2f, Bkgshape %.2f, FitVar %.2f --> Sum %.2f\n",yield_1030,bkgshape_1030,fitvar_1030,sum_GROUP1_1030);
    printf("30100 - Yield %.2f, Bkgshape %.2f, FitVar %.2f --> Sum %.2f\n",yield_30100,bkgshape_30100,fitvar_30100,sum_GROUP1_30100);
    printf("Ratio 0110/001  - Group1: %.2f\n",rel_GROUP1_Ratio1);
    printf("Ratio 1030/001  - Group1: %.2f\n",rel_GROUP1_Ratio2);
    printf("Ratio 30100/001 - Group1: %.2f\n",rel_GROUP1_Ratio3);

    //syst group 2: extraction
    if(ratiopanel < 4) { //NS yield
      c_001 = (TCanvas*)totalsyst_NSy_001[iAss]->Get("TotalSystematicSourcesNSYield");
      c_0110 = (TCanvas*)totalsyst_NSy_0110[iAss]->Get("TotalSystematicSourcesNSYield");
      c_1030 = (TCanvas*)totalsyst_NSy_1030[iAss]->Get("TotalSystematicSourcesNSYield");
      c_30100 = (TCanvas*)totalsyst_NSy_30100[iAss]->Get("TotalSystematicSourcesNSYield");
      gr_v2var_001 = (TGraphAsymmErrors*)c_001->FindObject("grv2NSYield");
      gr_v2var_0110 = (TGraphAsymmErrors*)c_0110->FindObject("grv2NSYield");
      gr_v2var_1030 = (TGraphAsymmErrors*)c_1030->FindObject("grv2NSYield");
      gr_v2var_30100 = (TGraphAsymmErrors*)c_30100->FindObject("grv2NSYield");
    } else { //NS width
      c_001 = (TCanvas*)totalsyst_NSw_001[iAss]->Get("TotalSystematicSourcesNSSigma");
      c_0110 = (TCanvas*)totalsyst_NSw_0110[iAss]->Get("TotalSystematicSourcesNSSigma");
      c_1030 = (TCanvas*)totalsyst_NSw_1030[iAss]->Get("TotalSystematicSourcesNSSigma");
      c_30100 = (TCanvas*)totalsyst_NSw_30100[iAss]->Get("TotalSystematicSourcesNSSigma");      
      gr_v2var_001 = (TGraphAsymmErrors*)c_001->FindObject("grv2NSSigma");
      gr_v2var_0110 = (TGraphAsymmErrors*)c_0110->FindObject("grv2NSSigma");
      gr_v2var_1030 = (TGraphAsymmErrors*)c_1030->FindObject("grv2NSSigma");
      gr_v2var_30100 = (TGraphAsymmErrors*)c_30100->FindObject("grv2NSSigma");
    }    

    v2relmin_001 = gr_v2var_001->GetEYlow()[iD]*(-1); //picks the downward baseline fit variations from point iD, and makes it negative
    v2relmax_001 = gr_v2var_001->GetEYhigh()[iD]; //picks the upward baseline fit variations from point iD
    v2relmin_0110 = gr_v2var_0110->GetEYlow()[iD]*(-1);
    v2relmax_0110 = gr_v2var_0110->GetEYhigh()[iD];
    v2relmin_1030 = gr_v2var_1030->GetEYlow()[iD]*(-1);
    v2relmax_1030 = gr_v2var_1030->GetEYhigh()[iD];
    v2relmin_30100 = gr_v2var_30100->GetEYlow()[iD]*(-1);
    v2relmax_30100 = gr_v2var_30100->GetEYhigh()[iD];

    v2absmin_001 = (v2relmin_001 + 1)*y001[iD]; //ok! it's necessary to add the +1: you have sth like -0.10, +0.07, etc. and you want the corresponding point position
    v2absmax_001 = (v2relmax_001 + 1)*y001[iD];
    v2absmin_0110 = (v2relmin_0110 + 1)*y0110[iD];
    v2absmax_0110 = (v2relmax_0110 + 1)*y0110[iD];
    v2absmin_1030 = (v2relmin_1030 + 1)*y1030[iD];
    v2absmax_1030 = (v2relmax_1030 + 1)*y1030[iD];
    v2absmin_30100 = (v2relmin_30100 + 1)*y30100[iD];
    v2absmax_30100 = (v2relmax_30100 + 1)*y30100[iD];        

    v2absmin_Ratio1 = v2absmin_0110/v2absmin_001;
    v2absmax_Ratio1 = v2absmax_0110/v2absmax_001;
    v2absmin_Ratio2 = v2absmin_1030/v2absmin_001;
    v2absmax_Ratio2 = v2absmax_1030/v2absmax_001;
    v2absmin_Ratio3 = v2absmin_30100/v2absmin_001;
    v2absmax_Ratio3 = v2absmax_30100/v2absmax_001;

    v2_ENVLOW_abs_Ratio1 = TMath::Min(v2absmin_Ratio1,TMath::Min(v2absmax_Ratio1,ratio1)); //min between std ratio, v2-min ratio, v2-max ratio (one of the two v2 coincides with std, but the direction can change, so this is more generic)
    v2_ENVHIG_abs_Ratio1 = TMath::Max(v2absmin_Ratio1,TMath::Max(v2absmax_Ratio1,ratio1));
    v2_ENVLOW_abs_Ratio2 = TMath::Min(v2absmin_Ratio2,TMath::Min(v2absmax_Ratio2,ratio2));
    v2_ENVHIG_abs_Ratio2 = TMath::Max(v2absmin_Ratio2,TMath::Max(v2absmax_Ratio2,ratio2));
    v2_ENVLOW_abs_Ratio3 = TMath::Min(v2absmin_Ratio3,TMath::Min(v2absmax_Ratio3,ratio3));
    v2_ENVHIG_abs_Ratio3 = TMath::Max(v2absmin_Ratio3,TMath::Max(v2absmax_Ratio3,ratio3));

    v2_ENVLOW_rel_Ratio1 = TMath::Abs(v2_ENVLOW_abs_Ratio1/ratio1 - 1); //back from absolute position of point to |relative difference| (note: the sign is positive, the direction is already in the name LOW/HIG)
    v2_ENVHIG_rel_Ratio1 = TMath::Abs(v2_ENVHIG_abs_Ratio1/ratio1 - 1);
    v2_ENVLOW_rel_Ratio2 = TMath::Abs(v2_ENVLOW_abs_Ratio2/ratio2 - 1);
    v2_ENVHIG_rel_Ratio2 = TMath::Abs(v2_ENVHIG_abs_Ratio2/ratio2 - 1);
    v2_ENVLOW_rel_Ratio3 = TMath::Abs(v2_ENVLOW_abs_Ratio3/ratio3 - 1);
    v2_ENVHIG_rel_Ratio3 = TMath::Abs(v2_ENVHIG_abs_Ratio3/ratio3 - 1);

    printf("001   - v2min rel %.2f, abs %.2f - v2max rel %.2f, abs %.2f\n",v2relmin_001,v2absmin_001,v2relmax_001,v2absmax_001);
    printf("0110  - v2min rel %.2f, abs %.2f - v2max rel %.2f, abs %.2f\n",v2relmin_0110,v2absmin_0110,v2relmax_0110,v2absmax_0110);
    printf("1030  - v2min rel %.2f, abs %.2f - v2max rel %.2f, abs %.2f\n",v2relmin_1030,v2absmin_1030,v2relmax_1030,v2absmax_1030);
    printf("30100 - v2min rel %.2f, abs %.2f - v2max rel %.2f, abs %.2f\n",v2relmin_30100,v2absmin_30100,v2relmax_30100,v2absmax_30100);
    printf("Ratio 0110/001  - v2min/max: %.2f, %.2f\n",v2absmin_Ratio1,v2absmax_Ratio1);
    printf("Ratio 1030/001  - v2min/max: %.2f, %.2f\n",v2absmin_Ratio2,v2absmax_Ratio2);
    printf("Ratio 30100/001 - v2min/max: %.2f, %.2f\n",v2absmin_Ratio3,v2absmax_Ratio3);
    printf("Ratio 0110/001  - ENV low/high (abs): %.2f, %.2f, (rel) %.2f, %.2f\n",v2_ENVLOW_abs_Ratio1,v2_ENVHIG_abs_Ratio1,v2_ENVLOW_rel_Ratio1,v2_ENVHIG_rel_Ratio1);
    printf("Ratio 1030/001  - ENV low/high (abs): %.2f, %.2f, (rel) %.2f, %.2f\n",v2_ENVLOW_abs_Ratio2,v2_ENVHIG_abs_Ratio2,v2_ENVLOW_rel_Ratio2,v2_ENVHIG_rel_Ratio2);
    printf("Ratio 30100/001 - ENV low/high (abs): %.2f, %.2f, (rel) %.2f, %.2f\n",v2_ENVLOW_abs_Ratio3,v2_ENVHIG_abs_Ratio3,v2_ENVLOW_rel_Ratio3,v2_ENVHIG_rel_Ratio3);

    //syst group 3: extraction
    if(ratiopanel < 4) { //NS yield
      hist_FDmin_BDflat_001 = (TH1D*)syst_BDflat_001[iAss]->Get("ValueHistoNSYield_10");
      hist_FDmax_BDflat_001 = (TH1D*)syst_BDflat_001[iAss]->Get("ValueHistoNSYield_11");
      hist_FDmin_BDextr_001 = (TH1D*)syst_BDextr_001[iAss]->Get("ValueHistoNSYield_10");
      hist_FDmax_BDextr_001 = (TH1D*)syst_BDextr_001[iAss]->Get("ValueHistoNSYield_11");
      hist_FDmin_BDflat_0110 = (TH1D*)syst_BDflat_0110[iAss]->Get("ValueHistoNSYield_10");
      hist_FDmax_BDflat_0110 = (TH1D*)syst_BDflat_0110[iAss]->Get("ValueHistoNSYield_11");
      hist_FDmin_BDextr_0110 = (TH1D*)syst_BDextr_0110[iAss]->Get("ValueHistoNSYield_10");
      hist_FDmax_BDextr_0110 = (TH1D*)syst_BDextr_0110[iAss]->Get("ValueHistoNSYield_11");
      hist_FDmin_BDflat_1030 = (TH1D*)syst_BDflat_1030[iAss]->Get("ValueHistoNSYield_10");
      hist_FDmax_BDflat_1030 = (TH1D*)syst_BDflat_1030[iAss]->Get("ValueHistoNSYield_11");
      hist_FDmin_BDextr_1030 = (TH1D*)syst_BDextr_1030[iAss]->Get("ValueHistoNSYield_10");
      hist_FDmax_BDextr_1030 = (TH1D*)syst_BDextr_1030[iAss]->Get("ValueHistoNSYield_11");
      hist_FDmin_BDflat_30100 = (TH1D*)syst_BDflat_30100[iAss]->Get("ValueHistoNSYield_10");
      hist_FDmax_BDflat_30100 = (TH1D*)syst_BDflat_30100[iAss]->Get("ValueHistoNSYield_11");
      hist_FDmin_BDextr_30100 = (TH1D*)syst_BDextr_30100[iAss]->Get("ValueHistoNSYield_10");
      hist_FDmax_BDextr_30100 = (TH1D*)syst_BDextr_30100[iAss]->Get("ValueHistoNSYield_11");                  
    } else { //NS width
      hist_FDmin_BDflat_001 = (TH1D*)syst_BDflat_001[iAss]->Get("ValueHistoNSSigma_10");
      hist_FDmax_BDflat_001 = (TH1D*)syst_BDflat_001[iAss]->Get("ValueHistoNSSigma_11");
      hist_FDmin_BDextr_001 = (TH1D*)syst_BDextr_001[iAss]->Get("ValueHistoNSSigma_10");
      hist_FDmax_BDextr_001 = (TH1D*)syst_BDextr_001[iAss]->Get("ValueHistoNSSigma_11");
      hist_FDmin_BDflat_0110 = (TH1D*)syst_BDflat_0110[iAss]->Get("ValueHistoNSSigma_10");
      hist_FDmax_BDflat_0110 = (TH1D*)syst_BDflat_0110[iAss]->Get("ValueHistoNSSigma_11");
      hist_FDmin_BDextr_0110 = (TH1D*)syst_BDextr_0110[iAss]->Get("ValueHistoNSSigma_10");
      hist_FDmax_BDextr_0110 = (TH1D*)syst_BDextr_0110[iAss]->Get("ValueHistoNSSigma_11");
      hist_FDmin_BDflat_1030 = (TH1D*)syst_BDflat_1030[iAss]->Get("ValueHistoNSSigma_10");
      hist_FDmax_BDflat_1030 = (TH1D*)syst_BDflat_1030[iAss]->Get("ValueHistoNSSigma_11");
      hist_FDmin_BDextr_1030 = (TH1D*)syst_BDextr_1030[iAss]->Get("ValueHistoNSSigma_10");
      hist_FDmax_BDextr_1030 = (TH1D*)syst_BDextr_1030[iAss]->Get("ValueHistoNSSigma_11");
      hist_FDmin_BDflat_30100 = (TH1D*)syst_BDflat_30100[iAss]->Get("ValueHistoNSSigma_10");
      hist_FDmax_BDflat_30100 = (TH1D*)syst_BDflat_30100[iAss]->Get("ValueHistoNSSigma_11");
      hist_FDmin_BDextr_30100 = (TH1D*)syst_BDextr_30100[iAss]->Get("ValueHistoNSSigma_10");
      hist_FDmax_BDextr_30100 = (TH1D*)syst_BDextr_30100[iAss]->Get("ValueHistoNSSigma_11"); 
    }  

    FDmin_BDflat_abs_001 = hist_FDmin_BDflat_001->GetBinContent(iD+1); //loop starts from 0, histo bins index from 1
    FDmax_BDflat_abs_001 = hist_FDmax_BDflat_001->GetBinContent(iD+1);
    FDmin_BDextr_abs_001 = hist_FDmin_BDextr_001->GetBinContent(iD+1);
    FDmax_BDextr_abs_001 = hist_FDmax_BDextr_001->GetBinContent(iD+1);
    FDmin_BDflat_abs_0110 = hist_FDmin_BDflat_0110->GetBinContent(iD+1);
    FDmax_BDflat_abs_0110 = hist_FDmax_BDflat_0110->GetBinContent(iD+1);
    FDmin_BDextr_abs_0110 = hist_FDmin_BDextr_0110->GetBinContent(iD+1);
    FDmax_BDextr_abs_0110 = hist_FDmax_BDextr_0110->GetBinContent(iD+1);
    FDmin_BDflat_abs_1030 = hist_FDmin_BDflat_1030->GetBinContent(iD+1);
    FDmax_BDflat_abs_1030 = hist_FDmax_BDflat_1030->GetBinContent(iD+1);
    FDmin_BDextr_abs_1030 = hist_FDmin_BDextr_1030->GetBinContent(iD+1);
    FDmax_BDextr_abs_1030 = hist_FDmax_BDextr_1030->GetBinContent(iD+1);
    FDmin_BDflat_abs_30100 = hist_FDmin_BDflat_30100->GetBinContent(iD+1);
    FDmax_BDflat_abs_30100 = hist_FDmax_BDflat_30100->GetBinContent(iD+1);
    FDmin_BDextr_abs_30100 = hist_FDmin_BDextr_30100->GetBinContent(iD+1);
    FDmax_BDextr_abs_30100 = hist_FDmax_BDextr_30100->GetBinContent(iD+1);

    FDmin_BDflat_abs_Ratio1 = FDmin_BDflat_abs_0110/FDmin_BDflat_abs_001;
    FDmax_BDflat_abs_Ratio1 = FDmax_BDflat_abs_0110/FDmax_BDflat_abs_001;
    FDmin_BDextr_abs_Ratio1 = FDmin_BDextr_abs_0110/FDmin_BDextr_abs_001;
    FDmax_BDextr_abs_Ratio1 = FDmax_BDextr_abs_0110/FDmax_BDextr_abs_001;
    FDmin_BDflat_abs_Ratio2 = FDmin_BDflat_abs_1030/FDmin_BDflat_abs_001;
    FDmax_BDflat_abs_Ratio2 = FDmax_BDflat_abs_1030/FDmax_BDflat_abs_001;
    FDmin_BDextr_abs_Ratio2 = FDmin_BDextr_abs_1030/FDmin_BDextr_abs_001;
    FDmax_BDextr_abs_Ratio2 = FDmax_BDextr_abs_1030/FDmax_BDextr_abs_001;
    FDmin_BDflat_abs_Ratio3 = FDmin_BDflat_abs_30100/FDmin_BDflat_abs_001;
    FDmax_BDflat_abs_Ratio3 = FDmax_BDflat_abs_30100/FDmax_BDflat_abs_001;
    FDmin_BDextr_abs_Ratio3 = FDmin_BDextr_abs_30100/FDmin_BDextr_abs_001;
    FDmax_BDextr_abs_Ratio3 = FDmax_BDextr_abs_30100/FDmax_BDextr_abs_001;    

    FD_ENVLOW_abs_Ratio1 = TMath::Min(FDmin_BDflat_abs_Ratio1,TMath::Min(FDmax_BDflat_abs_Ratio1,TMath::Min(FDmin_BDextr_abs_Ratio1,TMath::Min(FDmax_BDextr_abs_Ratio1,ratio1)))); //min between std ratio, FD-min and FD-max ratio each for BDflat and BDextr
    FD_ENVHIG_abs_Ratio1 = TMath::Max(FDmin_BDflat_abs_Ratio1,TMath::Max(FDmax_BDflat_abs_Ratio1,TMath::Max(FDmin_BDextr_abs_Ratio1,TMath::Max(FDmax_BDextr_abs_Ratio1,ratio1)))); //min between std ratio, FD-min and FD-max ratio each for BDflat and BDextr
    FD_ENVLOW_abs_Ratio2 = TMath::Min(FDmin_BDflat_abs_Ratio2,TMath::Min(FDmax_BDflat_abs_Ratio2,TMath::Min(FDmin_BDextr_abs_Ratio2,TMath::Min(FDmax_BDextr_abs_Ratio2,ratio2))));
    FD_ENVHIG_abs_Ratio2 = TMath::Max(FDmin_BDflat_abs_Ratio2,TMath::Max(FDmax_BDflat_abs_Ratio2,TMath::Max(FDmin_BDextr_abs_Ratio2,TMath::Max(FDmax_BDextr_abs_Ratio2,ratio2))));
    FD_ENVLOW_abs_Ratio3 = TMath::Min(FDmin_BDflat_abs_Ratio3,TMath::Min(FDmax_BDflat_abs_Ratio3,TMath::Min(FDmin_BDextr_abs_Ratio3,TMath::Min(FDmax_BDextr_abs_Ratio3,ratio3))));
    FD_ENVHIG_abs_Ratio3 = TMath::Max(FDmin_BDflat_abs_Ratio3,TMath::Max(FDmax_BDflat_abs_Ratio3,TMath::Max(FDmin_BDextr_abs_Ratio3,TMath::Max(FDmax_BDextr_abs_Ratio3,ratio3))));

    FD_ENVLOW_rel_Ratio1 = TMath::Abs(FD_ENVLOW_abs_Ratio1/ratio1 - 1); //go from absolute position of point to |relative difference| (note: the sign is positive, the direction is already in the name LOW/HIG)
    FD_ENVHIG_rel_Ratio1 = TMath::Abs(FD_ENVHIG_abs_Ratio1/ratio1 - 1);
    FD_ENVLOW_rel_Ratio2 = TMath::Abs(FD_ENVLOW_abs_Ratio2/ratio2 - 1);
    FD_ENVHIG_rel_Ratio2 = TMath::Abs(FD_ENVHIG_abs_Ratio2/ratio2 - 1);
    FD_ENVLOW_rel_Ratio3 = TMath::Abs(FD_ENVLOW_abs_Ratio3/ratio3 - 1);
    FD_ENVHIG_rel_Ratio3 = TMath::Abs(FD_ENVHIG_abs_Ratio3/ratio3 - 1);

    printf("001   - FDmin_BDflat %.2f, FDmax_BDflat %.2f - FDmin_BDextr %.2f, FDmax_BDextr %.2f\n",FDmin_BDflat_abs_001,FDmax_BDflat_abs_001,FDmin_BDextr_abs_001,FDmax_BDextr_abs_001);
    printf("0110  - FDmin_BDflat %.2f, FDmax_BDflat %.2f - FDmin_BDextr %.2f, FDmax_BDextr %.2f\n",FDmin_BDflat_abs_0110,FDmax_BDflat_abs_0110,FDmin_BDextr_abs_0110,FDmax_BDextr_abs_0110);
    printf("1030  - FDmin_BDflat %.2f, FDmax_BDflat %.2f - FDmin_BDextr %.2f, FDmax_BDextr %.2f\n",FDmin_BDflat_abs_1030,FDmax_BDflat_abs_1030,FDmin_BDextr_abs_1030,FDmax_BDextr_abs_1030);
    printf("30100 - FDmin_BDflat %.2f, FDmax_BDflat %.2f - FDmin_BDextr %.2f, FDmax_BDextr %.2f\n",FDmin_BDflat_abs_30100,FDmax_BDflat_abs_30100,FDmin_BDextr_abs_30100,FDmax_BDextr_abs_30100);
    printf("Ratio 0110/001  - FDmin_BDflat %.2f, FDmax_BDflat %.2f - FDmin_BDextr %.2f, FDmax_BDextr %.2f\n",FDmin_BDflat_abs_Ratio1,FDmax_BDflat_abs_Ratio1,FDmin_BDextr_abs_Ratio1,FDmax_BDextr_abs_Ratio1);
    printf("Ratio 1030/001  - FDmin_BDflat %.2f, FDmax_BDflat %.2f - FDmin_BDextr %.2f, FDmax_BDextr %.2f\n",FDmin_BDflat_abs_Ratio2,FDmax_BDflat_abs_Ratio2,FDmin_BDextr_abs_Ratio2,FDmax_BDextr_abs_Ratio2);
    printf("Ratio 30100/001 - FDmin_BDflat %.2f, FDmax_BDflat %.2f - FDmin_BDextr %.2f, FDmax_BDextr %.2f\n",FDmin_BDflat_abs_Ratio3,FDmax_BDflat_abs_Ratio3,FDmin_BDextr_abs_Ratio3,FDmax_BDextr_abs_Ratio3);
    printf("Ratio 0110/001  - ENV low/high (abs): %.2f, %.2f, (rel) %.2f, %.2f\n",FD_ENVLOW_abs_Ratio1,FD_ENVHIG_abs_Ratio1,FD_ENVLOW_rel_Ratio1,FD_ENVHIG_rel_Ratio1);
    printf("Ratio 1030/001  - ENV low/high (abs): %.2f, %.2f, (rel) %.2f, %.2f\n",FD_ENVLOW_abs_Ratio2,FD_ENVHIG_abs_Ratio2,FD_ENVLOW_rel_Ratio2,FD_ENVHIG_rel_Ratio2);
    printf("Ratio 30100/001 - ENV low/high (abs): %.2f, %.2f, (rel) %.2f, %.2f\n",FD_ENVLOW_abs_Ratio3,FD_ENVHIG_abs_Ratio3,FD_ENVLOW_rel_Ratio3,FD_ENVHIG_rel_Ratio3);

    
    //Sum up the systematic uncertainties!! (Simple sum in quadrature of the three components, up and down separately) -> go back to absolute errors, to insert in tgraph bars
    ratiosystUp1 = ratio1*TMath::Sqrt(rel_GROUP1_Ratio1*rel_GROUP1_Ratio1 + v2_ENVHIG_rel_Ratio1*v2_ENVHIG_rel_Ratio1 + FD_ENVHIG_rel_Ratio1*FD_ENVHIG_rel_Ratio1); //absolute errors, to be put in tgraph bars
    ratiosystDown1 = ratio1*TMath::Sqrt(rel_GROUP1_Ratio1*rel_GROUP1_Ratio1 + v2_ENVLOW_rel_Ratio1*v2_ENVLOW_rel_Ratio1 + FD_ENVLOW_rel_Ratio1*FD_ENVLOW_rel_Ratio1);
    ratiosystUp2 = ratio2*TMath::Sqrt(rel_GROUP1_Ratio2*rel_GROUP1_Ratio2 + v2_ENVHIG_rel_Ratio2*v2_ENVHIG_rel_Ratio2 + FD_ENVHIG_rel_Ratio2*FD_ENVHIG_rel_Ratio2);
    ratiosystDown2 = ratio2*TMath::Sqrt(rel_GROUP1_Ratio2*rel_GROUP1_Ratio2 + v2_ENVLOW_rel_Ratio2*v2_ENVLOW_rel_Ratio2 + FD_ENVLOW_rel_Ratio2*FD_ENVLOW_rel_Ratio2);
    ratiosystUp3 = ratio3*TMath::Sqrt(rel_GROUP1_Ratio3*rel_GROUP1_Ratio3 + v2_ENVHIG_rel_Ratio3*v2_ENVHIG_rel_Ratio3 + FD_ENVHIG_rel_Ratio3*FD_ENVHIG_rel_Ratio3);
    ratiosystDown3 = ratio3*TMath::Sqrt(rel_GROUP1_Ratio3*rel_GROUP1_Ratio3 + v2_ENVLOW_rel_Ratio3*v2_ENVLOW_rel_Ratio3 + FD_ENVLOW_rel_Ratio3*FD_ENVLOW_rel_Ratio3);
    printf("FINAL SYST: Ratio 0110/001  - Syst unc. low/high (rel) %.2f, %.2f (abs) %.2f, %.2f (ingr. Gr1, V2, FD low/high %.2f, %.2f, %.2f / %.2f, %.2f, %.2f)\n",ratiosystDown1/ratio1,ratiosystUp1/ratio1,ratiosystDown1,ratiosystUp1,rel_GROUP1_Ratio1,v2_ENVLOW_rel_Ratio1,FD_ENVLOW_rel_Ratio1,rel_GROUP1_Ratio1,v2_ENVHIG_rel_Ratio1,FD_ENVHIG_rel_Ratio1);
    printf("FINAL SYST: Ratio 1030/001  - Syst unc. low/high (rel) %.2f, %.2f (abs) %.2f, %.2f (ingr. Gr1, V2, FD low/high %.2f, %.2f, %.2f / %.2f, %.2f, %.2f)\n",ratiosystDown2/ratio2,ratiosystUp2/ratio2,ratiosystDown2,ratiosystUp2,rel_GROUP1_Ratio2,v2_ENVLOW_rel_Ratio2,FD_ENVLOW_rel_Ratio2,rel_GROUP1_Ratio2,v2_ENVHIG_rel_Ratio2,FD_ENVHIG_rel_Ratio2);
    printf("FINAL SYST: Ratio 30100/001 - Syst unc. low/high (rel) %.2f, %.2f (abs) %.2f, %.2f (ingr. Gr1, V2, FD low/high %.2f, %.2f, %.2f / %.2f, %.2f, %.2f)\n",ratiosystDown3/ratio3,ratiosystUp3/ratio3,ratiosystDown3,ratiosystUp3,rel_GROUP1_Ratio3,v2_ENVLOW_rel_Ratio3,FD_ENVLOW_rel_Ratio3,rel_GROUP1_Ratio3,v2_ENVHIG_rel_Ratio3,FD_ENVHIG_rel_Ratio3);

    //Finally, insert the values in the TGraphs!
    grRatio1->SetPoint(iD,x0110[iD],ratio1);   //i, x, y
    grRatio1->SetPointError(iD,ex[iD],ex[iD],ratiostat1,ratiostat1);  //i, exl, exh, eyl, eyh
    grRatiosyst1->SetPoint(iD,x0110[iD],ratio1);   //i, x, y
    grRatiosyst1->SetPointError(iD,0.65,0.65,ratiosystDown1,ratiosystUp1);  //i, exl, exh, eyl, eyh
    printf("ratio1 bin %d (x %f) = %f, stat %f, systl %f, systh %f, all absolute\n",iD,x0110[iD],ratio1,ratiostat1,ratiosystDown1,ratiosystUp1);

    grRatio2->SetPoint(iD,x1030[iD],ratio2);   //i, x, y
    grRatio2->SetPointError(iD,ex[iD],ex[iD],ratiostat2,ratiostat2);  //i, exl, exh, eyl, eyh
    grRatiosyst2->SetPoint(iD,x1030[iD],ratio2);   //i, x, y
    grRatiosyst2->SetPointError(iD,0.65,0.65,ratiosystDown2,ratiosystUp2);  //i, exl, exh, eyl, eyh
    printf("ratio2 bin %d (x %f) = %f, stat %f, systl %f, systh %f, all absolute\n",iD,x1030[iD],ratio2,ratiostat2,ratiosystDown2,ratiosystUp2);

    grRatio3->SetPoint(iD,x30100[iD],ratio3);   //i, x, y
    grRatio3->SetPointError(iD,ex[iD],ex[iD],ratiostat3,ratiostat3);  //i, exl, exh, eyl, eyh
    grRatiosyst3->SetPoint(iD,x30100[iD],ratio3);   //i, x, y
    grRatiosyst3->SetPointError(iD,0.65,0.65,ratiosystDown3,ratiosystUp3);  //i, exl, exh, eyl, eyh
    printf("ratio3 bin %d (x %f) = %f, stat %f, systl %f, systh %f, all absolute\n",iD,x30100[iD],ratio3,ratiostat3,ratiosystDown3,ratiosystUp3);
  }

 return;
}

void LoadFilesForSystPropagation() {

 syst_std_001[0] = TFile::Open("001/FitOutput_pp/FitSystematics_pp_ArithmeticAverage_0.3_99.0.root");
 syst_std_001[1] = TFile::Open("001/FitOutput_pp/FitSystematics_pp_ArithmeticAverage_0.3_1.0.root");
 syst_std_001[2] = TFile::Open("001/FitOutput_pp/FitSystematics_pp_ArithmeticAverage_1.0_2.0.root");
 syst_std_001[3] = TFile::Open("001/FitOutput_pp/FitSystematics_pp_ArithmeticAverage_2.0_3.0.root");
 syst_BDflat_001[0] = TFile::Open("001/FitOutput_pp_ForSyst_BDflat/FitSystematics_pp_ArithmeticAverage_0.3_99.0.root");
 syst_BDflat_001[1] = TFile::Open("001/FitOutput_pp_ForSyst_BDflat/FitSystematics_pp_ArithmeticAverage_0.3_1.0.root");
 syst_BDflat_001[2] = TFile::Open("001/FitOutput_pp_ForSyst_BDflat/FitSystematics_pp_ArithmeticAverage_1.0_2.0.root");
 syst_BDflat_001[3] = TFile::Open("001/FitOutput_pp_ForSyst_BDflat/FitSystematics_pp_ArithmeticAverage_2.0_3.0.root"); 
 syst_BDextr_001[0] = TFile::Open("001/FitOutput_pp_ForSyst_BDextr/FitSystematics_pp_ArithmeticAverage_0.3_99.0.root");
 syst_BDextr_001[1] = TFile::Open("001/FitOutput_pp_ForSyst_BDextr/FitSystematics_pp_ArithmeticAverage_0.3_1.0.root");
 syst_BDextr_001[2] = TFile::Open("001/FitOutput_pp_ForSyst_BDextr/FitSystematics_pp_ArithmeticAverage_1.0_2.0.root");
 syst_BDextr_001[3] = TFile::Open("001/FitOutput_pp_ForSyst_BDextr/FitSystematics_pp_ArithmeticAverage_2.0_3.0.root");  

 syst_std_0110[0] = TFile::Open("0110/FitOutput_pp/FitSystematics_pp_ArithmeticAverage_0.3_99.0.root");
 syst_std_0110[1] = TFile::Open("0110/FitOutput_pp/FitSystematics_pp_ArithmeticAverage_0.3_1.0.root");
 syst_std_0110[2] = TFile::Open("0110/FitOutput_pp/FitSystematics_pp_ArithmeticAverage_1.0_2.0.root");
 syst_std_0110[3] = TFile::Open("0110/FitOutput_pp/FitSystematics_pp_ArithmeticAverage_2.0_3.0.root");
 syst_BDflat_0110[0] = TFile::Open("0110/FitOutput_pp_ForSyst_BDflat/FitSystematics_pp_ArithmeticAverage_0.3_99.0.root");
 syst_BDflat_0110[1] = TFile::Open("0110/FitOutput_pp_ForSyst_BDflat/FitSystematics_pp_ArithmeticAverage_0.3_1.0.root");
 syst_BDflat_0110[2] = TFile::Open("0110/FitOutput_pp_ForSyst_BDflat/FitSystematics_pp_ArithmeticAverage_1.0_2.0.root");
 syst_BDflat_0110[3] = TFile::Open("0110/FitOutput_pp_ForSyst_BDflat/FitSystematics_pp_ArithmeticAverage_2.0_3.0.root"); 
 syst_BDextr_0110[0] = TFile::Open("0110/FitOutput_pp_ForSyst_BDextr/FitSystematics_pp_ArithmeticAverage_0.3_99.0.root");
 syst_BDextr_0110[1] = TFile::Open("0110/FitOutput_pp_ForSyst_BDextr/FitSystematics_pp_ArithmeticAverage_0.3_1.0.root");
 syst_BDextr_0110[2] = TFile::Open("0110/FitOutput_pp_ForSyst_BDextr/FitSystematics_pp_ArithmeticAverage_1.0_2.0.root");
 syst_BDextr_0110[3] = TFile::Open("0110/FitOutput_pp_ForSyst_BDextr/FitSystematics_pp_ArithmeticAverage_2.0_3.0.root");  

 syst_std_1030[0] = TFile::Open("1030/FitOutput_pp/FitSystematics_pp_ArithmeticAverage_0.3_99.0.root");
 syst_std_1030[1] = TFile::Open("1030/FitOutput_pp/FitSystematics_pp_ArithmeticAverage_0.3_1.0.root");
 syst_std_1030[2] = TFile::Open("1030/FitOutput_pp/FitSystematics_pp_ArithmeticAverage_1.0_2.0.root");
 syst_std_1030[3] = TFile::Open("1030/FitOutput_pp/FitSystematics_pp_ArithmeticAverage_2.0_3.0.root");
 syst_BDflat_1030[0] = TFile::Open("1030/FitOutput_pp_ForSyst_BDflat/FitSystematics_pp_ArithmeticAverage_0.3_99.0.root");
 syst_BDflat_1030[1] = TFile::Open("1030/FitOutput_pp_ForSyst_BDflat/FitSystematics_pp_ArithmeticAverage_0.3_1.0.root");
 syst_BDflat_1030[2] = TFile::Open("1030/FitOutput_pp_ForSyst_BDflat/FitSystematics_pp_ArithmeticAverage_1.0_2.0.root");
 syst_BDflat_1030[3] = TFile::Open("1030/FitOutput_pp_ForSyst_BDflat/FitSystematics_pp_ArithmeticAverage_2.0_3.0.root"); 
 syst_BDextr_1030[0] = TFile::Open("1030/FitOutput_pp_ForSyst_BDextr/FitSystematics_pp_ArithmeticAverage_0.3_99.0.root");
 syst_BDextr_1030[1] = TFile::Open("1030/FitOutput_pp_ForSyst_BDextr/FitSystematics_pp_ArithmeticAverage_0.3_1.0.root");
 syst_BDextr_1030[2] = TFile::Open("1030/FitOutput_pp_ForSyst_BDextr/FitSystematics_pp_ArithmeticAverage_1.0_2.0.root");
 syst_BDextr_1030[3] = TFile::Open("1030/FitOutput_pp_ForSyst_BDextr/FitSystematics_pp_ArithmeticAverage_2.0_3.0.root");  

 syst_std_30100[0] = TFile::Open("30100/FitOutput_pp/FitSystematics_pp_ArithmeticAverage_0.3_99.0.root");
 syst_std_30100[1] = TFile::Open("30100/FitOutput_pp/FitSystematics_pp_ArithmeticAverage_0.3_1.0.root");
 syst_std_30100[2] = TFile::Open("30100/FitOutput_pp/FitSystematics_pp_ArithmeticAverage_1.0_2.0.root");
 syst_std_30100[3] = TFile::Open("30100/FitOutput_pp/FitSystematics_pp_ArithmeticAverage_2.0_3.0.root");
 syst_BDflat_30100[0] = TFile::Open("30100/FitOutput_pp_ForSyst_BDflat/FitSystematics_pp_ArithmeticAverage_0.3_99.0.root");
 syst_BDflat_30100[1] = TFile::Open("30100/FitOutput_pp_ForSyst_BDflat/FitSystematics_pp_ArithmeticAverage_0.3_1.0.root");
 syst_BDflat_30100[2] = TFile::Open("30100/FitOutput_pp_ForSyst_BDflat/FitSystematics_pp_ArithmeticAverage_1.0_2.0.root");
 syst_BDflat_30100[3] = TFile::Open("30100/FitOutput_pp_ForSyst_BDflat/FitSystematics_pp_ArithmeticAverage_2.0_3.0.root"); 
 syst_BDextr_30100[0] = TFile::Open("30100/FitOutput_pp_ForSyst_BDextr/FitSystematics_pp_ArithmeticAverage_0.3_99.0.root");
 syst_BDextr_30100[1] = TFile::Open("30100/FitOutput_pp_ForSyst_BDextr/FitSystematics_pp_ArithmeticAverage_0.3_1.0.root");
 syst_BDextr_30100[2] = TFile::Open("30100/FitOutput_pp_ForSyst_BDextr/FitSystematics_pp_ArithmeticAverage_1.0_2.0.root");
 syst_BDextr_30100[3] = TFile::Open("30100/FitOutput_pp_ForSyst_BDextr/FitSystematics_pp_ArithmeticAverage_2.0_3.0.root");  

 totalsyst_NSy_001[0] = TFile::Open("001/Trends_pp/TotalSystematicSourcesNSYield_pthad0.3to99.0.root");
 totalsyst_NSw_001[0] = TFile::Open("001/Trends_pp/TotalSystematicSourcesNSSigma_pthad0.3to99.0.root");
 totalsyst_NSy_001[1] = TFile::Open("001/Trends_pp/TotalSystematicSourcesNSYield_pthad0.3to1.0.root");
 totalsyst_NSw_001[1] = TFile::Open("001/Trends_pp/TotalSystematicSourcesNSSigma_pthad0.3to1.0.root");
 totalsyst_NSy_001[2] = TFile::Open("001/Trends_pp/TotalSystematicSourcesNSYield_pthad1.0to2.0.root");
 totalsyst_NSw_001[2] = TFile::Open("001/Trends_pp/TotalSystematicSourcesNSSigma_pthad1.0to2.0.root");
 totalsyst_NSy_001[3] = TFile::Open("001/Trends_pp/TotalSystematicSourcesNSYield_pthad2.0to3.0.root");
 totalsyst_NSw_001[3] = TFile::Open("001/Trends_pp/TotalSystematicSourcesNSSigma_pthad2.0to3.0.root");
 
 totalsyst_NSy_0110[0] = TFile::Open("0110/Trends_pp/TotalSystematicSourcesNSYield_pthad0.3to99.0.root");
 totalsyst_NSw_0110[0] = TFile::Open("0110/Trends_pp/TotalSystematicSourcesNSSigma_pthad0.3to99.0.root");
 totalsyst_NSy_0110[1] = TFile::Open("0110/Trends_pp/TotalSystematicSourcesNSYield_pthad0.3to1.0.root");
 totalsyst_NSw_0110[1] = TFile::Open("0110/Trends_pp/TotalSystematicSourcesNSSigma_pthad0.3to1.0.root");
 totalsyst_NSy_0110[2] = TFile::Open("0110/Trends_pp/TotalSystematicSourcesNSYield_pthad1.0to2.0.root");
 totalsyst_NSw_0110[2] = TFile::Open("0110/Trends_pp/TotalSystematicSourcesNSSigma_pthad1.0to2.0.root");
 totalsyst_NSy_0110[3] = TFile::Open("0110/Trends_pp/TotalSystematicSourcesNSYield_pthad2.0to3.0.root");
 totalsyst_NSw_0110[3] = TFile::Open("0110/Trends_pp/TotalSystematicSourcesNSSigma_pthad2.0to3.0.root");  

 totalsyst_NSy_1030[0] = TFile::Open("1030/Trends_pp/TotalSystematicSourcesNSYield_pthad0.3to99.0.root");
 totalsyst_NSw_1030[0] = TFile::Open("1030/Trends_pp/TotalSystematicSourcesNSSigma_pthad0.3to99.0.root");
 totalsyst_NSy_1030[1] = TFile::Open("1030/Trends_pp/TotalSystematicSourcesNSYield_pthad0.3to1.0.root");
 totalsyst_NSw_1030[1] = TFile::Open("1030/Trends_pp/TotalSystematicSourcesNSSigma_pthad0.3to1.0.root");
 totalsyst_NSy_1030[2] = TFile::Open("1030/Trends_pp/TotalSystematicSourcesNSYield_pthad1.0to2.0.root");
 totalsyst_NSw_1030[2] = TFile::Open("1030/Trends_pp/TotalSystematicSourcesNSSigma_pthad1.0to2.0.root");
 totalsyst_NSy_1030[3] = TFile::Open("1030/Trends_pp/TotalSystematicSourcesNSYield_pthad2.0to3.0.root");
 totalsyst_NSw_1030[3] = TFile::Open("1030/Trends_pp/TotalSystematicSourcesNSSigma_pthad2.0to3.0.root");

 totalsyst_NSy_30100[0] = TFile::Open("30100/Trends_pp/TotalSystematicSourcesNSYield_pthad0.3to99.0.root");      
 totalsyst_NSw_30100[0] = TFile::Open("30100/Trends_pp/TotalSystematicSourcesNSSigma_pthad0.3to99.0.root");
 totalsyst_NSy_30100[1] = TFile::Open("30100/Trends_pp/TotalSystematicSourcesNSYield_pthad0.3to1.0.root");      
 totalsyst_NSw_30100[1] = TFile::Open("30100/Trends_pp/TotalSystematicSourcesNSSigma_pthad0.3to1.0.root");
 totalsyst_NSy_30100[2] = TFile::Open("30100/Trends_pp/TotalSystematicSourcesNSYield_pthad1.0to2.0.root");      
 totalsyst_NSw_30100[2] = TFile::Open("30100/Trends_pp/TotalSystematicSourcesNSSigma_pthad1.0to2.0.root");
 totalsyst_NSy_30100[3] = TFile::Open("30100/Trends_pp/TotalSystematicSourcesNSYield_pthad2.0to3.0.root");      
 totalsyst_NSw_30100[3] = TFile::Open("30100/Trends_pp/TotalSystematicSourcesNSSigma_pthad2.0to3.0.root");

 dPhisyst_001[0][0] = TFile::Open("001/ArithmeticAverageppDzeroDstarDplus3to5_assoc0.3to99.0.root");
 dPhisyst_001[0][1] = TFile::Open("001/ArithmeticAverageppDzeroDstarDplus5to8_assoc0.3to99.0.root");
 dPhisyst_001[0][2] = TFile::Open("001/ArithmeticAverageppDzeroDstarDplus8to16_assoc0.3to99.0.root");
 dPhisyst_001[0][3] = TFile::Open("001/ArithmeticAverageppDzeroDstarDplus16to24_assoc0.3to99.0.root");
 dPhisyst_001[1][0] = TFile::Open("001/ArithmeticAverageppDzeroDstarDplus3to5_assoc0.3to1.0.root");
 dPhisyst_001[1][1] = TFile::Open("001/ArithmeticAverageppDzeroDstarDplus5to8_assoc0.3to1.0.root");
 dPhisyst_001[1][2] = TFile::Open("001/ArithmeticAverageppDzeroDstarDplus8to16_assoc0.3to1.0.root");
 dPhisyst_001[1][3] = TFile::Open("001/ArithmeticAverageppDzeroDstarDplus16to24_assoc0.3to1.0.root"); 
 dPhisyst_001[2][0] = TFile::Open("001/ArithmeticAverageppDzeroDstarDplus3to5_assoc1.0to2.0.root");
 dPhisyst_001[2][1] = TFile::Open("001/ArithmeticAverageppDzeroDstarDplus5to8_assoc1.0to2.0.root");
 dPhisyst_001[2][2] = TFile::Open("001/ArithmeticAverageppDzeroDstarDplus8to16_assoc1.0to2.0.root");
 dPhisyst_001[2][3] = TFile::Open("001/ArithmeticAverageppDzeroDstarDplus16to24_assoc1.0to2.0.root");
 dPhisyst_001[3][0] = TFile::Open("001/ArithmeticAverageppDzeroDstarDplus3to5_assoc2.0to3.0.root");
 dPhisyst_001[3][1] = TFile::Open("001/ArithmeticAverageppDzeroDstarDplus5to8_assoc2.0to3.0.root");
 dPhisyst_001[3][2] = TFile::Open("001/ArithmeticAverageppDzeroDstarDplus8to16_assoc2.0to3.0.root");
 dPhisyst_001[3][3] = TFile::Open("001/ArithmeticAverageppDzeroDstarDplus16to24_assoc2.0to3.0.root"); 

 dPhisyst_0110[0][0] = TFile::Open("0110/ArithmeticAverageppDzeroDstarDplus3to5_assoc0.3to99.0.root");
 dPhisyst_0110[0][1] = TFile::Open("0110/ArithmeticAverageppDzeroDstarDplus5to8_assoc0.3to99.0.root");
 dPhisyst_0110[0][2] = TFile::Open("0110/ArithmeticAverageppDzeroDstarDplus8to16_assoc0.3to99.0.root");
 dPhisyst_0110[0][3] = TFile::Open("0110/ArithmeticAverageppDzeroDstarDplus16to24_assoc0.3to99.0.root");
 dPhisyst_0110[1][0] = TFile::Open("0110/ArithmeticAverageppDzeroDstarDplus3to5_assoc0.3to1.0.root");
 dPhisyst_0110[1][1] = TFile::Open("0110/ArithmeticAverageppDzeroDstarDplus5to8_assoc0.3to1.0.root");
 dPhisyst_0110[1][2] = TFile::Open("0110/ArithmeticAverageppDzeroDstarDplus8to16_assoc0.3to1.0.root");
 dPhisyst_0110[1][3] = TFile::Open("0110/ArithmeticAverageppDzeroDstarDplus16to24_assoc0.3to1.0.root"); 
 dPhisyst_0110[2][0] = TFile::Open("0110/ArithmeticAverageppDzeroDstarDplus3to5_assoc1.0to2.0.root");
 dPhisyst_0110[2][1] = TFile::Open("0110/ArithmeticAverageppDzeroDstarDplus5to8_assoc1.0to2.0.root");
 dPhisyst_0110[2][2] = TFile::Open("0110/ArithmeticAverageppDzeroDstarDplus8to16_assoc1.0to2.0.root");
 dPhisyst_0110[2][3] = TFile::Open("0110/ArithmeticAverageppDzeroDstarDplus16to24_assoc1.0to2.0.root");
 dPhisyst_0110[3][0] = TFile::Open("0110/ArithmeticAverageppDzeroDstarDplus3to5_assoc2.0to3.0.root");
 dPhisyst_0110[3][1] = TFile::Open("0110/ArithmeticAverageppDzeroDstarDplus5to8_assoc2.0to3.0.root");
 dPhisyst_0110[3][2] = TFile::Open("0110/ArithmeticAverageppDzeroDstarDplus8to16_assoc2.0to3.0.root");
 dPhisyst_0110[3][3] = TFile::Open("0110/ArithmeticAverageppDzeroDstarDplus16to24_assoc2.0to3.0.root"); 

 dPhisyst_1030[0][0] = TFile::Open("1030/ArithmeticAverageppDzeroDstarDplus3to5_assoc0.3to99.0.root");
 dPhisyst_1030[0][1] = TFile::Open("1030/ArithmeticAverageppDzeroDstarDplus5to8_assoc0.3to99.0.root");
 dPhisyst_1030[0][2] = TFile::Open("1030/ArithmeticAverageppDzeroDstarDplus8to16_assoc0.3to99.0.root");
 dPhisyst_1030[0][3] = TFile::Open("1030/ArithmeticAverageppDzeroDstarDplus16to24_assoc0.3to99.0.root");
 dPhisyst_1030[1][0] = TFile::Open("1030/ArithmeticAverageppDzeroDstarDplus3to5_assoc0.3to1.0.root");
 dPhisyst_1030[1][1] = TFile::Open("1030/ArithmeticAverageppDzeroDstarDplus5to8_assoc0.3to1.0.root");
 dPhisyst_1030[1][2] = TFile::Open("1030/ArithmeticAverageppDzeroDstarDplus8to16_assoc0.3to1.0.root");
 dPhisyst_1030[1][3] = TFile::Open("1030/ArithmeticAverageppDzeroDstarDplus16to24_assoc0.3to1.0.root"); 
 dPhisyst_1030[2][0] = TFile::Open("1030/ArithmeticAverageppDzeroDstarDplus3to5_assoc1.0to2.0.root");
 dPhisyst_1030[2][1] = TFile::Open("1030/ArithmeticAverageppDzeroDstarDplus5to8_assoc1.0to2.0.root");
 dPhisyst_1030[2][2] = TFile::Open("1030/ArithmeticAverageppDzeroDstarDplus8to16_assoc1.0to2.0.root");
 dPhisyst_1030[2][3] = TFile::Open("1030/ArithmeticAverageppDzeroDstarDplus16to24_assoc1.0to2.0.root");
 dPhisyst_1030[3][0] = TFile::Open("1030/ArithmeticAverageppDzeroDstarDplus3to5_assoc2.0to3.0.root");
 dPhisyst_1030[3][1] = TFile::Open("1030/ArithmeticAverageppDzeroDstarDplus5to8_assoc2.0to3.0.root");
 dPhisyst_1030[3][2] = TFile::Open("1030/ArithmeticAverageppDzeroDstarDplus8to16_assoc2.0to3.0.root");
 dPhisyst_1030[3][3] = TFile::Open("1030/ArithmeticAverageppDzeroDstarDplus16to24_assoc2.0to3.0.root"); 

 dPhisyst_30100[0][0] = TFile::Open("30100/ArithmeticAverageppDzeroDstarDplus3to5_assoc0.3to99.0.root");
 dPhisyst_30100[0][1] = TFile::Open("30100/ArithmeticAverageppDzeroDstarDplus5to8_assoc0.3to99.0.root");
 dPhisyst_30100[0][2] = TFile::Open("30100/ArithmeticAverageppDzeroDstarDplus8to16_assoc0.3to99.0.root");
 dPhisyst_30100[0][3] = TFile::Open("30100/ArithmeticAverageppDzeroDstarDplus16to24_assoc0.3to99.0.root");
 dPhisyst_30100[1][0] = TFile::Open("30100/ArithmeticAverageppDzeroDstarDplus3to5_assoc0.3to1.0.root");
 dPhisyst_30100[1][1] = TFile::Open("30100/ArithmeticAverageppDzeroDstarDplus5to8_assoc0.3to1.0.root");
 dPhisyst_30100[1][2] = TFile::Open("30100/ArithmeticAverageppDzeroDstarDplus8to16_assoc0.3to1.0.root");
 dPhisyst_30100[1][3] = TFile::Open("30100/ArithmeticAverageppDzeroDstarDplus16to24_assoc0.3to1.0.root"); 
 dPhisyst_30100[2][0] = TFile::Open("30100/ArithmeticAverageppDzeroDstarDplus3to5_assoc1.0to2.0.root");
 dPhisyst_30100[2][1] = TFile::Open("30100/ArithmeticAverageppDzeroDstarDplus5to8_assoc1.0to2.0.root");
 dPhisyst_30100[2][2] = TFile::Open("30100/ArithmeticAverageppDzeroDstarDplus8to16_assoc1.0to2.0.root");
 dPhisyst_30100[2][3] = TFile::Open("30100/ArithmeticAverageppDzeroDstarDplus16to24_assoc1.0to2.0.root");
 dPhisyst_30100[3][0] = TFile::Open("30100/ArithmeticAverageppDzeroDstarDplus3to5_assoc2.0to3.0.root");
 dPhisyst_30100[3][1] = TFile::Open("30100/ArithmeticAverageppDzeroDstarDplus5to8_assoc2.0to3.0.root");
 dPhisyst_30100[3][2] = TFile::Open("30100/ArithmeticAverageppDzeroDstarDplus8to16_assoc2.0to3.0.root");
 dPhisyst_30100[3][3] = TFile::Open("30100/ArithmeticAverageppDzeroDstarDplus16to24_assoc2.0to3.0.root");    
/*
 for(int i=0; i<4; i++) {
  printf("syst_std_001[%d] = %p\n",i,syst_std_001[i]); //for the syst uncertainties
  printf("syst_BDflat_001[%d] = %p\n",i,syst_BDflat_001[i]); //for the syst uncertainties
  printf("syst_BDextr_001[%d] = %p\n",i,syst_BDextr_001[i]); //for the syst uncertainties
  printf("syst_std_0110[%d] = %p\n",i,syst_std_0110[i]); //for the syst uncertainties  
  printf("syst_BDflat_0110[%d] = %p\n",i,syst_BDflat_0110[i]); //for the syst uncertainties
  printf("syst_BDextr_0110[%d] = %p\n",i,syst_BDextr_0110[i]); //for the syst uncertainties
  printf("syst_std_1030[%d] = %p\n",i,syst_std_1030[i]); //for the syst uncertainties
  printf("syst_BDflat_1030[%d] = %p\n",i,syst_BDflat_1030[i]); //for the syst uncertainties
  printf("syst_BDextr_1030[%d] = %p\n",i,syst_BDextr_1030[i]); //for the syst uncertainties
  printf("syst_std_30100[%d] = %p\n",i,syst_std_30100[i]); //for the syst uncertainties
  printf("syst_BDflat_30100[%d] = %p\n",i,syst_BDflat_30100[i]); //for the syst uncertainties
  printf("syst_BDextr_30100[%d] = %p\n",i,syst_BDextr_30100[i]); //for the syst uncertainties
  printf("totalsyst_NSy_001[%d] = %p\n",i,totalsyst_NSy_001[i]); //for the syst uncertainties
  printf("totalsyst_NSw_001[%d] = %p\n",i,totalsyst_NSw_001[i]); //for the syst uncertainties
  printf("totalsyst_NSy_0110[%d] = %p\n",i,totalsyst_NSy_0110[i]); //for the syst uncertainties
  printf("totalsyst_NSw_0110[%d] = %p\n",i,totalsyst_NSw_0110[i]); //for the syst uncertainties
  printf("totalsyst_NSy_1030[%d] = %p\n",i,totalsyst_NSy_1030[i]); //for the syst uncertainties
  printf("totalsyst_NSw_1030[%d] = %p\n",i,totalsyst_NSw_1030[i]); //for the syst uncertainties
  printf("totalsyst_NSy_30100[%d] = %p\n",i,totalsyst_NSy_30100[i]); //for the syst uncertainties
  printf("totalsyst_NSw_30100[%d] = %p\n",i,totalsyst_NSw_30100[i]); //for the syst uncertainties
   for(int j=0; j<4; j++) {
    printf("dPhisyst_001[%d][%d] = %p\n",i,j,dPhisyst_001[i][j]); //for the syst uncertainties
    printf("dPhisyst_0110[%d][%d] = %p\n",i,j,dPhisyst_0110[i][j]); //for the syst uncertainties
    printf("dPhisyst_1030[%d][%d] = %p\n",i,j,dPhisyst_1030[i][j]); //for the syst uncertainties
    printf("dPhisyst_30100[%d][%d] = %p\n",i,j,dPhisyst_30100[i][j]); //for the syst uncertainties
   }
 }
 */
}