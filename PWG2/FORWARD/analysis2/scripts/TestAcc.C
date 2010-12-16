#include <TMath.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TMarker.h>
#include <TLine.h>
#include <TArc.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <THStack.h>

//_____________________________________________________________________
void AcceptanceCorrection(Char_t r, UShort_t t, Float_t& oldm, Float_t& newm)
{
  // AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  //AliFMDRing fmdring(ring);
  //fmdring.Init();
  // Float_t   rad       = pars->GetMaxR(r)-pars->GetMinR(r);
  // Float_t   slen      = pars->GetStripLength(r,t);
  // Float_t   sblen     = pars->GetBaseStripLength(r,t);
  const Double_t ic1[] = { 4.9895, 15.3560 };
  const Double_t ic2[] = { 1.8007, 17.2000 };
  const Double_t oc1[] = { 4.2231, 26.6638 };
  const Double_t oc2[] = { 1.8357, 27.9500 };

  Double_t  minR      = (r == 'I' ?  4.5213 : 15.4);
  Double_t  maxR      = (r == 'I' ? 17.2    : 28.0);
  Double_t  rad       = maxR - minR;
  
  Int_t     nStrips   = (r == 'I' ? 512 : 256);
  Int_t     nSec      = (r == 'I' ?  20 :  40);
  Float_t   segment   = rad / nStrips;
  Float_t   basearc   = 2 * TMath::Pi() / (0.5 * nSec);
  Float_t   radius    = minR + t * segment;
  Float_t   baselen   = 0.5 * basearc * radius;

  const Double_t* c1  = (r == 'I' ? ic1 : oc1);
  const Double_t* c2  = (r == 'I' ? ic2 : oc2);

  // If the radius of the strip is smaller than the radius corresponding 
  // to the first corner we have a full strip length 
  Float_t   cr        = TMath::Sqrt(c1[0]*c1[0]+c1[1]*c1[1]);
  if (radius <= cr) newm = 1;
  else {
    // Next, we should find the end-point of the strip - that is, 
    // the coordinates where the line from c1 to c2 intersects a circle 
    // with radius given by the strip. 
    // (See http://mathworld.wolfram.com/Circle-LineIntersection.html)
    Float_t D   = c1[0] * c2[1] - c1[1] * c2[0];
    Float_t dx  = c2[0] - c1[0];
    Float_t dy  = c2[1] - c1[1];
    Float_t dr  = TMath::Sqrt(dx*dx+dy*dy);
    Float_t det = radius * radius * dr * dr - D*D;
    
    if      (det <  0) newm = 1; // No intersection
    else if (det == 0) newm = 1; // Exactly tangent
    else {
      Float_t x   = (+D * dy + dx * TMath::Sqrt(det)) / dr / dr;
      Float_t y   = (-D * dx + dy * TMath::Sqrt(det)) / dr / dr;
      Float_t th  = TMath::ATan2(x, y);

      newm = th / (basearc/2);
    }
  }

  Float_t   slope     = (c1[1] - c2[1]) / (c1[0] - c2[0]);
  Float_t   constant  = (c2[1] * c1[0] - c2[0] * c1[1]) / (c1[0]-c2[0]);
  
  Float_t   d         = (TMath::Power(TMath::Abs(radius*slope),2) + 
			 TMath::Power(radius,2) - TMath::Power(constant,2));
  
  // If below corners return 1
  Float_t arclen = baselen;
  if (d > 0) {
    Float_t   x         = ((-TMath::Sqrt(d) - slope * constant) / 
			   (1+TMath::Power(slope, 2)));
    Float_t   y         = slope*x + constant;
    
    // If x is larger than corner x or y less than corner y, we have full
    // length strip
    if(x < c1[0] && y> c1[1]) {
      //One sector since theta is by definition half-hybrid
      Float_t   theta     = TMath::ATan2(x,y);
      arclen              = radius * theta;
    }
  }
  // Calculate the area of a strip with no cut
  Float_t   basearea1 = 0.5 * baselen * TMath::Power(radius,2);
  Float_t   basearea2 = 0.5 * baselen * TMath::Power((radius-segment),2);
  Float_t   basearea  = basearea1 - basearea2;
  
  // Calculate the area of a strip with cut
  Float_t   area1     = 0.5 * arclen * TMath::Power(radius,2);
  Float_t   area2     = 0.5 * arclen * TMath::Power((radius-segment),2);
  Float_t   area      = area1 - area2;
  
  // Acceptance is ratio 
  oldm = area/basearea;
}

void DrawSolution(Char_t r, UShort_t dt=16, UShort_t offT=128)
{
  TCanvas* c = new TCanvas(Form("c%c", r), r == 'I' ? 
			   "Inner" : "Outer", 800, 500);
  
  const Double_t ic1[] = { 4.9895, 15.3560 };
  const Double_t ic2[] = { 1.8007, 17.2000 };
  const Double_t oc1[] = { 4.2231, 26.6638 };
  const Double_t oc2[] = { 1.8357, 27.9500 };

  const Double_t* c1   = (r == 'I' ? ic1 : oc1);
  const Double_t* c2   = (r == 'I' ? ic2 : oc2);
  Double_t  minR       = (r == 'I' ?  4.5213 : 15.4);
  Double_t  maxR       = (r == 'I' ? 17.2    : 28.0);
  Double_t  rad        = maxR - minR;
  Float_t   cr         = TMath::Sqrt(c1[0]*c1[0]+c1[1]*c1[1]);
  Double_t  theta      = (r == 'I' ? 18 : 9);
  TH2*   frame = new TH2F("frame", "Frame", 100, -10, 10, 100, 0, 30);
  frame->Draw();

  TLine* line = new TLine(c1[0], c1[1], c2[0], c2[1]);
  line->SetLineColor(kBlue+1);
  line->Draw();

  UShort_t nT = (r == 'I' ? 512 : 256);
  for (Int_t t = offT; t < nT; t += dt) {
    Double_t radius = minR + t * rad / nT;

    TArc*    circle = new TArc(0, 0, radius, 90-theta, 90+theta);
    circle->SetFillColor(0);
    circle->SetFillStyle(0);
    circle->Draw();

    // Now find the intersection 
    if (radius <= cr) continue;

    // Next, we should find the end-point of the strip - that is, 
    // the coordinates where the line from c1 to c2 intersects a circle 
    // with radius given by the strip. 
    // (See http://mathworld.wolfram.com/Circle-LineIntersection.html)
    Float_t D   = c1[0] * c2[1] - c1[1] * c2[0];
    Float_t dx  = c2[0] - c1[0];
    Float_t dy  = c2[1] - c1[1];
    Float_t dr  = TMath::Sqrt(dx*dx+dy*dy);
    Float_t det = radius * radius * dr * dr - D*D;
    
    if (det <  0) continue;
    if (det == 0) continue;


    Float_t x   = (+D * dy - dx * TMath::Sqrt(det)) / dr / dr;
    Float_t y   = (-D * dx - dy * TMath::Sqrt(det)) / dr / dr;
    
    TMarker* m = new TMarker(x, y, 20);
    m->SetMarkerColor(kRed+1);
    m->Draw();

    x   = (+D * dy + dx * TMath::Sqrt(det)) / dr / dr;
    y   = (-D * dx + dy * TMath::Sqrt(det)) / dr / dr;

    // Float_t th = TMath::ATan2(x,y);
    // Printf("theta of strip %3d: %f degrees", 180. / TMath::Pi() * th);

    m = new TMarker(x, y, 20);
    m->SetMarkerColor(kGreen+1);
    m->Draw();
  }
  c->cd();
}
    
  
void TestAcc()
{
  TCanvas* c =  new TCanvas("c", "C");
  TGraph*      innerOld = new TGraph(512);
  TGraph*      outerOld = new TGraph(256);
  TGraph*      innerNew = new TGraph(512);
  TGraph*      outerNew = new TGraph(256);
  innerOld->SetName("innerOld");
  innerOld->SetName("Inner type, old");
  innerOld->SetMarkerStyle(21);
  innerOld->SetMarkerColor(kRed+1);
  outerOld->SetName("outerOld");
  outerOld->SetName("Outer type, old");
  outerOld->SetMarkerStyle(21);
  outerOld->SetMarkerColor(kBlue+1);
  innerNew->SetName("innerNew");
  innerNew->SetName("Inner type, new");
  innerNew->SetMarkerStyle(20);
  innerNew->SetMarkerColor(kGreen+1);
  outerNew->SetName("outerNew");
  outerNew->SetName("Outer type, new");
  outerNew->SetMarkerStyle(20);
  outerNew->SetMarkerColor(kMagenta+1);

  TMultiGraph* all      = new TMultiGraph("all", "Acceptances");
  all->Add(innerOld);
  all->Add(outerOld);
  all->Add(innerNew);
  all->Add(outerNew);

  TH2* innerCor = new TH2F("innerCor", "Inner correlation", 100,0,1,100,0,1);
  TH2* outerCor = new TH2F("outerCor", "Outer correlation", 100,0,1,100,0,1);
  innerCor->SetMarkerStyle(20);
  outerCor->SetMarkerStyle(20);
  innerCor->SetMarkerColor(kRed+1);
  outerCor->SetMarkerColor(kBlue+1);
  // innerCor->SetMarkerSize(1.2);
  // outerCor->SetMarkerSize(1.2);
  THStack* stack = new THStack("corr", "Correlations");
  stack->Add(innerCor);
  stack->Add(outerCor);

  for (Int_t i = 0; i < 512; i++) { 
    Float_t oldAcc, newAcc;
    
    AcceptanceCorrection('I', i, oldAcc, newAcc);
    innerOld->SetPoint(i, i, oldAcc);
    innerNew->SetPoint(i, i, newAcc);
    // Printf("Inner strip # %3d:  old=%8.6f, new=%8.6f", i, oldAcc, newAcc);

    innerCor->Fill(oldAcc, newAcc);
    if (i >= 256) continue;
    
    AcceptanceCorrection('O', i, oldAcc, newAcc);
    outerOld->SetPoint(i, i, oldAcc);
    outerNew->SetPoint(i, i, newAcc);
    // Printf("Outer strip # %3d:  old=%8.6f, new=%8.6f", i, oldAcc, newAcc);
    outerCor->Fill(oldAcc, newAcc);
  }

  all->Draw("ap");

  DrawSolution('I',4,256);
  DrawSolution('O',4,128);

  TCanvas* c2 = new TCanvas("cc", "cc");
  stack->Draw("nostack p");
  TLine* l = new TLine(0,0,1,1);
  l->SetLineStyle(2);
  l->Draw();

  c2->cd();
}
