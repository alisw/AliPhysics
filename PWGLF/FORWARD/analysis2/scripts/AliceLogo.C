#ifndef ALICELOGO_C
#define ALICELOGO_C
#include <TGraph.h>
#include <TBox.h>
#include <TEllipse.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TImage.h>
#include <TH1.h>
#include <TLine.h>
#include <TColor.h>

/** 
 * Structure to draw the ALICE logo on a pad (canvas).  The logo and
 * text is drawn using ROOT primitives and the result is somewhat
 * scalable.
 *
 * Note, you can save the logo as a PNG, PDF, or SVG to include it in
 * some other document or the like.  
 *
 * How to use
 * 
 * @code 
 * {
 *   TH1* h = new TH1F("h", "h", 60, -3, 3);
 *   h->SetFillColor(kRed+1);
 *   h->SetFillStyle(3001);
 *   h->FillRandom("gaus");
 *   h->Draw();
 *   
 *   gROOT->LoadMacro("AliceLogo.C");
 *   AliceLogo* al = new AliceLogo();
 *   // al->Draw(0, .75, .4, .25, "preliminary");
 *   // al->Draw(0, .75, .4, .25, "performance");
 *   // al->Draw(0, .75, .4, .25, "","text l3 many particles");
 *   al->Draw(0, .75, .4, .25, "preliminary", "text l3 many particles");
 * }
 * @endcode
 *
 * There are many options for colours and options for preliminary,
 * performance and just the logo.
 *
 * The logo produced does not exactly match the logos found on the
 * ALICE web-site 
 *
 *   https://aliceinfo.cern.ch/Figure/alice_logo
 *
 * In particular, the 'C' was hard to get right.  However, the logo
 * does follow the charter as described in 
 * 
 *   http://aliweb.cern.ch/system/files/documents/charter/charte_A5_120628.pdf
 *
 * Use the static member function(s) AliceLogo::Check to see how well
 * this code reproduces the logos.  Note, you need the official logos
 * in PNG format for this.  Here's a list of places to get them:
 *
 * - https://aliceinfo.cern.ch/Figure/sites/aliceinfo.cern.ch.Figure/files/Figures/General/elopez/2012-Jul-04-Preliminary_Logo.png
 * - https://aliceinfo.cern.ch/Figure/sites/aliceinfo.cern.ch.Figure/files/Figures/General/elopez/2012-Jul-04-Performance_Logo.png
 * - https://aliceinfo.cern.ch/Figure/sites/aliceinfo.cern.ch.Figure/files/Figures/General/elopez/2012-Jul-04-4_Color_Logo_CB.png
 *
 * @author Christian Holm Christensen <cholm@nbi.dk>
 */
struct AliceLogo 
{
  /** 
   * Tags
   */
  enum ETag {
    kNone = 0, 
    kPreliminary, 
    kPerformance
  };
  /** 
   * Colour options 
   */
  enum {
    kBW          = 0,    // Black and white only 
    kText        = 0x1,  // Text is coloured dark-blue
    kL3          = 0x2,  // Colour L3 red
    kParticles   = 0x4,  // Particles coloured red
    kMany        = 0x8,  // Many particles in many colours
    kPad         = 0x10, // BAckground from pad 
    kReverse     = 0x20  // Text is white 
  };
  /**
   * Constructor 
   */
  AliceLogo(Bool_t debug=false)
    : fPad(0), 
      fAspect(1),
      fL3Thick(.05),
      fL3Height(.63), 
      fL3Width(0),
      fDebug(debug)
  {
    fRed    = TColor::GetColor(226,   0, 26);
    fYellow = TColor::GetColor(238, 125, 17);
    fPink   = TColor::GetColor(202,  71, 67);
    fBlue   = TColor::GetColor(40,   58, 68);
  }
  /** 
   * Parse a string with colour options 
   * 
   * @param what 
   * 
   * @return Colour flags 
   */
  static UInt_t ParseColors(const TString& what)
  {
    UInt_t col = AliceLogo::kBW;
    if (what.Contains("l3",    TString::kIgnoreCase)) col |= kL3;
    if (what.Contains("text",  TString::kIgnoreCase)) col |= kText;
    if (what.Contains("part",  TString::kIgnoreCase)) col |= kParticles;
    if (what.Contains("many",  TString::kIgnoreCase)) col |= kMany;
    if (what.Contains("pad",   TString::kIgnoreCase)) col |= kPad;
    if (what.Contains("rever", TString::kIgnoreCase)) col |= kReverse;
    return col;
  }
  /** 
   * Parse a tag string 
   * 
   * @param what String to parse 
   * 
   * @return Tag 
   */
  static UShort_t ParseTag(const TString& what)
  {
    UShort_t tag = AliceLogo::kNone;
    if      (what.BeginsWith("pre", TString::kIgnoreCase)) tag = kPreliminary;
    else if (what.BeginsWith("per", TString::kIgnoreCase)) tag = kPerformance;
    return tag;
  }
  /**
   * Draw everything 
   * 
   * @param p      Pad to draw in
   * @param x      Starting X position
   * @param y      Starting Y position
   * @param h      Relative height 
   * @param tag    Tag 
   * @param colors Color options
   */
  void Draw(TVirtualPad* p, Double_t x=.77, Double_t y=.35, Double_t h=.25, 
	    const TString& tag="", const TString& colors="")
  {
    Draw(p, x, y, h, ParseTag(tag), ParseColors(colors));
  }
  /** 
   * Draw everything 
   * 
   * @param p      Pad to draw in 
   * @param x      Starting X position
   * @param y      Starting Y position
   * @param h      Relative height 
   * @param tag    Tag 
   * @param color  Color options
   */
  void Draw(TVirtualPad* p, Double_t x, Double_t y, Double_t h,
	    UShort_t tag, UInt_t color)
  {
    if (!p) p = gPad;
    if (!p) { 
      Error("AliceLogo::Draw", "No pad to draw in");
      return;
    }
    if (fDebug) 
      Printf("Pad to draw logo in: %s", p->GetName());
    // p->SetFixedAspectRatio(true);
    Int_t    pH  = p->YtoPixel(0);
    if (fDebug) {
      Int_t    pW  = p->XtoPixel(1);
      Double_t pA  = Double_t(pW) / pH;
      Printf("pA=%d/%d=%lf", pW, pH, pA);
    }


    if (tag == 0) { 
      fAspect   =  1437. / 1933. ;
      fL3Height = .75;
      fL3Thick  *= 1.1;
    }
    else { 
      fAspect   = 1835. / 2272. ;
      fL3Height = .63;
    }
    fL3Width   = fL3Height / fAspect;
    Int_t    hP   = Int_t(h * pH+.5); // Height in pixels 
    Int_t    wP   = Int_t(fAspect * hP + .5); // Width in pixels 
    Int_t    xP   = p->XtoPixel(x); // X pos in pixels 
    Int_t    rP   = xP + wP; // Right in pixels; 
    Double_t topX = p->PixeltoX(rP);
    Double_t topY = y + h; // p->PixeltoY(tP); // y + h;
    if (fDebug){
      Int_t    yP   = p->YtoPixel(y); // X pos in pixels 
      Int_t    tP   = yP - hP; // Top in pixels
      Printf("(x,y)=(%d,%d) hP=%d wP=%d xP=%d yP=%d rP=%d tP=%d", 
	     xP, yP, hP, wP, xP, yP, rP, tP);
    }
    if (topX > 1 || topY > 1) { 
      Warning("AliceLogo::Draw", "with (x,y)=(%f,%f) and h=%f -> top=(%f,%f)",
	      x, y, h, topX, topY);
      topX = TMath::Min(topX, 1.);
      topY = TMath::Min(topY, 1.);
    }

    p->cd();
    fPad = new TPad("logo", "logo", x, y, topX, topY);
    fPad->SetFixedAspectRatio(true);
    fPad->SetFillColor(fDebug ? kGreen-4 : 0);
    fPad->SetFillStyle(fDebug ? 4030     : 0);
    fPad->SetBorderSize(fDebug ? 1 : 0);
    fPad->SetBorderMode(fDebug ? 1 : 0);
    fPad->SetLineColor(fDebug ? kBlue+1 : 0);
    fPad->SetLineWidth(fDebug ? 1 : 0);
    fPad->SetTopMargin(0);
    fPad->SetBottomMargin(0);
    fPad->SetLeftMargin(0);
    fPad->SetRightMargin(0);
    fPad->Draw();
    if (fDebug)
      Printf("(x,y)=(%f,%f) (tx,ty)=(%f,%f) (x1,y1)=(%f,%f) (x2,y2)=(%f,%f)", 
	     x, y, topX, topY, GetX1(), GetY1(), GetX2(), GetY2());
    
    if (fDebug) {
      Int_t    ppH  = fPad->YtoPixel(0);
      Int_t    ppW  = fPad->XtoPixel(1);
      Double_t ppA  = Double_t(ppW) / ppH;
      Printf("ppA=%d/%d=%lf A=%lf",ppW, ppH, ppA, fAspect);      
      // Printf("ppA/A=%lf/%lf=%lf", ppA, fAspect, ppA/fAspect);
    }

    Int_t bg = (color & kPad) ? p->GetFillColor() : Int_t(kWhite);

    DrawL3(color, (color & kReverse && !(color & kL3) ? bg : kWhite));

    Int_t txtFg = kBlack;
    if (color & kText)
      txtFg = (color & kReverse ? kWhite : fBlue);
    DrawALICE(txtFg, bg);
    if      (tag == kPreliminary) DrawLabel(txtFg, bg, true);
    else if (tag == kPerformance) DrawLabel(txtFg, bg, false);
   
    p->cd();
  }
  /** 
   * Draw the L3 outline.  Note, this is put in a sub-pad of it's own. 
   * 
   * @param color Colour flags
   * @param bg    Background color 
   */
  void DrawL3(UInt_t color, Int_t bg)
  {
    printf("Drawing L3 image ");
    fPad->cd();
    Double_t       o = TMath::Max((1-fL3Width)/2,0.);
    TVirtualPad* pad = new TPad("l3","L3",o, 1-fL3Height, 
				TMath::Min(fL3Width+o,1.), 1);
    pad->SetFixedAspectRatio(true);
    pad->SetFillColor(fDebug ? kGreen-4 : 0);
    pad->SetFillStyle(fDebug ? 4030     : 0);
    pad->SetBorderSize(fDebug ? 1 : 0);
    pad->SetBorderMode(fDebug ? 1 : 0);
    pad->SetLineColor(fDebug ? kRed+1 : 0);
    pad->SetLineWidth(fDebug ? 1 : 0);
    pad->SetTopMargin(0);
    pad->SetBottomMargin(0);
    pad->SetLeftMargin(0);
    pad->SetRightMargin(0);
    pad->Draw();
    
    Int_t lcol = (color & kL3 ? fRed : kBlack);
    if (color & kReverse && lcol == kBlack) lcol = kWhite;

    Double_t dA = 2 * TMath::Pi() / 8;
    Double_t ro = .5 + .5 * (1 -  TMath::Cos(dA/2));
    Double_t ri = ro - 1. / fL3Height * fL3Thick;
    Double_t x0 = .5;
    Double_t y0 = .5;
    TGraph*  gi = new TGraph(9);
    TGraph*  go = new TGraph(9);
    gi->SetName("L3_inner");
    go->SetName("L3_outer");
    SetAttr(gi, gi, bg);
    SetAttr(go, go, lcol);
    for (Int_t i = 0; i < 9; i++) { 
      Double_t a = (i + .5) * dA;
      go->SetPoint(i, x0 + ro * TMath::Cos(a), y0 + ro * TMath::Sin(a));
      gi->SetPoint(i, x0 + ri * TMath::Cos(a), y0 + ri * TMath::Sin(a));
    }
    pad->cd();
    printf(".");
    go->Draw("lf same");
    printf(".");
    gi->Draw("lf same");

    Double_t len = .26;
    if (color & kMany) {
      Double_t w = 0.02;
      for (Int_t i = 0; i < 36; i++) {
	Int_t col = fYellow;
	switch (i) { 
	case 1: case 5: case 10: case 20: case 24: case 27: case 33: case 34:
	  col = fBlue; break;
	case 3: case 4: case 8: case 16: case 19: case 23: case 29: case 30: 
	case 32: 
	  col = fPink; break;
	case 9: case 12: case 14: case 18: case 25: case 28: case 35:
	  col = fRed; break;
	default:
	  break;
	}
	double l = len;
	switch (i) { 
	case 1: case 9: case 8: case 10: case 19: case 21: case 33:
	  l += 0.005; break; 
	case 2: case 7:  case 11: case 12: case 15: case 16: case 20: case 34:
	  l += 0.02;  break;
	case 3:  case 6: case 13: case 17: case 25: case 29:
	  l += 0.01;  break;
	case 23: case 27: case 28: 
	  l -= 0.005;
	default:      break;
	}

	printf(".");
	DrawSpoke(i * 360/36, w, l, col);
      }
    }
    else {
      Double_t w  = .015;
      Int_t col = (color & kParticles ? fRed : kBlack);
      if (bg != kWhite) col = kWhite;
      for (Int_t i = 0; i < 18; i++) {
	double l = len;
	switch (i) { 
	case 0: case 2:                     
	  break;
	case 1:  case 6:  case 8:         
	case 10: case 17:
	  l += .02; break;
	case 3: case 4: case 5: case 7: case 12: 
	  l += .01; break;
	case 15: 
	  l += .005; break;
	case 16:
	  l -= .005; break;
	default:
	  break;
	}
	printf(".");
	DrawSpoke(i * 360/18, w, l, col);
      }
    }
    fPad->cd();
    printf("\n");
  }
  /** 
   * Draw a spoke in the L3 outline. 
   * 
   * @param ang Rotation angle in degrees. 
   * @param w   Width
   * @param len Length of spoke 
   * @param fg  Foreground colour 
   */
  void DrawSpoke(Int_t ang, Double_t w, Double_t len, Int_t fg)
  {
    Double_t r1 = .12;

    Double_t xx = r1;
    Double_t yy = w/2;
    Double_t xl = r1+len;

    Double_t rad = TMath::Pi() / 180 * ang;
    Double_t xtl  = .5 + xx * TMath::Cos(rad) - yy * TMath::Sin(rad);
    Double_t ytl  = .5 + xx * TMath::Sin(rad) + yy * TMath::Cos(rad);
    Double_t xbl  = .5 + xx * TMath::Cos(rad) + yy * TMath::Sin(rad);
    Double_t ybl  = .5 + xx * TMath::Sin(rad) - yy * TMath::Cos(rad);
    Double_t xtr  = .5 + xl * TMath::Cos(rad) - yy * TMath::Sin(rad);
    Double_t ytr  = .5 + xl * TMath::Sin(rad) + yy * TMath::Cos(rad);
    Double_t xbr  = .5 + xl * TMath::Cos(rad) + yy * TMath::Sin(rad);
    Double_t ybr  = .5 + xl * TMath::Sin(rad) - yy * TMath::Cos(rad);

    TGraph* tri = new TGraph(4);
    SetAttr(tri, tri, fg);
    tri->SetName(Form("spoke_tr_%03d", ang));
    tri->SetPoint(0, 0.5,0.5);
    tri->SetPoint(3, 0.5,0.5);
    tri->SetPoint(1, xtl, ytl);
    tri->SetPoint(2, xbl, ybl);
    tri->Draw("same lf");

    TGraph* lne = new TGraph(4);
    SetAttr(lne, lne, fg);
    lne->SetName(Form("spoke_ln_%03d", ang));
    lne->SetPoint(0, xtl, ytl);
    lne->SetPoint(1, xbl, ybl);
    lne->SetPoint(2, xbr, ybr);
    lne->SetPoint(3, xtr, ytr);
    lne->SetPoint(4, xtl, ytl);
    lne->Draw("same lf");
    
  }
  /** 
   * Draw the text ALICE in a separate pad 
   * 
   */
  void DrawALICE(Int_t fg, Int_t bg)
  {
    printf("Drawing ALICE title ");
    fPad->cd();
    Double_t o  = TMath::Max((1-fL3Width)/2,0.);
    Double_t r  = fL3Width/2;
    Double_t t  = r*TMath::Sin(TMath::Pi()/8-.05);
    Double_t t2 = r*TMath::Sin(TMath::Pi()/8+.05);
    // Middle-bottom-left point is 
    Double_t mblX = o;
    Double_t mblY = (1-fL3Height/2) - t;
    // Bottom-left point is 
    Double_t blY = (1-fL3Height);
    Double_t blX = .5 - t2;
    // Slope
    Double_t a = (blY - mblY) / (blX - mblX);
    // Y is then, using y  = y1 + a (x - x1)
    Double_t b = mblY + a * (.5 - mblX);

    Double_t top = 1-fL3Height-fL3Thick;
    Double_t hh  = 2 * (top - b);
    Double_t bot = TMath::Max(top - hh,0.);
    TVirtualPad* pad = new TPad("alice", "alice", o, bot+.005, 
				TMath::Min(o+fL3Width,1.), top+.005);
    pad->SetFixedAspectRatio(true);
    pad->SetFillColor(fDebug ? kGreen-4 : 0);
    pad->SetFillStyle(fDebug ? 4030     : 0);
    pad->SetBorderSize(fDebug ? 1 : 0);
    pad->SetBorderMode(fDebug ? 1 : 0);
    pad->SetLineColor(fDebug ? kMagenta+1 : 0);
    pad->SetLineWidth(fDebug ? 1 : 0);
    pad->SetTopMargin(0);
    pad->SetBottomMargin(0);
    pad->SetLeftMargin(0);
    pad->SetRightMargin(0);
    pad->Draw();
    pad->cd();

    Double_t wp = pad->XtoPixel(1);
    Double_t hp = pad->YtoPixel(0);
    a           = hp / wp;
    
    Double_t w = .15;
    printf(".");
    DrawA(0,            0.01, w, a, fg, bg);
    printf(".");
    DrawL(.5 - 5.8*a*w, 0.01, w, a, fg);
    printf(".");
    DrawI(.5 - 1.5*a*w, 0.01, w, a, fg);
    printf(".");
    DrawC(.5 + 2.6*a*w, 0.00, /*w,*/ a, fg, bg);
    printf(".");
    DrawE(1  - 3.5*a*w, 0.01, w, a, fg);
    fPad->cd();

    printf("\n");
  }
  /** 
   * Draw the letter A 
   * 
   * @param x    Base X coordinate 
   * @param y    Base Y coordinate 
   * @param w    Line width 
   * @param a    Aspect ratio of pad 
   * @param fg   Foreground colour
   * @param bg   Background colour
   */
  void DrawA(Double_t x, Double_t y, Double_t w, Double_t a,
	     Int_t fg, Int_t bg)
  {
    Double_t u = 4.5 * w;
    Double_t r = u / 2;
    TBox* b1 = new TBox(x, y, x + a * w, y + 1 - r);
    SetAttr(b1, b1, fg);
    b1->Draw();
    b1->DrawBox(x + (u-w) * a, y,         x + u * a, y + 1 - r);
    b1->DrawBox(x,             y+1-r-2*w, x + u * a,     y + 1 - r -w);
    TEllipse* a1 = new TEllipse(x+u*a/2, y+1-r-.01, (r-.005)*a, r, 0, 180);
    SetAttr(a1, a1, fg);
    a1->Draw("only");
    TEllipse* a2 = new TEllipse(x+u*a/2, y+1-r-.01, (r-w-.005)*a, r-w, 0, 180);
    SetAttr(a2, a2, bg);
    a2->Draw("only");

#if 0
    // To draw as lines - note line width does not scale 
    Int_t lw = w*gPad->YtoPixel(0);
    Double_t x1 = x+a*w/2;
    Double_t x2 = x+a*(u-w/2);
    Double_t y1 = y + 1 - r - 1.5*w;
    TLine*    l1 = new TLine(x1, y, x1, y+1-r);
    l1->SetLineWidth(lw);
    l1->SetLineColor(kBlue);
    l1->Draw();
    l1->DrawLine(x2, y,x2, y+1-r);
    l1->DrawLine(x1, y1,x2, y1);
    
    TEllipse* a3 = new TEllipse(x+u*a/2, y+1-r-.01, 
				(r-w/2-.005)*a, r-w/2, 0, 180);
    a3->SetFillStyle(0);
    a3->SetFillColor(0);
    // Int_t lw = w*gPad->YtoPixel(0);
    Printf("Setting line width to %d", lw);
    a3->SetLineWidth(lw);
    a3->SetLineColor(kBlue);
    a3->Draw("only");
#endif
  }
  /** 
   * Draw the letter L
   * 
   * @param x    Base X coordinate 
   * @param y    Base Y coordinate 
   * @param w    Line width 
   * @param a    Aspect ratio of pad 
   * @param fg   Foreground colour
   */
  void DrawL(Double_t x, Double_t y, Double_t w, Double_t a, Int_t fg)
  {
    TBox* b1 = new TBox(x, y, x+3.1 * a * w, y + w);
    SetAttr(b1, b1, fg);
    b1->Draw();
    b1->DrawBox(x, y , x + a * w,y+.99);
  }
  /** 
   * Draw the letter I
   * 
   * @param x    Base X coordinate 
   * @param y    Base Y coordinate 
   * @param w    Line width 
   * @param a    Aspect ratio of pad 
   * @param fg   Foreground colour
   */
  void DrawI(Double_t x, Double_t y, Double_t w, Double_t a, Int_t fg)
  {    
    TBox* b1 = new TBox(x, y, x+3 * a * w, y + w);
    SetAttr(b1, b1, fg);
    b1->Draw();
    b1->DrawBox(x+a*w,y,x+2*a*w,y+.99);
    b1->DrawBox(x,y+1-w,x+3*a*w,y+.99);
  }
  /** 
   * Draw the letter C
   * 
   * @param x    Base X coordinate 
   * @param y    Base Y coordinate 
   * @param a    Aspect ratio of pad 
   * @param fg   Foreground colour
   * @param bg   Background colour
   */
  void DrawC(Double_t x, Double_t y, Double_t a, Int_t fg, Int_t bg)
  {
    Int_t    nO   = 13;
    Double_t dx   = 0; 
    Double_t pO[] = { 
      /* 0 */  0.698245  +dx, 0.271355,
      /* 1 */  0.590204  +dx, 0.0755058,
      /* 2 */  0.424947  +dx, 0.00,
      /* 3 */  0.220343  +dx, 0.015668+.005,
      /* 4 */  0.0629551 +dx, 0.15668,
      /* 5 */  0.00786939+dx, 0.329027,
      /* 6 */  -0.005,        0.564046,
      /* 7 */  0.00786939+dx, 1-0.329027,
      /* 8 */  0.0629551 +dx, 1-0.15668,
      /* 9 */  0.220343  +dx, 1-0.015668-.005,
      /* 10 */ 0.424947  +dx, .99,
      /* 11 */ 0.590204  +dx, 1-0.0755058,
      /* 12 */ 0.6982450 +dx, 1-0.271355, 0 };
    TGraph* gO = new TGraph(nO);
    gO->SetName("c_ext");
    gO->SetLineColor(fDebug ? kGreen+1 : fg);
    gO->SetLineWidth(fDebug ? 1        : 0);
    gO->SetFillColor(fDebug ? kGreen+1 : fg);
    gO->SetFillStyle(fDebug ? 0 : 1001);
    for (Int_t i = 0; i < nO; i++) 
      gO->SetPoint(i, x+pO[2*i]*a,y+pO[2*i+1]);
    gO->Draw("fc same");
    Int_t    nI   = 11;
    Double_t pI[] = { 
      /* 0 */ 0.551869, 1-0.713617,
      /* 1 */ 0.481045, 1-0.821127,
      /* 2 */ 0.363004, 1-0.853799,
      /* 3 */ 0.213486, 1-0.783625,
      /* 4 */ 0.158508, 1-0.635277,
      /* 4 */ 0.158,    0.5, // 1-0.635277,
      /* 5 */ 0.158508, 0.635277,
      /* 6 */ 0.213486, 0.773625,
      /* 7 */ 0.363004, 0.853799,
      /* 8 */ 0.481045, 0.821127,
      /* 9 */ 0.551869, 0.713617, 0 };
    TGraph* gI = new TGraph(nI);
    gI->SetName("c_int");
    gI->SetLineColor(fDebug ? kGreen+1 : bg);
    gI->SetLineWidth(fDebug ? 1        : 0);
    gI->SetFillColor(fDebug ? kGreen+1 : bg);
    gI->SetFillStyle(fDebug ? 0 : 1001);
    for (Int_t i = 0; i < nI; i++) 
      gI->SetPoint(i, x+pI[2*i]*a,y+pI[2*i+1]);
    gI->Draw("fc same");


    TGraph* gC = new TGraph(5);
    gC->SetName("c_cut");
    SetAttr(gC, gC, bg);
    gC->SetPoint(0, gO->GetX()[0],    gO->GetY()[0]);
    gC->SetPoint(1, gI->GetX()[0],    gI->GetY()[0]);
    gC->SetPoint(2, gI->GetX()[nI-1], gI->GetY()[nI-1]);
    gC->SetPoint(3, gO->GetX()[nO-1], gO->GetY()[nO-1]);
    gC->SetPoint(4, gO->GetX()[0],    gO->GetY()[0]);
    gC->Draw("fl same");
  }
  /** 
   * Draw the letter E 
   * 
   * @param x    Base X coordinate 
   * @param y    Base Y coordinate 
   * @param w    Line width 
   * @param a    Aspect ratio of pad 
   * @param fg   Foreground colour
   */
  void DrawE(Double_t x, Double_t y, Double_t w, Double_t a, Int_t fg)
  {
    Double_t u = 3.5;
    TBox* b1 = new TBox(x, y, x + a * w, y + .99);
    SetAttr(b1, b1, fg);
    b1->Draw();
    b1->DrawBox(x, y,            x + u*a*w,      y+w);
    b1->DrawBox(x, y+.99-w,      x + u*a*w,      y+.99);
    b1->DrawBox(x, y+.5-w/2+.02, x + (u-.1)*a*w, y+.5+w/2+.02);
  }
  /** 
   * Draw the letter F 
   * 
   * @param x    Base X coordinate 
   * @param y    Base Y coordinate 
   * @param w    Line width 
   * @param a    Aspect ratio of pad 
   * @param fg   Foreground colour
   */
  void DrawF(Double_t x, Double_t y, Double_t w, Double_t a, Int_t fg)
  {
    Double_t u = 3.5;
    TBox* b1 = new TBox(x, y, x + a * w, y + .99);
    SetAttr(b1, b1, fg);
    b1->Draw();
    b1->DrawBox(x, y+.99-w,      x + u*a*w,      y+.99);
    b1->DrawBox(x, y+.5-w/2+.02, x + (u-.1)*a*w, y+.5+w/2+.02);
  }
  /** 
   * Draw the label text ("Preliminary" or "Performance") in it's own pad 
   * 
   * @param fg Foreground colour 
   * @param bg Background colour 
   * @param preliminary If true, draw preliminary, otherwise performance
   */
  void DrawLabel(Int_t fg, Int_t bg, bool preliminary=true)
  {
    printf("Drawing stamp ");
    fPad->cd();

    Double_t top = .08;
    Double_t bot = 0;
    TVirtualPad* pad = new TPad("label", "label", 0, bot, 1, top);
    pad->SetFixedAspectRatio(true);
    pad->SetFillColor(fDebug ? kGreen-4 : 0);
    pad->SetFillStyle(fDebug ? 4030     : 0);
    pad->SetBorderSize(fDebug ? 1 : 0);
    pad->SetBorderMode(fDebug ? 1 : 0);
    pad->SetLineColor(fDebug ? kCyan+2 : 0);
    pad->SetLineWidth(fDebug ? 1 : 0);
    pad->SetTopMargin(0);
    pad->SetBottomMargin(0);
    pad->SetLeftMargin(0);
    pad->SetRightMargin(0);
    pad->Draw();
    pad->cd();

    Double_t wp = pad->XtoPixel(1);
    Double_t hp = pad->YtoPixel(0);
    Double_t a  = hp / wp;
    
    Double_t w = .15;
    if (preliminary) {
      printf("."); DrawP(0,            0.01, w, a, fg, bg);
      printf("."); DrawR(6.3  * w * a, 0.01, w, a, fg, bg);
      printf("."); DrawE(13   * w * a, 0.01, w, a, fg);
      printf("."); DrawL(18.7 * w * a, 0.01, w, a, fg);
      printf("."); DrawI(23.3 * w * a, 0.01, w, a, fg);
      printf("."); DrawM(29   * w * a, 0.01, w, a, fg);
      printf("."); DrawI(37.2 * w * a, 0.01, w, a, fg);
      printf("."); DrawN(42.5 * w * a, 0.01, w, a, fg);
      printf("."); DrawA(49.4 * w * a, 0.01, w, a, fg, bg);
      printf("."); DrawR(56.4 * w * a, 0.01, w, a, fg, bg);
      printf("."); DrawY(62.3 * w * a, 0.01, w, a, fg);
    }
    else { 
      printf("."); DrawP(0,            0.01, w, a, fg, bg);
      printf("."); DrawE(6.3  * w * a, 0.01, w, a, fg);
      printf("."); DrawR(12   * w * a, 0.01, w, a, fg, bg);
      printf("."); DrawF(18.7 * w * a, 0.01, w, a, fg);
      printf("."); DrawO(24   * w * a, 0.01, w, a, fg, bg);
      printf("."); DrawR(30.7 * w * a, 0.01, w, a, fg, bg);
      printf("."); DrawM(37.5 * w * a, 0.01, w, a, fg);
      printf("."); DrawA(45.5 * w * a, 0.01, w, a, fg, bg);
      printf("."); DrawN(52.4 * w * a, 0.01, w, a, fg);
      printf("."); DrawC(58.7 * w * a, 0.01, /*w,*/ a, fg, bg);
      printf("."); DrawE(64.8 * w * a, 0.01, w, a, fg);
    }

    fPad->cd();
    printf("\n"); 
  }    
  /** 
   * Draw the letter O
   * 
   * @param x    Base X coordinate 
   * @param y    Base Y coordinate 
   * @param w    Line width 
   * @param a    Aspect ratio of pad 
   * @param fg   Foreground colour
   * @param bg   Background colour
   */
  void DrawO(Double_t x, Double_t y, Double_t w, Double_t a, Int_t fg, Int_t bg)
  {
    Double_t fw = 4.8;

    TEllipse* a1 = new TEllipse(x+fw/2*a*w, y + .495, a*w*fw/2, .495, 0, 360);
    SetAttr(a1, a1, fg);
    a1->Draw("only");
    TEllipse* a2 = new TEllipse(x+fw/2*a*w, y + .495, a*w*(fw/2-1), .495-w);
    SetAttr(a2, a2, bg);
    a2->Draw("only");
  }
  /** 
   * Draw the letter P 
   * 
   * @param x    Base X coordinate 
   * @param y    Base Y coordinate 
   * @param w    Line width 
   * @param a    Aspect ratio of pad 
   * @param fg   Foreground colour
   * @param bg   Background colour
   */
  void DrawP(Double_t x, Double_t y, Double_t w, Double_t a, Int_t fg, Int_t bg)
  {
    Double_t u  = 4.3 * w;
    Double_t r  = 4.3 * w / 2;

    TBox* b1 = new TBox(x, y, x + a * w, y + .99);
    SetAttr(b1, b1, fg);
    b1->Draw();
    b1->DrawBox(x + a * w, y + .99 - w, x + 2.5 * a * w, y + .99);
    b1->DrawBox(x + a * w, y + .99 - u, x + 2.5 * a * w, y + .99 - u + w);
    TEllipse* a1 = new TEllipse(x+2.5*a*w, y + .99 - r, a*r, r-.01, -90, 90);
    SetAttr(a1, a1, fg);
    a1->Draw("only");
    TEllipse* a2 = new TEllipse(x+2.5*a*w, y + .99 - r, a*(r-w), r-w-.01, 
				-90, 90);
    SetAttr(a2, a2, bg);
    a2->Draw("only");
  }
  /** 
   * Draw the letter R
   * 
   * @param x    Base X coordinate 
   * @param y    Base Y coordinate 
   * @param w    Line width 
   * @param a    Aspect ratio of pad 
   * @param fg   Foreground colour
   * @param bg   Background colour
   */
  void DrawR(Double_t x, Double_t y, Double_t w, Double_t a, Int_t fg, Int_t bg)
  {
    DrawP(x, y, w, a, fg, bg);

    Double_t u   = 4.3 * w;
    Double_t fw  = 2.5 + 4.3 / 2;
    Double_t ang = 60 * TMath::Pi() / 180;
    Double_t o   = w / TMath::Sin(ang);

    TGraph* g = new TGraph(5);
    SetAttr(g, g, fg);
    // Top right
    g->SetPoint(0, x+a*(2*w+o), y+1-u+.01);
    // Top left 
    g->SetPoint(1, x+a*2*w,     y+1-u+.01);
    // Bottom left 
    g->SetPoint(2, x+a*(fw*w-o+.06), y);
    // Bottom right
    g->SetPoint(3, x+a*(fw*w+.06),     y);
    g->SetPoint(4, g->GetX()[0], g->GetY()[0]);
    g->Draw("same lf");
  }
  /** 
   * Draw the letter M
   * 
   * @param x    Base X coordinate 
   * @param y    Base Y coordinate 
   * @param w    Line width 
   * @param a    Aspect ratio of pad 
   * @param fg   Foreground colour
   */
  void DrawM(Double_t x, Double_t y, Double_t w, Double_t a, Int_t fg) 
  {
    TBox* b1 = new TBox(x, y, x + a * w, y + 1);
    SetAttr(b1, b1, fg);
    b1->Draw();
    b1->DrawBox(x+4.8*a*w, y, x+5.8*a*w, y+1);
    
    Double_t ang = 28 * TMath::Pi() / 180;
    Double_t o = w / TMath::Sin(ang);

    TGraph* gl = new TGraph(7);
    SetAttr(gl, gl, fg);
    gl->SetPoint(0, x+a*w, y+1);
    gl->SetPoint(1, x+a*w, y+1-o);
    gl->SetPoint(2, x+a*w*(1+3.8/2), y+w);
    gl->SetPoint(3, x+a*w*4.8, y+1-o);
    gl->SetPoint(4, x+a*w*4.8, y+1);
    gl->SetPoint(5, x+a*w*(1+3.8/2), y+w+o);
    gl->SetPoint(6, x+a*w, y+1);
    gl->Draw("same lf");

  }
  /** 
   * Draw the letter N 
   * 
   * @param x    Base X coordinate 
   * @param y    Base Y coordinate 
   * @param w    Line width 
   * @param a    Aspect ratio of pad 
   * @param fg   Foreground colour
   */
  void DrawN(Double_t x, Double_t y, Double_t w, Double_t a, Int_t fg) 
  {
    Double_t fw = 4.6;
    TBox* b1 = new TBox(x, y, x + a * w, y + 1);
    SetAttr(b1, b1, fg);
    b1->Draw();
    b1->DrawBox(x+a*w*(fw-1), y, x+a*w*(fw), y+1);
    
    Double_t ang = 28 * TMath::Pi() / 180;
    Double_t o = w / TMath::Sin(ang);

    TGraph* gl = new TGraph(5);
    gl->SetLineColor(fDebug ? kGreen+1 : fg);
    gl->SetLineWidth(fDebug ? 1        : 0);
    gl->SetFillColor(fDebug ? kGreen+1 : fg);
    gl->SetFillStyle(fDebug ? 0 : 1001);
    gl->SetPoint(0, x+a*w, y+1);
    gl->SetPoint(1, x+a*w, y+1-o);
    gl->SetPoint(2, x+a*w*(fw-1), y);
    gl->SetPoint(3, x+a*w*(fw-1), y+o);
    gl->SetPoint(4, x+a*w, y+1);
    gl->Draw("same lf");
  }
  /** 
   * Draw the letter Y
   * 
   * @param x    Base X coordinate 
   * @param y    Base Y coordinate 
   * @param w    Line width 
   * @param a    Aspect ratio of pad 
   * @param fg   Foreground colour
   */
  void DrawY(Double_t x, Double_t y, Double_t w, Double_t a, Int_t fg)
  {
    Double_t u  = 4 * w;
    Double_t fw = 5.5;

    TBox* b1 = new TBox(x+a*w*(fw/2-.5), y, x + a * w *(fw/2+.5), y + 1 - u);
    b1->SetLineColor(fDebug ? kGreen+1 : fg);
    b1->SetFillColor(fDebug ? 0 : fg);
    b1->SetFillStyle(fDebug ? 0 : 1001);
    b1->Draw();

    Double_t ang = 40 * TMath::Pi() / 180;
    Double_t o  = w / TMath::Sin(ang);
    Double_t o2 = w / TMath::Cos(ang);


    TGraph* g = new TGraph(8);
    SetAttr(g, g, fg);
    g->SetPoint(0, x+a*w*(fw/2-.5), y + 1-u);
    g->SetPoint(1, x,             y + 1);
    g->SetPoint(2, x+a*o2,        y + 1);
    g->SetPoint(3, x+a*w*fw/2,    y + 1 - u + o);
    g->SetPoint(4, x+a*w*fw-a*o2, y+1);
    g->SetPoint(5, x+a*w*fw,      y+1);
    g->SetPoint(6, x+a*w*(fw/2+.5), y + 1-u);
    g->SetPoint(7, x+a*w*(fw/2-.5), y + 1-u);

    
    g->Draw("same lf");
  }
  /** 
   * Set Line and Fill attributes
   * 
   * @param lne Line attribute interface 
   * @param fll Fill attribute interface 
   * @param col Colour 
   */
  void SetAttr(TAttLine* lne, TAttFill* fll, Int_t col) 
  {
    if (lne) {
      lne->SetLineColor(fDebug ? kGreen+1 : col);
      lne->SetLineWidth(fDebug ? 1        : 0);
    }
    if (fll) { 
      fll->SetFillColor(fDebug ? kGreen+1 : col);
      fll->SetFillStyle(fDebug ? 0        : 1001);
    }
  }
  Double_t GetX1() const { return fPad ? fPad->GetXlowNDC() : -1; }
  Double_t GetX2() const { return fPad ? GetX1()+fPad->GetWNDC() : -1; }
  Double_t GetY1() const { return fPad ? fPad->GetYlowNDC() : -1; }
  Double_t GetY2() const { return fPad ? GetY1()+fPad->GetHNDC() : -1; }
  Double_t GetCX() const { return fPad ? (GetX1()+GetX2())/2 : -1; }
  Double_t GetCY() const { return fPad ? (GetY1()+GetY2())/2 : -1; }
  /** 
   * run a simple test
   * 
   * @param what 
   */
  static void Test(const TString& what)
  {
    Test(AliceLogo::ParseTag(what), AliceLogo::ParseColors(what), 
	 what.Contains("debug", TString::kIgnoreCase));
  }
  /** 
   * Run a simple test
   * 
   * @param tag 
   * @param col 
   * @param debug 
   */
  static void Test(UShort_t tag, UInt_t col, Bool_t debug)
  {
    TCanvas* c = new TCanvas("test", "Test of AliceLogo");
    c->SetFillColor(kYellow-10);
    TH1* h = new TH1F("h", "h", 60, -3, 3);
    h->SetFillColor(kRed+1);
    h->SetFillStyle(3001);
    h->FillRandom("gaus");
    h->Draw();
    
    AliceLogo* dut = new AliceLogo(debug);
    dut->Draw(c, .75, .4, .25, tag, col);

    c->SaveAs("test.png");
    c->SaveAs("test.pdf");
    c->SaveAs("test.svg");
  }
  /** 
   * Check compatiblity by overlaying this logo on 'official' images
   * 
   * @param what 
   */
  static void Check(const TString& what)
  {
    Check(AliceLogo::ParseTag(what), AliceLogo::ParseColors(what), 
	 what.Contains("debug", TString::kIgnoreCase));
  }
  /** 
   * Check compatiblity by overlaying this logo on 'official' images
   * 
   * @param tag 
   * @param col 
   * @param debug 
   */
  static void Check(UShort_t tag, UInt_t col, Bool_t debug=true)
  {
    
    TImage* img = TImage::Create();
    TString imgName = "2012-Jul-04-Base_Logo.png";
    // Int_t   logoW   = 0; // img->GetWidth();
    // Int_t   logoH   = 0; // img->GetHeight();
    switch (tag) { 
    case AliceLogo::kPreliminary: 
      // logoW   = 1849;
      // logoH   = 2277;
      imgName = "2012-Jul-04-Preliminary_Logo.png"; 
      break;
    case AliceLogo::kPerformance: 
      imgName = "2012-Jul-04-Performance_Logo.png";
      break;
    default:
      if (col & kMany) 
	imgName = "2012-Jul-04-4_Color_Logo_CB.png";
      else
	imgName = "2012-Jul-04-Base_Logo.png";
    }
    Printf("Reading %s", imgName.Data());
    img->ReadImage(imgName);
    Int_t    logoW = img->GetWidth();
    Int_t    logoH = img->GetHeight();
    Double_t logoA = Float_t(logoH) / logoW; 

    Printf("Logo dimensions (wxh): (%d x %d)", logoW, logoH);
    Int_t myHeight = 800;
    TCanvas* c = new TCanvas("check", "Check of AliceLogo", 
			     (myHeight-22)/logoA, myHeight);
    c->SetTopMargin(0);
    c->SetRightMargin(0);
    c->SetBottomMargin(0);
    c->SetLeftMargin(0);
    c->SetFillColor(debug ? kRed-4 : kWhite);
    c->SetBorderSize(0);
    c->SetBorderMode(0);
    if (col & AliceLogo::kReverse /*&& !(col & AliceLogo::kL3)*/ && !debug) 
      c->SetFillColor(kBlack);
    
    if (debug) img->Draw();

    AliceLogo* dut = new AliceLogo(debug);
    dut->Draw(c, 0, 0, 1, tag, col);

    c->SaveAs("check.svg");
  }

  TVirtualPad* fPad;      // Our pad 
  Double_t     fAspect;   // Aspect ratio of our pad 
  Double_t     fL3Thick;  // Fraction of canvas height
  Double_t     fL3Height; // Fraction of canvas height 
  Double_t     fL3Width;  // Width of L3 sub-pad
  Bool_t       fDebug;    // Debug (make outline only)
  Int_t        fRed;      // The red color 
  Int_t        fYellow;   // The yellow-ish color
  Int_t        fPink;     // The dark-red pink-ish color 
  Int_t        fBlue;     // The blue gray-ish color 
};

#endif
//
// EOF
//
