//
// This macro is a part of "Alice PPR fast flow analysis package"
// 
// The macro Draw event plane resolution and V2 resolution
// for a given multiplicity in function of V2
//
// INPUT PARAMETERS
// type: 
// 0 - event plane resolution
// 1 - V2 reslution
// 
// Sylwester Radomski, GSI
// mail: S.Radomski@gsi.de
// 31 Oct 2002
//

AliFlowDrawV2(int type, const char *fileName = 0) {

  gROOT->SetStyle("Plain");

  if (type > 2) {
    ::Error("AliFLowDrawV2", "Wrong Type [0-1]: %d ", type);
    return;
  }
    
  const char *dataName = "flowPicoEvent.root";
  TFile *dataFile = new TFile(dataName,"UPDATE");
  TNtuple *data = (TNtuple *) dataFile->Get("flowData");

  TGraphErrors *g = new TGraphErrors();
  TH1D *h;

  const char* patt = "Mult == %d && trueV2 > %f && trueV2 < %f";
  const char* varName[2] = {"Psi - truePsi >> htemp(90, -90, 90)",
			    "100 * (trueV2-V2)/trueV2 >> htemp(100,-200,200)"};

  const char* xTitle = "V_{2}";
  const char* yTitle[2] = {"Event Plane Resolution [deg]",  
			   "V_{2} Resolution [%]"};

  char buff[100];

  Int_t mult = 1000;
  const Int_t nV2Point = 5;  
  const Double_t v2Points[nV2Point] = {0.02, 0.04, 0.06, 0.08, 0.1};

  TCanvas *c = new TCanvas();
  c->SetFillColor(10);
  c->SetTicks(1,1);

  TCanvas *d = new TCanvas();
  
  for(Int_t i=0; i<nV2Point; i++) {
    
    Double_t totalV2 = v2Points[i];
    
    sprintf(buff, patt, mult, totalV2-0.001, totalV2+0.001);
    data->Draw(varName[type], buff);
    h = (TH1D*)gPad->GetPrimitive("htemp"); 
    
    ::Info("AliFlowDrawV2",buff);

    g->SetPoint(i, v2Points[i], h->GetRMS());
    g->SetPointError(i, 0, h->GetRMS()/20. );
  }
  
  g->SetMinimum(0);
  g->SetMarkerSize(1);
  g->SetMarkerStyle(21);
  
  c->cd();
  g->Draw("APL");
  g->GetHistogram()->SetXTitle(xTitle);
  g->GetHistogram()->SetYTitle(yTitle[type]);
  
  c->Modified();
  c->Update();

  if (fileName) {
    
    char buffer[60];
    sprintf(buffer, "plots/%s.eps", fileName);
    c->Print(buffer, "eps");
    
    sprintf(buffer, "plots/%s.gif", fileName);
    c->Print(buffer, "gif");
  }  
}
