//
// This macro is a part of "Alice PPR fast flow analysis package"
//
// The macro draws event plane resolution and V2 resolution in function
// of V2 and multiplicity
// 
// IMPUT PARAMETERS
// type:
// 0 - event plane resolution
// 1 - V2 resolution
//
// desc 0/1 - draw Graph descriptions
//
// Sylwester Radomski, GSI
// mail: S.Radomski@gsi.de
// Oct 31, 2002
//

AliFlowDrawSummary(int type = 0, int desc = 0, const char* fileName = 0) {

  gROOT->SetStyle("Plain");
  
  if (type > 1) {
    ::Error("AliFlowDrawSummary","Wrong Type [0-1] : %d", type);
    return;
  }

  const Int_t startMarker = 21;
  const char *dataName = "flowPicoEvent.root";
  TFile *dataFile = new TFile(dataName,"UPDATE");
  TNtuple *data = (TNtuple *) dataFile->Get("flowData");

  TGraphErrors *g;
  TH1D *h;

  const char* patt = "Mult == %d && trueV2 > %f && trueV2 < %f";
  
  const char* varName[2] = {"Psi - truePsi >> htemp(90, -90, 90) ",
                            "100 * (trueV2-V2)/trueV2 >> htemp(100,-200,200)"};

  const char* yTitle[2] = {"Event Plane Resolution [deg]",
			   "V_{2} Resolution [%]"};

  char buff[100];

  Double_t v2 = 0.02;
  
  const Int_t nMultPoint = 11;
  const Int_t nV2Point = 4;
  
  const Int_t multPoints[nMultPoint] = {200, 400, 600, 800, 1000, 1500, 2000, 2500, 3000,
					  4000, 5000};// 6000, 8000, 10000};
  
  const Double_t v2Points[nV2Point] = {0.04, 0.06, 0.08, 0.1};
  //const Double_t v2Points[nV2Point] = {0.06};

  TCanvas *c = new TCanvas();
  c->SetFillColor(10);
  c->SetTicks(1,1);

  TCanvas *d = new TCanvas();

  for (Int_t j=0; j<nV2Point; j++) {
    
    g = new TGraphErrors();
    g->SetMarkerStyle(18+j); 
    
    for(Int_t i=0; i<nMultPoint; i++) {
      
      Int_t totalMult = multPoints[i];
      Double_t totalV2 = v2Points[j];

      sprintf(buff, patt, totalMult, totalV2-0.001, totalV2+0.001);
      data->Draw(varName[type], buff);
      h = (TH1D*)gPad->GetPrimitive("htemp"); 
      
      ::Info("AliFlowDrawSummary", buff);

      g->SetPoint(i, multPoints[i], h->GetRMS());
      g->SetPointError(i, 0, h->GetRMS() / TMath::Sqrt(h->GetEntries()));
    }
    
    g->SetMinimum(0);
    g->SetMarkerSize(1);
    g->SetMarkerStyle(startMarker+j);

    c->cd();
    if (j==0) {
       g->Draw("APL");
       g->GetHistogram()->SetXTitle("Number of Good Tracks");
       g->GetHistogram()->SetYTitle(yTitle[type]);
    }
    else g->Draw("PL");
    d->cd();
  }

  if (desc) {

    TText *text;
    TMarker *marker;
    c->cd();
    
    const char* info = "v2 = %.2f";
    
    for (Int_t j=0; j<nV2Point; j++) {
      
      Double_t y;
      if (type == 0) y = 30 - 3*y; 
      else y = 70 - 5*j;
      
      marker = new TMarker(3000, y+1, startMarker+j);
      marker->Draw("a");
      marker->ls();
      
      sprintf(buff, info, v2Points[j]);
      text = new TText(3100, y, buff);
      //text->SetTextColor(j-1);
      text->Draw();
    }
  }

  if (fileName) {
    
    char buffer[60];
    sprintf(buffer, "plots/%s.eps", fileName);
    c->Print(buffer, "eps");
    
    sprintf(buffer, "plots/%s.gif", fileName);
    c->Print(buffer, "gif");
  }
}
