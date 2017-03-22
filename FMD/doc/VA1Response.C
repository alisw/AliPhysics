//____________________________________________________________________
// 
// $Id: VA1Response.C 13249 2006-03-24 16:09:36Z cholm $
//
// Script to try to fit the reponse function of the VA1 signals, based
// on a finite number of ALTRO samples. 
//
/** Make Va1 response
    @param n 
    @param B 
    @param dc 
    @param errors 
    @ingroup simple_script
*/
void 
VA1Response(Int_t n=4, Float_t B=6, Float_t dc=.01, Bool_t errors=kFALSE,
	    Bool_t doFit=kFALSE) 
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetLabelFont(132, "xyz");
  gStyle->SetTitleFont(132, "xyz");
  gStyle->SetTitleSize(0.08, "y");
  gStyle->SetTitleOffset(0.5, "y");
  gStyle->SetTitleSize(0.06, "x");
  gStyle->SetTitleOffset(0.7, "x");

  TCanvas* c = new TCanvas("c", "c", 800, 500);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(0);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  c->SetGridx();
  c->SetGridy();

  TF1* response = new TF1("response", "[0] * (1 - exp(-[1] * x))", 0, 1.4);
  response->SetParameters(1, B);
  response->SetParNames("A", "B");
  response->SetLineColor(2);
  
  TGraph* graph = 0;
  if (n >= 2) {
    if (errors) graph = new TGraphErrors(n);
    else        graph = new TGraph(n);
    for (Int_t i = 0; i < n; i++) {
      Float_t t = Float_t(i + 1) / n;
      Float_t q = gRandom->Gaus(response->Eval(t), dc);
      graph->SetPoint(i, t, q);
      if (errors) ((TGraphErrors*)graph)->SetPointError(i, 0, dc);
    }
  }

  response->Draw();
  response->GetHistogram()->GetYaxis()->SetRangeUser(0, 1.05);
  response->GetHistogram()->GetXaxis()->SetRangeUser(0, 1.1);
  response->GetHistogram()->GetXaxis()->SetNdivisions(6, kTRUE);
  response->GetHistogram()->GetYaxis()->SetNdivisions(10, kTRUE);
  response->GetHistogram()->SetXTitle("t");
  response->GetHistogram()->SetYTitle(Form("1-e^{-%3.1f t}", B));
  
  if (graph) {
    
    
    graph->Draw("P*");
    TString fitOpt("E");
    if (!errors) fitOpt.Append("W");

    if (doFit) {
      TF1* fit = new TF1("fit",  "[0] * (1 - exp(-[1] * x))",  0, 1);
      fit->SetParameters(.5, B/2);
      fit->SetParNames("A", "B");
      fit->SetLineColor(3);
      graph->Fit("fit", fitOpt.Data());
      graph->Fit("fit", fitOpt.Data());
      
      std::cout << "Chi^2/NDF = " << fit->GetChisquare() << "/" << fit->GetNDF()
		<< " = " << std::flush;
      if (fit->GetNDF() == 0) 
	std::cout << " undefined!" << std::endl;
      else
	std::cout << (fit->GetChisquare() / fit->GetNDF()) << std::endl;
      std::cout << "f(t) = " 
		<< fit->GetParameter(0) << "+/-" << fit->GetParError(0) 
		<< " * (1 - exp(" 
		<< fit->GetParameter(1) << "+/-" << fit->GetParError(1) 
		<< " * t))" << std::endl;
    }
  }

  c->Modified();
  c->Update();
  c->cd();
  c->SaveAs("va1_response.png");
}
//____________________________________________________________________
//
// EOF
//
