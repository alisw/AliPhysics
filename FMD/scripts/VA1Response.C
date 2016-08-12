//____________________________________________________________________
// 
// $Id$
//
// Script to try to fit the reponse function of the VA1 signals, based
// on a finite number of ALTRO samples. 
//
/** Make Va1 response
    @param n 
    @param B 
    @param dc 
    @param errors 
    @ingroup FMD_simple_script
*/
void 
VA1Response(Int_t n=2, Float_t B=5, Float_t dc=.01, Bool_t errors=kFALSE) 
{

  TF1* response = new TF1("response", "[0] * (1 - exp(-[1] * x))", 0, 1.4);
  response->SetParameters(1, B);
  response->SetParNames("A", "B");
  response->SetLineColor(2);
  
  TF1* fit = new TF1("fit",  "[0] * (1 - exp(-[1] * x))",  0, 1);
  fit->SetParameters(.5, B/2);
  fit->SetParNames("A", "B");
  fit->SetLineColor(3);
  
  TGraph* graph = 0;
  if (errors) graph = new TGraphErrors(n);
  else        graph = new TGraph(n);
  for (Int_t i = 0; i < n; i++) {
    Float_t t = Float_t(i + 1) / n;
    Float_t c = gRandom->Gaus(response->Eval(t), dc);
    graph->SetPoint(i, t, c);
    if (errors) ((TGraphErrors*)graph)->SetPointError(i, 0, dc);
  }
  
  response->Draw();
  response->GetHistogram()->GetYaxis()->SetRangeUser(0, 1.4);
  response->GetHistogram()->GetXaxis()->SetRangeUser(0, 1.4);
  graph->Draw("P*");
  TString fitOpt("E");
  if (!errors) fitOpt.Append("W");
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
//____________________________________________________________________
//
// EOF
//
