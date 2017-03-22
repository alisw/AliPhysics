//____________________________________________________________________
//
// $Id: VA1Train.C 13249 2006-03-24 16:09:36Z cholm $
//
// Small script that shows a signal train from a VA1 pre-amp. 
// 
/** Make VA1 sample train
    @ingroup simple_script
 */
void 
VA1Train() 
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetLabelFont(132, "xyz");
  gStyle->SetTitleFont(132, "xyz");
  gStyle->SetTitleSize(0.08, "y");
  gStyle->SetTitleOffset(0.5, "y");
  gStyle->SetTitleSize(0.06, "x");
  gStyle->SetTitleOffset(0.7, "x");

  TCanvas* c = new TCanvas("c", "C", 800, 500);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(0);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  c->SetGridx();
  c->SetGridy();

  TArrayF measurements(6);
  std::cout << "Measurements are: " << std::flush;
  for (Int_t i = 0; i < measurements.fN; i++) {
    measurements[i] = gRandom->Uniform(0,1);
    std::cout << measurements[i] << " " << std::flush;
  }
  std::cout << std::endl;

  Float_t last = 0;
  Float_t B    = 5;
  TH2* frame = new TH2F("frame", "Frame", measurements.fN, -.1, 
			measurements.fN+.1, 10, 0, 1.05);
  frame->Draw();
  // frame->GetYaxis()->SetRangeUser(0, 1.05);
  // frame->GetXaxis()->SetRangeUser(-.1, measurements.fN+.1);
  frame->GetXaxis()->SetNdivisions(measurements.fN+2, kTRUE);
  frame->GetYaxis()->SetNdivisions(10, kTRUE);
  frame->SetXTitle("t");
  frame->SetYTitle(Form("1-e^{-%3.1f t}", B));

  for (Int_t i = 0; i < measurements.fN; i++) {
    TF1* f = new TF1("f", "[2] + exp(-[1] * (x - [3])) * ([0] - [2])", 
		     i, i + 1);
    f->SetParameter(3, i);
    f->SetParameter(1, B);
    f->SetParameter(2, measurements[i]);
    f->SetParameter(0, last);
    
    f->Draw("same");  
    f->SetLineColor(kRed+i/2);
    last = measurements[i]; 
  }

  TLatex* l = new TLatex(measurements.fN/2, .9, 
			 "#Delta_{i} + (#Delta_{i-1}-#Delta_{i}) e^{-b (t-i)}");
  l->SetTextAlign(22);
  l->SetTextFont(132);
  l->Draw();

  c->SaveAs("va1_train.png");
}
//____________________________________________________________________
//
// EOF
//
