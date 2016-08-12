//____________________________________________________________________
//
// $Id$
//
// Small script that shows a signal train from a VA1 pre-amp. 
// 
/** Make VA1 sample train
    @ingroup FMD_simple_script
 */
void 
VA1Train() 
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas("c", "C", 800, 400);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(0);
  

  TArrayF measurements(6);
  std::cout << "Measurements are: " << std::flush;
  for (Int_t i = 0; i < measurements.fN; i++) {
    measurements[i] = gRandom->Uniform(0,1);
    std::cout << measurements[i] << " " << std::flush;
  }
  std::cout << std::endl;

  Float_t last = 0;
  Float_t B    = 5;
  TH2* frame = new TH2F("frame", "Frame", measurements.fN, 0, 
			measurements.fN, 10, 0, 1.1);
  frame->Draw();
  for (Int_t i = 0; i < measurements.fN; i++) {
    TF1* f = new TF1("f", "[2] + exp(-[1] * (x - [3])) * ([0] - [2])", 
		     i, i + 1);
    f->SetParameter(3, i);
    f->SetParameter(1, B);
    f->SetParameter(2, measurements[i]);
    f->SetParameter(0, last);
    
    if (measurements[i] > last) {      
      // f = new TF1("f", "[0] * (1 - exp(-[1] * (x - [2]))) + [3]", i, i+1);
      // f->SetParameters(measurements[i] - last, B, i, last);
      // f->SetParameter(0, measurements[i]);
      // f->SetParameter(2, last);
    }
    else {
      // f->SetParameter(2, measurements[i]);
      // f->SetParameter(0, last);
      // f = new TF1("f", "[0] * (exp(-[1] * (x - [2]))) + [3]", i, i+1);
      // f->SetParameters(last - measurements[i], B, i, measurements[i]);
    }
    f->Draw("same");  
    last = measurements[i]; 
  }
}
//____________________________________________________________________
//
// EOF
//
