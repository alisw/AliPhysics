#include <iomanip>

void 
PrintOne(Char_t r)
{
  AliFMDRing* ring = AliFMDGeometry::Instance()->GetRing(r);
  if (!ring) { 
    std::cerr << "Ring " << r << " not found" << std::endl;
    return;
  }
  const TObjArray& vertices = ring->GetVerticies();
  Double_t x[] = { -1, -1, -1, -1, -1, -1 };
  Double_t y[] = { -1, -1, -1, -1, -1, -1 };
  for (int i = 1; i < 4; i++) { 
    TVector2* v = static_cast<TVector2*>(vertices.At(i));
    int       j = i-1;
    x[j]        = v->X();
    y[j]        = v->Y();
  }
  for (int i = 3; i > 0; i--) { 
    TVector2* v= static_cast<TVector2*>(vertices.At(i));
    int       j = 3+(3-i);
    x[j]        = -v->X();
    y[j]        = v->Y();
  }
  std::cout << "  Double_t r" << r << "[][2] = { ";
  TGraph* g = new TGraph(7);
  g->SetFillColor(kGray);
  g->SetFillStyle(3001);
  for (int i = 0; i < 6; i++)  {
    std::cout << "{" << std::setw(12) << x[i] 
	      << "," << std::setw(12) << y[i] << "}";
    if (i != 5) 
      std::cout << ",   // " << i << "\n                       ";
    else 
      std::cout << " }; // " << i << std::endl;
    g->SetPoint(i, x[i], y[i]);
  }
  g->SetPoint(6, x[0], y[0]);

  // TCanvas* c = new TCanvas(Form("c%c", r));
  // c->cd();
  // g->Draw("afl");
}

void
PrintSensorVertices()
{
  AliFMDGeometry::Instance()->Init();
  
  PrintOne('I');
  PrintOne('O');
}
