#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TError.h"

namespace Reve{
class TTriangleSet;
}

Reve::TriangleSet *ts1=0, *ts2=0, *ts3=0;

void triangleset()
{
  {
    ts1 = Reve::TriangleSet::ReadTrivialFile("broken_torus.tring");
    ts1->SetName("RandomColors");
    ts1->GenerateTriangleNormals();
    ts1->GenerateRandomColors();
    ts1->SetColor(0);
    gReve->AddRenderElement(ts1);
  }
  {
    ts2 = Reve::TriangleSet::ReadTrivialFile("broken_torus.tring");
    ts2->SetName("SmallBlue");
    ts2->GenerateTriangleNormals();
    ts2->SetColor(4);
    TGeoHMatrix m;
    m.RotateY(90);
    Double_t scale[3] = { 0.8, 0.8, 1.2 };
    m.SetScale(scale);
    ts2->SetTransMatrix(m);
    gReve->AddRenderElement(ts2);
  }
  {
    ts3 = Reve::TriangleSet::ReadTrivialFile("broken_torus.tring");
    ts3->SetName("Spectrum");
    ts3->GenerateTriangleNormals();
    gStyle->SetPalette(1, 0);
    ts3->GenerateZNormalColors(50, -50, 50, kTRUE, kTRUE);
    ts3->SetColor(0);
    TGeoHMatrix m;
    m.RotateZ(90);
    Double_t scale[3] = { 1.3, 1.0, 1.6 };
    m.SetScale(scale);
    ts3->SetTransMatrix(m);
    gReve->AddRenderElement(ts3);
  }

  gReve->Redraw3D(kTRUE);
}
