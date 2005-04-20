// $Id$
//
// Test macro for iterating over the whole plane

#include <iomanip>

class AliMpVPadIterator;

void MarkPads(AliMpVPadIterator& it, Double_t xmax, Double_t ymax, 
              Bool_t print = kTRUE)
{
// Marks pads according their position.
// Fills histogram with pad indices.
// Measures time that takes processing of full plane. 
// ---

  Int_t num=0;

  TH2C* histo = new TH2C("pads", "pads", 401, -200, 200, 501, -250, 250); 

  TStopwatch timer;
  timer.Start();  

  for (it.First(); ! it.IsDone(); it.Next()){
   
    if (print) cout << endl 
                    << setw(5) << ++num 
                    << " " << it.CurrentItem() << endl;   
    
    // mark pads positions
    TVector2 posi = it.CurrentItem().Position();
    TMarker* marker = new TMarker( posi.X()/xmax, posi.Y()/ymax, 2);
    marker->Draw();

    // fill pads indices in the histogram
    histo->Fill(it.CurrentItem().GetIndices().GetFirst(),
                it.CurrentItem().GetIndices().GetSecond());   		
  }
  
  TCanvas *canv2 = new TCanvas("canv2");
  canv2->cd();
  //histo->SetMinimum(1.5);
  histo->Draw("box");

  timer.Stop();
  //timer.Print();
}

void testPlaneAreaIterator(AliMpStationType station = kStation1,
                  AliMpPlaneType planeType = kBendingPlane)
{
  AliMpPlane* plane = AliMpPlane::Create(station, planeType);
  AliMpPlaneSegmentation planeSeg(plane);

  AliMpArea area;
  if ( station == kStation1 )
    area = AliMpArea(TVector2(0.,0.),TVector2(900.,900.));
  else   
    area = AliMpArea(TVector2(0.,0.),TVector2(1200.,1200.));
  AliMpVPadIterator* iter = planeSeg.CreateIterator(area);

  TCanvas* graph = new TCanvas("Graph");
  graph->Divide(2);
  graph->cd(1);
  AliMpVPainter::CreatePainter(plane->GetFrontSector())->Draw("ZSSMP");
  graph->cd(2);
  AliMpVPainter::CreatePainter(plane->GetBackSector())->Draw("ZSSMP");
  TCanvas *canv = new TCanvas("canv");
  canv->Range(-1,-1,1,1);
  
  MarkPads(*iter, TMath::Abs(area.Position().X())+area.Dimensions().X(),
                  TMath::Abs(area.Position().Y())+area.Dimensions().Y(), kTRUE);
}
