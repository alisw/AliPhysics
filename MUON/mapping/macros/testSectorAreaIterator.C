// $Id$
// $MpId: testSectorAreaIterator.C,v 1.5 2005/10/28 15:37:12 ivana Exp $
//
// Test macro for iterating over the whole sector

#include <iomanip>

class AliMpVPadIterator;

void MarkPads(AliMpVPadIterator& it, Double_t xmax, Double_t ymax, 
              Bool_t print = kTRUE)
{
// Marks pads according their position.
// Fills histogram with pad indices.
// Measures time that takes processing of full quadrant 
// ---

  Int_t num=0;

  TH2C* histo = new TH2C("pads", "pads", 201, 0, 200, 251, 0, 250); 

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

void testSectorAreaIterator(AliMp::StationType station = AliMp::kStation1,
                            AliMp::PlaneType plane = AliMp::kBendingPlane,
	     	            Bool_t rootInput = false)
{
  AliMpSector *sector = 0;
  if (!rootInput) {
    AliMpSectorReader r(station, plane);
    sector=r.BuildSector();
  }
  else  {
    TString filePath = AliMpFiles::SectorFilePath(station,plane);
    filePath.ReplaceAll("zones.dat", "sector.root"); 

    TFile f(filePath.Data(), "READ");
    sector = (AliMpSector*)f.Get("Sector");
  }  

  AliMpSectorSegmentation segmentation(sector);

  AliMpArea area;
  if ( station == AliMp::kStation1 )
    area = AliMpArea(TVector2(45.,45.),TVector2(45.,45.));
  else   
    area = AliMpArea(TVector2(60.,60.),TVector2(60.,60.));
  AliMpVPadIterator* iter = segmentation.CreateIterator(area);

  TCanvas* graph = new TCanvas("Graph");
  AliMpVPainter::CreatePainter(sector)->Draw("ZSSMP");

  TCanvas *canv = new TCanvas("canv");
  canv->Range(-1,-1,1,1);
  MarkPads(*iter, TMath::Abs(area.Position().X())+area.Dimensions().X(),
                  TMath::Abs(area.Position().Y())+area.Dimensions().Y(), kTRUE);
}
