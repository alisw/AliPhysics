// $Id$
// $MpId: testSectorAreaIterator.C,v 1.5 2005/10/28 15:37:12 ivana Exp $
//
// Test macro for iterating over the whole sector

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpSector.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpSectorReader.h"
#include "AliMpArea.h"
#include "AliMpVPadIterator.h"
#include "AliMpVPainter.h"

#include <Riostream.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TH2.h>
#include <TStopwatch.h>

#endif

class AliMpVPadIterator;

TCanvas* CreateTCanvas(const TString& name, const TString& title,
                       AliMq::Station12Type station, AliMp::PlaneType plane)
{
  TString newName(name);
  TString newTitle(title);
  TString unique = AliMq::Station12TypeName(station) + AliMp::PlaneTypeName(plane);
  newName += unique;
  newTitle += unique;
  return new TCanvas(newName.Data(), newTitle.Data());
}                     

void MarkPads(AliMpVPadIterator& it, Double_t xmax, Double_t ymax, 
              AliMq::Station12Type station, AliMp::PlaneType plane,
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
    histo->Fill(it.CurrentItem().GetIx(), it.CurrentItem().GetIy());   		
  }
  
  TCanvas* canv2 = CreateTCanvas("canv2 ", " ", station, plane);
  canv2->cd();
  //histo->SetMinimum(1.5);
  histo->Draw("box");

  timer.Stop();
  //timer.Print();
}

void testSectorAreaIterator(AliMq::Station12Type station, AliMp::PlaneType plane)
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSectorReader r(dataStreams, station, plane);
  AliMpSector* sector = r.BuildSector();
  AliMpSectorSegmentation segmentation(sector);

  AliMpArea area;
  if ( station == AliMq::kStation1 )
    area = AliMpArea(TVector2(45.,45.),TVector2(45.,45.));
  else   
    area = AliMpArea(TVector2(60.,60.),TVector2(60.,60.));
  AliMpVPadIterator* iter = segmentation.CreateIterator(area);

  CreateTCanvas("Graph ", " ", station, plane);
  AliMpVPainter::CreatePainter(sector)->Draw("ZSSMP");

  TCanvas* canv = CreateTCanvas("canv ", " ", station, plane);
  canv->Range(-1,-1,1,1);
  MarkPads(*iter, TMath::Abs(area.Position().X())+area.Dimensions().X(),
                  TMath::Abs(area.Position().Y())+area.Dimensions().Y(), 
                  station, plane, kTRUE);
}
     
void testSt12SectorAreaIterator()
{
  AliMq::Station12Type  station[2] = { AliMq::kStation1, AliMq::kStation2 }; 
  AliMp::PlaneType      plane[2]   = { AliMp::kBendingPlane, AliMp::kNonBendingPlane };
  
  for ( Int_t is = 0; is < 2; is++ ) {
    for ( Int_t ip = 0; ip < 2; ip++ ) {
    
      cout << "Running testSectorAreaIterator for " 
           << AliMq::Station12TypeName(station[is]) << "  "
           << AliMp::PlaneTypeName(plane[ip])  << " ... " << endl;
       
      testSectorAreaIterator(station[is], plane[ip]);
    
      cout << "... end running " << endl << endl;
    }  
  }   
}  
  
