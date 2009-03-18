// $Id$
// $MpId: testNeighboursPadIterator.C,v 1.10 2005/10/28 15:36:08 ivana Exp $
//
// Test macro for reading  sector, and iterate over it

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpSector.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpSectorReader.h"
#include "AliMpPad.h"
#include "AliMpArea.h"
#include "AliMpNeighboursPadIterator.h"

#include <Riostream.h>
#include <TCanvas.h>
#include <TMarker.h>

#endif

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

void testNeighboursPadIterator(AliMq::Station12Type station, AliMp::PlaneType plane,
                               Int_t i=50, Int_t j=50)
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSectorReader r(dataStreams, station, plane);
  AliMpSector* sector = r.BuildSector();
  AliMpSectorSegmentation segm(sector);  
  
  CreateTCanvas("canv ", "", station, plane);

  const Double_t xmax=75;
  const Double_t ymax=120;
  
  Int_t num=0;

  AliMpPad pad = segm.PadByIndices(AliMpIntPair(i,j));
  AliMpNeighboursPadIterator it = AliMpNeighboursPadIterator(&segm, pad,kFALSE);
  
  for (it.First(); ! it.IsDone(); it.Next()) {
    AliMpIntPair indices = it.CurrentItem().GetIndices();
    cout<<"Iterator number "<< num++ << " at "<< indices <<endl;
    TMarker* marker = new TMarker( (Double_t)indices.GetFirst() /xmax,
                                   (Double_t)indices.GetSecond()/ymax,
                                   2);
    marker->Draw();
  }
  
  AliMpVPadIterator* iter2
    = segm.CreateIterator(AliMpArea(pad.Position(),2.*pad.Dimensions()*1.1));

  Int_t counter = 0;
  for( iter2->First(); !iter2->IsDone() && i<10; iter2->Next()) {
    //Int_t ix = iter2->CurrentItem().GetIndices().GetFirst();
    //Int_t iy = iter2->CurrentItem().GetIndices().GetSecond();
    cout<<"Iterator number "<< counter++ << " at "<< iter2->CurrentItem().GetIndices() <<endl;
  }
  
  delete iter2;
  delete sector;
}

void testSt12NeighboursPadIterator()
{
  AliMq::Station12Type  station[2] = { AliMq::kStation1, AliMq::kStation2 }; 
  AliMp::PlaneType      plane[2]   = { AliMp::kBendingPlane, AliMp::kNonBendingPlane };
  
  for ( Int_t is = 0; is < 2; is++ ) {
    for ( Int_t ip = 0; ip < 2; ip++ ) {
    
      cout << "Running testNeighboursPadIterator for " 
           << AliMq::Station12TypeName(station[is]) << "  "
           << AliMp::PlaneTypeName(plane[ip])  << " ... " << endl;
       
      testNeighboursPadIterator(station[is], plane[ip]);
    
      cout << "... end running " << endl << endl;
    }  
  }   
}  
  
