// $Id$
// $MpId: testSectorPadIterators.C,v 1.11 2006/03/15 13:07:07 ivana Exp $
//
// Test macro for reading  sector, and iterate over it

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpSector.h"
#include "AliMpSectorPadIterator.h"
#include "AliMpSectorReader.h"
#include "AliMpArea.h"
#include "AliMpVPadIterator.h"
#include "AliMpVPainter.h"
#include "AliMpEncodePair.h"

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

void testSectorPadIterators(AliMq::Station12Type station, AliMp::PlaneType plane)
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSectorReader r(dataStreams, station, plane);
  AliMpSector* sector = r.BuildSector();
  
  Int_t num=0;
  
  //new TCanvas("canv");
  CreateTCanvas("", "", station, plane);

  const Double_t xmax=150;
  const Double_t ymax=250;

  AliMpSectorPadIterator it = AliMpSectorPadIterator(sector);

  for (it.First(); ! it.IsDone(); it.Next()) {
    Long_t indices = it.CurrentItem().GetIndices();
    cout<<"Iterator number "<< num << " at ";
    AliMp::PairPut(cout, indices) <<  endl;
    num++;
    TMarker* marker = new TMarker( (Double_t)AliMp::PairFirst(indices) /xmax,
                                   (Double_t)AliMp::PairSecond(indices)/ymax,
                                   2);
    marker->Draw();
  }
  
  delete sector;
}

void testSt12SectorPadIterators()
{
  AliMq::Station12Type  station[2] = { AliMq::kStation1, AliMq::kStation2 }; 
  AliMp::PlaneType      plane[2]   = { AliMp::kBendingPlane, AliMp::kNonBendingPlane };
  
  for ( Int_t is = 0; is < 2; is++ ) {
    for ( Int_t ip = 0; ip < 2; ip++ ) {
    
      cout << "Running testSectorPadIterators for " 
           << AliMq::Station12TypeName(station[is]) << "  "
           << AliMp::PlaneTypeName(plane[ip])  << " ... " << endl;
       
      testSectorPadIterators(station[is], plane[ip]);
    
      cout << "... end running " << endl << endl;
    }  
  }   
}  
  
