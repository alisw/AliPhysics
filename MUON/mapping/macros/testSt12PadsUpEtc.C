// $Id$
// $MpId: testPadsUpEtc.C,v 1.8 2005/10/28 15:36:08 ivana Exp $
//
// Test macro that starts from a given pad and prints 
// all pads up, down, right, left from this pad
// (up to the plane border).

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpSector.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpSectorReader.h"

#include <Riostream.h>

#endif

void testPadsUpEtc(AliMq::Station12Type station, AliMp::PlaneType  plane)
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSectorReader r(dataStreams, station, plane);
  AliMpSector* sector = r.BuildSector();
  AliMpSectorSegmentation segmentation(sector);
  
  AliMpIntPair indices(85, 101);
  if( plane == AliMp::kNonBendingPlane) indices = AliMpIntPair(129, 10);
 
  AliMpPad pad;
  if (segmentation.HasPad(indices)) {
	  
    pad = segmentation.PadByIndices(indices);
    cout << "Pad:      "  << pad << endl << endl;
  
    cout << "######### GO UP  ############### " << endl;

    AliMpPadPair nextPads(pad, pad);
    while (nextPads.GetFirst().IsValid()) {
      nextPads = segmentation.PadsUp(nextPads.GetFirst());
      cout << " up    1: "  << nextPads.GetFirst() << endl;    
      cout << "       2: "  << nextPads.GetSecond() << endl;    
    }

    cout << "######### GO DOWN  ############### " << endl;

    nextPads = AliMpPadPair(pad, pad);
    while (nextPads.GetFirst().IsValid()) {
      nextPads = segmentation.PadsDown(nextPads.GetFirst());
      cout << " down  1: "  << nextPads.GetFirst() << endl;    
      cout << "       2: "  << nextPads.GetSecond() << endl;    
    }

    cout << "######### GO RIGHT  ############### " << endl;

    nextPads = AliMpPadPair(pad, pad);
    while (nextPads.GetFirst().IsValid()) {
      nextPads = segmentation.PadsRight(nextPads.GetFirst());
      cout << " right 1: "  << nextPads.GetFirst() << endl;    
      cout << "       2: "  << nextPads.GetSecond() << endl;    
    }

    cout << "######### GO LEFT  ############### " << endl;

    nextPads = AliMpPadPair(pad, pad);
    while (nextPads.GetFirst().IsValid()) {
      nextPads = segmentation.PadsLeft(nextPads.GetFirst());
      cout << " left  1: "  << nextPads.GetFirst() << endl;    
      cout << "       2: "  << nextPads.GetSecond() << endl;    
    }
  }  
}

void testPadsUpEtc()
{
  AliMq::Station12Type  station[2] = { AliMq::kStation1, AliMq::kStation2 }; 
  AliMp::PlaneType      plane[2]   = { AliMp::kBendingPlane, AliMp::kNonBendingPlane };
  
  for ( Int_t is = 0; is < 2; is++ ) {
    for ( Int_t ip = 0; ip < 2; ip++ ) {
    
      cout << "Running testPadsUpEtc for " 
           << AliMq::Station12TypeName(station[is]) << "  "
           << AliMp::PlaneTypeName(plane[ip])  << " ... " << endl;
       
      testPadsUpEtc(station[is], plane[ip]);
    
      cout << "... end running " << endl << endl;
    }  
  }   
}  
  
