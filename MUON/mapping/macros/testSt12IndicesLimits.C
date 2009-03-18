// $Id$
// $MpId: testIndicesLimits.C,v 1.5 2005/10/28 15:36:08 ivana Exp $
//
// Test macro for indices limits.

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpSector.h"
#include "AliMpSectorReader.h"
#include "AliMpSectorSegmentation.h" 
#include "AliMpRow.h"
#include "AliMpVRowSegment.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotifType.h"

#include <Riostream.h>
#include <TCanvas.h>
#include <TH2.h>

#endif

void testIndicesLimits(AliMq::Station12Type station,AliMp::PlaneType plane)
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSectorReader r(dataStreams, station, plane);

  AliMpSector *sector=r.BuildSector();
  AliMpSectorSegmentation segmentation(sector);

  // Loop over rows
  for (Int_t i=0; i<sector->GetNofRows(); i++) {
    AliMpRow* row = sector->GetRow(i);
    cout << i
         << "th row limits: "
	 << row->GetLowIndicesLimit() << "  "
	 << row->GetHighIndicesLimit() << endl;
    
    // Loop over row segments
    for (Int_t j=0; j<row->GetNofRowSegments(); j++) {
      AliMpVRowSegment* rowSegment = row->GetRowSegment(j);
      cout << "   "
           << j
           << "th row segment limits: "
	   << rowSegment->GetLowIndicesLimit() << "  "
	   << rowSegment->GetHighIndicesLimit() << endl;
      
      // Loop over motif positions
      for (Int_t k=0; k<rowSegment->GetNofMotifs(); k++) {
        Int_t mposID = rowSegment->GetMotifPositionId(k);
        AliMpMotifPosition* mPos = 
	  sector->GetMotifMap()->FindMotifPosition(mposID);     
        if (mPos) {
          cout << "      "
               << mPos->GetID()
               << " motif position limits: "
	       << mPos->GetLowIndicesLimit() << "  "
	       << mPos->GetHighIndicesLimit() << endl;
	}
	else {
	  cerr << "Motif position "
	       <<  mposID << " not found in the map" << endl; 
	}           
      }
    }  	   
  }   
  
  delete sector;
}			       
      
void testSt12IndicesLimits()
{
  AliMq::Station12Type  station[2] = { AliMq::kStation1, AliMq::kStation2 }; 
  AliMp::PlaneType      plane[2]   = { AliMp::kBendingPlane, AliMp::kNonBendingPlane };
  
  for ( Int_t is = 0; is < 2; is++ ) {
    for ( Int_t ip = 0; ip < 2; ip++ ) {
    
      cout << "Running testIndicesLimits for " 
           << AliMq::Station12TypeName(station[is]) << "  "
           << AliMp::PlaneTypeName(plane[ip])  << " ... " << endl;
       
      testIndicesLimits(station[is], plane[ip]);
    
      cout << "... end running " << endl << endl;
    }  
  }   
}  
  
 

 
