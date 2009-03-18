// $Id$
// $MpId: testExistingPads2.C,v 1.3 2005/08/24 08:53:27 ivana Exp $
//
// Extended testExistingPads macro for testing which pad is seen as 
// "existing" by AliMpSectorSegmentation or AliMpFastSegmentation.
// The loop over pads is performed in two ways:
// 1) via sector area pad iterator
// 2) via indices
// To run macro:
// root [0] .L testExistingPads2.C+
// root [1] testExistingPads2();


#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpSector.h"
#include "AliMpSectorReader.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpFastSegmentation.h"
#include "AliMpArea.h"
#include "AliMpVPadIterator.h"

#include <Riostream.h>
#include <TString.h>

#endif

void testExistingPads2(AliMq::Station12Type station, AliMp::PlaneType plane,
                       Bool_t useFastSegmentation = kTRUE) 
{
  // Read data
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);
  AliMpSectorReader r(dataStreams, station, plane);
  AliMpSector* sector = r.BuildSector();

  // Segmentation
  AliMpVSegmentation* segmentation;
  if ( useFastSegmentation)
    segmentation = new AliMpFastSegmentation(new AliMpSectorSegmentation(sector));
  else
    segmentation = new AliMpSectorSegmentation(sector);
    
  // Output files 
  TString outName = "testExistingPads.";
  outName += AliMq::Station12TypeName(station);
  outName += AliMp::PlaneTypeName(plane);
  TString out1Name = outName + ".ixiy.out";
  TString out2Name = outName + ".iter.out";

  ofstream out1(out1Name.Data(), ios::out);
  ofstream out2(out2Name.Data(), ios::out);
  
  // First loop over indices
  cout << "Iterating via indices ..." << endl;
  Int_t counter1 = 0;
  if ( sector->GetDirection() == AliMp::kX )
    for (Int_t i=1; i<segmentation->MaxPadIndexX()+1; i++) {
      for (Int_t j=1; j<segmentation->MaxPadIndexY()+1; j++) {
        AliMpIntPair indices(i,j);
        if ( segmentation->HasPad(indices) ) 
          out1 << std::setw(4) << ++counter1 << "  "
	       << segmentation->PadByIndices(indices) << endl;;
      }
    }

  if ( sector->GetDirection() == AliMp::kY )
    for (Int_t j=1; j<segmentation->MaxPadIndexY()+1; j++) {
      for (Int_t i=1; i<segmentation->MaxPadIndexX()+1; i++) {
        AliMpIntPair indices(i,j);
        if ( segmentation->HasPad(indices) ) 
          out1 << std::setw(4) << ++counter1 << "  "
	       << segmentation->PadByIndices(indices) << endl;;
      }
    }
  
  // Second loop via pad area iterator
  cout << "Iterating via iterator ..." << endl;
  AliMpArea area(TVector2(60.,60.), TVector2(60.,60.) );
  AliMpVPadIterator* it = segmentation->CreateIterator(area);
  Int_t counter2 = 0;
  for (it->First(); ! it->IsDone(); it->Next()){
   
    out2 << std::setw(4) << ++counter2 << "  "
	 << it->CurrentItem() << endl;
  }  
  
  delete segmentation;
}

void testSt12ExistingPads2()
{
  AliMq::Station12Type  station[2] = { AliMq::kStation1, AliMq::kStation2 }; 
  AliMp::PlaneType      plane[2]   = { AliMp::kBendingPlane, AliMp::kNonBendingPlane };
  
  for ( Int_t is = 0; is < 2; is++ ) {
    for ( Int_t ip = 0; ip < 2; ip++ ) {
    
      cout << "Running testExistingPads2 for " 
           << AliMq::Station12TypeName(station[is]) << "  "
           << AliMp::PlaneTypeName(plane[ip])  << " ... " << endl;
       
      testExistingPads2(station[is], plane[ip]);
    
      cout << "... end running " << endl << endl;
    }  
  }   
}  
  
