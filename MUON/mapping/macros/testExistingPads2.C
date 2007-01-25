// $Id$
// $MpId: testExistingPads2.C,v 1.3 2005/08/24 08:53:27 ivana Exp $
//
// Extended testExistingPads macro for testing which pad is seen as 
// "existing" by AliMpSector.
// The loop over pads is performed in two ways:
// 1) via sector area pad iterator
// 2) via indices

#include <iomanip>

void testExistingPads2(AliMp::StationType station = AliMp::kStation1,
                       AliMp::PlaneType plane = AliMp::kBendingPlane) 
{
  AliMpSectorReader r(station, plane);
  AliMpSector* sector = r.BuildSector();
  AliMpSectorSegmentation segmentation(sector);
  
  ofstream out1("testExistingPads2.ixiy.out", ios::out);
  ofstream out2("testExistingPads2.iter.out", ios::out);
  
  // First loop over indices
  cout << "Iterating via indices ..." << endl;
  Int_t counter1 = 0;
  if ( sector->GetDirection() == AliMp::kX )
    for (Int_t i=1; i<segmentation.MaxPadIndexX()+1; i++) {
      for (Int_t j=1; j<segmentation.MaxPadIndexY()+1; j++) {
        AliMpIntPair indices(i,j);
        if ( segmentation.HasPad(indices) ) 
          out1 << std::setw(4) << ++counter1 << "  "
	       << segmentation.PadByIndices(indices) << endl;;
      }
    }

  if ( sector->GetDirection() == AliMp::kY )
    for (Int_t j=1; j<segmentation.MaxPadIndexY()+1; j++) {
      for (Int_t i=1; i<segmentation.MaxPadIndexX()+1; i++) {
        AliMpIntPair indices(i,j);
        if ( segmentation.HasPad(indices) ) 
          out1 << std::setw(4) << ++counter1 << "  "
	       << segmentation.PadByIndices(indices) << endl;;
      }
    }
  
  // Second loop via pad area iterator
  cout << "Iterating via iterator ..." << endl;
  AliMpArea area(TVector2(60.,60.), TVector2(60.,60.) );
  AliMpVPadIterator* it = segmentation.CreateIterator(area);
  Int_t counter2 = 0;
  for (it->First(); ! it->IsDone(); it->Next()){
   
    out2 << std::setw(4) << ++counter2 << "  "
	 << it->CurrentItem() << endl;
  }  
}
