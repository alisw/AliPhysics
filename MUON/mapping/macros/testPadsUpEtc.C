// $Id$
//
// Test macro that starts from a given pad and prints 
// all pads up, down, right, left from this pad
// (up to the plane border).

void testPadsUpEtc(AliMpStationType station = kStation1,
                   AliMpPlaneType  planeType = kBendingPlane)
{
  AliMpReader reader(station, planeType);  
  AliMpSector* sector = reader.BuildSector();
  AliMpSectorSegmentation segmentation(sector);
  
  AliMpIntPair indices(85, 101);
  if( planeType == kNonBendingPlane) indices = AliMpIntPair(129, 10);
 
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
