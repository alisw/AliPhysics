// $Id$
//
// Test macro that starts from a given pad and prints 
// all pads up, down, right, left from this pad
// (up to the plane border).

void testPadsUpEtc(AliMpStationType station = kStation1,
                   AliMpPlaneType  planeType = kBendingPlane)
{
  AliMpPlane* plane = AliMpPlane::Create(station, planeType);
  AliMpPlaneSegmentation planeSegmentation(plane);
  
  //AliMpIntPair indices(85, 101);
  AliMpIntPair indices(-129, 10);
 
  AliMpPad pad;
  if (planeSegmentation.HasPad(indices)) {
	  
    pad = planeSegmentation.PadByIndices(indices);
    cout << "Pad:      "  << pad << endl << endl;
  
    cout << "######### GO UP  ############### " << endl;

    AliMpPadPair nextPads(pad, pad);
    while (nextPads.GetFirst().IsValid()) {
      nextPads = planeSegmentation.PadsUp(nextPads.GetFirst());
      cout << " up    1: "  << nextPads.GetFirst() << endl;    
      cout << "       2: "  << nextPads.GetSecond() << endl;    
    }

    cout << "######### GO DOWN  ############### " << endl;

    nextPads = AliMpPadPair(pad, pad);
    while (nextPads.GetFirst().IsValid()) {
      nextPads = planeSegmentation.PadsDown(nextPads.GetFirst());
      cout << " down  1: "  << nextPads.GetFirst() << endl;    
      cout << "       2: "  << nextPads.GetSecond() << endl;    
    }

    cout << "######### GO RIGHT  ############### " << endl;

    nextPads = AliMpPadPair(pad, pad);
    while (nextPads.GetFirst().IsValid()) {
      nextPads = planeSegmentation.PadsRight(nextPads.GetFirst());
      cout << " right 1: "  << nextPads.GetFirst() << endl;    
      cout << "       2: "  << nextPads.GetSecond() << endl;    
    }

    cout << "######### GO LEFT  ############### " << endl;

    nextPads = AliMpPadPair(pad, pad);
    while (nextPads.GetFirst().IsValid()) {
      nextPads = planeSegmentation.PadsLeft(nextPads.GetFirst());
      cout << " left  1: "  << nextPads.GetFirst() << endl;    
      cout << "       2: "  << nextPads.GetSecond() << endl;    
    }
  }  
}
