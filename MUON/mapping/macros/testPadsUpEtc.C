// $Id$
// $MpId: testPadsUpEtc.C,v 1.7 2005/09/26 16:05:25 ivana Exp $
//
// Test macro that starts from a given pad and prints 
// all pads up, down, right, left from this pad
// (up to the plane border).

void testPadsUpEtc(AliMpStationType station = kStation1,
                   AliMpPlaneType  plane = kBendingPlane,
		   Bool_t rootInput = false)
{
  AliMpSector *sector = 0;
  if (!rootInput) {
    AliMpSectorReader r(station, plane);
    sector=r.BuildSector();
  }
  else  {
    TString filePath = AliMpFiles::Instance()->SectorFilePath(station,plane);
    filePath.ReplaceAll("zones.dat", "sector.root"); 

    TFile f(filePath.Data(), "READ");
    sector = (AliMpSector*)f.Get("Sector");
  }  
    
  AliMpSectorSegmentation segmentation(sector);
  
  AliMpIntPair indices(85, 101);
  if( plane == kNonBendingPlane) indices = AliMpIntPair(129, 10);
 
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
