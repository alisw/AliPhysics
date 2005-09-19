// $Id$
// $MpId: testPadDimensions.C,v 1.4 2005/08/24 08:53:27 ivana Exp $
//
// Test macro for testing retrieving of pad dimensions from
// the map in AliMpSectorSegmentation.

void testPadDimensions(AliMpStationType station = kStation1,
                       AliMpPlaneType plane = kBendingPlane) 
{
  AliMpSectorReader r(station, plane);
  AliMpSector* sector=r.BuildSector();
  AliMpSectorSegmentation segmentation(sector);  
  segmentation.PrintZones(); 
  
  TVector2 previousDimensions;
  for (Int_t i=1; i<segmentation.MaxPadIndexX()+1;i++) 
    for (Int_t j=1;j<segmentation.MaxPadIndexY()+1;++j) {

      AliMpIntPair indices(i,j);
      if (segmentation.HasPad(indices)) {

        // Check pad dimensions
	AliMpPad pad = segmentation.PadByIndices(indices);
	TVector2 dimensions = segmentation.PadDimensions(segmentation.Zone(pad));
	
	if ( dimensions.X() != previousDimensions.X() || 
	     dimensions.Y() != previousDimensions.Y() ) {

          // Print dimensions
	  cout << "Pad: " << indices;
	  cout << "  dimensions: (" << dimensions.X() << ", " << dimensions.Y() << ")" 
	       << endl;
          
	  previousDimensions = dimensions;
        }	     
     }
   }
}
