// $Id$
//
// Test macro for testing retrieving of pad dimensions from
// the map in AliMpSectorSegmentation.

void testPadDimensions(AliMpStationType station = kStation1,
                       AliMpPlaneType plane = kBendingPlane) 
{
  AliMpReader r(station, plane);
  AliMpSector* sector=r.BuildSector();
  AliMpSectorSegmentation segmentation(sector);  

  TVector2 previousDimensions;
  for (Int_t i=0; i<150;i++) 
    for (Int_t j=0;j<200;++j) {

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
