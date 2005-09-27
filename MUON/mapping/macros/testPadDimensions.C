// $Id$
// $MpId: testPadDimensions.C,v 1.5 2005/09/26 16:05:25 ivana Exp $
//
// Test macro for testing retrieving of pad dimensions from
// the map in AliMpSectorSegmentation.

void testPadDimensions(AliMpStationType station = kStation1,
                       AliMpPlaneType plane = kBendingPlane,
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
