// $Id$
//
// Test macro for which verify that all FindPosition, FindIndices
// and FindLocation methods are consistents between them.

void testSectorFind(AliMpStationType station = kStation1,
                    AliMpPlaneType plane = kBendingPlane)
{
  AliMpReader r(station, plane);

  AliMpSector *sector=r.BuildSector();
  AliMpSectorSegmentation segmentation(sector);
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(sector);
  painter->Draw("ZSSMP");

  for (Int_t i=0; i<90;i++){
    cout<<"Verifying column "<<i<<"....."<<endl;

    for (Int_t j=0;j<230;++j)
      segmentation.CircleTest(AliMpIntPair(i,j));    
  }
}
