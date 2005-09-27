// $Id$
// $MpId: testSectorFind.C,v 1.9 2005/09/26 16:05:25 ivana Exp $
//
// Test macro for which verify that all FindPosition, FindIndices
// and FindLocation methods are consistents between them.

void testSectorFind(AliMpStationType station = kStation1,
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
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(sector);
  painter->Draw("ZSSMP");

  for (Int_t i=0; i<90;i++){
    cout<<"Verifying column "<<i<<"....."<<endl;

    for (Int_t j=0;j<230;++j)
      segmentation.CircleTest(AliMpIntPair(i,j));    
  }
}
