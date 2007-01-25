// $Id$
// $MpId: testSectorFind.C,v 1.10 2005/10/28 15:36:08 ivana Exp $
//
// Test macro for which verify that all FindPosition, FindIndices
// and FindLocation methods are consistents between them.

void testSectorFind(AliMp::StationType station = AliMp::kStation1,
                    AliMp::PlaneType plane = AliMp::kBendingPlane,
		    Bool_t rootInput = false) 
{
  AliMpSector *sector = 0;
  if (!rootInput) {
    AliMpSectorReader r(station, plane);
    sector=r.BuildSector();
  }
  else  {
    TString filePath = AliMpFiles::SectorFilePath(station,plane);
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
