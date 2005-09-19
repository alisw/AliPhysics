// $Id$
// $MpId: testExistingPads.C,v 1.10 2005/08/24 08:53:27 ivana Exp $
//
// Test macro for testing which pad is seen as "existing" by AliMpSector.

void testExistingPads(AliMpStationType station = kStation1,
                      AliMpPlaneType plane = kBendingPlane) 
{
  AliMpSectorReader r(station, plane);
  AliMpSector *sector=r.BuildSector();
  AliMpSectorSegmentation segmentation(sector);
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(sector);

  TCanvas* c1 = new TCanvas("view",
                            "AliMpSectorPainter::Draw() output (view per pad)");
  painter->Draw("ZSSMP");
  c1->Update();

  Int_t maxPadIndexX = segmentation.MaxPadIndexX();
  Int_t maxPadIndexY = segmentation.MaxPadIndexY();
  
  // Define histogram limits
  Int_t nx = (maxPadIndexX/10 + 1)*10;
  Int_t ny = (maxPadIndexY/10 + 1)*10;
  TH2C* histo = new TH2C("histo","Existing pads", 
                          nx, -0.5, nx-0.5, ny, -0.5, ny-0.5);

  TCanvas* c2 = new TCanvas("c2","Only existing pads are filled");

  AliMpSectorSegmentation segmentation(sector);
  for (Int_t i=0; i<maxPadIndexX+1;i++){
    for (Int_t j=0;j<maxPadIndexY+1;++j){

      AliMpIntPair indices(i,j);
      if (segmentation.HasPad(indices)) histo->Fill(i,j);
    }
  }

  c2->cd();
  histo->Draw("box");
}
