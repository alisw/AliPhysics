// $Id$
//
// Test macro for testing which pad is seen as "existing" by AliMpSector.

void testExistingPads(AliMpStationType station = kStation1,
                      AliMpPlaneType plane = kBendingPlane) 
{
  AliMpReader r(station, plane);

  AliMpSector *sector=r.BuildSector();
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(sector);

  TCanvas* c1 = new TCanvas("view",
                            "AliMpSectorPainter::Draw() output (view per pad)");
  painter->Draw("ZSSMP");
  c1->Update();

  TH2C* histo = new TH2C("histo","Existing pads",150,-0.5,149.5,
                                               200,-0.5,199.5);
  TCanvas* c2 = new TCanvas("c2","Only existing pads are filled");

  AliMpSectorSegmentation segmentation(sector);
  for (Int_t i=0; i<150;i++){
    for (Int_t j=0;j<200;++j){

      AliMpIntPair indices(i,j);
      if (segmentation.HasPad(indices)) histo->Fill(i,j);
    }
  }

  c2->cd();
  histo->Draw("box");
}
