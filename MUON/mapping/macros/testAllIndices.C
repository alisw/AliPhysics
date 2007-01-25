// $Id$
// $MpId: testAllIndices.C,v 1.7 2005/08/24 08:53:27 ivana Exp $
//
// Test macro for testing which pad is seen as "existing" by AliMpSector.

void testAllIndices(AliMp::StationType station = AliMp::kStation1,
                    AliMp::PlaneType plane = AliMp::kBendingPlane) 
{
  AliMpSectorReader r(station, plane);

  AliMpSector *sector=r.BuildSector();
  AliMpSectorSegmentation segmentation(sector);
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(sector);

  TCanvas* c1 = new TCanvas("view",
                            "MSectorPainter::Draw() output (view per pad)");
  painter->Draw("ZSSMP");
  c1->Update();

  Int_t maxPadIndexX = segmentation.MaxPadIndexX();
  Int_t maxPadIndexY = segmentation.MaxPadIndexY();
  
  // Define histogram limits
  Int_t nx = (maxPadIndexX/10 + 1)*10;
  Int_t ny = (maxPadIndexY/10 + 1)*10;
  TH2C* histo  = new TH2C("histo","Existing pads", 
                          nx, -0.5, nx-0.5, ny, -0.5, ny-0.5);

  Int_t nx2 = 950/2;
  Int_t ny2 = 950/2;
  if (station == AliMp::kStation2) {
    nx2 = 1200/2;
    ny2 = 1200/2;
  }
  TH2F* histo2 = new TH2F("histo2","Existing positions",
                          nx2, 0, nx2*2, ny2, 0, ny2*2);

  // Define canvas
  TCanvas* c2 = new TCanvas("c2","Only existing pads are filled");
  TCanvas* c3 = new TCanvas("c3","Positions");
  
  for (Int_t irow=0;irow<sector->GetNofRows();irow++){
    AliMpRow* row = sector->GetRow(irow);
    
    for (Int_t  iseg=0;iseg<row->GetNofRowSegments();iseg++){
      AliMpVRowSegment* seg = row->GetRowSegment(iseg);
      
      for (Int_t imot=0;imot<seg->GetNofMotifs();imot++){
        AliMpMotifPosition* motifPos 
         = sector->GetMotifMap()->FindMotifPosition(seg->GetMotifPositionId(imot));
         
        for (Int_t gassNum=0;gassNum<64;gassNum++){
          if (motifPos->GetMotif()->GetMotifType()->FindConnectionByGassiNum(gassNum)){
          
            AliMpPad pad = segmentation.PadByLocation(AliMpIntPair(motifPos->GetID(),gassNum));
            if (pad != AliMpPad::Invalid()) {
              histo->Fill (pad.GetIndices().GetFirst(),
                          pad.GetIndices().GetSecond());
              histo2->Fill(pad.Position().X(),
                           pad.Position().Y());
            }
          }
        }
      }
    }
  }
  c2->cd();
  histo->Draw("col");
  c3->cd();
  histo2->Draw("box");
}
