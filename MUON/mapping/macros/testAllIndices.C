// $Id$
//
// Test macro for testing which pad is seen as "existing" by AliMpSector.

void testAllIndices() 
{
  if (!gInterpreter->IsLoaded("mlibs.C")){ 
    gROOT->LoadMacro("mlibs.C");
    gInterpreter->ProcessLine("mlibs()");
  }  

  AliMpReader r(kStation1, kNonBendingPlane);

  AliMpSector *sector=r.BuildSector();
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(sector);

  TCanvas* c1 = new TCanvas("view",
                            "MSectorPainter::Draw() output (view per pad)");
  painter->Draw("ZSSMP");
  c1->Update();

  TH2C* histo = new TH2C("histo","Existing pads",90,-0.5,89.5,
                                                 230,-0.5,229.5);
  TH2F* histo2 = new TH2F("histo2","Existing positions",950/2,0,950,
                                                        950/2,0,950);
  TCanvas* c2 = new TCanvas("c2","Only existing pads are filled");
  TCanvas* c3 = new TCanvas("c3","Positions");

  AliMpSectorSegmentation segmentation(sector);
  
  
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
