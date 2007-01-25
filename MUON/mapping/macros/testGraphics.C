// $Id$
// $MpId: testGraphics.C,v 1.13 2005/10/28 15:36:08 ivana Exp $
//
// Test macro for drawing sector data.

void testGraphics(AliMp::StationType station = AliMp::kStation1,
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
    
  AliMpVPainter *painter=AliMpVPainter::CreatePainter(sector);

  TCanvas* canvas  = new TCanvas();
  TCanvas* canvas2 = new TCanvas();
  TCanvas* c[4];
  for (int i=0;i<4;++i) {
    c[i] = new TCanvas();
    c[i]->Divide(2,2);
  }  

 //first, paint the whole sector
  canvas->cd();
  painter->Draw("");
  //now paint each rows
  c[0]->cd(1);
  painter->Draw("R");
  //paint each row segments in each row
  c[0]->cd(2);
  painter->Draw("RS");
  //paint each motifs, in each row segments in each row
  c[0]->cd(3);
  painter->Draw("RSMP");
  //paint each pads, in each motifs, in each row segments in each row
  c[0]->cd(4);
  painter->Draw("RSMT");

  ///////////////////////////////
  //now paint each rows, wwith its name
  c[1]->cd(1);
  painter->Draw("RT");
  //paint each row segments in each row, and its name
  c[1]->cd(2);
  painter->Draw("RST");
  //paint each motifs, in each row segments in each row
  c[1]->cd(3);
  painter->Draw("RSMX");
  c[1]->cd(4);
  painter->Draw("RSMI");

  ///////////////////////////////
  //now paint each zones
  c[2]->cd(1);
  painter->Draw("Z");
  //paint each sub-zones, in each zones
  c[2]->cd(2);
  painter->Draw("ZS");
  //paint each row segments, in each sub-zone, ...
  c[2]->cd(3);
  painter->Draw("ZSS");
  // each motifs, in each row segments, ...
  c[2]->cd(4);
  painter->Draw("ZSSM");

  ///////////////////////////////
  //now paint each zones with its name
  c[3]->cd(1);
  painter->Draw("ZT");
  //paint each sub-zones, in each zones with its name
  c[3]->cd(2);
  painter->Draw("ZST");
  //paint each row segments, in each sub-zone, ... with its name
  c[3]->cd(3);
  painter->Draw("ZSST");
  // each motifs, in each row segments, ... with its name
  c[3]->cd(4);
  painter->Draw("ZSSMT");
  // now, draw a specific motif, in a whole canvas, and
  // print all its pad names
  Int_t id = sector->GetRow(5)->GetRowSegment(0)->GetMotifPositionId(0);
  AliMpMotifPosition* motifPos = sector->GetMotifMap()->FindMotifPosition(id);
  motifPainter = AliMpVPainter::CreatePainter(motifPos);
  TCanvas* onepad = new TCanvas("onepad","One motif");
  motifPainter->Draw("PT");
  
  //////////////////////////////  
  //now paint motifs with their real contours and mani Ids
  canvas2->cd();
  painter->Draw("RSMCI");
}
