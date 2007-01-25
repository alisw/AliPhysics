// $Id$
// $MpId: testUpdateGlobalIndices.C,v 1.7 2005/08/24 08:53:27 ivana Exp $
//
// Tests updating global indices of motif positions from file.

void testUpdateGlobalIndices()
{
  AliMpSectorReader reader(AliMp::kStation1, AliMp::kNonBendingPlane);  
  //reader.SetVerboseLevel(1);
  
  // Read data 
  AliMpSector* sector = reader.BuildSector();

  sector->GetMotifMap()->UpdateGlobalIndices("motif_map.dat");  
  
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(sector);
  TCanvas* canvas = new TCanvas();
  canvas->cd();
  painter->Draw("ZSSMI");
}          
 
