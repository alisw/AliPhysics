// $Id$
//
// Tests updating global indices of motif positions from file.

void testUpdateGlobalIndices()
{
  AliMpReader reader(kStation1, kNonBendingPlane);  
  //reader.SetVerboseLevel(1);
  
  // Read data 
  AliMpSector* sector = reader.BuildSector();

  sector->GetMotifMap()->UpdateGlobalIndices("motif_map.dat");  
  
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(sector);
  TCanvas* canvas = new TCanvas();
  canvas->cd();
  painter->Draw("ZSSMI");
}          
 
