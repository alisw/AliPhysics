// $Id$
// $MpId: testSectorPadIterators.C,v 1.11 2006/03/15 13:07:07 ivana Exp $
//
// Test macro for reading  sector, and iterate over it

void testSectorPadIterators(AliMp::StationType station = AliMp::kStation1,
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
  
  Int_t num=0;
  
  TCanvas *can = new TCanvas("canv");

  const Double_t xmax=150;
  const Double_t ymax=250;

  AliMpSectorPadIterator it = AliMpSectorPadIterator(sector);

  for (it.First(); ! it.IsDone(); it.Next()) {
    AliMpIntPair indices = it.CurrentItem().GetIndices();
    cout<<"Iterator number "<< num << " at "<< indices <<endl;
    num++;
    TMarker* marker = new TMarker( (Double_t)indices.GetFirst() /xmax,
                                   (Double_t)indices.GetSecond()/ymax,
                                   2);
    marker->Draw();
  }
  
  delete sector;
}
