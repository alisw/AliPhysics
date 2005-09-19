// $Id$
// $MpId: testSectorPadIterators.C,v 1.7 2005/08/24 08:53:27 ivana Exp $
//
// Test macro for reading  sector, and iterate over it

void testSectorPadIterators(AliMpStationType station = kStation1,
                            AliMpPlaneType plane = kBendingPlane)
{
  AliMpSectorReader r(station, plane);
  AliMpSector* sect = r.BuildSector();
  
  Int_t num=0;
  
  TCanvas *can = new TCanvas("canv");

  const Double_t xmax=150;
  const Double_t ymax=250;

  AliMpSectorPadIterator it = AliMpSectorPadIterator(sect);

  for (it.First(); ! it.IsDone(); it.Next()) {
    AliMpIntPair indices = it.CurrentItem().GetIndices();
    cout<<"Iterator number "<< num++ << " at "<< indices <<endl;
    TMarker* marker = new TMarker( (Double_t)indices.GetFirst() /xmax,
                                   (Double_t)indices.GetSecond()/ymax,
                                   2);
    marker->Draw();
  }
  
  delete sect;
}
