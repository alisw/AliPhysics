// $Id$
//
// Test macro for reading  sector, and iterate over it

void testNeighboursPadIterator(AliMpStationType station = kStation1,
                               AliMpPlaneType plane = kBendingPlane, 
                               Int_t i=50, Int_t j=50)
{
  if (!gInterpreter->IsLoaded("mlibs.C")){ 
    gROOT->LoadMacro("mlibs.C");
    gInterpreter->ProcessLine("mlibs()");
  }  

  AliMpReader r(station, plane);
  AliMpSector* sect = r.BuildSector();
  AliMpSectorSegmentation segm(sect);  
  
  TCanvas *can = new TCanvas("canv");

  const Double_t xmax=75;
  const Double_t ymax=120;
  
  Int_t num=0;

  AliMpPad pad = segm.PadByIndices(AliMpIntPair(i,j));
  AliMpNeighboursPadIterator it = AliMpNeighboursPadIterator(&segm, pad,kFALSE);
  
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
