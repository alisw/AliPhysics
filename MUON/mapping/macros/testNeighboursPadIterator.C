// $Id$
// $MpId: testNeighboursPadIterator.C,v 1.10 2005/10/28 15:36:08 ivana Exp $
//
// Test macro for reading  sector, and iterate over it

void testNeighboursPadIterator(AliMp::StationType station = AliMp::kStation1,
                               AliMp::PlaneType plane = AliMp::kBendingPlane, 
		               Bool_t rootInput = false,
                               Int_t i=50, Int_t j=50)
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

  AliMpSectorSegmentation segm(sector);  
  
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
  
  AliMpVPadIterator* iter2
    = segm.CreateIterator(AliMpArea(pad.Position(),2.*pad.Dimensions()*1.1));

  Int_t i=0;
  for( iter2->First(); !iter2->IsDone() && i<10; iter2->Next()) {
    Int_t ix = iter2->CurrentItem().GetIndices().GetFirst();
    Int_t iy = iter2->CurrentItem().GetIndices().GetSecond();
    cout<<"Iterator number "<< i << " at "<< iter2->CurrentItem().GetIndices() <<endl;
    i++;
  }
  
  delete iter2;
  delete sector;
}
