// $Id$
// $MpId: testAnyPadIterators.C,v 1.17 2006/03/15 13:07:07 ivana Exp $
//
// Test macro for reading  sector, and iterate over it

class AliMpVPadIterator;
void MarkPads(AliMpVPadIterator& it,Int_t xmax,Int_t ymax,Bool_t print=kTRUE)
{
// This function works with pad iterator, no matter which kind of iterator
// it is. So it can be used for drawing all pad of the sector (AliMpSectorPadIterator)
// or all pad around a given pad (AliMpNeighboursPadIterator), as with pads
// of a given motif type (AliMpMotifTypePadIterator)

  Int_t num=0;

  for (it.First(); ! it.IsDone(); it.Next()){
    AliMpIntPair indices = it.CurrentItem().GetIndices();
    if (print) cout<<"Iterator number "<< num << " at "<< indices <<endl;
    num++;
    TMarker* marker = new TMarker( (Double_t)indices.GetFirst() /xmax,
                                   (Double_t)indices.GetSecond()/ymax,
                                   2);
    marker->Draw();
  }
}

void testAnyPadIterators(AliMp::StationType station = AliMp::kStation1,
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
    
  TCanvas *canv = new TCanvas("canv");
  canv->Divide(2,2);
  //canv_1->Divide(2);
  
  canv->cd(1);
  MarkPads(AliMpSectorPadIterator(sector), 150, 250, kFALSE);
  canv->cd(2);
  AliMpVMotif* motif = sector->FindMotif(TVector2(30,3));

  if (motif) {
    AliMpMotifType* motifType = motif->GetMotifType();
    MarkPads(AliMpMotifTypePadIterator(motifType),15,15);
    cout<<"______________ MotifType " << motifType->GetID() 
        <<"__________________________"<<endl;
  } else cout<<"No motif found at given position..."<<endl;
  
  canv->cd(3);
  //MarkPads(*AliMpPadIteratorPtr(AliMpSectorSegmentation(sector)->CreateIterator(AliMpIntPair(i,j)))
  AliMpSectorSegmentation segm(sector);
  AliMpPad pad = segm.PadByIndices(AliMpIntPair(i,j));
  MarkPads(AliMpNeighboursPadIterator(&AliMpSectorSegmentation(sector),pad)
           ,i+8,j+8);
  cout<<"________________ Neighbours __________________________"<<endl;
  canv->cd(4);
  Int_t motifPosId = 20 | AliMpConstants::ManuMask(plane); 
  if (plane == AliMp::kNonBendingPlane) motifPosId = 19;
  AliMpMotifPosition* motifPos = sector->GetMotifMap()->FindMotifPosition(motifPosId);
  if (motifPos){
    //MarkPads(*AliMpPadIteratorPtr(motifPos->CreateIterator()),15,15);
    MarkPads(AliMpMotifPositionPadIterator(motifPos),15,15);
    cout<<"_________________ MotifPosition _________________________"<<endl;
  }
}
