// $Id$
// $MpId: testAnyPadIterators.C,v 1.12 2005/08/24 08:53:27 ivana Exp $
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
    if (print) cout<<"Iterator number "<< num++ << " at "<< indices <<endl;
    TMarker* marker = new TMarker( (Double_t)indices.GetFirst() /xmax,
                                   (Double_t)indices.GetSecond()/ymax,
                                   2);
    marker->Draw();
  }
}

void testAnyPadIterators(AliMpStationType station = kStation1,
                         AliMpPlaneType plane = kBendingPlane, 
                         Int_t i=50, Int_t j=50)
{
  AliMpSectorReader r(station, plane);
  AliMpSector* sect = r.BuildSector();
    
  TCanvas *canv = new TCanvas("canv");
  canv->Divide(2,2);
  //canv_1->Divide(2);
  
  canv->cd(1);
  MarkPads(AliMpSectorPadIterator(sect), 150, 250, kFALSE);
  canv->cd(2);
  AliMpVMotif* motif = sect->FindMotif(TVector2(300,30));

  if (motif) {
    AliMpMotifType* motifType = motif->GetMotifType();
    MarkPads(AliMpMotifTypePadIterator(motifType),15,15);
    cout<<"______________ MotifType " << motifType->GetID() 
        <<"__________________________"<<endl;
  } else cout<<"No motif found at given position..."<<endl;
  
  canv->cd(3);
  //MarkPads(*AliMpPadIteratorPtr(AliMpSectorSegmentation(sect)->CreateIterator(AliMpIntPair(i,j)))
  AliMpSectorSegmentation segm(sect);
  AliMpPad pad = segm.PadByIndices(AliMpIntPair(i,j));
  MarkPads(AliMpNeighboursPadIterator(&AliMpSectorSegmentation(sect),pad)
           ,i+8,j+8);
  cout<<"________________ Neighbours __________________________"<<endl;
  canv->cd(4);
  Int_t motifPosId = 20; 
  if (plane == kNonBendingPlane) motifPosId = 19;
  AliMpMotifPosition* motifPos = sect->GetMotifMap()->FindMotifPosition(motifPosId);
  if (motifPos){
    //MarkPads(*AliMpPadIteratorPtr(motifPos->CreateIterator()),15,15);
    MarkPads(AliMpMotifPositionPadIterator(motifPos),15,15);
    cout<<"_________________ MotifPosition _________________________"<<endl;
  }
}
