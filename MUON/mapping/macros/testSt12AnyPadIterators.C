// $Id$
// $MpId: testAnyPadIterators.C,v 1.17 2006/03/15 13:07:07 ivana Exp $
//
// Test macro for reading  sector, and iterate over it

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpSector.h"
#include "AliMpSectorReader.h"
#include "AliMpSectorSegmentation.h" 
#include "AliMpMotifType.h"
#include "AliMpMotifMap.h"
#include "AliMpVMotif.h"
#include "AliMpVPadIterator.h"
#include "AliMpSectorPadIterator.h"
#include "AliMpNeighboursPadIterator.h"
#include "AliMpMotifPositionPadIterator.h"
#include "AliMpConstants.h"

#include <Riostream.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TH2.h>

#endif

TCanvas* CreateTCanvas(const TString& name, const TString& title,
                       AliMq::Station12Type station, AliMp::PlaneType plane)
{
  TString newName(name);
  TString newTitle(title);
  TString unique = AliMq::Station12TypeName(station) + AliMp::PlaneTypeName(plane);
  newName += unique;
  newTitle += unique;
  return new TCanvas(newName.Data(), newTitle.Data());
}                     

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

void testAnyPadIterators(AliMq::Station12Type station, AliMp::PlaneType plane,
                         Int_t i=50, Int_t j=50)
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSectorReader r(dataStreams, station, plane);
  AliMpSector* sector = r.BuildSector();
  AliMpSectorSegmentation segmentation(sector);
    
  TCanvas* canv = CreateTCanvas("canv ", "", station, plane);
  canv->Divide(2,2);
  
  canv->cd(1);
  AliMpSectorPadIterator its(sector);
  MarkPads(its, 150, 250, kFALSE);

  canv->cd(2);
  AliMpVMotif* motif = sector->FindMotif(TVector2(30,3));
  if ( motif ) {
    AliMpMotifType* motifType = motif->GetMotifType();
    AliMpMotifTypePadIterator itm(motifType);
    MarkPads(itm,15,15);
    cout << "______________ MotifType "  << motifType->GetID() 
         << "__________________________" << endl;
  } 
  else 
    cout << "No motif found at given position..." << endl;
  
  canv->cd(3);
  AliMpPad pad = segmentation.PadByIndices(AliMpIntPair(i,j));
  AliMpNeighboursPadIterator itn(&segmentation,pad);
  MarkPads(itn,i+8,j+8);
  cout<<"________________ Neighbours __________________________"<<endl;
  
  canv->cd(4);
  Int_t motifPosId = 20 | AliMpConstants::ManuMask(plane); 
  if ( plane == AliMp::kNonBendingPlane ) motifPosId = 19;
  AliMpMotifPosition* motifPos = sector->GetMotifMap()->FindMotifPosition(motifPosId);
  if ( motifPos ){
    AliMpMotifPositionPadIterator itmp(motifPos);
    MarkPads(itmp,15,15);
    cout<<"_________________ MotifPosition _________________________"<<endl;
  }
}

void testSt12AnyPadIterators()
{
  AliMq::Station12Type  station[2] = { AliMq::kStation1, AliMq::kStation2 }; 
  AliMp::PlaneType      plane[2]   = { AliMp::kBendingPlane, AliMp::kNonBendingPlane };
  
  for ( Int_t is = 0; is < 2; is++ ) {
    for ( Int_t ip = 0; ip < 2; ip++ ) {
    
      cout << "Running testAnyPadIterators for " 
           << AliMq::Station12TypeName(station[is]) << "  "
           << AliMp::PlaneTypeName(plane[ip])  << " ... " << endl;
       
      testAnyPadIterators(station[is], plane[ip]);
    
      cout << "... end running " << endl << endl;
    }  
  }   
}  
  
