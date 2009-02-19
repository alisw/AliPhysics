// $Id$
// $MpId: testExistingPads.C,v 1.12 2005/10/28 15:36:07 ivana Exp $
//
// Test macro for testing which pad is seen as "existing" by 
// by AliMpSectorSegmentation and AliMpFastSegmentation
// To run macro:
// root [0] .L testExistingPads.C+
// root [1] testExistingPads();

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpSector.h"
#include "AliMpSectorReader.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpFastSegmentation.h"
#include "AliMpArea.h"
#include "AliMpVPadIterator.h"
#include "AliMpVPainter.h"

#include <Riostream.h>
#include <TCanvas.h>
#include <TH2.h>

#include <iomanip>
#include <fstream>

#endif

void testExistingPads(AliMq::Station12Type station = AliMq::kStation1,
                      AliMp::PlaneType plane = AliMp::kBendingPlane)
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSectorReader r(dataStreams, station, plane);
  AliMpSector* sector = r.BuildSector();
  AliMpSectorSegmentation* segmentation = new AliMpSectorSegmentation(sector);
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(sector);

  TCanvas* c1 = new TCanvas("view",
                            "AliMpSectorPainter::Draw() output (view per pad)");
  painter->Draw("ZSSMP");
  c1->Update();

  Int_t maxPadIndexX = segmentation->MaxPadIndexX();
  Int_t maxPadIndexY = segmentation->MaxPadIndexY();
  
  // Define histogram limits
  Int_t nx = (maxPadIndexX/10 + 1)*10;
  Int_t ny = (maxPadIndexY/10 + 1)*10;
  TH2C* histo = new TH2C("histo","Existing pads", 
                          nx, -0.5, nx-0.5, ny, -0.5, ny-0.5);

  TCanvas* c2 = new TCanvas("c2","Only existing pads are filled");

  for (Int_t i=0; i<maxPadIndexX+1;i++){
    for (Int_t j=0;j<maxPadIndexY+1;++j){

      AliMpIntPair indices(i,j);
      if (segmentation->HasPad(indices)) histo->Fill(i,j);
    }
  }

  c2->cd();
  histo->Draw("box");

  // the same plot with fast segmentation
  TH2C* histo2 = new TH2C("histo2","Existing pads2", 
                          nx, -0.5, nx-0.5, ny, -0.5, ny-0.5);

  TCanvas* c3 = new TCanvas("c3","Only existing pads are filled");

  AliMpFastSegmentation* fast = new AliMpFastSegmentation(segmentation);
  for (Int_t i=0; i<maxPadIndexX+1;i++){
    for (Int_t j=0;j<maxPadIndexY+1;++j){

      AliMpIntPair indices(i,j);
      if (fast->HasPad(indices)) histo2->Fill(i,j);
    }
  }

  c3->cd();
  histo2->Draw("box");
  
  delete fast;
}
