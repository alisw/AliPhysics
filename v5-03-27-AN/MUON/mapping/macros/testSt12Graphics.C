// $Id$
// $MpId: testGraphics.C,v 1.13 2005/10/28 15:36:08 ivana Exp $
//
// Test macro for drawing sector data.

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpSector.h"
#include "AliMpSectorReader.h"
#include "AliMpRow.h"
#include "AliMpVRowSegment.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpVPainter.h"

#include <Riostream.h>
#include <TCanvas.h>

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

void testGraphics(AliMq::Station12Type station, AliMp::PlaneType plane) 
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSectorReader r(dataStreams, station, plane);
  AliMpSector* sector = r.BuildSector();
    
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(sector);

  TCanvas* canvas  = new TCanvas();
  TCanvas* canvas2 = new TCanvas();
  TCanvas* c[4];
  for (int i=0;i<4;++i) {
    c[i] = new TCanvas();
    c[i]->Divide(2,2);
  }  

  //first, paint the whole sector
  canvas->cd();
  painter->Draw("");
  //now paint each rows
  c[0]->cd(1);
  painter->Draw("R");
  //paint each row segments in each row
  c[0]->cd(2);
  painter->Draw("RS");
  //paint each motifs, in each row segments in each row
  c[0]->cd(3);
  painter->Draw("RSMP");
  //paint each pads, in each motifs, in each row segments in each row
  c[0]->cd(4);
  painter->Draw("RSMT");

  ///////////////////////////////
  //now paint each rows, wwith its name
  c[1]->cd(1);
  painter->Draw("RT");
  //paint each row segments in each row, and its name
  c[1]->cd(2);
  painter->Draw("RST");
  //paint each motifs, in each row segments in each row
  c[1]->cd(3);
  painter->Draw("RSMX");
  c[1]->cd(4);
  painter->Draw("RSMI");

  ///////////////////////////////
  //now paint each zones
  c[2]->cd(1);
  painter->Draw("Z");
  //paint each sub-zones, in each zones
  c[2]->cd(2);
  painter->Draw("ZS");
  //paint each row segments, in each sub-zone, ...
  c[2]->cd(3);
  painter->Draw("ZSS");
  // each motifs, in each row segments, ...
  c[2]->cd(4);
  painter->Draw("ZSSM");

  ///////////////////////////////
  //now paint each zones with its name
  c[3]->cd(1);
  painter->Draw("ZT");
  //paint each sub-zones, in each zones with its name
  c[3]->cd(2);
  painter->Draw("ZST");
  //paint each row segments, in each sub-zone, ... with its name
  c[3]->cd(3);
  painter->Draw("ZSST");
  // each motifs, in each row segments, ... with its name
  c[3]->cd(4);
  painter->Draw("ZSSMT");

  // now, draw a specific motif, in a whole canvas, and
  // print all its pad names
  Int_t id = sector->GetRow(5)->GetRowSegment(0)->GetMotifPositionId(0);
  AliMpMotifPosition* motifPos = sector->GetMotifMap()->FindMotifPosition(id);
  AliMpVPainter* motifPainter = AliMpVPainter::CreatePainter(motifPos);
  CreateTCanvas("onepad ","One motif ", station, plane);
  motifPainter->Draw("PT");
  
  //////////////////////////////  
  //now paint motifs with their real contours and mani Ids
  canvas2->cd();
  painter->Draw("RSMCI");
}

void testSt12Graphics()
{
  AliMq::Station12Type  station[2] = { AliMq::kStation1, AliMq::kStation2 }; 
  AliMp::PlaneType      plane[2]   = { AliMp::kBendingPlane, AliMp::kNonBendingPlane };
  
  for ( Int_t is = 0; is < 2; is++ ) {
    for ( Int_t ip = 0; ip < 2; ip++ ) {
    
      cout << "Running testGraphics for " 
           << AliMq::Station12TypeName(station[is]) << "  "
           << AliMp::PlaneTypeName(plane[ip])  << " ... " << endl;
       
      testGraphics(station[is], plane[ip]);
    
      cout << "... end running " << endl << endl;
    }  
  }   
}  
  
