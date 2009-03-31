// $Id$
// $MpId: testAllIndices.C,v 1.7 2005/08/24 08:53:27 ivana Exp $
//
// Test macro for testing which pad is seen as "existing" by AliMpSector.

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpSector.h"
#include "AliMpSectorReader.h"
#include "AliMpSectorSegmentation.h" 
#include "AliMpVPainter.h" 
#include "AliMpRow.h"
#include "AliMpVRowSegment.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotifType.h"

#include <Riostream.h>
#include <TCanvas.h>
#include <TString.h>
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

void testAllIndices(AliMq::Station12Type station, AliMp::PlaneType plane) 
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSectorReader r(dataStreams, station, plane);

  AliMpSector *sector=r.BuildSector();
  AliMpSectorSegmentation segmentation(sector);
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(sector);

  TCanvas* c1 = CreateTCanvas("view ", 
                              "MSectorPainter::Draw() output (view per pad) ", 
                              station, plane);
  painter->Draw("ZSSMP");
  c1->Update();

  Int_t maxPadIndexX = segmentation.MaxPadIndexX();
  Int_t maxPadIndexY = segmentation.MaxPadIndexY();
  
  // Define histogram limits
  Int_t nx = (maxPadIndexX/10 + 1)*10;
  Int_t ny = (maxPadIndexY/10 + 1)*10;
  TH2C* histo  = new TH2C("histo","Existing pads", 
                          nx, -0.5, nx-0.5, ny, -0.5, ny-0.5);

  Int_t nx2 = 95/2;
  Int_t ny2 = 95/2;
  if (station == AliMq::kStation2) {
    nx2 = 120/2;
    ny2 = 120/2;
  }
  TH2F* histo2 = new TH2F("histo2","Existing positions",
                          nx2, 0, nx2*2, ny2, 0, ny2*2);

  // Define canvas
  TCanvas* c2 = CreateTCanvas("c2 ", "Only existing pads are filled ", station, plane);
  TCanvas* c3 = CreateTCanvas("c3 ", "Positions ", station, plane);
  
  for ( Int_t irow=0; irow<sector->GetNofRows(); irow++ ) {
    AliMpRow* row = sector->GetRow(irow);
    
    for ( Int_t  iseg=0; iseg<row->GetNofRowSegments(); iseg++ ) {
      AliMpVRowSegment* seg = row->GetRowSegment(iseg);
      
      for ( Int_t imot=0; imot<seg->GetNofMotifs(); imot++) {
        AliMpMotifPosition* motifPos 
         = sector->GetMotifMap()->FindMotifPosition(seg->GetMotifPositionId(imot));
         
        for ( Int_t gassNum=0; gassNum<64; gassNum++ ) {
          if (motifPos->GetMotif()->GetMotifType()->FindConnectionByGassiNum(gassNum)){
          
            AliMpPad pad = segmentation.PadByLocation(motifPos->GetID(),gassNum);
            if (pad != AliMpPad::Invalid()) {
              histo->Fill(pad.GetIx(), pad.GetIy());
              histo2->Fill(pad.Position().X(),
                           pad.Position().Y());
            }
          }
        }
      }
    }
  }
  c2->cd();
  histo->Draw("col");
  c3->cd();
  histo2->Draw("box");
}

void testSt12AllIndices()
{
  AliMq::Station12Type  station[2] = { AliMq::kStation1, AliMq::kStation2 }; 
  AliMp::PlaneType      plane[2]   = { AliMp::kBendingPlane, AliMp::kNonBendingPlane };
  
  for ( Int_t is = 0; is < 2; is++ ) {
    for ( Int_t ip = 0; ip < 2; ip++ ) {
    
      cout << "Running testAllIndices for " 
           << AliMq::Station12TypeName(station[is]) << "  "
           << AliMp::PlaneTypeName(plane[ip])  << " ... " << endl;
       
      testAllIndices(station[is], plane[ip]);
    
      cout << "... end running " << endl << endl;
    }  
  }   
}  
  
