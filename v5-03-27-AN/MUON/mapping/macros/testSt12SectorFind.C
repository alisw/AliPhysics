// $Id$
// $MpId: testSectorFind.C,v 1.10 2005/10/28 15:36:08 ivana Exp $
//
// Test macro for which verify that all FindPosition, FindIndices
// and FindLocation methods are consistents between them.

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpSector.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpSectorReader.h"
#include "AliMpArea.h"
#include "AliMpVPadIterator.h"
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

void testSectorFind(AliMq::Station12Type station, AliMp::PlaneType plane) 
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSectorReader r(dataStreams, station, plane);
  AliMpSector* sector = r.BuildSector();
  AliMpSectorSegmentation segmentation(sector);
  
  CreateTCanvas("", "", station, plane);
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(sector);
  painter->Draw("ZSSMP");

  for (Int_t i=0; i<90;i++){
    cout<<"Verifying column "<<i<<"....."<<endl;

    for (Int_t j=0;j<230;++j)
      segmentation.CircleTest(i,j);    
  }
}

void testSt12SectorFind()
{
  AliMq::Station12Type  station[2] = { AliMq::kStation1, AliMq::kStation2 }; 
  AliMp::PlaneType      plane[2]   = { AliMp::kBendingPlane, AliMp::kNonBendingPlane };
  
  for ( Int_t is = 0; is < 2; is++ ) {
    for ( Int_t ip = 0; ip < 2; ip++ ) {
    
      cout << "Running testSectorFind for " 
           << AliMq::Station12TypeName(station[is]) << "  "
           << AliMp::PlaneTypeName(plane[ip])  << " ... " << endl;
       
      testSectorFind(station[is], plane[ip]);
    
      cout << "... end running " << endl << endl;
    }  
  }   
}  
  
