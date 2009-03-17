// $Id$
// $MpId: testPrintLimits.C,v 1.9 2005/10/28 15:36:08 ivana Exp $
//
// Test macro for making an output file, where all mapping elements
// indices & positions are written.

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

void testPrintLimits(AliMq::Station12Type station, AliMp::PlaneType  plane,
		     ostream& out=cout)
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSectorReader r(dataStreams, station, plane);
  AliMpSector* sector = r.BuildSector();
  
  new TCanvas();

  AliMpVPainter* painter = AliMpVPainter::CreatePainter(sector);
  painter->Draw("ZSSMP");

  AliMpIntPair low,high;
  TVector2 rlow,rhigh;

  for (Int_t irow=0;irow<sector->GetNofRows();irow++){
    // For each row

    AliMpRow* row = sector->GetRow(irow);
    low  = row->GetLowIndicesLimit();
    high = row->GetHighIndicesLimit();
    rlow = TVector2(row->Position().X()-row->Dimensions().X(),
                    row->Position().Y()-row->Dimensions().Y());
    rhigh = TVector2(row->Position().X()+row->Dimensions().X(),
                     row->Position().Y()+row->Dimensions().Y());
    out<<"_______________________________________________________________"<<endl;
    out<<"Row "<<irow<<" between "<<low<<" and "<<high
       <<"-->("<<rlow.X()<<','<<rlow.Y()<<") and ("
       <<rhigh.X()<<','<<rhigh.Y()<<')'<<endl;
    out<<"_______________________________________________________________"<<endl;
  
    for (Int_t iseg=0;iseg<row->GetNofRowSegments();iseg++){
      // For each row segment  
  
      AliMpVRowSegment* seg = row->GetRowSegment(iseg);
      low  = seg->GetLowIndicesLimit();
      high = seg->GetHighIndicesLimit();
      rlow = TVector2(seg->Position().X()-seg->Dimensions().X(),
                      seg->Position().Y()-seg->Dimensions().Y());
      rhigh = TVector2(seg->Position().X()+seg->Dimensions().X(),
                       seg->Position().Y()+seg->Dimensions().Y());
      out<<"-----------------------------------------------------------"<<endl;
      out<<"     Segment "<<iseg<<" between "<<low<<" and "<<high
         <<"-->("<<rlow.X()<<','<<rlow.Y()<<") and ("
         <<rhigh.X()<<','<<rhigh.Y()<<')'<<endl;
      out<<"-----------------------------------------------------------"<<endl;

      for (Int_t imotif=0;imotif<seg->GetNofMotifs();imotif++){
         // For each motif pos
  
        AliMpMotifPosition* motifPos 
          = sector->GetMotifMap()
	      ->FindMotifPosition(seg->GetMotifPositionId(imotif));
        AliMpVMotif* motif = motifPos->GetMotif();
      
        low  = motifPos->GetLowIndicesLimit();
        high = motifPos->GetHighIndicesLimit();
        rlow = TVector2(motifPos->Position().X()-motif->Dimensions().X(),
                  motifPos->Position().Y()-motif->Dimensions().Y());
        rhigh = TVector2(motifPos->Position().X()+motif->Dimensions().X(),
                   motifPos->Position().Y()+motif->Dimensions().Y());
        out<<"          Motif "<<imotif<<" between "<<low<<" and "<<high
           <<"-->("<<rlow.X()<<','<<rlow.Y()<<") and ("
           <<rhigh.X()<<','<<rhigh.Y()<<')'<<endl;
      }
    }
  }
}

void testPrintLimits()
{
  AliMq::Station12Type  station[2] = { AliMq::kStation1, AliMq::kStation2 }; 
  AliMp::PlaneType      plane[2]   = { AliMp::kBendingPlane, AliMp::kNonBendingPlane };
  
  for ( Int_t is = 0; is < 2; is++ ) {
    for ( Int_t ip = 0; ip < 2; ip++ ) {
    
      cout << "Running testPrintLimits for " 
           << AliMq::Station12TypeName(station[is]) << "  "
           << AliMp::PlaneTypeName(plane[ip])  << " ... " << endl;
       
      testPrintLimits(station[is], plane[ip]);
    
      cout << "... end running " << endl << endl;
    }  
  }   
}  
  
