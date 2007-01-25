// $Id$
// $MpId: testPrintLimits.C,v 1.9 2005/10/28 15:36:08 ivana Exp $
//
// Test macro for making an output file, where all mapping elements
// indices & positions are written.

void testPrintLimits(AliMp::StationType station = AliMp::kStation1,
                     AliMp::PlaneType plane = AliMp::kBendingPlane, 
		     Bool_t rootInput = false,
		     ostream& out=cout)
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
