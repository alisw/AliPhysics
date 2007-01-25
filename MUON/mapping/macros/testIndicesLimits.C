// $Id$
// $MpId: testIndicesLimits.C,v 1.5 2005/10/28 15:36:08 ivana Exp $
//
// Test macro for indices limits.

void testIndicesLimits(AliMp::StationType station = AliMp::kStation1,
                       AliMp::PlaneType plane = AliMp::kBendingPlane, 
		       Bool_t rootInput = false)
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

  // Loop over rows
  for (Int_t i=0; i<sector->GetNofRows(); i++) {
    AliMpRow* row = sector->GetRow(i);
    cout << i
         << "th row limits: "
	 << row->GetLowIndicesLimit() << "  "
	 << row->GetHighIndicesLimit() << endl;
    
    // Loop over row segments
    for (Int_t j=0; j<row->GetNofRowSegments(); j++) {
      AliMpVRowSegment* rowSegment = row->GetRowSegment(j);
      cout << "   "
           << j
           << "th row segment limits: "
	   << rowSegment->GetLowIndicesLimit() << "  "
	   << rowSegment->GetHighIndicesLimit() << endl;
      
      // Loop over motif positions
      for (Int_t k=0; k<rowSegment->GetNofMotifs(); k++) {
        Int_t mposID = rowSegment->GetMotifPositionId(k);
        AliMpMotifPosition* mPos = 
	  sector->GetMotifMap()->FindMotifPosition(mposID);     
        if (mPos) {
          cout << "      "
               << mPos->GetID()
               << " motif position limits: "
	       << mPos->GetLowIndicesLimit() << "  "
	       << mPos->GetHighIndicesLimit() << endl;
	}
	else {
	  cerr << "Motif position "
	       <<  mposID << " not found in the map" << endl; 
	}           
      }
    }  	   
  }   
  
  delete sector;
}			       
      

 
