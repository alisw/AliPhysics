// $Id$
//
// Test macro for indices limits.

void testIndicesLimits(AliMpStationType station = kStation1,
                       AliMpPlaneType plane = kBendingPlane) 
{
  AliMpReader reader(station, plane);  
  //reader.SetVerboseLevel(1);
  
  // Read data 
  AliMpSector* sector = reader.BuildSector();

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
      

 
