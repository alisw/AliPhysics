// $Id$
// $MpId: testReadSector.C,v 1.13 2005/09/02 10:57:47 ivana Exp $
//
// Test macro for reading sector data.

#include <iomanip>

void testReadSector(AliMpStationType station = kStation1,
                    AliMpPlaneType plane = kBendingPlane) 
{
  AliMpSectorReader reader(station, plane);  
  //reader.SetVerboseLevel(1);
  
  // Read data 
  AliMpSector* sector = reader.BuildSector();
  
  cout << endl;

  // Sector geometry
  sector->PrintGeometry();

  cout << endl;

  // Find row test
  if (plane == kBendingPlane)
    cout << "0th row low border " << sector->FindRow(TVector2(0., 0.))->GetID() << endl;
  if (plane == kNonBendingPlane)
    cout << "0th row low border " << sector->FindRow(TVector2(0., 0.215))->GetID() << endl;
  cout << "in 0th row         " << sector->FindRow(TVector2(0., 25.))->GetID() << endl;
  cout << "0th row up border  " << sector->FindRow(TVector2(0., 67.2))->GetID() << endl;
  cout << "in 4th row         " << sector->FindRow(TVector2(0., 300.))->GetID() << endl;
  if (plane == kBendingPlane)
    cout << "12th row up border " << sector->FindRow(TVector2(0., 894.6))->GetID() << endl;
  if (plane == kNonBendingPlane)
    cout << "12th row up border " << sector->FindRow(TVector2(0., 848.4))->GetID() << endl;
  cout << endl;
  
  // Find motif position test
  Int_t ids[15] = { 19, 14, 9, 32, 36, 136, 187, 212, 207, 220, 1, 131, 239, 243, 250 };  
  for (Int_t i=0; i<15 ; i++) {
    Int_t id = ids[i];
    cout << "Motif pos " << std::setw(3) << id;
    if (!sector->FindRowSegment(id)) {
       cout << " not found." << endl;
    }
    else {	 
       cout << " found in : "
            << sector->FindRow(id)->GetID() << " row, "
            << " motif id: "
            << sector->FindRowSegment(id)->GetMotif(0)->GetID().Data()
	    << endl;
    }      
  }
  cout << endl;

  // Find motif by coordinates test
  for (Int_t i=0; i<2 ; i++) {
    TVector2 pos(5., 186. - i*20.);  // i=0 in motif 1001,
                                     // i=1 outside (below) motif 1001
    AliMpMotif* motif = sector->FindMotif(pos);
    cout << "In the position " << pos.X() << " " << pos.Y();
    
    if (motif)
      cout << " found motif " << motif->GetID() << endl;
    else  
      cout << " motif not found " << endl;
  }

  // Find special motif test
  if (plane == kNonBendingPlane) {
  
    Int_t ids[6] = { 20, 46, 47, 74, 75, 76 };
    for (Int_t i=0; i<6 ; i++) {
      
      Int_t id = ids[i];
      cout << "Motif pos " << id;
      if (!sector->FindRowSegment(id)) {
         cout << " not found." << endl;
      }
      else {	 
         cout << " found in : "
              << sector->FindRow(id)->GetID() << " row, "
	      << " position : "
	      << sector->FindPosition(id).X() << "  " <<  sector->FindPosition(id).Y()
	      << endl;
      }
    }
  }             
  cout << endl;

  // Motif map
  sector->GetMotifMap()->Print();    
  
  delete sector;
}			       
      

 
