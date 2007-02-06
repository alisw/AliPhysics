// $Id$
// $MpId: testReadSector.C,v 1.16 2006/01/11 10:00:50 ivana Exp $
//
// Test macro for reading sector data.

#include <iomanip>

void testReadSector(AliMp::StationType station = AliMp::kStation1,
                    AliMp::PlaneType plane = AliMp::kBendingPlane, 
	 	    Bool_t rootInput = false)
{
  AliMpSector *sector = 0;
  if (!rootInput) {
    AliMpSectorReader r(station, plane);
    //reader.SetVerboseLevel(1);
    sector=r.BuildSector();

    // Write sector on Root file
    TString filePath = AliMpFiles::SectorFilePath(station,plane);
    filePath.ReplaceAll("zones.dat", "sector.root"); 
  
    TFile f(filePath.Data(), "RECREATE");
    sector->Write();
    f.Close();
  }
  else  {
    TString filePath = AliMpFiles::SectorFilePath(station,plane);
    filePath.ReplaceAll("zones.dat", "sector.root"); 

    TFile f(filePath.Data(), "READ");
    sector = (AliMpSector*)f.Get("Sector");
  }  

  cout << endl;

  // Sector geometry
  sector->PrintGeometry();

  cout << endl;

  // Find row test
  if (plane == AliMp::kBendingPlane)
    cout << "0th row low border " << sector->FindRow(TVector2(0., 0.))->GetID() << endl;
  if (plane == AliMp::kNonBendingPlane)
    cout << "0th row low border " << sector->FindRow(TVector2(0., 0.215))->GetID() << endl;
  cout << "in 0th row         " << sector->FindRow(TVector2(0., 2.5))->GetID() << endl;
  cout << "0th row up border  " << sector->FindRow(TVector2(0., 6.72))->GetID() << endl;
  cout << "in 4th row         " << sector->FindRow(TVector2(0., 30.))->GetID() << endl;
  if (plane == AliMp::kBendingPlane)
    cout << "12th row up border " << sector->FindRow(TVector2(0., 89.46))->GetID() << endl;
  if (plane == AliMp::kNonBendingPlane)
    cout << "12th row up border " << sector->FindRow(TVector2(0., 84.84))->GetID() << endl;
  cout << endl;
  
  // Find motif position test
  Int_t ids[15] = { 19, 14, 9, 32, 36, 136, 187, 212, 207, 220, 1, 131, 239, 243, 250 };  
  for (Int_t i=0; i<15 ; i++) {
    Int_t id = ids[i];
    id |= AliMpConstants::ManuMask(plane);
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
    TVector2 pos(0.5, 18.6 - i*2.);  // i=0 in motif 1001,
                                     // i=1 outside (below) motif 1001
    AliMpMotif* motif = sector->FindMotif(pos);
    cout << "In the position " << pos.X() << " " << pos.Y();
    
    if (motif)
      cout << " found motif " << motif->GetID() << endl;
    else  
      cout << " motif not found " << endl;
  }

  // Find special motif test
  if (plane == AliMp::kNonBendingPlane) {
  
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
      

 
