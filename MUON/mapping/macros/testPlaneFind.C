// $Id$
//
// Test macro that runs the segmentation circle test
// (finds pad by location, position and indices in a circle)
// over the whole plane.

void testPlane(AliMpPlane* plane)
{
  AliMpPlaneSegmentation planeSegmentation(plane);
  
  //AliMpVPainter* painter = AliMpVPainter::CreatePainter(sector);
  //painter->Draw("ZSSMP");

  for (Int_t isec=0; isec<4; isec++) {
             // 0  x>=0, y>=0
             // 1  x>=0, y <0
             // 2  x <0, y>=0
             // 3  x <0, y <0

    cout << "Verifying sector " << isec << "....." << endl;
    
    for (Int_t i = AliMpConstants::StartPadIndex(); i<AliMpConstants::StartPadIndex()+90; i++){
     
      Int_t iscale = 1;
      if (isec >1) iscale = -1;
    
      cout << "Verifying column " << i << "....." << endl;

      for (Int_t j = AliMpConstants::StartPadIndex(); j<AliMpConstants::StartPadIndex()+20; j++) {

        Int_t jscale = 1;  
        if (isec == 1 || isec == 3) jscale = -1;
	
	if (planeSegmentation.HasPad(AliMpIntPair(i*iscale, j*jscale))) {
	
	  planeSegmentation.PadByIndices(AliMpIntPair(i*iscale, j*jscale)).Print();
          cout << "test result " 
	       << planeSegmentation.CircleTest(AliMpIntPair(i*iscale, j*jscale))
	       << endl;    
	}
	else {
	  cout << " has not indices  " << AliMpIntPair(i*iscale, j*jscale) << endl;
	}  
      }
    }
  }  
}

void testPlaneFind(AliMpStationType station = kStation1,
                   AliMpPlaneType planeType = kBendingPlane)
{
  cout << endl;

  cout << "Testing plane1 ..." << endl;
  AliMpPlane* plane1 = AliMpPlane::Create(station, planeType);
  testPlane(plane1);
  delete plane1;
  cout << endl; 

  cout << "Testing plane2 ..." << endl; 
  AliMpPlane* plane2 
    = AliMpPlane::Create(station, planeType, 
                         TVector2(), TVector2(-10., 10.), TVector2(), TVector2());
  testPlane(plane2);
  delete plane2;  
}
