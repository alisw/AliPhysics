//==========================================================================================================================================

// Macro to create the config file (AliMFTGeometry.root) for the geometry of the ALICE Muon Forward Tracker
//
// Contact author: antonio.uras@cern.ch

//==========================================================================================================================================

void SetMFTGeometry() {

  const Int_t nPlanes = 5;
  
  const Float_t zCenter[nPlanes]          = {  -46.0,   -49.3,   -53.1,   -68.7,   -76.8 };   // expressed in cm
				          
  const Float_t rMin[nPlanes]             = {   2.30,    2.30,    2.30,    3.30,    3.60 };   // expressed in cm  
  const Float_t rMax[nPlanes]             = {  10.00,   10.50,   11.10,   13.90,   15.30 };   // expressed in cm

  const Float_t pixelSizeX[nPlanes]       = { 28.e-4,  28.e-4,  28.e-4,  28.e-4,  28.e-4 };   // expressed in cm
  const Float_t pixelSizeY[nPlanes]       = { 28.e-4,  28.e-4,  28.e-4,  28.e-4,  28.e-4 };   // expressed in cm

  const Float_t thicknessActive[nPlanes]  = {  50.e-4,   50.e-4,   50.e-4,   50.e-4,   50.e-4 };   // expressed in cm
  const Float_t thicknessSupport[nPlanes] = {   1.4  ,    1.4  ,    1.4  ,    1.4  ,    1.4   };   // expressed in cm
  const Float_t thicknessReadout[nPlanes] = {  50.e-4,   50.e-4,   50.e-4,   50.e-4,   50.e-4 };   // expressed in cm

  const Float_t equivalentSilicon[nPlanes]            = { 600.e-4, 600.e-4, 600.e-4, 600.e-4, 600.e-4 };    // expressed in cm
  const Float_t equivalentSiliconBeforeFront[nPlanes] = {   0.e-4,   0.e-4,   0.e-4,   0.e-4,   0.e-4 };    // expressed in cm
  const Float_t equivalentSiliconBeforeBack[nPlanes]  = { 550.e-4, 550.e-4, 550.e-4, 550.e-4, 550.e-4 };    // expressed in cm

  const Float_t hasPixelRectangularPatternAlongY[nPlanes] = {0., 0., 0., 0., 0.};
				         
  TNtuple *geomMFT = new TNtuple("AliMFTGeometry", "ALICE MFT Geometry", "zCenter:rMin:rMax:pixelSizeX:pixelSizeY:thicknessActive:thicknessSupport:thicknessReadout:equivalentSilicon:equivalentSiliconBeforeFront:equivalentSiliconBeforeBack:hasPixelRectangularPatternAlongY");

  for (Int_t iPlane=0; iPlane<nPlanes; iPlane++) geomMFT -> Fill(zCenter[iPlane],
								 rMin[iPlane],
								 rMax[iPlane],
								 pixelSizeX[iPlane],
								 pixelSizeY[iPlane],
								 thicknessActive[iPlane],
								 thicknessSupport[iPlane],
								 thicknessReadout[iPlane],
								 equivalentSilicon[iPlane],
								 equivalentSiliconBeforeFront[iPlane],
								 equivalentSiliconBeforeBack[iPlane],
								 hasPixelRectangularPatternAlongY[iPlane]);

  TFile *fileGeomMFT = new TFile("AliMFTGeometry.root", "recreate");
  geomMFT -> Write();
  fileGeomMFT -> Close();

}

//==========================================================================================================================================
