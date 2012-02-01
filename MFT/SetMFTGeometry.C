//==========================================================================================================================================

// Macro to create the config file (AliMFTGeometry.root) for the geometry of the ALICE Muon Forward Tracker
//
// Contact author: antonio.uras@cern.ch

//==========================================================================================================================================

void SetMFTGeometry() {

  const Int_t nPlanes = 5;
  
  const Float_t zCenter[nPlanes]          = {  -50.0,   -58.0,   -66.0,   -72.0,   -76.0 };   // expressed in cm
				          
  const Float_t rMin[nPlanes]             = {   2.20,    2.60,    3.00,    3.30,    3.60 };   // expressed in cm  
  const Float_t rMax[nPlanes]             = {   9.70,   11.00,   12.40,   13.40,   14.00 };   // expressed in cm
				          
  const Float_t pixelSizeX[nPlanes]       = { 20.e-4,  20.e-4,  20.e-4,  20.e-4,  20.e-4 };   // expressed in cm
  const Float_t pixelSizeY[nPlanes]       = { 20.e-4,  20.e-4,  20.e-4,  20.e-4,  20.e-4 };   // expressed in cm

  const Float_t thicknessActive[nPlanes]  = {  50.e-4,   50.e-4,   50.e-4,   50.e-4,   50.e-4 };   // expressed in cm
  const Float_t thicknessSupport[nPlanes] = {2000.e-4, 2000.e-4, 2000.e-4, 2000.e-4, 2000.e-4 };   // expressed in cm
  const Float_t thicknessReadout[nPlanes] = {  50.e-4,   50.e-4,   50.e-4,   50.e-4,   50.e-4 };   // expressed in cm

  const Float_t equivalentSilicon[nPlanes]            = { 300.e-4, 300.e-4, 300.e-4, 300.e-4, 300.e-4};    // expressed in cm
  const Float_t equivalentSiliconBeforeFront[nPlanes] = {   0.e-4,   0.e-4,   0.e-4,   0.e-4,   0.e-4};    // expressed in cm
  const Float_t equivalentSiliconBeforeBack[nPlanes]  = { 250.e-4, 250.e-4, 250.e-4, 250.e-4, 250.e-4};    // expressed in cm

  TNtuple *geomMFT = new TNtuple("AliMFTGeometry", "ALICE MFT Geometry", "zCenter:rMin:rMax:pixelSizeX:pixelSizeY:thicknessActive:thicknessSupport:thicknessReadout:equivalentSilicon:equivalentSiliconBeforeFront:equivalentSiliconBeforeBack");

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
								 equivalentSiliconBeforeBack[iPlane]);

  TFile *fileGeomMFT = new TFile("AliMFTGeometry.root", "recreate");
  geomMFT -> Write();
  fileGeomMFT -> Close();

}

//==========================================================================================================================================
