#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#endif

//---------------------------------------
void CreateITSU()
{
  //
  gSystem->Load("libITSUpgradeBase.so");
  gSystem->Load("libITSUpgradeSim.so");


  // build ITS upgrade detector
  // sensitive area 13x15mm (X,Z) with 20x20 micron pitch, 2mm dead zone on readout side and 50 micron guardring
  const double kSensThick = 18e-4;
  const double kPitchX = 20e-4;
  const double kPitchZ = 20e-4;
  const int    kNRow   = 650; 
  const int    kNCol   = 750;
  const int    kNChips = 2;
  const double kLrTick03 = 195e-4;   // -> effective thickness for ~0.3%X layers
  const double kLrTick08 = 600e-4;   // -> effective thickness for ~0.8%X layers
  //
  const double kReadOutEdge = 0.2;   // width of the readout edge (passive bottom)
  const double kGuardRing   = 50e-4; // width of passive area on left/right/top of the sensor
  // create segmentations:
  AliITSUSegmentationPix* seg0 = new AliITSUSegmentationPix(0,        // segID (0:9)
							    kNChips,  // chips per module
							    kNChips*kNCol,    // ncols (total for module)
							    kNRow,    // nrows
							    kPitchX,  // default row pitch in cm
							    kPitchZ,  // default col pitch in cm
							    kSensThick,  // sensor thickness in cm
							    -1,     // no special left col between chips
							    -1,     // no special right col between chips
							    kGuardRing, // left
							    kGuardRing, // right
							    kGuardRing, // top
							    kReadOutEdge  // bottom
							    );    // see AliITSUSegmentationPix.h for extra options
  seg0->Store(AliITSUGeomTGeo::GetITSsegmentationFileName());
  //
  AliITSUSegmentationPix* seg1 = new AliITSUSegmentationPix(1,        // segID (0:9)
							    kNChips,  // chips per module
							    kNChips*kNCol,    // ncols (total for module)
							    2*kNRow,    // nrows for oute layers
							    kPitchX,  // default row pitch in cm
							    kPitchZ,  // default col pitch in cm
							    kSensThick,  // sensor thickness in cm
							    -1,     // no special left col between chips
							    -1,     // no special right col between chips
							    kGuardRing, // left
							    kGuardRing, // right
							    kReadOutEdge, // top   !!! readout from both sides
							    kReadOutEdge  // bottom
							    );    // see AliITSUSegmentationPix.h for extra options
  seg1->Store(AliITSUGeomTGeo::GetITSsegmentationFileName());
  //
  seg0->Print();
  seg1->Print();
  //
  const double kMinOvl = 0.005; // require active zones overlap
  const double kPhi0 = 0.;  // az.angle of 1st stave
  const double kTilt = 10.; // tilt in degrees
  double dzLr,rLr,ovlA,xActProj;
  AliITSUSegmentationPix* seg=0;
  int nStaveLr,nModPerStaveLr,idLr;
  //      virtual void   DefineLayerTurbo(const Int_t nlay, const Double_t r,  const Double_t zlen, const Int_t nladd,   const Int_t nmod, const Double_t width,
  //				  const Double_t tilt,   const Double_t lthick = 0.,    const Double_t dthick = 0.,   const UInt_t detType=0);
  AliITSUv0 *ITS  = new AliITSUv0("ITS Upgrade",7);
  ITS->SetStaveModel(AliITSUv0::kModel22);
  //
  // INNER LAYERS
  idLr = 0;
  rLr = 2.2;
  dzLr = 2*11.2;   // min Z to cover
  seg = seg0;
  nModPerStaveLr = 1+dzLr/seg->Dz();
  ovlA = -1;
  xActProj = seg->DxActive()*TMath::Cos(kTilt*TMath::DegToRad()); // effective r-phi coverage by single stave
  nStaveLr = 1 + rLr*TMath::Pi()*2/xActProj;
  do { ovlA = 1.-rLr*TMath::Pi()*2/nStaveLr/xActProj; } while ( kMinOvl>=0 && ovlA<kMinOvl && nStaveLr++ );		
  ITS->DefineLayerTurbo(idLr, kPhi0, rLr, nModPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, seg->Dx(), kTilt, kLrTick03, seg->Dy(), seg->GetDetTypeID());
  printf("Add Lr%d: R=%.1f DZ:%.1f Staves:%3d NMod/Stave:%3d -> Active Overlap:%.1f (%d micron)\n",
	 idLr,rLr,nModPerStaveLr*seg->Dz()/2,nStaveLr,nModPerStaveLr,ovlA*100,int(ovlA*xActProj*1.e4));
  //
  idLr = 1;
  rLr = 2.8;
  dzLr = 2*12.1;
  seg = seg0;
  nModPerStaveLr = 1+dzLr/seg->Dz();
  ovlA = -1;
  xActProj = seg->DxActive()*TMath::Cos(kTilt*TMath::DegToRad()); // effective r-phi coverage by single stave
  nStaveLr = 1 + rLr*TMath::Pi()*2/xActProj;
  do { ovlA = 1.-rLr*TMath::Pi()*2/nStaveLr/xActProj; } while ( kMinOvl>=0 && ovlA<kMinOvl && nStaveLr++ );		
  ITS->DefineLayerTurbo(idLr, kPhi0, rLr, nModPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, seg->Dx(), kTilt, kLrTick03, seg->Dy(), seg->GetDetTypeID());
  printf("Add Lr%d: R=%.1f DZ:%.1f Staves:%3d NMod/Stave:%3d -> Active Overlap:%.1f (%d micron)\n",
	 idLr,rLr,nModPerStaveLr*seg->Dz()/2,nStaveLr,nModPerStaveLr,ovlA*100,int(ovlA*xActProj*1e4));
  //
  idLr = 2;
  rLr = 3.6;
  dzLr = 2*13.4;
  seg = seg0;
  nModPerStaveLr = 1+dzLr/seg->Dz();
  ovlA = -1;
  xActProj = seg->DxActive()*TMath::Cos(kTilt*TMath::DegToRad()); // effective r-phi coverage by single stave
  nStaveLr = 1 + rLr*TMath::Pi()*2/xActProj;
  do { ovlA = 1.-rLr*TMath::Pi()*2/nStaveLr/xActProj; } while ( kMinOvl>=0 && ovlA<kMinOvl && nStaveLr++ );		
  ITS->DefineLayerTurbo(idLr, kPhi0, rLr, nModPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, seg->Dx(), kTilt, kLrTick03, seg->Dy(), seg->GetDetTypeID());
  printf("Add Lr%d: R=%.1f DZ:%.1f Staves:%3d NMod/Stave:%3d -> Active Overlap:%.1f (%d micron)\n",
	 idLr,rLr,nModPerStaveLr*seg->Dz()/2,nStaveLr,nModPerStaveLr,ovlA*100,int(ovlA*xActProj*1e4));
  //
  // 
  // MIDDLE LAYERS (double side readout sensors)
  idLr = 3;
  rLr = 20.0;
  dzLr = 2*39.0;
  seg = seg1;
  nModPerStaveLr = 1+dzLr/seg->Dz();
  ovlA = -1;
  xActProj = seg->DxActive()*TMath::Cos(kTilt*TMath::DegToRad()); // effective r-phi coverage by single stave
  nStaveLr = 1 + rLr*TMath::Pi()*2/xActProj;
  do { ovlA = 1.-rLr*TMath::Pi()*2/nStaveLr/xActProj; } while ( kMinOvl>=0 && ovlA<kMinOvl && nStaveLr++ );		
  ITS->DefineLayerTurbo(idLr, kPhi0, rLr, nModPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, seg->Dx(), kTilt, kLrTick08, seg->Dy(), seg->GetDetTypeID());
  printf("Add Lr%d: R=%.1f DZ:%.1f Staves:%3d NMod/Stave:%3d -> Active Overlap:%.1f (%d micron)\n",
	 idLr,rLr,nModPerStaveLr*seg->Dz()/2,nStaveLr,nModPerStaveLr,ovlA*100,int(ovlA*xActProj*1e4));
  //
  idLr = 4;
  rLr = 22.0;
  dzLr = 2*41.8;
  seg = seg1;
  nModPerStaveLr = 1+dzLr/seg->Dz();
  ovlA = -1;
  xActProj = seg->DxActive()*TMath::Cos(kTilt*TMath::DegToRad()); // effective r-phi coverage by single stave
  nStaveLr = 1 + rLr*TMath::Pi()*2/xActProj;
  do { ovlA = 1.-rLr*TMath::Pi()*2/nStaveLr/xActProj; } while ( kMinOvl>=0 && ovlA<kMinOvl && nStaveLr++ );		
  ITS->DefineLayerTurbo(idLr, kPhi0, rLr, nModPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, seg->Dx(), kTilt, kLrTick08, seg->Dy(), seg->GetDetTypeID());
  printf("Add Lr%d: R=%.1f DZ:%.1f Staves:%3d NMod/Stave:%3d -> Active Overlap:%.1f (%d micron)\n",
	 idLr,rLr,nModPerStaveLr*seg->Dz()/2,nStaveLr,nModPerStaveLr,ovlA*100,int(ovlA*xActProj*1e4));
  //
  // 
  // OUTER LAYERS (double side readout sensors)
  idLr = 5;
  rLr = 40.0;
  dzLr = 2*71.2;
  seg = seg1;
  nModPerStaveLr = 1+dzLr/seg->Dz();
  ovlA = -1;
  xActProj = seg->DxActive()*TMath::Cos(kTilt*TMath::DegToRad()); // effective r-phi coverage by single stave
  nStaveLr = 1 + rLr*TMath::Pi()*2/xActProj;
  do { ovlA = 1.-rLr*TMath::Pi()*2/nStaveLr/xActProj; } while ( kMinOvl>=0 && ovlA<kMinOvl && nStaveLr++ );		
  ITS->DefineLayerTurbo(idLr, kPhi0, rLr, nModPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, seg->Dx(), kTilt, kLrTick08, seg->Dy(), seg->GetDetTypeID());
  printf("Add Lr%d: R=%.1f DZ:%.1f Staves:%3d NMod/Stave:%3d -> Active Overlap:%.1f (%d micron)\n",
	 idLr,rLr,nModPerStaveLr*seg->Dz()/2,nStaveLr,nModPerStaveLr,ovlA*100,int(ovlA*xActProj*1e4));
  //
  idLr = 6;
  rLr = 43.0;
  dzLr = 2*74.3;
  seg = seg1;
  nModPerStaveLr = 1+dzLr/seg->Dz();
  ovlA = -1;
  xActProj = seg->DxActive()*TMath::Cos(kTilt*TMath::DegToRad()); // effective r-phi coverage by single stave
  nStaveLr = 1 + rLr*TMath::Pi()*2/xActProj;
  do { ovlA = 1.-rLr*TMath::Pi()*2/nStaveLr/xActProj; } while ( kMinOvl>=0 && ovlA<kMinOvl && nStaveLr++ );		
  ITS->DefineLayerTurbo(idLr, kPhi0, rLr, nModPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, seg->Dx(), kTilt, kLrTick08, seg->Dy(), seg->GetDetTypeID());
  printf("Add Lr%d: R=%.1f DZ:%.1f Staves:%3d NMod/Stave:%3d -> Active Overlap:%.1f (%d micron)\n",
	 idLr,rLr,nModPerStaveLr*seg->Dz()/2,nStaveLr,nModPerStaveLr,ovlA*100,int(ovlA*xActProj*1e4));
  //
  
}
