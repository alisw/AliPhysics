#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#endif

//---------------------------------------
Int_t getNStaves(AliITSMFTSegmentationPix* seg, double tilt, double r0, double minOvl);

void CreateITSU()
{
  //
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");


  // build ITS upgrade detector
  // sensitive area 13x15mm (X,Z) with 20x20 micron pitch, 2mm dead zone on readout side and 50 micron guardring
  const double kSensThick = 18e-4;
  const double kPitchX = 20e-4;
  const double kPitchZ = 20e-4;
  const int    kNRow   = 650; 
  const int    kNCol   = 750;
  const int    kNChips = 2;
  const double kLrTick03 = 120e-4;   // -> effective thickness for ~0.3%X layers
  const double kLrTick08 = 600e-4;   // -> effective thickness for ~0.8%X layers
  //
  const double kReadOutEdge = 0.2;   // width of the readout edge (passive bottom)
  const double kGuardRing   = 50e-4; // width of passive area on left/right/top of the sensor
  // create segmentations:
  AliITSMFTSegmentationPix* seg0 = new AliITSMFTSegmentationPix(0,        // segID (0:9)
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
							    );    // see AliITSMFTSegmentationPix.h for extra options
  seg0->Store(AliITSUGeomTGeo::GetITSsegmentationFileName());
  //
  AliITSMFTSegmentationPix* seg1 = new AliITSMFTSegmentationPix(1,        // segID (0:9)
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
							    );    // see AliITSMFTSegmentationPix.h for extra options
  seg1->Store(AliITSUGeomTGeo::GetITSsegmentationFileName());
  //
  seg0->Print();
  seg1->Print();
  //
  const double kMinOvl = 0.005; // require active zones overlap
  const double kPhi0 = 0.;  // az.angle of 1st stave
  const double kTilt = -10.; // tilt in degrees
  double dzLr,rLr,ovlA,xActProj;
  AliITSMFTSegmentationPix* seg=0;
  int nStaveLr,nModPerStaveLr,idLr;
  //      virtual void   DefineLayerTurbo(const Int_t nlay, const Double_t r,  const Double_t zlen, const Int_t nladd,   const Int_t nmod, const Double_t width,
  //				  const Double_t tilt,   const Double_t lthick = 0.,    const Double_t dthick = 0.,   const UInt_t detType=0);
  AliITSUv0 *ITS  = new AliITSUv0("ITS Upgrade",7);
  ITS->SetStaveModel(AliITSUv0::kModel22);
  //
  const int kNWrapVol = 3;
  const double wrpRMin[kNWrapVol]  = { 2.05, 15.0, 32.0};
  const double wrpRMax[kNWrapVol]  = { 8.0, 27.0, 45.0};
  const double wrpZSpan[kNWrapVol] = {28.0, 86.0, 152.0};
  ITS->SetNWrapVolumes(kNWrapVol); // define wrapper volumes for layers
  for (int iw=0;iw<kNWrapVol;iw++) ITS->DefineWrapVolume(iw,wrpRMin[iw],wrpRMax[iw],wrpZSpan[iw]);
  // 
  // INNER LAYERS
  idLr = 0;
  rLr = 2.2; // 2.165?
  dzLr = 2*11.2;   // min Z to cover
  seg = seg0;
  nModPerStaveLr = 1+dzLr/seg->Dz();
  nStaveLr = getNStaves(seg,kTilt,rLr,kMinOvl);
  ITS->DefineLayerTurbo(idLr, kPhi0, rLr, nModPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, seg->Dx(), kTilt, kLrTick03, seg->Dy(), seg->GetChipTypeID());
  printf("Add Lr%d: R=%.1f DZ:%.1f Staves:%3d NMod/Stave:%3d\n",idLr,rLr,nModPerStaveLr*seg->Dz()/2,nStaveLr,nModPerStaveLr);
  //
  idLr = 1;
  rLr = 2.8; // 2.77 ?
  dzLr = 2*12.1;
  seg = seg0;
  nModPerStaveLr = 1+dzLr/seg->Dz();
  nStaveLr = getNStaves(seg,kTilt,rLr,kMinOvl);
  ITS->DefineLayerTurbo(idLr, kPhi0, rLr, nModPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, seg->Dx(), kTilt, kLrTick03, seg->Dy(), seg->GetChipTypeID());
  printf("Add Lr%d: R=%.1f DZ:%.1f Staves:%3d NMod/Stave:%3d\n",idLr,rLr,nModPerStaveLr*seg->Dz()/2,nStaveLr,nModPerStaveLr);
  //
  idLr = 2;
  rLr = 3.6; // 3.58 ?
  dzLr = 2*13.4;
  seg = seg0;
  nModPerStaveLr = 1+dzLr/seg->Dz();
  nStaveLr = getNStaves(seg,kTilt,rLr,kMinOvl);
  ITS->DefineLayerTurbo(idLr, kPhi0, rLr, nModPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, seg->Dx(), kTilt, kLrTick03, seg->Dy(), seg->GetChipTypeID());
  printf("Add Lr%d: R=%.1f DZ:%.1f Staves:%3d NMod/Stave:%3d\n",idLr,rLr,nModPerStaveLr*seg->Dz()/2,nStaveLr,nModPerStaveLr);
  //
  // 
  // MIDDLE LAYERS (double side readout sensors)
  idLr = 3;
  rLr = 20.0;
  dzLr = 2*39.0;
  seg = seg1;
  nModPerStaveLr = 1+dzLr/seg->Dz();
  nStaveLr = getNStaves(seg,kTilt,rLr,kMinOvl);
  ITS->DefineLayerTurbo(idLr, kPhi0, rLr, nModPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, seg->Dx(), kTilt, kLrTick08, seg->Dy(), seg->GetChipTypeID());
  printf("Add Lr%d: R=%.1f DZ:%.1f Staves:%3d NMod/Stave:%3d\n",idLr,rLr,nModPerStaveLr*seg->Dz()/2,nStaveLr,nModPerStaveLr);
  //
  idLr = 4;
  rLr = 22.0;
  dzLr = 2*41.8;
  seg = seg1;
  nModPerStaveLr = 1+dzLr/seg->Dz();
  nStaveLr = getNStaves(seg,kTilt,rLr,kMinOvl);
  ITS->DefineLayerTurbo(idLr, kPhi0, rLr, nModPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, seg->Dx(), kTilt, kLrTick08, seg->Dy(), seg->GetChipTypeID());
  printf("Add Lr%d: R=%.1f DZ:%.1f Staves:%3d NMod/Stave:%3d\n",idLr,rLr,nModPerStaveLr*seg->Dz()/2,nStaveLr,nModPerStaveLr);
  //
  // 
  // OUTER LAYERS (double side readout sensors)
  idLr = 5;
  rLr = 40.0;
  dzLr = 2*71.2;
  seg = seg1;
  nModPerStaveLr = 1+dzLr/seg->Dz();
  nStaveLr = getNStaves(seg,kTilt,rLr,kMinOvl);
  ITS->DefineLayerTurbo(idLr, kPhi0, rLr, nModPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, seg->Dx(), kTilt, kLrTick08, seg->Dy(), seg->GetChipTypeID());
  printf("Add Lr%d: R=%.1f DZ:%.1f Staves:%3d NMod/Stave:%3d\n",idLr,rLr,nModPerStaveLr*seg->Dz()/2,nStaveLr,nModPerStaveLr);
  //
  idLr = 6;
  rLr = 43.0;
  dzLr = 2*74.3;
  seg = seg1;
  nModPerStaveLr = 1+dzLr/seg->Dz();
  nStaveLr = getNStaves(seg,kTilt,rLr,kMinOvl);
  ITS->DefineLayerTurbo(idLr, kPhi0, rLr, nModPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, seg->Dx(), kTilt, kLrTick08, seg->Dy(), seg->GetChipTypeID());
  printf("Add Lr%d: R=%.1f DZ:%.1f Staves:%3d NMod/Stave:%3d\n",idLr,rLr,nModPerStaveLr*seg->Dz()/2,nStaveLr,nModPerStaveLr);
  //  
}

Int_t getNStaves(AliITSMFTSegmentationPix* seg, double tilt, double r0, double minOvl)
{
  double dphi = (90.-tilt)*TMath::DegToRad();
  double cs = TMath::Cos(dphi);
  double sn = TMath::Sin(dphi);  
  double dx = seg->Dx();
  double tL = -dx/2 + seg->GetGuardBot();
  double tU =  dx/2 - seg->GetGuardTop();
  //
  double xL = r0 + cs*tL;
  double yL =      sn*tL;
  double xU = r0 + cs*tU;
  double yU =      sn*tU;
  double phiL = TMath::ATan2(yL,xL);
  double phiU = TMath::ATan2(yU,xU);
  double dphi = TMath::Abs(phiL-phiU);
  if (dphi>TMath::Pi()) dphi = TMath::Abs(dphi-TMath::Pi()*2);
  double span = dphi*r0;
  //
  double ov = -1;
  int nStaveLr = 1 + r0*TMath::Pi()*2/span;
  do { ov = 1.-r0*TMath::Pi()*2/nStaveLr/span; } while ( minOvl>=0 && ov<minOvl && nStaveLr++ );
  printf("Recomend %2d staves for R=%5.2f, ActiveOvl=%5.2f\% (%6.f micron)\n",nStaveLr,r0,ov*100,ov*span*1e4);
  return nStaveLr;
}
