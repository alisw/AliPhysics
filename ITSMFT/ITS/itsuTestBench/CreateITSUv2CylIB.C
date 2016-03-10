#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TMath.h>
#endif

//---------------------------------------
double radii2Turbo(double rMin,double rMid,double rMax, double sensW)
{
  // compute turbo angle from radii and sensor width
  return TMath::ASin((rMax*rMax-rMin*rMin)/(2*rMid*sensW))*TMath::RadToDeg();
}

double radii2Phi(double rMin,double rMid,double rMax, double sensW)
{
  // compute phi coverage
  return 2*TMath::ACos((rMax+rMin)*
		       (rMid*rMid+rMin*rMax-sensW*sensW/4.)/
		       (4.*rMid*rMax*rMin));
}

Int_t getNStaves(AliITSMFTSegmentationPix* seg, double tilt, double r0, double minOvl);


void CreateITSUv2CylIB()
{
  //
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");
  //
  // build ITS upgrade detector
  // sensitive area 13x15mm (X,Z) with 20x20 micron pitch, 2mm dead zone on readout side and 50 micron guardring
  const double kSensThick = 18e-4;
  const double kPitchX = 20e-4;
  const double kPitchZ = 20e-4;
  const int    kNRowOB = 650; 
  const int    kNColOB = 1500;
  //
  const int    kNRowIB = 100; 
  const int    kNColIB = 3000;
  //
  const double kMinOvlIB = 20e-4;
  const double kTurboIB = 0.3;
  const double kSiThickIB = 50e-4;
  const double kSiThickOB = 50e-4;
  //  const double kSensThick = 120e-4;   // -> sensor Si thickness
  //
  const double kReadOutEdgeOB = 0.2;   // width of the readout edge (passive bottom)
  const double kGuardRingOB   = 50e-4; // width of passive area on left/right/top of the sensor
  const double kReadOutEdgeIB = 0.0;   // width of the readout edge (passive bottom)
  const double kGuardRingIB   = 0.0;   // width of passive area on left/right/top of the sensor
  //
  const int kNLr = 7;
  const int kNLrInner = 3;
  const int kBuildLevel = 0;
  enum {kRmn,kRmd,kRmx,kNModPerStave,kPhi0,kNStave,kNPar};
  // Radii are from last TDR (ALICE-TDR-017.pdf Tab. 1.1, rMid is mean value)
  const double tdr5dat[kNLr][kNPar] = { 
    { -1, 1.8, -1,  5.,   0., -1}, // for each inner layer: rMin,rMid,rMax,NChip/Stave, phi0, nStaves
    { -1, 2.4, -1,  5.,   0., -1},
    { -1, 3.0, -1,  5.,   0., -1},
    {-1,  19.6 ,   -1,  4.,  0.  , 24},  // for others: -, rMid, -, NMod/HStave, phi0, nStaves // 24 was 49
    {-1,  24.55, -1,    4.,  0.  , 30},  // 30 was 61
    {-1,  34.39, -1,    7.,  0.  , 42},  // 42 was 88
    {-1,  39.34, -1,    7.,  0.  , 48}   // 48 was 100
  };
  const int nChipsPerModule = 7; // For OB: how many chips in a row

  // create segmentations:
  AliITSMFTSegmentationPix* seg0 = new AliITSMFTSegmentationPix(0,        // segID (0:9)
							    1,  // chips per module
							    kNColOB,    // ncols (total for module)
							    kNRowOB,  // nrows
							    kPitchX,  // default row pitch in cm
							    kPitchZ,  // default col pitch in cm
							    kSensThick,  // sensor thickness in cm
							    -1,     // no special left col between chips
							    -1,     // no special right col between chips
							    kGuardRingOB, // left
							    kGuardRingOB, // right
							    kGuardRingOB, // top
							    kReadOutEdgeOB  // bottom
							    );    // see AliITSMFTSegmentationPix.h for extra options
  seg0->Store(AliITSUGeomTGeo::GetITSsegmentationFileName());
  //
  seg0->Print();
  //
  // special narrow sensors for cyl. layer
  AliITSMFTSegmentationPix* seg1 = new AliITSMFTSegmentationPix(1,        // segID (0:9)
							    1,  // chips per module
							    kNColIB,    // ncols (total for module)
							    kNRowIB,  // nrows
							    kPitchX,  // default row pitch in cm
							    kPitchZ,  // default col pitch in cm
							    kSensThick,  // sensor thickness in cm
							    -1,     // no special left col between chips
							    -1,     // no special right col between chips
							    kGuardRingIB, // left
							    kGuardRingIB, // right
							    kGuardRingIB, // top
							    kReadOutEdgeIB  // bottom
							    );    // see AliITSMFTSegmentationPix.h for extra options
  seg1->Store(AliITSUGeomTGeo::GetITSsegmentationFileName());
  //
  seg1->Print();
  //
  double dzLr,rLr,phi0,turbo;
  int nStaveLr,nModPerStaveLr,idLr;
  //
  AliITSUv2 *ITS  = new AliITSUv2("ITS Upgrade",kNLr);
  ITS->SetStaveModelIB(AliITSUv2::kIBModelDummy);
  ITS->SetStaveModelOB(AliITSUv2::kOBModel2);
  //
  const int kNWrapVol = 3;
  const double wrpRMin[kNWrapVol]  = { 1.7, 15.0, 32.0};
  const double wrpRMax[kNWrapVol]  = { 7.0, 27.0+2.5, 43.0+1.9};
  const double wrpZSpan[kNWrapVol] = {32.0, 90.1, 152.7};
  //
  ITS->SetNWrapVolumes(kNWrapVol); // define wrapper volumes for layers
  for (int iw=0;iw<kNWrapVol;iw++) ITS->DefineWrapVolume(iw,wrpRMin[iw],wrpRMax[iw],wrpZSpan[iw]);
  //
  for (int idLr=0;idLr<kNLr;idLr++) {
    rLr   = tdr5dat[idLr][kRmd];
    phi0  = tdr5dat[idLr][kPhi0]; 
    //
    nStaveLr = TMath::Nint(tdr5dat[idLr][kNStave]);
    nModPerStaveLr =  TMath::Nint(tdr5dat[idLr][kNModPerStave]);
    int nChipsPerStaveLr = nModPerStaveLr;
    //
    if (idLr>=kNLrInner) {
      nChipsPerStaveLr *= nChipsPerModule;
      ITS->DefineLayer(idLr, phi0, rLr, nChipsPerStaveLr*seg0->Dz(), nStaveLr, nModPerStaveLr, 
		       kSiThickOB, seg0->Dy(), seg0->GetChipTypeID(),kBuildLevel);
    } else {
      if (nStaveLr<0) nStaveLr = getNStaves(seg1,kTurboIB,rLr,kMinOvlIB); // calculate automatically
      ITS->DefineLayerTurbo(idLr, phi0, rLr, nChipsPerStaveLr*seg1->Dz(), nStaveLr, nChipsPerStaveLr, 
			    seg1->Dx(), kTurboIB, kSiThickIB, seg1->Dy(), seg1->GetChipTypeID(),kBuildLevel);
    }

  }

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
  printf("Reccommend %2d staves for R=%5.2f, ActiveOvl=%5.2f\% (%6.f micron)\n",nStaveLr,r0,ov*100,ov*span*1e4);
  return nStaveLr;
}
