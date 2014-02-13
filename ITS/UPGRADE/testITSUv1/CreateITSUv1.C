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

void CreateITSUv1()
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
  const double kLrThick03 = 120e-4;   // -> effective thickness for ~0.3%X layers
  const double kLrThick08 = 600e-4;   // -> effective thickness for ~0.8%X layers
  //
  const double kReadOutEdge = 0.2;   // width of the readout edge (passive bottom)
  const double kGuardRing   = 50e-4; // width of passive area on left/right/top of the sensor
  //
  const int kNLr = 7;
  const int kNLrInner = 3;
  const int kBuildLevel = 3;
  enum {kRmn,kRmd,kRmx,kNModPerStave,kPhi0,kNStave,kNPar};
  // Radii are from last TDR (ALICE-TDR-017.pdf Tab. 1.1, rMid is mean value)
  const double tdr5dat[kNLr][kNPar] = { 
    {2.24, 2.34, 2.67,  9., 16.37, 12}, // for each inner layer: rMin,rMid,rMax,NMod/Stave,phi0, nStave
    {3.01, 3.15, 3.46,  9., 12.03, 16},
    {3.78, 3.93, 4.21,  9., 10.02, 20},
    {-1,  19.6 ,   -1,  4.,  0.  , 24},  // for others: -, rMid, -, NMod/Stave, phi0, nStave // 24 was 49
    {-1,  24.55, -1,    4.,  0.  , 30},  // 30 was 61
    {-1,  34.39, -1,    7.,  0.  , 42},  // 42 was 88
    {-1,  39.34, -1,    7.,  0.  , 48}   // 48 was 100
  };
  const int nChipsPerModule = 7; // For OB we have to specify how many chips

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
  const double kTilt = -10.; // tilt in degrees for outer layers
  double dzLr,rLr,phi0,turbo,thick;
  AliITSUSegmentationPix* seg=0;
  int nStaveLr,nModPerStaveLr,idLr;
  //      virtual void   DefineLayerTurbo(const Int_t nlay, const Double_t r,  const Double_t zlen, const Int_t nladd,   const Int_t nmod, const Double_t width,
  //				  const Double_t tilt,   const Double_t lthick = 0.,    const Double_t dthick = 0.,   const UInt_t detType=0);
  AliITSUv1 *ITS  = new AliITSUv1("ITS Upgrade",7);
  ITS->SetStaveModelIB(AliITSUv1::kIBModel22);
  ITS->SetStaveModelOB(AliITSUv1::kOBModel1);
  //
  const int kNWrapVol = 3;
  const double wrpRMin[kNWrapVol]  = { 2.1, 15.0, 32.0};
  const double wrpRMax[kNWrapVol]  = { 7.0, 27.0, 43.0};
  const double wrpZSpan[kNWrapVol] = {28.0, 96.0, 158.0};
  ITS->SetNWrapVolumes(kNWrapVol); // define wrapper volumes for layers
  for (int iw=0;iw<kNWrapVol;iw++) ITS->DefineWrapVolume(iw,wrpRMin[iw],wrpRMax[iw],wrpZSpan[iw]);
  //
  for (int idLr=0;idLr<kNLr;idLr++) {
    rLr   = tdr5dat[idLr][kRmd];
    phi0  = tdr5dat[idLr][kPhi0]; 
    if (idLr<kNLrInner) {
      seg = seg0;
      turbo = -radii2Turbo(tdr5dat[idLr][kRmn],rLr,tdr5dat[idLr][kRmx],seg->Dx());	
      thick = kLrThick03;
    }
    else {
      seg   = seg1;
      turbo = kTilt;
      thick = kLrThick08;
    }
    nStaveLr = TMath::Nint(tdr5dat[idLr][kNStave]);
    if (nStaveLr<0) nStaveLr = getNStaves(seg,kTilt,rLr,kMinOvl); // calculate automatically
    nModPerStaveLr =  TMath::Nint(tdr5dat[idLr][kNModPerStave]);
    int nChipsPerStaveLr = nModPerStaveLr;
    if (idLr>=kNLrInner) {
      nChipsPerStaveLr *= nChipsPerModule;
      ITS->DefineLayer(idLr, phi0, rLr, nChipsPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, 
		       thick, seg->Dy(), seg->GetChipTypeID(),kBuildLevel);
      printf("Add Lr%d: R=%6.2f DZ:%6.2f Staves:%3d NMod/Stave:%3d\n",idLr,rLr,nChipsPerStaveLr*seg->Dz(),nStaveLr,nModPerStaveLr);
    } else {
      ITS->DefineLayerTurbo(idLr, phi0, rLr, nChipsPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, 
			  seg->Dx(), turbo, thick, seg->Dy(), seg->GetChipTypeID());
      printf("Add Lr%d: R=%6.2f DZ:%6.2f Turbo:%+6.2f Staves:%3d NMod/Stave:%3d\n",idLr,rLr,nChipsPerStaveLr*seg->Dz(),turbo,nStaveLr,nModPerStaveLr);
    }
    //
  }
  //  
}

Int_t getNStaves(AliITSUSegmentationPix* seg, double tilt, double r0, double minOvl)
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
