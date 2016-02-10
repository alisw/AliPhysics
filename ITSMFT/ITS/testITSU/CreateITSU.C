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
  const double kLrThick03 = 120e-4;   // -> effective thickness for ~0.3%X layers
  const double kLrThick08 = 600e-4;   // -> effective thickness for ~0.8%X layers
  //
  const double kReadOutEdge = 0.2;   // width of the readout edge (passive bottom)
  const double kGuardRing   = 50e-4; // width of passive area on left/right/top of the sensor
  //
  const int kNLr = 7;
  const int kNLrInner = 3;
  enum {kRmn,kRmd,kRmx,kNModPerStave,kPhi0,kNStave,kNPar};
  const double tdr5dat[kNLr][kNPar] = { 
    {2.24, 2.34, 2.67,  9., 16.37, 12}, // for each inner layer: rMin,rMid,rMax,NMod/Stave,phi0, nStave
    {3.01, 3.15, 3.46,  9., 12.03, 16},
    {3.78, 3.93, 4.21,  9., 10.02, 20},
    {-1,   19.4, -1,    28., 0.  , 49},  // for others: -, rMid, -, NMod/Stave, phi0, nStave
    {-1,   24.7, -1,    28., 0.  , 61},
    {-1,   35.32,-1,    49., 0.  , 88},
    {-1,   40.52,-1,    49., 0.  , 100}
  };

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
  const double kTilt = -10.; // tilt in degrees for outer layers
  double dzLr,rLr,phi0,turbo,thick;
  AliITSMFTSegmentationPix* seg=0;
  int nStaveLr,nModPerStaveLr,idLr;
  //      virtual void   DefineLayerTurbo(const Int_t nlay, const Double_t r,  const Double_t zlen, const Int_t nladd,   const Int_t nmod, const Double_t width,
  //				  const Double_t tilt,   const Double_t lthick = 0.,    const Double_t dthick = 0.,   const UInt_t detType=0);
  AliITSUv0 *ITS  = new AliITSUv0("ITS Upgrade",7);
  ITS->SetStaveModel(AliITSUv0::kModel22);
  //
  const int kNWrapVol = 3;
  const double wrpRMin[kNWrapVol]  = { 2.1, 15.0, 32.0};
  const double wrpRMax[kNWrapVol]  = { 7.0, 27.0, 43.0};
  const double wrpZSpan[kNWrapVol] = {28.0, 86.0, 150.0};
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
    ITS->DefineLayerTurbo(idLr, phi0, rLr, nModPerStaveLr*seg->Dz(), nStaveLr, nModPerStaveLr, 
			  seg->Dx(), turbo, thick, seg->Dy(), seg->GetChipTypeID());
    printf("Add Lr%d: R=%6.2f DZ:%6.2f Turbo:%+6.2f Staves:%3d NMod/Stave:%3d\n",idLr,rLr,nModPerStaveLr*seg->Dz(),turbo,nStaveLr,nModPerStaveLr);
    //
  }
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
  printf("Reccommend %2d staves for R=%5.2f, ActiveOvl=%5.2f\% (%6.f micron)\n",nStaveLr,r0,ov*100,ov*span*1e4);
  return nStaveLr;
}
