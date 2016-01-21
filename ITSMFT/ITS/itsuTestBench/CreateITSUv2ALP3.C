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

void CreateITSUv2ALP3()
{
  //
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");
  //
  // build ITS upgrade detector
  // pALPIDE3 15x30 mm^2  (X,Z) with 26.88 x 29.24 micron pitch
  const double kSensThick = 18e-4;
  const double kPitchZ = 29.24e-4;
  const double kPitchX = 26.88e-4;
  const int    kNRow   = 512; 
  const int    kNCol   = 1024;
  const double kSiThickIB = 50e-4;
  const double kSiThickOB = 50e-4;
  //  const double kSensThick = 120e-4;   // -> sensor Si thickness
  //
  const double kReadOutEdge = 0.12;   // width of the readout edge (passive bottom)
  const double kTopEdge = 37.44e-4;   // dead area on top
  const double kLeftRightEdge   = 29.12e-4; // width of passive area on left/right of the sensor
  //
  const int kNLr = 7;
  const int kNLrInner = 3;
  const int kBuildLevel = 0;
  enum {kRmn,kRmd,kRmx,kNModPerStave,kPhi0,kNStave,kNPar};
  // Radii are from last TDR (ALICE-TDR-017.pdf Tab. 1.1, rMid is mean value)
  const double tdr5dat[kNLr][kNPar] = { 
    {2.24, 2.34, 2.67,  9., 16.42, 12}, // for each inner layer: rMin,rMid,rMax,NChip/Stave, phi0, nStaves
    {3.01, 3.15, 3.46,  9., 12.18, 16},
    {3.78, 3.93, 4.21,  9.,  9.55, 20},
    {-1,  19.6 ,   -1,  4.,  0.  , 24},  // for others: -, rMid, -, NMod/HStave, phi0, nStaves // 24 was 49
    {-1,  24.55, -1,    4.,  0.  , 30},  // 30 was 61
    {-1,  34.39, -1,    7.,  0.  , 42},  // 42 was 88
    {-1,  39.34, -1,    7.,  0.  , 48}   // 48 was 100
  };
  const int nChipsPerModule = 7; // For OB: how many chips in a row

  // create segmentations:
  AliITSMFTSegmentationPix* seg0 = new AliITSMFTSegmentationPix(0,        // segID (0:9)
							    1,  // chips per module
							    kNCol,    // ncols (total for module)
							    kNRow,    // nrows
							    kPitchX,  // default row pitch in cm
							    kPitchZ,  // default col pitch in cm
							    kSensThick,  // sensor thickness in cm
							    -1,     // no special left col between chips
							    -1,     // no special right col between chips
							    kLeftRightEdge, // left
							    kLeftRightEdge, // right
							    kTopEdge, // top
							    kReadOutEdge  // bottom
							    );    // see AliITSMFTSegmentationPix.h for extra options
  seg0->Store(AliITSUGeomTGeo::GetITSsegmentationFileName());
  //
  seg0->Print();
  //
  double dzLr,rLr,phi0,turbo;
  int nStaveLr,nModPerStaveLr,idLr;
  //
  AliITSUv2 *ITS  = new AliITSUv2("ITS Upgrade",kNLr);
  ITS->SetStaveModelIB(AliITSUv2::kIBModel4);
  ITS->SetStaveModelOB(AliITSUv2::kOBModel2);
  //
  const int kNWrapVol = 3;
  const double wrpRMin[kNWrapVol]  = { 2.1, 15.0, 32.0};
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
      turbo = radii2Turbo(tdr5dat[idLr][kRmn],rLr,tdr5dat[idLr][kRmx],seg0->Dx());	
      ITS->DefineLayerTurbo(idLr, phi0, rLr, nChipsPerStaveLr*seg0->Dz(), nStaveLr, nChipsPerStaveLr, 
			    seg0->Dx(), turbo, kSiThickIB, seg0->Dy(), seg0->GetChipTypeID(),kBuildLevel);
    }

  }

}
