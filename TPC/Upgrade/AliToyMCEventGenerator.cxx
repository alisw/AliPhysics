#include <iostream>

#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TH2F.h>
#include <TGeoGlobalMagField.h>
#include <TSpline.h>
#include <TObjString.h>

#include <AliLog.h>
#include <AliTPCROC.h>
#include <AliTrackPointArray.h>
#include <AliTrackerBase.h>
#include <AliCDBManager.h>
#include <AliTPCParam.h>
#include <AliGeomManager.h>
#include <AliTPCcalibDB.h>
#include <AliTPCclusterMI.h>
#include <AliTPCSpaceCharge3D.h>
#include <AliTPCROC.h>

#include "AliToyMCEvent.h"
#include "AliToyMCTrack.h"

#include "AliToyMCEventGenerator.h"

ClassImp(AliToyMCEventGenerator);


AliToyMCEventGenerator::AliToyMCEventGenerator()
  :TObject()
  ,fTPCParam(0x0)
  ,fEvent(0x0)
  ,fSpaceCharge(0x0)
  ,fSpaceChargeFile("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2_eps5_50kHz_precal.root")
  ,fOutputFileName("toyMC.root")
  ,fOutFile(0x0)
  ,fOutTree(0x0)
{
  fTPCParam = AliTPCcalibDB::Instance()->GetParameters();
  fTPCParam->ReadGeoMatrices();
}
//________________________________________________________________
AliToyMCEventGenerator::AliToyMCEventGenerator(const AliToyMCEventGenerator &gen)
  :TObject(gen)
  ,fTPCParam(gen.fTPCParam)
  ,fEvent(0x0)
  ,fSpaceCharge(gen.fSpaceCharge)
  ,fSpaceChargeFile(gen.fSpaceChargeFile)
  ,fOutputFileName(gen.fOutputFileName)
  ,fOutFile(0x0)
  ,fOutTree(0x0)
{
  //
}
//________________________________________________________________
AliToyMCEventGenerator::~AliToyMCEventGenerator() 
{
  delete fSpaceCharge;
}

//________________________________________________________________
Bool_t AliToyMCEventGenerator::DistortTrack(AliToyMCTrack &trackIn, Double_t t0)
{
  //
  //
  //

  if(!fTPCParam) {
    fTPCParam = AliTPCcalibDB::Instance()->GetParameters();
    fTPCParam->ReadGeoMatrices();
  }

  // make it big enough to hold all points
  // store real number of generated points in the unique id
  const Int_t nMaxPoints=3000;
  AliTrackPointArray pointArray0(nMaxPoints);  //undistorted
  AliTrackPointArray pointArray1(nMaxPoints);  //distorted

  //Create space point of undistorted and distorted clusters along the propagated track trajectory
  CreateSpacePoints(trackIn,pointArray0,pointArray1);
  //Convert the space points into clusters in the local frame
  //for undistorted and distorted clusters using the same function
  ConvertTrackPointsToLocalClusters(pointArray0,trackIn,t0,0);
  ConvertTrackPointsToLocalClusters(pointArray1,trackIn,t0,1);
  
  return 1;
}

//________________________________________________________________
void AliToyMCEventGenerator::CreateSpacePoints(AliToyMCTrack &trackIn,
                                               AliTrackPointArray &arrUdist,
                                               AliTrackPointArray &arrDist)
{
  //
  // sample the track from the inner to the outer wall of the TPC
  // a graph is filled in local coordinates for later
  //
  
  const Double_t kMaxSnp = 0.85;
  const Double_t kMaxZ0  = fTPCParam->GetZLength();
  const Double_t kMass   = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  
  const Double_t iFCRadius =  83.5; //radius constants found in AliTPCCorrection.cxx
  const Double_t oFCRadius = 254.5;
  
  AliToyMCTrack track(trackIn);
  //!!! TODO: make this adjustable perhaps
  const Double_t stepSize=0.1;
  Double_t xyz[3]  = {0.,0.,0.};
  Float_t  xyzf[3] = {0.,0.,0.};
  
  //!!! when does the propagation not work, how often does it happen?
  if (!AliTrackerBase::PropagateTrackTo(&track,iFCRadius,kMass,5,kTRUE,kMaxSnp)) {
    AliError(Form("Propagation to IFC: %.2f failed\n",iFCRadius));
    return;
  }
  
  Int_t npoints=0;
  
  for (Double_t radius=iFCRadius; radius<oFCRadius; radius+=stepSize){
    
    //!!! changed from return 0 to continue -> Please check
    if (!AliTrackerBase::PropagateTrackTo(&track,radius,kMass,1,kTRUE,kMaxSnp)) {
      AliError(Form("Propagation to %.2f failed\n",radius));
      continue;
    }
    track.GetXYZ(xyz);
    
    //!!! Why is this smeared
    xyz[0]+=gRandom->Gaus(0,0.000005);
    xyz[1]+=gRandom->Gaus(0,0.000005);
    xyz[2]+=gRandom->Gaus(0,0.000005);
    
    xyzf[0]=Float_t(xyz[0]);
    xyzf[1]=Float_t(xyz[1]);
    xyzf[2]=Float_t(xyz[2]);
    
    if (TMath::Abs(track.GetZ())>kMaxZ0) continue;
    //       if (TMath::Abs(track.GetX())<iFCRadius) continue;
    //       if (TMath::Abs(track.GetX())>oFCRadius) continue;
    
    
    AliTrackPoint pUdist;                               // undistorted space point
    AliTrackPoint pDist;                                // distorted space point
    // Set undistorted point
    SetPoint(xyzf,pUdist);
    arrUdist.AddPoint(npoints, &pUdist);
    Int_t sector=pUdist.GetVolumeID();    
    // abuse volume ID for the sector number
    pUdist.SetVolumeID(sector);
    
    // set distorted point
    Float_t distPoint[3]={xyz[0],xyz[1],xyz[2]};
    fSpaceCharge->DistortPoint(distPoint, sector);
    SetPoint(distPoint, pDist);
    arrDist.AddPoint(npoints, &pDist);
    // abuse volume ID for the sector number
    sector=pDist.GetVolumeID();
    pDist.SetVolumeID(sector);
    
    ++npoints;
  }
  
  arrUdist.SetUniqueID(npoints);
  arrDist.SetUniqueID(npoints);
}

//________________________________________________________________
void AliToyMCEventGenerator::SetPoint(Float_t xyz[3], AliTrackPoint &point)
{
  //
  // make AliTrackPoint out of AliTPCclusterMI
  //
  
  //covariance at the local frame
  //assume 1mm distortion in y and z
  Int_t i[3]={0,0,0};
  const Double_t kSigmaY=0.1;
  const Double_t kSigmaZ=0.1;
  Float_t cov[6]={0,0,0, kSigmaY*kSigmaY,0,kSigmaZ*kSigmaZ};
  
  const Float_t alpha = -TMath::ATan2(xyz[1],xyz[0]);
  const Float_t sin   = TMath::Sin(alpha), cos = TMath::Cos(alpha);
  
  Float_t newcov[6];
  newcov[0] = cov[0]*cos*cos           + 2*cov[1]*sin*cos         + cov[3]*sin*sin;
  newcov[1] = cov[1]*(cos*cos-sin*sin) + (cov[3]-cov[0])*sin*cos;
  newcov[2] = cov[2]*cos               + cov[4]*sin;
  newcov[3] = cov[0]*sin*sin           - 2*cov[1]*sin*cos          +  cov[3]*cos*cos;
  newcov[4] = cov[4]*cos               - cov[2]*sin;
  newcov[5] = cov[5];
  
  // voluem ID to add later ....
  point.SetXYZ(xyz);
  point.SetCov(newcov);
  point.SetVolumeID(fTPCParam->Transform0to1(xyz,i));

  // TODO: Add sampled dE/dx  (use SetCharge)
}

//________________________________________________________________
void AliToyMCEventGenerator::ConvertTrackPointsToLocalClusters(AliTrackPointArray &arrPoints,
                                                               AliToyMCTrack &tr, Double_t t0, Int_t type)
{
  //
  //
  //


  const Int_t npoints=Int_t(arrPoints.GetUniqueID());
  Int_t secOld=-1;

  // create an array for graphs which are used for local interpolation
  // we need a new graph if the sector changed since then the local frame will also change
  TObjArray arrGraphsXY(72);
  arrGraphsXY.SetOwner();
  TObjArray arrGraphsXZ(72);
  arrGraphsXZ.SetOwner();
  AliTrackPoint p;
  //create initial graph
  TGraph *grXY=0x0;
  TGraph *grXZ=0x0;
  //row -> sector mapping
  Int_t rowMap[159];
  for (Int_t irow=0; irow<159; ++irow) rowMap[irow]=-1;

  // 1. Step
  // Make from the list of global space points graphs in the local frame
  // one graph per sector is needed
  for (Int_t ipoint=0; ipoint<npoints; ++ipoint){
    arrPoints.GetPoint(p,ipoint);
    Float_t xyz[3] = {p.GetX(),p.GetY(),p.GetZ()};
    Int_t index[3] = {0,0,0};
    Int_t row = fTPCParam->GetPadRow(xyz,index);
    // rotate space point to local frame
    // the angle is given by the VolumeID which was set in CrateSpacePoints
    const Int_t sec=p.GetVolumeID();
    if (row<0 || (sec<36 && row>62) || row>95 ) continue;
    Double_t angle=((sec%18)*20.+10.)/TMath::RadToDeg();
    AliTrackPoint pRot=p.Rotate(angle);
    const Int_t secrow=row+(sec>35)*63;
    if (rowMap[secrow]==-1) rowMap[secrow]=sec;
    // check if we need a new graph (sector change)
    if (secOld!=sec){
      grXY=new TGraph;
      grXZ=new TGraph;
      arrGraphsXY.AddAt(grXY,sec);
      arrGraphsXZ.AddAt(grXZ,sec);
    }
    
    //add coordinates in local frame for later interpolation
    grXY->SetPoint(grXY->GetN(), pRot.GetX(), pRot.GetY());
    grXZ->SetPoint(grXZ->GetN(), pRot.GetX(), pRot.GetZ());
    secOld=sec;
  }

  // 2. Step
  // create in the center of each row a space point by using the graph to interpolate
  // the the center of the row. This is done in xy and xz
  TSpline3 *splXY=0x0;
  TSpline3 *splXZ=0x0;
  AliTPCclusterMI tempCl;
  secOld=-1;
  for (Int_t irow=0; irow<159; ++irow ){
    const Int_t sec     = rowMap[irow];
    if (sec==-1) continue;
    const Int_t secrow = irow<63?irow:irow-63;
    Double_t localX = fTPCParam->GetPadRowRadii(sec,secrow);
    // get graph for the current row
    if (sec!=secOld){
      delete splXY;
      splXY=0x0;
      delete splXZ;
      splXZ=0x0;
      
      grXY=(TGraph*)arrGraphsXY.At(sec);
      grXZ=(TGraph*)arrGraphsXZ.At(sec);
      if (!grXY) continue;

      splXY=new TSpline3("splXY",grXY);
      splXZ=new TSpline3("splXZ",grXZ);
    }
    secOld=sec;

    // check we are in an active area
    if (splXY->FindX(localX)<1 || splXZ->FindX(localX)<1) continue;
    
    //get interpolated value at the center for the pad row
    //  using splines
    const Double_t localY=splXY->Eval(localX/*,0x0,"S"*/);
    const Double_t localZ=splXZ->Eval(localX/*,0x0,"S"*/);
    Float_t xyz[3]={localX,localY,localZ};

    if (!SetupCluster(tempCl,xyz,sec,t0)) continue;
    
    if (type==0) tr.AddSpacePoint(tempCl);
    else tr.AddDistortedSpacePoint(tempCl);
//     printf("SetupCluster %3d: (%.2f, %.2f, %.2f), %d, %.2f\n",irow,xyz[0],xyz[1],xyz[2],sec,t0);
  }

  delete splXY;
  splXY=0x0;
  delete splXZ;
  splXZ=0x0;
  
}

//________________________________________________________________
Bool_t AliToyMCEventGenerator::SetupCluster(AliTPCclusterMI &tempCl, Float_t xyz[3], Int_t sec, Double_t t0)
{
  //
  //
  //

  const Double_t kSigmaY   = 0.1;
  const Double_t kSigmaZ   = 0.1;
  const Double_t kMaxZ0    = fTPCParam->GetZLength();
  //TODO: Get this from the OCDB at some point?
  const Double_t kDriftVel = fTPCParam->GetDriftV();
  
  tempCl.SetX(xyz[0]);
  tempCl.SetY(xyz[1]);
  tempCl.SetZ(xyz[2]);
  
  tempCl.SetSigmaY2(kSigmaY*kSigmaY);
  tempCl.SetSigmaZ2(kSigmaZ*kSigmaZ);

  // transform from the local coordinates to the coordinates expressed in pad coordinates
  Int_t index[3] = {0,sec,0};
  fTPCParam->Transform2to3(xyz,index);
  fTPCParam->Transform3to4(xyz,index);
  fTPCParam->Transform4to8(xyz,index);

  const Int_t   row   = index[2];
  const Int_t   nPads = fTPCParam->GetNPads(sec, row);
  // pad is fractional, but it needs to be shifted from the center
  // to the edge of the row
  const Float_t pad   = xyz[1] + nPads/2;

  tempCl.SetRow(row);
  tempCl.SetPad(pad);
  Float_t timeBin=Float_t(t0 + (kMaxZ0-TMath::Abs(tempCl.GetZ()))/kDriftVel);
  tempCl.SetTimeBin(timeBin); // set time as t0  + drift time from dist z
  tempCl.SetDetector(sec);

  //check if we are in the active area
  if (pad<0 || pad>=nPads) return kFALSE;

  return kTRUE;
}

//________________________________________________________________
Bool_t AliToyMCEventGenerator::ConnectOutputFile()
{
  //
  // Create the output file name and tree and connect the event
  //

  fOutFile = new TFile(fOutputFileName.Data(),"recreate");

  if (!fOutFile || !fOutFile->IsOpen()){
    delete fOutFile;
    fOutFile=0x0;
    return kFALSE;
  }
  
  fOutTree = new TTree("toyMCtree","Tree with toyMC simulation");
  fOutTree->Branch("event","AliToyMCEvent",&fEvent);

  return kTRUE;
}

//________________________________________________________________
Bool_t AliToyMCEventGenerator::CloseOutputFile()
{
  //
  // close the output file
  //
  if (!fOutFile) return kFALSE;
  fOutFile->Write();
  fOutFile->Close();
  delete fOutFile;
  fOutFile=0x0;

  return kTRUE;
}

//________________________________________________________________
void AliToyMCEventGenerator::FillTree()
{
  // fill the tree
  if (fOutTree) fOutTree->Fill();
}

//________________________________________________________________
void AliToyMCEventGenerator::SetSpaceCharge(EEpsilon epsilon, EGasType gasType/*=kNeCO2_9010*/, ECollRate collRate/*=k50kHz*/)
{
  //
  // Set the space charge conditions
  //
  fSpaceChargeFile="$ALICE_ROOT/TPC/Calib/maps/SC";
  switch (gasType) {
    case kNeCO2_9010:
      fSpaceChargeFile.Append("_NeCO2");
      break;
  }
  switch (epsilon) {
    case kEps5:
      fSpaceChargeFile.Append("_eps5");
      break;
    case kEps10:
      fSpaceChargeFile.Append("_eps10");
      break;
    case kEps20:
      fSpaceChargeFile.Append("_eps20");
      break;
  }
  switch (collRate) {
    case k50kHz:
      fSpaceChargeFile.Append("_50kHz");
      break;
  }
  fSpaceChargeFile.Append("_precal.root");
}

//________________________________________________________________
void AliToyMCEventGenerator::InitSpaceCharge()
{
  //
  // init the space charge conditions
  // this should be called after the tree was connected
  //

  AliInfo(Form("Using space charge map file: '%s'",fSpaceChargeFile.Data()));
  
  TFile f(fSpaceChargeFile.Data());
  fSpaceCharge=(AliTPCSpaceCharge3D*)f.Get("map");

  if (fOutTree){
    AliInfo("Attaching space charge map file name to the tree");
    fOutTree->GetUserInfo()->Add(new TObjString(fSpaceChargeFile.Data()));
  }
}
