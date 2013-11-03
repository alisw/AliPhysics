#include <iostream>
#include <fstream>

#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TH2F.h>
#include <TGeoGlobalMagField.h>
#include <TSpline.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TLinearFitter.h>

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
#include <AliExternalTrackParam.h>

#include "AliToyMCEvent.h"
#include "AliToyMCTrack.h"

#include "AliToyMCEventGenerator.h"

ClassImp(AliToyMCEventGenerator);


AliToyMCEventGenerator::AliToyMCEventGenerator()
  :TObject()
  ,fTPCParam(0x0)
  ,fEvent(0x0)
  ,fCurrentTrack(0)
  ,fTPCCorrection(0x0)
  ,fTPCCorrectionAv(0x0)
  ,fSCList(0x0)
  ,fSCListFile()
  ,fCorrectionFile("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2_eps5_50kHz_precal.lookup.root")
  ,fOutputFileName("toyMC.root")
  ,fOutFile(0x0)
  ,fOutTree(0x0)
  ,fUseStepCorrection(kFALSE)
  ,fUseMaterialBudget(kFALSE)
  ,fIsLaser(kTRUE)
  ,fPrereadSCList(kFALSE)
  ,fCalculateScaling(kTRUE)
{
  fTPCParam = AliTPCcalibDB::Instance()->GetParameters();
  fTPCParam->ReadGeoMatrices();
  gRandom->SetSeed();
}
//________________________________________________________________
AliToyMCEventGenerator::AliToyMCEventGenerator(const AliToyMCEventGenerator &gen)
  :TObject(gen)
  ,fTPCParam(gen.fTPCParam)
  ,fEvent(0x0)
  ,fCurrentTrack(0)
  ,fTPCCorrection(gen.fTPCCorrection)
  ,fTPCCorrectionAv(0x0)
  ,fSCList(gen.fSCList)
  ,fSCListFile(gen.fSCListFile)
  ,fCorrectionFile(gen.fCorrectionFile)
  ,fOutputFileName(gen.fOutputFileName)
  ,fOutFile(0x0)
  ,fOutTree(0x0)
  ,fUseStepCorrection(gen.fUseStepCorrection)
  ,fUseMaterialBudget(gen.fUseMaterialBudget)
  ,fIsLaser(gen.fIsLaser)
  ,fPrereadSCList(gen.fPrereadSCList)
  ,fCalculateScaling(gen.fCalculateScaling)
{
  //
  gRandom->SetSeed();
}
//________________________________________________________________
AliToyMCEventGenerator::~AliToyMCEventGenerator() 
{
  if (HasSCList() &&!fPrereadSCList) delete fTPCCorrection;
  delete fSCList;
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

  MakeITSClusters(trackIn/*,t0*/);
  MakeTPCClusters(trackIn, t0);
  MakeTRDClusters(trackIn/*,t0*/);

  return kTRUE;
}
//________________________________________________________________
void AliToyMCEventGenerator::MakeITSClusters(AliToyMCTrack &trackIn/*, Double_t t0*/)
{
  //Upgrade ITS parameters
  const Int_t nITSLayers = 7;
  const Double_t ITSRadii[nITSLayers] = {2.2, 2.8, 3.6, 20.0, 22.0, 41.0, 43.0};
  const Double_t lengthITS[nITSLayers] = {22.4, 24.2, 26.8, 78.0, 83.6, 142.4, 148.6};

  //resolution of the point is 10um
  const Float_t sigmaY = 0.001;
  const Float_t sigmaZ = 0.001;

  AliTrackPoint point;
  
  const Double_t kMaxSnp = 0.85;
//   const Double_t kMaxZ0  = fTPCParam->GetZLength();
  const Double_t kMass   = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  
  AliExternalTrackParam track(trackIn);
  Double_t xyz[3]  = {0.,0.,0.};
      
  if (!AliTrackerBase::PropagateTrackTo(&track,ITSRadii[0],kMass,5,kTRUE,kMaxSnp,0,kFALSE,fUseMaterialBudget)) {
    AliError(Form("Propagation to %.2f failed\n",ITSRadii[0]));
    return;
  }

  for(Int_t iLayer = 0; iLayer<nITSLayers; iLayer++){

    if (!AliTrackerBase::PropagateTrackTo(&track,ITSRadii[iLayer],kMass,1,kTRUE,kMaxSnp,0,kFALSE,fUseMaterialBudget)) {
      AliError(Form("Propagation to %.2f failed\n",ITSRadii[iLayer]));
      continue;
    }
    // since I don't know how to extract the volumeid of the ITS layers I use the following strategy:
    //  - rotate the track to the next integer angle
    //  - use coordinates there and set as volumeid the integer angle
    //  - in the reco one can then rotate the cluster to the global frame using the angle stored in the volume id
    track.GetXYZ(xyz);
    
//     if (TMath::Abs(track.GetZ())>kMaxZ0) continue;
    if (TMath::Abs(track.GetZ())>lengthITS[iLayer]/2) continue;


    // smear the ideal positions with the cluster resolution
    xyz[1]+=gRandom->Gaus(0,sigmaY);
    xyz[2]+=gRandom->Gaus(0,sigmaZ);

    Float_t xyzf[3]={xyz[0],xyz[1],xyz[2]};

    SetPoint(xyzf,sigmaY,sigmaZ,point);
    
    trackIn.AddITSPoint(point)->SetUniqueID(trackIn.GetUniqueID());
  }

}
//________________________________________________________________
void AliToyMCEventGenerator::MakeTRDClusters(AliToyMCTrack &trackIn/*, Double_t t0*/)

{
  //Uses current TRD parameters
  const Int_t nTRDLayers = 6;
  const Double_t distToMid = 3.2 + 30./2; //dist to middle of drift region (radiator + half drift region)
  const Double_t TRDRadii[nTRDLayers] = {294.5 + distToMid, 307.1 + distToMid, 319.7 + distToMid, 332.3 + distToMid, 344.9 + distToMid, 357.5 + distToMid};
  const Double_t lengthTRD[nTRDLayers] = {604.0, 634.0, 656.0, 686.0, 700.0, 700.0};
  
  const Float_t sigmaY = 0.06;
  const Float_t sigmaZ = 0.2;

  AliTrackPoint point;
  
  const Double_t kMaxSnp = 0.85;
  const Double_t kMaxZ0  = fTPCParam->GetZLength();
  const Double_t kMass   = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  
  AliExternalTrackParam track(trackIn);
  Double_t xyz[3]  = {0.,0.,0.};
      
  if (!AliTrackerBase::PropagateTrackTo(&track,TRDRadii[0],kMass,5,kTRUE,kMaxSnp,0,kFALSE,fUseMaterialBudget)) {
    AliError(Form("Propagation to %.2f failed\n",TRDRadii[0]));
    return;
  }

  for(Int_t iLayer = 0; iLayer<nTRDLayers; iLayer++){

    if (!AliTrackerBase::PropagateTrackTo(&track,TRDRadii[iLayer],kMass,1,kTRUE,kMaxSnp,0,kFALSE,fUseMaterialBudget)) {
      AliError(Form("Propagation to %.2f failed\n",TRDRadii[iLayer]));
      continue;
    }
    track.GetXYZ(xyz);
    
    if (TMath::Abs(track.GetZ())>kMaxZ0) continue;
    if (TMath::Abs(track.GetZ())>lengthTRD[iLayer]/2) continue;
    

    // smear the ideal positions with the cluster resolution
    xyz[1]+=gRandom->Gaus(0,sigmaY);
    xyz[2]+=gRandom->Gaus(0,sigmaZ);
    
    Float_t xyzf[3]={xyz[0],xyz[1],xyz[2]};

    SetPoint(xyzf,sigmaY,sigmaZ,point);

    trackIn.AddTRDPoint(point)->SetUniqueID(trackIn.GetUniqueID());
  }

}
//________________________________________________________________
void AliToyMCEventGenerator::MakeTPCClusters(AliToyMCTrack &trackIn, Double_t t0)
{

  // make it big enough to hold all points
  // store real number of generated points in the unique id
  const Int_t nMaxPoints=3000;
  static AliTrackPointArray pointArray0(nMaxPoints);  //undistorted
  static AliTrackPointArray pointArray1(nMaxPoints);  //distorted
  
  //Create space point of undistorted and distorted clusters along the propagated track trajectory
  CreateSpacePoints(trackIn,pointArray0,pointArray1);
  //Convert the space points into clusters in the local frame
  //for undistorted and distorted clusters using the same function
  ConvertTrackPointsToLocalClusters(pointArray0,trackIn,t0,0);
  ConvertTrackPointsToLocalClusters(pointArray1,trackIn,t0,1);
  

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
  
  Double_t kMaxSnp = 0.85;
  if (fIsLaser) kMaxSnp=0.99;
  const Double_t kMaxZ0  = fTPCParam->GetZLength();
  const Double_t kMass   = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  
  const Double_t iFCRadius =  83.5; //radius constants found in AliTPCCorrection.cxx
  const Double_t oFCRadius = 254.5;

  const Float_t kSigmaY=0.1;
  const Float_t kSigmaZ=0.1;
  
  AliExternalTrackParam track(trackIn);
  //!!! TODO: make this adjustable perhaps
  const Double_t stepSize=0.1;
  Double_t xyz[3]  = {0.,0.,0.};
  Float_t  xyzf[3] = {0.,0.,0.};
  
  //!!! when does the propagation not work, how often does it happen?
  if (!AliTrackerBase::PropagateTrackTo(&track,iFCRadius,kMass,5,kTRUE,kMaxSnp,0,kFALSE,fUseMaterialBudget) && !fIsLaser) {
    AliError(Form("Propagation to IFC: %.2f failed\n",iFCRadius));
    return;
  }
  
  Int_t npoints=0;
  
  for (Double_t radius=iFCRadius; radius<oFCRadius; radius+=stepSize){
    //!!! changed from return 0 to continue -> Please check
    if (!AliTrackerBase::PropagateTrackTo(&track,radius,kMass,1,kTRUE,kMaxSnp,0,kFALSE,fUseMaterialBudget)) {
      AliError(Form("Propagation to r=%.2f (snp=%.2f) failed\n",radius,track.GetSnp()));
      continue;
    }
    track.GetXYZ(xyz);
    
    //!!! Why is this smeared
    
    xyzf[0]=Float_t(xyz[0]);
    xyzf[1]=Float_t(xyz[1]);
    xyzf[2]=Float_t(xyz[2]);
    
    if (TMath::Abs(track.GetZ())>kMaxZ0) continue;
    //       if (TMath::Abs(track.GetX())<iFCRadius) continue;
    //       if (TMath::Abs(track.GetX())>oFCRadius) continue;
    
    
    AliTrackPoint pUdist;                               // undistorted space point
    AliTrackPoint pDist;                                // distorted space point
    // Set undistorted point
    SetPoint(xyzf,kSigmaY,kSigmaZ,pUdist);
    arrUdist.AddPoint(npoints, &pUdist);
    Int_t sector=pUdist.GetVolumeID();    
    
    // set distorted point
    Float_t distPoint[3]={xyz[0],xyz[1],xyz[2]};
    Float_t dxyz[3]={0.,0.,0.};
    if (!fUseStepCorrection){
      fTPCCorrection->DistortPoint(distPoint, sector);
    } else {
      fTPCCorrection->GetCorrectionIntegralDz(distPoint,sector,dxyz,5);
      distPoint[0]-=dxyz[0];
      distPoint[1]-=dxyz[1];
      distPoint[2]-=dxyz[2];
    }
    SetPoint(distPoint,kSigmaY,kSigmaZ,pDist);
    arrDist.AddPoint(npoints, &pDist);
    
    ++npoints;
  }
  
  arrUdist.SetUniqueID(npoints);
  arrDist.SetUniqueID(npoints);
}

//________________________________________________________________
void AliToyMCEventGenerator::SetPoint(Float_t xyz[3], Float_t sigmaY, Float_t sigmaZ, AliTrackPoint &point)
{
  //
  // make AliTrackPoint out of AliTPCclusterMI
  //
  
  //covariance at the local frame
  //assume 1mm distortion in y and z
  Int_t i[3]={0,0,0};
  Float_t cov[6]={0,0,0, sigmaY*sigmaY,0,sigmaZ*sigmaZ};
  
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
  // abuse volume ID for the sector number
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

//       if(grXY->GetN()>1 && grXZ->GetN()>1) { //causes segmentation violation if N==1
//        	splXY=new TSpline3("splXY",grXY);
// 	splXZ=new TSpline3("splXZ",grXZ);
//       }
//       else {
	//TODO: make a cluster also in the sector w only one space point?
// 	continue;
	// Double_t tempX=0., tempY = 0., tempZ = 0.;
	
	// grXY->GetPoint(0,tempX,localY);
	// grXZ->GetPoint(0,tempX,localZ);
//       }

    }
    secOld=sec;

    // check we are in an active area
//     if (splXY->FindX(localX)<1 || splXZ->FindX(localX)<1) continue;
    if ( localX<grXY->GetX()[0] || localX>grXY->GetX()[grXY->GetN()-1] || localX<grXZ->GetX()[0] || localX>grXZ->GetX()[grXZ->GetN()-1]) continue;
    
    //get interpolated value at the center for the pad row
    //  using splines
//     const Double_t localY=splXY->Eval(localX/*,0x0,"S"*/);
//     const Double_t localZ=splXZ->Eval(localX/*,0x0,"S"*/);
    const Double_t localY=grXY->Eval(localX/*,0x0,"S"*/);
    const Double_t localZ=grXZ->Eval(localX/*,0x0,"S"*/);
    Float_t xyz[3]={localX,localY,localZ};

    if (!SetupCluster(tempCl,xyz,sec,t0)) continue;
    tempCl.SetLabel(tr.GetUniqueID(), 0);

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

  // intrinsic cluster resolution is 1mm
  const Double_t kSigmaY   = 0.1;
  const Double_t kSigmaZ   = 0.1;
  const Double_t kMaxZ0    = fTPCParam->GetZLength();
  //TODO: Get this from the OCDB at some point?
  const Double_t kDriftVel = fTPCParam->GetDriftV();

  // smear the ideal positions with the cluster resolution
  xyz[1]+=gRandom->Gaus(0,kSigmaY);
  xyz[2]+=gRandom->Gaus(0,kSigmaZ);
  
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

  gROOT->cd();

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
  if (fOutTree&&fEvent) fOutTree->Fill();
}

//________________________________________________________________
void AliToyMCEventGenerator::SetSpaceCharge(EEpsilon epsilon, EGasType gasType/*=kNeCO2_9010*/,
                                            ECollRate collRate/*=k50kHz*/, ECorrection corrType/*=kLookup*/)
{
  //
  // Set the space charge conditions
  //
  fCorrectionFile="$ALICE_ROOT/TPC/Calib/maps/SC";
  switch (gasType) {
    case kNeCO2_9010:
      fCorrectionFile.Append("_NeCO2");
      break;
    case kNeCO2N2_90105:
      fCorrectionFile.Append("_NeCO2N2");
      break;
  }
  switch (epsilon) {
    case kEps5:
      fCorrectionFile.Append("_eps5");
      break;
    case kEps10:
      fCorrectionFile.Append("_eps10");
      break;
    case kEps20:
      fCorrectionFile.Append("_eps20");
      break;
    case kEps25:
      fCorrectionFile.Append("_eps25");
      break;
    case kEps30:
      fCorrectionFile.Append("_eps30");
      break;
    case kEps35:
      fCorrectionFile.Append("_eps35");
      break;
    case kEps40:
      fCorrectionFile.Append("_eps40");
      break;
  }
  switch (collRate) {
    case k50kHz:
      fCorrectionFile.Append("_50kHz");
      break;
  }
  switch (corrType) {
    case kLookup:
      fCorrectionFile.Append("_precal.lookup.root");
      break;
    case kSpaceChargeFile:
      fCorrectionFile.Append("_precal.root");
      break;
  }
}

//________________________________________________________________
void AliToyMCEventGenerator::InitSpaceCharge()
{
  //
  // init the space charge conditions
  // this should be called after the tree was connected
  //

  if (HasSCList()) {
    InitSpaceChargeList();
    return;
  }
  
  AliInfo(Form("Using space charge map file: '%s'",fCorrectionFile.Data()));

  SetCorrectionFromFile(fCorrectionFile, fTPCCorrection);

  if (fOutTree){
    AliInfo("Attaching space charge map file name to the tree");
    fOutTree->GetUserInfo()->Add(new TObjString(fCorrectionFile.Data()));
  }
}

//________________________________________________________________
void AliToyMCEventGenerator::InitSpaceChargeList()
{
  //
  // init space charge conditions from a list of files
  // this should be called after the tree was connected
  //

  std::ifstream file(fSCListFile.Data());
  TString list;
  list.ReadFile(file);
  
  TObjArray *arr=list.Tokenize("\n");
  if (!arr->GetEntriesFast()) {
    delete arr;
    AliFatal(Form("No SC file initialised. SC list '%s' seems empty",fSCListFile.Data()));
    return;
  }

  // it is assumed that in case of an input list
  // fCorrectionFile contains the name of the average correction
  // to be then used in the reconstruction
  if (fOutTree){
    if (fCorrectionFile.IsNull()) {
      AliFatal("List of SC files set, but no average map is specified. Use 'SetSpaceChargeFile' to do so");
      return;
    }
    AliInfo("Attaching average space charge map file name to the tree");
    fOutTree->GetUserInfo()->Add(new TObjString(fCorrectionFile.Data()));
  }
  
  // check for an average map
  if (!fCorrectionFile.IsNull()) {
    SetCorrectionFromFile(fCorrectionFile, fTPCCorrectionAv);
  }
  
  // case of non preread
  // store the names of the files
  if (!fPrereadSCList) {
    if (fSCList) delete fSCList;
    fSCList=arr;
    return;
  }

  // case of preread
  // load all SC files and set them to the list
  if (!fSCList) fSCList=new TObjArray;
  fSCList->SetOwner();
  fSCList->Delete();
  
  for (Int_t ifile=0; ifile<arr->GetEntriesFast(); ++ifile) {
    TString scname=arr->At(ifile)->GetName();
    AliTPCCorrection *sc=0x0;
    SetCorrectionFromFile(scname, sc);
    if (!sc) continue;
    fSCList->Add(sc);
  }
}

//________________________________________________________________
void AliToyMCEventGenerator::IterateSC(Int_t ipos)
{
  //
  // In case a list of SC files is set iterate over them
  //

  if (!HasSCList()) return;

  if (ipos<0) {
    if (fOutTree) ipos=fOutTree->GetEntries();
    else ipos=0;
  }

  TObject *sc=fSCList->At(ipos%fSCList->GetEntriesFast());
  AliInfo(Form("Event: %d - SC: %s",ipos,sc->GetName()));
  // case SC files have been preread
  if (fPrereadSCList) {
    fTPCCorrection=(AliTPCCorrection*)sc;
    return;
  }

  // case no preread was done
  TString &file=((TObjString*)sc)->String();
  delete fTPCCorrection;
  fTPCCorrection=0x0;

  SetCorrectionFromFile(file, fTPCCorrection);
  
  if (!fTPCCorrection) {
    AliFatal(Form("Could not read SC map from SC file '%s'",file.Data()));
    return;
  }

}

//________________________________________________________________
Float_t AliToyMCEventGenerator::GetSCScalingFactor(AliTPCCorrection *corr, AliTPCCorrection *averageCorr, Float_t &chi2)
{
  //
  // calculate the scaling factor
  // between the fluctuation map and the average map
  //

  TLinearFitter fitter(2,"pol1");
  
  for (Float_t iz=-245; iz<=245; iz+=10) {
    Short_t roc=(iz>=0)?0:18;
    for (Float_t ir=86; ir<250; ir+=10) {
      for (Float_t iphi=0; iphi<TMath::TwoPi(); iphi+=10*TMath::DegToRad()){
        Float_t x=ir*(Float_t)TMath::Cos(iphi);
        Float_t y=ir*(Float_t)TMath::Sin(iphi);
        Float_t x3[3]    = {x,y,iz};
        Float_t dx3[3]   = {0.,0.,0.};
        Float_t dx3av[3] = {0.,0.,0.};
        
        corr->GetDistortion(x3,roc,dx3);
        averageCorr->GetDistortion(x3,roc,dx3av);

        Double_t dr   = TMath::Sqrt(dx3[0]*dx3[0]     +  dx3[1]*dx3[1]);
        Double_t drav = TMath::Sqrt(dx3av[0]*dx3av[0] +  dx3av[1]*dx3av[1]);

        fitter.AddPoint(&drav,dr);
      }
    }
  }
  
  fitter.Eval();
  chi2 = fitter.GetChisquare()/fitter.GetNpoints();
  return fitter.GetParameter(1);
}

//________________________________________________________________
void AliToyMCEventGenerator::SetSCScalingFactor()
{
  //
  // if running with a SC list calculate the scaling factor
  // between the fluctuation map and the average map
  //

  if ( !(fCalculateScaling && HasSCList() && fTPCCorrection && fTPCCorrectionAv && fEvent) ) return;

  // loop over several z, r and phi bins and find a factor that minimises
  // the distortions between the current map and the average map

  Float_t chi2   = 0.;
  Float_t factor = GetSCScalingFactor(fTPCCorrection, fTPCCorrectionAv, chi2);

  fEvent->SetSCscale(factor);
  fEvent->SetSCscaleChi2(chi2);
}

//________________________________________________________________
void AliToyMCEventGenerator::SetCorrectionFromFile(TString file, AliTPCCorrection* &corr)
{
  //
  // read the correction from file and set it to corr
  //

  corr=0x0;
  TString corrName("map");
  
  // allow for specifying an object name for the AliTPCCorrection in the file name
  // separated by a ':'
  TObjArray *arr=file.Tokenize(":");
  if (arr->GetEntriesFast()>1) {
    file=arr->At(0)->GetName();
    corrName=arr->At(1)->GetName();
  }
  delete arr;
  
  
  TFile f(file.Data());
  if (!f.IsOpen()||f.IsZombie()) {
    printf("E-AliToyMCEventGenerator::Could not open SC file '%s'",file.Data());
    return;
  }
  
  corr=(AliTPCSpaceCharge3D*)f.Get(corrName.Data());
  if (corr) corr->SetName(f.GetName());
}


