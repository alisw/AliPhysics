
#include <TDatabasePDG.h>
#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TFile.h>
#include <TPRegexp.h>

#include <AliExternalTrackParam.h>
#include <AliTPCcalibDB.h>
#include <AliTPCclusterMI.h>
#include <AliTPCSpaceCharge3D.h>
#include <AliTrackerBase.h>
#include <AliTrackPointArray.h>
#include <AliLog.h>
#include <AliTPCParam.h>
#include <AliTPCROC.h>
#include <TTreeStream.h>

#include "AliToyMCTrack.h"
#include "AliToyMCEvent.h"

#include "AliToyMCReconstruction.h"


//____________________________________________________________________________________
AliToyMCReconstruction::AliToyMCReconstruction() : TObject()
, fSeedingRow(140)
, fSeedingDist(10)
, fClusterType(0)
, fCorrectionType(kNoCorrection)
, fDoTrackFit(kTRUE)
, fUseMaterial(kFALSE)
, fIdealTracking(kFALSE)
, fTime0(-1)
, fCreateT0seed(kFALSE)
, fStreamer(0x0)
, fTree(0x0)
, fEvent(0x0)
, fTPCParam(0x0)
, fSpaceCharge(0x0)
{
  //
  //  ctor
  //
  fTPCParam=AliTPCcalibDB::Instance()->GetParameters();
}

//____________________________________________________________________________________
AliToyMCReconstruction::~AliToyMCReconstruction()
{
  //
  //  dtor
  //

  delete fStreamer;
//   delete fTree;
}

//____________________________________________________________________________________
void AliToyMCReconstruction::RunReco(const char* file, Int_t nmaxEv)
{
  //
  //
  //

  TFile f(file);
  if (!f.IsOpen() || f.IsZombie()) {
    printf("ERROR: couldn't open the file '%s'\n", file);
    return;
  }
  
 fTree=(TTree*)f.Get("toyMCtree");
  if (!fTree) {
    printf("ERROR: couldn't read the 'toyMCtree' from file '%s'\n", file);
    return;
  }

  fEvent=0x0;
  fTree->SetBranchAddress("event",&fEvent);
  
  // read spacecharge from the Userinfo ot the tree
  InitSpaceCharge();
  
  TString debugName=file;
  debugName.ReplaceAll(".root","");
  debugName.Append(Form(".%1d.%1d_%1d_%1d_%03d_%02d",
                        fUseMaterial,fIdealTracking,fClusterType,
                        Int_t(fCorrectionType),fSeedingRow,fSeedingDist));
  debugName.Append(".debug.root");
  
  gSystem->Exec(Form("test -f %s && rm %s", debugName.Data(), debugName.Data()));
  if (!fStreamer) fStreamer=new TTreeSRedirector(debugName.Data());
  
  gROOT->cd();

  static AliExternalTrackParam dummySeedT0;
  static AliExternalTrackParam dummySeed;
  static AliExternalTrackParam dummyTrack;

  AliExternalTrackParam t0seed;
  AliExternalTrackParam seed;
  AliExternalTrackParam track;
  AliExternalTrackParam tOrig;

  AliExternalTrackParam *dummy;
  
  Int_t maxev=fTree->GetEntries();
  if (nmaxEv>0&&nmaxEv<maxev) maxev=nmaxEv;
  
  for (Int_t iev=0; iev<maxev; ++iev){
    printf("==============  Processing Event %6d =================\n",iev);
    fTree->GetEvent(iev);
    for (Int_t itr=0; itr<fEvent->GetNumberOfTracks(); ++itr){
      printf(" > ======  Processing Track %6d ========  \n",itr);
      const AliToyMCTrack *tr=fEvent->GetTrack(itr);
      tOrig = *tr;

      
      // set dummy 
      t0seed    = dummySeedT0;
      seed      = dummySeed;
      track     = dummyTrack;
      
      Float_t z0=fEvent->GetZ();
      Float_t t0=fEvent->GetT0();

      Float_t vdrift=GetVDrift();
      Float_t zLength=GetZLength(0);

      // crate time0 seed, steered by fCreateT0seed
      printf("t0 seed\n");
      fTime0=-1.;
      fCreateT0seed=kTRUE;
      dummy = GetSeedFromTrack(tr);
      
      if (dummy) {
        t0seed = *dummy;
        delete dummy;

        // crate real seed using the time 0 from the first seed
        // set fCreateT0seed now to false to get the seed in z coordinates
        fTime0 = t0seed.GetZ()-zLength/vdrift;
        fCreateT0seed = kFALSE;
        printf("seed (%.2g)\n",fTime0);
        dummy  = GetSeedFromTrack(tr);
        if (dummy) {
          seed = *dummy;
          delete dummy;

          // create fitted track
          if (fDoTrackFit){
            printf("track\n");
            dummy = GetFittedTrackFromSeed(tr, &seed);
            track = *dummy;
            delete dummy;
          }

          // propagate seed to 0
          const Double_t kMaxSnp = 0.85;
          const Double_t kMass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
//           AliTrackerBase::PropagateTrackTo(&seed,0,kMass,5,kTRUE,kMaxSnp,0,kFALSE,kFALSE);
          
        }
      }

      Int_t ctype(fCorrectionType);
      
      if (fStreamer) {
        (*fStreamer) << "Tracks" <<
        "iev="         << iev             <<
        "z0="          << z0              <<
        "t0="          << t0              <<
        "fTime0="      << fTime0          <<
        "itr="         << itr             <<
        "clsType="     << fClusterType    <<
        "corrType="    << ctype           <<
        "seedRow="     << fSeedingRow     <<
        "seedDist="    << fSeedingDist    <<
        "vdrift="      << vdrift          <<
        "zLength="     << zLength         <<
        "t0seed.="     << &t0seed         <<
        "seed.="       << &seed           <<
        "track.="      << &track          <<
        "tOrig.="      << &tOrig          <<
        "\n";
      }
      
      
    }
  }

  delete fStreamer;
  fStreamer=0x0;

  delete fEvent;
  fEvent = 0x0;
  
  delete fTree;
  fTree=0x0;
  f.Close();
}

//____________________________________________________________________________________
AliExternalTrackParam* AliToyMCReconstruction::GetSeedFromTrack(const AliToyMCTrack * const tr)
{
  //
  // if we don't have a valid time0 informaion (fTime0) available yet
  // assume we create a seed for the time0 estimate
  //

  // seed point informaion
  AliTrackPoint    seedPoint[3];
  const AliTPCclusterMI *seedCluster[3]={0x0,0x0,0x0};
  
  // number of clusters to loop over
  const Int_t ncls=(fClusterType==0)?tr->GetNumberOfSpacePoints():tr->GetNumberOfDistSpacePoints();
  
  UChar_t nextSeedRow=fSeedingRow;
  Int_t   nseeds=0;
  
  //assumes sorted clusters
  for (Int_t icl=0;icl<ncls;++icl) {
    const AliTPCclusterMI *cl=tr->GetSpacePoint(icl);
    if (fClusterType==1) cl=tr->GetDistortedSpacePoint(icl);
    if (!cl) continue;
    // use row in sector
    const UChar_t row=cl->GetRow() + 63*(cl->GetDetector()>35);
    // skip clusters without proper pad row
    if (row>200) continue;
    
    //check seeding row
    // if we are in the last row and still miss a seed we use the last row
    //   even if the row spacing will not be equal
    if (row>=nextSeedRow || icl==ncls-1){
      seedCluster[nseeds]=cl;
      SetTrackPointFromCluster(cl, seedPoint[nseeds]);
      ++nseeds;
      nextSeedRow+=fSeedingDist;
      
      if (nseeds==3) break;
    }
  }
  
  // check we really have 3 seeds
  if (nseeds!=3) {
    AliError(Form("Seeding failed for parameters %d, %d\n",fSeedingDist,fSeedingRow));
    return 0x0;
  }
  
  // do cluster correction for fCorrectionType:
  //   0 - no correction
  //   1 - TPC center
  //   2 - average eta
  //   3 - ideal
  // assign the cluster abs time as z component to all seeds
  for (Int_t iseed=0; iseed<3; ++iseed) {
    Float_t xyz[3]={0,0,0};
    seedPoint[iseed].GetXYZ(xyz);
    const Float_t r=TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
    
    const Int_t sector=seedCluster[iseed]->GetDetector();
    const Int_t sign=1-2*((sector/18)%2);
    
    if ( (fClusterType == 1) && (fCorrectionType != kNoCorrection) ) {
      printf("correction type: %d\n",(Int_t)fCorrectionType);

      // the settings below are for the T0 seed
      // for known T0 the z position is already calculated in SetTrackPointFromCluster
      if ( fCreateT0seed ){
        if ( fCorrectionType == kTPCCenter  ) xyz[2] = 125.*sign;
        //!!! TODO: is this the correct association?
        if ( fCorrectionType == kAverageEta ) xyz[2] = TMath::Tan(45./2.*TMath::DegToRad())*r*sign;
      }
      
      if ( fCorrectionType == kIdeal      ) xyz[2] = seedCluster[iseed]->GetZ();
      
      //!!! TODO: to be replaced with the proper correction
      fSpaceCharge->CorrectPoint(xyz, seedCluster[iseed]->GetDetector());
    }

    // after the correction set the time bin as z-Position in case of a T0 seed
    if ( fCreateT0seed )
      xyz[2]=seedCluster[iseed]->GetTimeBin();
    
    seedPoint[iseed].SetXYZ(xyz);
  }
  
  const Double_t kMaxSnp = 0.85;
  const Double_t kMass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  
  AliExternalTrackParam *seed = AliTrackerBase::MakeSeed(seedPoint[0], seedPoint[1], seedPoint[2]);
  seed->ResetCovariance(10);

  if (fCreateT0seed){
    // if fTime0 < 0 we assume that we create a seed for the T0 estimate
    AliTrackerBase::PropagateTrackTo(seed,0,kMass,5,kTRUE,kMaxSnp,0,kFALSE,kFALSE);
    if (TMath::Abs(seed->GetX())>3) {
      printf("Could not propagate track to 0, %.2f, %.2f, %.2f\n",seed->GetX(),seed->GetAlpha(),seed->Pt());
    }
  }
  
  return seed;
  
}

//____________________________________________________________________________________
void AliToyMCReconstruction::SetTrackPointFromCluster(const AliTPCclusterMI *cl, AliTrackPoint &p )
{
  //
  // make AliTrackPoint out of AliTPCclusterMI
  //
  
  if (!cl) return;
    Float_t xyz[3]={0.,0.,0.};
  //   ClusterToSpacePoint(cl,xyz);
  //   cl->GetGlobalCov(cov);
  //TODO: what to do with the covariance matrix???
  //TODO: the problem is that it is used in GetAngle in AliTrackPoint
  //TODO: which is used by AliTrackerBase::MakeSeed to get alpha correct ...
  //TODO: for the moment simply assign 1 permill squared
  // in AliTrackPoint the cov is xx, xy, xz, yy, yz, zz
  //   Float_t cov[6]={xyz[0]*xyz[0]*1e-6,xyz[0]*xyz[1]*1e-6,xyz[0]*xyz[2]*1e-6,
  //                   xyz[1]*xyz[1]*1e-6,xyz[1]*xyz[2]*1e-6,xyz[2]*xyz[2]*1e-6};
  //   cl->GetGlobalXYZ(xyz);
  //   cl->GetGlobalCov(cov);
  // voluem ID to add later ....
  //   p.SetXYZ(xyz);
  //   p.SetCov(cov);
  AliTrackPoint *tp=const_cast<AliTPCclusterMI*>(cl)->MakePoint();
  p=*tp;
  delete tp;
  //   cl->Print();
  //   p.Print();
  p.SetVolumeID(cl->GetDetector());
  
  
  if ( !fCreateT0seed && !fIdealTracking ) {
    p.GetXYZ(xyz);
    const Int_t sector=cl->GetDetector();
    const Int_t sign=1-2*((sector/18)%2);
    const Float_t zT0=( GetZLength(sector) - (cl->GetTimeBin()-fTime0)*GetVDrift() )*sign;
    printf(" z:  %.2f  %.2f\n",xyz[2],zT0);
    xyz[2]=zT0;
    p.SetXYZ(xyz);
  }
  
  
  //   p.Rotate(p.GetAngle()).Print();
}

//____________________________________________________________________________________
void AliToyMCReconstruction::ClusterToSpacePoint(const AliTPCclusterMI *cl, Float_t xyz[3])
{
  //
  // convert the cluster to a space point in global coordinates
  //
  if (!cl) return;
  xyz[0]=cl->GetRow();
  xyz[1]=cl->GetPad();
  xyz[2]=cl->GetTimeBin(); // this will not be correct at all
  Int_t i[3]={0,cl->GetDetector(),cl->GetRow()};
  //   printf("%.2f, %.2f, %.2f - %d, %d, %d\n",xyz[0],xyz[1],xyz[2],i[0],i[1],i[2]);
  fTPCParam->Transform8to4(xyz,i);
  //   printf("%.2f, %.2f, %.2f - %d, %d, %d\n",xyz[0],xyz[1],xyz[2],i[0],i[1],i[2]);
  fTPCParam->Transform4to3(xyz,i);
  //   printf("%.2f, %.2f, %.2f - %d, %d, %d\n",xyz[0],xyz[1],xyz[2],i[0],i[1],i[2]);
  fTPCParam->Transform2to1(xyz,i);
  //   printf("%.2f, %.2f, %.2f - %d, %d, %d\n",xyz[0],xyz[1],xyz[2],i[0],i[1],i[2]);
}

//____________________________________________________________________________________
AliExternalTrackParam* AliToyMCReconstruction::GetFittedTrackFromSeed(const AliToyMCTrack *tr, const AliExternalTrackParam *seed)
{
  //
  //
  //

  // create track
  AliExternalTrackParam *track = new AliExternalTrackParam(*seed);

  Int_t ncls=(fClusterType == 0)?tr->GetNumberOfSpacePoints():tr->GetNumberOfDistSpacePoints();

  const AliTPCROC * roc = AliTPCROC::Instance();
  
  const Double_t kRTPC0  = roc->GetPadRowRadii(0,0);
  const Double_t kRTPC1  = roc->GetPadRowRadii(36,roc->GetNRows(36)-1);
  const Double_t kMaxSnp = 0.85;
  const Double_t kMaxR   = 500.;
  const Double_t kMaxZ   = 500.;
  
  //   const Double_t kMaxZ0=220;
//   const Double_t kZcut=3;
  
  const Double_t refX = tr->GetX();
  
  const Double_t kMass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  
  // loop over all other points and add to the track
  for (Int_t ipoint=ncls-1; ipoint>=0; --ipoint){
    AliTrackPoint pIn;
    const AliTPCclusterMI *cl=tr->GetSpacePoint(ipoint);
    if (fClusterType == 1) cl=tr->GetDistortedSpacePoint(ipoint);
    SetTrackPointFromCluster(cl, pIn);
    if (fCorrectionType != kNoCorrection){
      Float_t xyz[3]={0,0,0};
      pIn.GetXYZ(xyz);
//       if ( fCorrectionType == kIdeal ) xyz[2] = cl->GetZ();
      fSpaceCharge->CorrectPoint(xyz, cl->GetDetector());
      pIn.SetXYZ(xyz);
    }
    // rotate the cluster to the local detector frame
    track->Rotate(((cl->GetDetector()%18)*20+10)*TMath::DegToRad());
    AliTrackPoint prot = pIn.Rotate(track->GetAlpha());   // rotate to the local frame - non distoted  point
    if (TMath::Abs(prot.GetX())<kRTPC0) continue;
    if (TMath::Abs(prot.GetX())>kRTPC1) continue;
    //
    Bool_t res=kTRUE;
    if (fUseMaterial) res=AliTrackerBase::PropagateTrackTo(track,prot.GetX(),kMass,5,kFALSE,kMaxSnp);
    else res=AliTrackerBase::PropagateTrackTo(track,prot.GetX(),kMass,5,kFALSE,kMaxSnp,0,kFALSE,kFALSE);

    if (!res) break;
    
    if (TMath::Abs(track->GetZ())>kMaxZ) break;
    if (TMath::Abs(track->GetX())>kMaxR) break;
//     if (TMath::Abs(track->GetZ())<kZcut)continue;
    //
    Double_t pointPos[2]={0,0};
    Double_t pointCov[3]={0,0,0};
    pointPos[0]=prot.GetY();//local y
    pointPos[1]=prot.GetZ();//local z
    pointCov[0]=prot.GetCov()[3];//simay^2
    pointCov[1]=prot.GetCov()[4];//sigmayz
    pointCov[2]=prot.GetCov()[5];//sigmaz^2
    
    if (!track->Update(pointPos,pointCov)) {printf("no update\n"); break;}
  }

  if (fUseMaterial) AliTrackerBase::PropagateTrackTo2(track,refX,kMass,5.,kTRUE,kMaxSnp);
  else AliTrackerBase::PropagateTrackTo2(track,refX,kMass,5.,kTRUE,kMaxSnp,0,kFALSE,kFALSE);

  // rotate fittet track to the frame of the original track and propagate to same reference
  track->Rotate(tr->GetAlpha());
  
  if (fUseMaterial) AliTrackerBase::PropagateTrackTo2(track,refX,kMass,1.,kFALSE,kMaxSnp);
  else AliTrackerBase::PropagateTrackTo2(track,refX,kMass,1.,kFALSE,kMaxSnp,0,kFALSE,kFALSE);
  
  return track;
}

//____________________________________________________________________________________
void AliToyMCReconstruction::InitSpaceCharge()
{
  //
  // Init the space charge map
  //

  TString filename="$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2_eps5_50kHz_precal.root";
  if (fTree) {
    TList *l=fTree->GetUserInfo();
    for (Int_t i=0; i<l->GetEntries(); ++i) {
      TString s(l->At(i)->GetName());
      if (s.Contains("SC_")) {
        filename=s;
        break;
      }
    }
  }

  printf("Initialising the space charge map using the file: '%s'\n",filename.Data());
  TFile f(filename.Data());
  fSpaceCharge=(AliTPCSpaceCharge3D*)f.Get("map");
  
  //   fSpaceCharge = new AliTPCSpaceCharge3D();
  //   fSpaceCharge->SetSCDataFileName("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2_eps10_50kHz.root");
  //   fSpaceCharge->SetOmegaTauT1T2(0.325,1,1); // Ne CO2
  // //   fSpaceCharge->SetOmegaTauT1T2(0.41,1,1.05); // Ar CO2
  //   fSpaceCharge->InitSpaceCharge3DDistortion();
  
}

//____________________________________________________________________________________
Double_t AliToyMCReconstruction::GetVDrift() const
{
  //
  //
  //
  return fTPCParam->GetDriftV();
}

//____________________________________________________________________________________
Double_t AliToyMCReconstruction::GetZLength(Int_t roc) const
{
  //
  //
  //
  if (roc<0 || roc>71) return -1;
  return fTPCParam->GetZLength(roc);
}

//____________________________________________________________________________________
TTree* AliToyMCReconstruction::ConnectTrees (const char* files) {
  TString s=gSystem->GetFromPipe(Form("ls %s",files));

  TTree *tFirst=0x0;
  TObjArray *arrFiles=s.Tokenize("\n");
  
  for (Int_t ifile=0; ifile<arrFiles->GetEntriesFast(); ++ifile){
    TString name(arrFiles->At(ifile)->GetName());
    
    TPRegexp reg(".*([0-9]_[0-9]_[0-9]_[0-9]{3}_[0-9]{2}).debug.root");
    TObjArray *arrMatch=0x0;
    arrMatch=reg.MatchS(name);
    
    if (!tFirst) {
      TFile *f=TFile::Open(name.Data());
      if (!f) continue;
      TTree *t=(TTree*)f->Get("Tracks");
      if (!t) {
        delete f;
        continue;
      }
      
      t->SetName(arrMatch->At(1)->GetName());
      tFirst=t;
    } else {
      tFirst->AddFriend(Form("t%s=Tracks",arrMatch->At(1)->GetName()), name.Data());
//       tFirst->AddFriend(Form("t%d=Tracks",ifile), name.Data());
    }
  }

  tFirst->GetListOfFriends()->Print();
  return tFirst;
}

