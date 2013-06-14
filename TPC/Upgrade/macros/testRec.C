Int_t GetTimeAtVertex(Float_t &tVtx, Float_t &x, AliToyMCTrack *tr, Int_t clsType=0, Int_t seedRow=140, Int_t seedDist=10, Int_t correctionType=0);
void SetTrackPointFromCluster(AliTPCclusterMI *cl, AliTrackPoint &p);
void ClusterToSpacePoint(AliTPCclusterMI *cl, Float_t xyz[3]);
void InitSpaceCharge();
/*

root.exe -l $ALICE_ROOT/TPC/Upgrade/macros/{loadlibs.C,ConfigOCDB.C}
.x $ALICE_ROOT/TPC/Upgrade/macros/testRec.C

*/

AliTPCParam *fTPCParam=AliTPCcalibDB::Instance()->GetParameters();
AliTPCSpaceCharge3D *fSpaceCharge=0x0;
TTreeSRedirector *fStreamer=0x0;

Int_t fEvent=-1;
Int_t fTrack=-1;
Float_t fT0event=-1;
Float_t fZevent=-1;

void testRec(Int_t nmaxEv=-1)
{
  //
  //
  //

  TFile f("toyMC.root");
  TTree *t=(TTree*)f.Get("toyMCtree");
  AliToyMCEvent *ev=0x0;
  t->SetBranchAddress("event",&ev);

  gSystem->Exec("rm debug.root");
  if (!fStreamer) fStreamer=new TTreeSRedirector("debug.root");
  
  gROOT->cd();

//   const Double_t kDriftVel = fTPCParam->GetDriftV()/1000000;
  const Double_t kDriftVel = fTPCParam->GetDriftV();
  const Double_t kMaxZ0=fTPCParam->GetZLength();

  TH1F *h0=new TH1F("h0","h0",1000,0,0);
  TH1F *hX=new TH1F("hX","hX",1000,0,0);
  TH1F *h1=new TH1F("h1","h1",1000,0,0);
  TH1I *hcount0=new TH1I("count0","Failed extrapolation1",5,0,5);
  TH1I *hcount1=new TH1I("count1","Failed extrapolation2",5,0,5);
  
  Int_t maxev=t->GetEntries();
  if (nmaxEv>0&&nmaxEv<maxev) maxev=nmaxEv;
  
  for (Int_t iev=0; iev<maxev; ++iev){
    t->GetEvent(iev);
    fEvent=iev;
    for (Int_t itr=0; itr<ev->GetNumberOfTracks(); ++itr){
      printf("==============  Processing Track %6d\n",itr);
      fTrack=itr;
      fT0event = ev->GetT0();
      fZevent  = ev->GetZ();

      //Float_t &tVtx, Float_t &x, AliToyMCTrack *tr,
      //  Int_t clsType=0, Int_t seedRow=140, Int_t seedDist=10, Int_t correctionType=0
      // correctionType: 0 none, 1 center, 2 mean tan,
      //                 3 full from seed (iterative), 4 ideal (real z-Position)
      AliToyMCTrack *tr=ev->GetTrack(itr);
      tr->SetUniqueID(itr);
      Float_t tVtx0=0;
      Float_t xmin=0;
      Int_t ret0=GetTimeAtVertex(tVtx0,xmin,tr);
      hX->Fill(xmin);
      Float_t tVtx1=0;
      Int_t ret1=GetTimeAtVertex(tVtx1,xmin,tr,1);
      //fully distorted
      GetTimeAtVertex(tVtx1,xmin,tr,1); // seeding at the outside
      GetTimeAtVertex(tVtx1,xmin,tr,1,70); // seeding in the center
      GetTimeAtVertex(tVtx1,xmin,tr,1,0); // seeding at the inside
      //correction at tpc center
      GetTimeAtVertex(tVtx1,xmin,tr,1,140, 10, 1);
      //correction with mean tan theta
      GetTimeAtVertex(tVtx1,xmin,tr,1,140, 10, 2);
      
      hcount0->Fill(ret0);
      hcount1->Fill(ret1);
      if (ret0==0) {
        h0->Fill(tVtx0);
      }

      if (ret1==0) {
        h1->Fill(tVtx1);
      }
//       printf("TVtx: %f, %f\n",tVtx0,0);
    }
  }

  TCanvas *c=(TCanvas*)gROOT->GetListOfCanvases()->FindObject("cOutput");
  if (!c) c=new TCanvas("cOutput","Results");
  c->Clear();
  c->Divide(2,2);
  c->cd(1);
  h0->Draw();
  h1->SetLineColor(kRed);
  h1->Draw("same");
  c->cd(2);
  hcount0->Draw();
  hcount1->SetLineColor(kRed);
  hcount1->Draw("same");
  c->cd(3);
  hX->Draw();

  delete fStreamer;
  fStreamer=0x0;
}

//____________________________________________________________________________
Float_t GetTimeAtVertex(Float_t &tVtx,  Float_t &x, AliToyMCTrack *tr, Int_t clsType, Int_t seedRow, Int_t seedDist, Int_t correctionType)
{
  //
  // clsType:    0 undistorted; 1 distorted
  // seedRow:    seeding row
  // seedDist:   distance of seeding points
  // correctionType: 0 none, 1 center, 2 mean tan,
  //                 3 full from seed (iterative), 4 ideal (real z-Position)
  //

  // seed point informaion
  AliTrackPoint    seedPoint[3];
  AliTPCclusterMI *seedCluster[3]={0x0,0x0,0x0};

  // number of clusters to loop over
  const Int_t ncls=(clsType==0)?tr->GetNumberOfSpacePoints():tr->GetNumberOfDistSpacePoints();

  UChar_t nextSeedRow=seedRow;
  Int_t   seed=0;

  //assumes sorted clusters
  for (Int_t icl=0;icl<ncls;++icl) {
    AliTPCclusterMI *cl=tr->GetSpacePoint(icl);
    if (!cl) continue;
    // use row in sector
    const UChar_t row=cl->GetRow() + 63*(cl->GetDetector()>35);
    // skip clusters without proper pad row
    if (row>200) continue;

    //check seeding row
    // if we are in the last row and still miss a seed we use the last row
    //   even if the row spacing will not be equal
    if (row>=nextSeedRow || icl==ncls-1){
      seedCluster[seed]=cl;
      SetTrackPointFromCluster(cl, seedPoint[seed]);
//       printf("\nSeed point %d: %d, %d, %.2f, %.2f, %.2f, %.2f, %.2f\n",seed, cl->GetDetector(), row, seedPoint[seed].GetX(),seedPoint[seed].GetY(),seedPoint[seed].GetZ(), seedPoint[seed].GetAngle(), ((cl->GetDetector()%18)*20.+10.)/180.*TMath::Pi());
      ++seed;
      nextSeedRow+=seedDist;

      if (seed==3) break;
    }
  }

  // check we really have 3 seeds
  if (seed!=3) {
    printf("Seeding failed for parameters %d, %d\n",seedRow,seedDist, seed);
    return 1;
  }
  
  // do cluster correction and 
  // assign the cluster abs time as z component to all seeds
  for (Int_t iseed=0; iseed<3; ++iseed) {
    Float_t xyz[3]={0,0,0};
    seedPoint[iseed].GetXYZ(xyz);
    
    Int_t sector=seedCluster[iseed]->GetDetector();
    Int_t sign=1-2*((sector/18)%2);
    
    if (clsType && correctionType) {
      if (correctionType==1) xyz[2]=125.;
      if (correctionType==2) xyz[2]=TMath::Tan(45./2.*TMath::DegToRad())*xyz[1]*sign;
//       if (correctionType==3) xyz[2]=125.;
      if (correctionType==4) xyz[2]=seedCluster->GetZ();

      if (!fSpaceCharge)   InitSpaceCharge();
      fSpaceCharge->CorrectPoint(xyz, seedCluster[iseed]->GetDetector());
    }

    //set different sign for c-Side
    xyz[2]=seedCluster[iseed]->GetTimeBin()/* * sign*/;
    seedPoint[iseed].SetXYZ(xyz);
  }

  // create seed and Propagate to r=0;

  const Double_t kMaxSnp = 0.85;
  const Double_t kMass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  
  AliExternalTrackParam *track = 0x0;
  track = AliTrackerBase::MakeSeed(seedPoint[0], seedPoint[1], seedPoint[2]);
  track->ResetCovariance(10);
//   printf("Track: %.2f, %.2f, %.2f, %.2f, %.2f\n",track->GetX(),track->GetY(),track->GetZ(), track->GetAlpha(),track->Phi());
  AliExternalTrackParam pInit(*track);

  // NOTE:
  // when propagating with the time binwe need to switch off the material correction
  // otherwise we will be quite off ...
  // 
  AliTrackerBase::PropagateTrackTo(track,0,kMass,5,kTRUE,kMaxSnp,0,kFALSE,kFALSE);
  if (TMath::Abs(track->GetX())>3) {
    printf("Could not propagate track to 0, %.2f, %.2f, %.2f\n",track->GetX(),track->GetAlpha(),track->Pt());
    return 2;
  }
//   printf("Track2: %.2f, %.2f, %.2f, %.2f\n",track->GetX(),track->GetY(),track->GetZ(), track->GetAlpha());

  // simple linear fit
  TGraph gr;
  for (Int_t i=0; i<3; ++i)
    gr.SetPoint(gr.GetN(),seedPoint[i].GetX(),seedPoint[i].GetZ());
//   gr.Print();
  TF1 fpol1("fpol1","pol1");
  gr.Fit(&fpol1,"QN");
  Float_t fitT0=fpol1.Eval(0);
//   fpol1.Print();
  AliExternalTrackParam pOrig(*tr);
  (*fStreamer) << "Tracks" <<
    "iev="      << fEvent <<
    "t0="       << fT0event <<
    "z0="       << fZevent <<
    "itrack="   << fTrack <<
    "clsType="  << clsType <<
    "seedRow="  << seedRow <<
    "seedDist=" << seedDist <<
    "corrType=" << correctionType <<
    "track.="   << tr    <<
    "seed.="    << track <<
    "seedI.="   << &pInit <<
//     "seedcl0.=" << seedCluster[0] <<
//     "seedcl1.=" << seedCluster[1] <<
//     "seedcl2.=" << seedCluster[2] <<
//     "seedp0.="  << &seedPoint[0] <<
//     "seedp1.="  << &seedPoint[1] <<
//     "seedp2.="  << &seedPoint[2] <<
    "fitT0="    << fitT0 <<
    "\n";
  
  tVtx=track->GetZ();
  x=track->GetX();
  delete track;
  return 0;
}

//____________________________________________________________________________
void SetTrackPointFromCluster(AliTPCclusterMI *cl, AliTrackPoint &p ) {
  //
  // make AliTrackPoint out of AliTPCclusterMI
  //

  if (!cl) return;
//   Float_t xyz[3]={0.,0.,0.};
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
  AliTrackPoint *tp=cl->MakePoint();
  p=*tp;
  delete tp;
//   cl->Print();
//   p.Print();
  p.SetVolumeID(cl->GetDetector());
//   p.Rotate(p.GetAngle()).Print();
}

//____________________________________________________________________________
void ClusterToSpacePoint(AliTPCclusterMI *cl, Float_t xyz[3])
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

//____________________________________________________________________________
void InitSpaceCharge()
{
  //
  // Init the space charge map
  //
  TFile f("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2_eps5_50kHz_precal.root");
  fSpaceCharge=(AliTPCSpaceCharge3D*)f.Get("map");
  
//   fSpaceCharge = new AliTPCSpaceCharge3D();
//   fSpaceCharge->SetSCDataFileName("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2_eps10_50kHz.root");
//   fSpaceCharge->SetOmegaTauT1T2(0.325,1,1); // Ne CO2
// //   fSpaceCharge->SetOmegaTauT1T2(0.41,1,1.05); // Ar CO2
//   fSpaceCharge->InitSpaceCharge3DDistortion();
  
}

//____________________________________________________________________________


