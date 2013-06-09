Float_t GetTimeAtVertex(AliToyMCTrack *tr, Int_t clsType=0, Int_t seedRow=140, Int_t seedDist=10, Int_t correctionType=0);
void SetTrackPointFromCluster(AliTPCclusterMI *cl, AliTrackPoint &p);
void ClusterToSpacePoint(AliTPCclusterMI *cl, Float_t xyz[3]);
void InitSpaceCharge();
/*

root.exe -l $ALICE_ROOT/TPC/Upgrade/macros/{loadlibs.C,ConfigOCDB.C}


*/

AliTPCParam *fTPCParam=AliTPCcalibDB::Instance()->GetParameters();
AliTPCSpaceCharge3D *fSpaceCharge=0x0;

void testRec()
{
  //
  //
  //

  TFile f("toyMC.root");
  TTree *t=(TTree*)f.Get("toyMCtree");
  AliToyMCEvent *ev=0x0;
  t->SetBranchAddress("event",&ev);
  t->GetEvent(0);

  gROOT->cd();

//   const Double_t kDriftVel = fTPCParam->GetDriftV()/1000000;
  const Double_t kDriftVel = fTPCParam->GetDriftV();
  const Double_t kMaxZ0=fTPCParam->GetZLength();

  TH1F *h0=new TH1F("h0","h0",1000,0,0);
  TH1F *h1=new TH1F("h1","h1",1000,0,0);
  
  // why are there bad clusters
  for (Int_t itr=0; itr<ev->GetNumberOfTracks(); ++itr){
    printf("Processing Track %6d\n",itr);
    AliToyMCTrack *tr=ev->GetTrack(itr);
    Float_t tVtx0=GetTimeAtVertex(tr);
//     Float_t tVtx1=GetTimeAtVertex(tr,1);
//     if (tVtx0<1e20) h0->Fill(tVtx0);
//     if (tVtx1<1e20) h1->Fill(tVtx1);
  }

  h0->Draw();
  h1->SetLineColor(kRed);
  h1->Draw("same");
}

//____________________________________________________________________________
Float_t GetTimeAtVertex(AliToyMCTrack *tr, Int_t clsType, Int_t seedRow, Int_t seedDist, Int_t correctionType)
{
  //
  // clsType:    0 undistorted; 1 distorted
  // seedRow:    seeding row
  // seedDist:   distance of seeding points
  // correctionType: 0 none, 1 center, 2 full from seed (iterative), 3 ideal (real z-Position)
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
      printf("Seed: %d, %d, %d, %.2f, %.2f, %.2f\n",seed, cl->GetDetector(), row, seedPoint[seed].GetY(),seedPoint[seed].GetX(),seedPoint[seed].GetZ());
      ++seed;
      nextSeedRow+=seedDist;

      if (seed==3) break;
    }
  }

  // check we really have 3 seeds
  if (seed!=3) {
    printf("Seeding failed for parameters %d, %d\n",seedRow,seedDist);
    return 1e30;
  }
  
  // do cluster correction and 
  // assign the cluster abs time as z component to all seeds
  for (Int_t iseed=0; iseed<3; ++iseed) {
    Float_t xyz[3]={0,0,0};
    seedPoint[iseed].GetXYZ(xyz);
    if (clsType && correctionType) {
      if (correctionType==1) xyz[2]=125.;
      if (correctionType==3) xyz[2]=seedCluster->GetZ();

      if (!fSpaceCharge)   InitSpaceCharge();
      fSpaceCharge->CorrectPoint(xyz, cl->GetDetector());
    }

    xyz[2]=cl->GetTimeBin();
    
    seedPoint[iseed].SetXYZ(xyz);
  }

  // create seed and Propagate to r=0;

  const Double_t kMaxSnp = 0.85;
  const Double_t kMass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  
  AliExternalTrackParam *track = 0x0;
  track = AliTrackerBase::MakeSeed(seedPoint[0], seedPoint[1], seedPoint[2]);
  track->ResetCovariance(10);
  
  if (!AliTrackerBase::PropagateTrackTo(track,0,kMass,5,kTRUE,kMaxSnp)) {
//     AliError("Could not propagate track to 0");
    return 1e30;
  }
  
  return track->GetZ();
}

//____________________________________________________________________________
void SetTrackPointFromCluster(AliTPCclusterMI *cl, AliTrackPoint &p ) {
  //
  // make AliTrackPoint out of AliTPCclusterMI
  //
  
  if (!cl) return;
  Float_t xyz[3]={0.,0.,0.};
  Float_t cov[6];
  ClusterToSpacePoint(cl,xyz);
  cl->GetGlobalCov(cov); //what to do with the covariance matrix???
  // voluem ID to add later ....
  p.SetXYZ(xyz);
  p.SetCov(cov);
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
  printf("%.2f, %.2f, %.2f - %d, %d, %d\n",xyz[0],xyz[1],xyz[2],i[0],i[1],i[2]);
  fTPCParam->Transform8to4(xyz,i);
  printf("%.2f, %.2f, %.2f - %d, %d, %d\n",xyz[0],xyz[1],xyz[2],i[0],i[1],i[2]);
  fTPCParam->Transform4to3(xyz,i);
  printf("%.2f, %.2f, %.2f - %d, %d, %d\n",xyz[0],xyz[1],xyz[2],i[0],i[1],i[2]);
  fTPCParam->Transform2to1(xyz,i);
  printf("%.2f, %.2f, %.2f - %d, %d, %d\n",xyz[0],xyz[1],xyz[2],i[0],i[1],i[2]);
}

//____________________________________________________________________________
void InitSpaceCharge()
{
  //
  // Init the space charge map
  //
  fSpaceCharge = new AliTPCSpaceCharge3D();
  fSpaceCharge->SetSCDataFileName("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2_eps5_50kHz.root");
  fSpaceCharge->SetOmegaTauT1T2(0.325,1,1); // Ne CO2
  //fSpaceCharge->SetOmegaTauT1T2(0.41,1,1.05); // Ar CO2
  fSpaceCharge->InitSpaceCharge3DDistortion();
  
}

//____________________________________________________________________________
