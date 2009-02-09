void trd_qaRec()
{
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libTRDqaRec.so");

  if(!TFile::Open("TRD.DebugInfoGen.root")){
    printf("No debug file for InfoGen task.\n");
    return;
  }
  TTree *t = (TTree*)gFile->Get("trackInfo");
  AliTRDtrackInfo *fTrackInfo = 0x0;
  t->SetBranchAddress("TrackInfo.", &fTrackInfo);
  gROOT->cd();

  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  AliEveEventManager::AssertGeometry();

  AliTRDReconstructor *reco = new AliTRDReconstructor();
  reco->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());

  AliEveTRDTrackList *tracks = new AliEveTRDTrackList("TRD QA Tracks");



  AliTRDtrackV1 *fTrack = 0x0, *track = 0x0;
  for (Int_t it=0; it<t->GetEntries(); it++){
    if(!t->GetEntry(it)) continue;
    if(!fTrackInfo) continue;
    if(!(fTrack = fTrackInfo->GetTrack())) continue;
    
    track = new AliTRDtrackV1(*fTrack);
    track->SetOwner();
    track->SetReconstructor(reco);
    tracks->AddElement(new AliEveTRDTrack(track));
    printf("Trk[%3d] ESD[%d] Ncls[%d]\n", it, fTrackInfo->GetESDinfo()->GetId(), fTrack->GetNumberOfClusters());
    if(it>= 100) break;
  }
  gEve->AddElement(tracks);
  gEve->Redraw3D();
}
