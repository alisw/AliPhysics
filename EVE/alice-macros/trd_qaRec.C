void trd_qaRec()
{
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libTRDqaRec.so");

  if(!TFile::Open("TRD.DebugInfoGen.root")){
    printf("No debug file for InfoGen task.\n");
    return;
  }
  TTree *t = (TTree*)gFile->Get("trackInfo");
  AliTRDtrackInfo *fTrackInfo = new AliTRDtrackInfo();
  t->SetBranchAddress("TrackInfo.", &fTrackInfo);
  gROOT->cd();

  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT");
  cdb->SetRun(0);
  AliEveEventManager::AssertGeometry();
  AliEveTRDTrackList *tracks = new AliEveTRDTrackList("TRD QA Tracks");



  AliTRDtrackV1 *fTrack = 0x0;
  AliTRDReconstructor *reco = new AliTRDReconstructor();
  for (Int_t it=0; it<t->GetEntries(); it++){
    if(!t->GetEntry(it)) continue;
    if(!fTrackInfo) continue;
    if(!(fTrack = fTrackInfo->GetTRDtrack())) continue;
    
    fTrack->SetReconstructor(reco);
    tracks->AddElement(new AliEveTRDTrack(fTrack));
    //printf("Trk[%3d] ESD[%d] Ncls[%d]\n", it, fTrackInfo->GetESDinfo()->GetId(), fTrack->GetNumberOfClusters());
  }
  gEve->AddElement(tracks);
  gEve->Redraw3D();
}