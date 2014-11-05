AliITSUReconstructor* rec = 0;
AliITSUTrackerGlo* trk=0;

void chkTracker(int run=0)
{
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetSpecificStorage("GRP/GRP/Data",
  			  Form("local://%s",gSystem->pwd()));
  man->SetSpecificStorage("ITS/Align/Data",
			  Form("local://%s",gSystem->pwd()));
  man->SetSpecificStorage("ITS/Calib/RecoParam",
			  Form("local://%s",gSystem->pwd()));
  man->SetRun(run);
  if ( !TGeoGlobalMagField::Instance()->GetField() ) {
    printf("Loading field map...\n");
    AliGRPManager grpMan;
    if( !grpMan.ReadGRPEntry() ) { 
      printf("Cannot get GRP entry\n"); 
    }
    if( !grpMan.SetMagField() ) { 
      printf("Problem with magnetic field setup\n"); 
    }
  }
  AliGeomManager::LoadGeometry("geometry.root");
  //
  //
  AliCDBEntry* ent = man->Get("ITS/Calib/RecoParam");
  AliITSURecoParam* par = (AliITSURecoParam*)((TObjArray*)ent->GetObject())->At(1);
  //
  rec = new AliITSUReconstructor();
  rec->SetRecoParam(par);
  //
  rec->Init();
  trk = (AliITSUTrackerGlo*) rec->CreateTracker();
}
