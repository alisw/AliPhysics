
/*
aliroot -b -q $ALICE_ROOT/HLT/TPCLib/macros/testTPCFastTransformOCDBObject.C'(NULL, 127941 )'
*/

void testTPCFastTransformOCDBObject( const Char_t* cdbUri=NULL,
				     Int_t runMin=0, Int_t runMax=AliCDBRunRange::Infinity()) {

  // --------------------------------------
  // -- Setup CDB
  // --------------------------------------

  AliCDBManager* man = AliCDBManager::Instance();
  if (!man) {
    cerr << "Error : Can not get AliCDBManager" << end;
    exit;
  }

  TString storage;
  if (!man->IsDefaultStorageSet()) {
    if ( cdbUri ) {
      storage = cdbUri;
      if ( storage.Contains("://") == 0 ) {
	storage = "local://"; 
	storage += cdbUri;
      }
    } 
    else {
      storage="local://$ALICE_ROOT/OCDB";
    }
    man->SetDefaultStorage(storage);
  } 
  else {
    storage = man->GetDefaultStorage()->GetURI();
  }

  man->SetRun( runMin );
  AliGRPManager grp;
  grp.ReadGRPEntry();
  grp.SetMagField();

  const AliGRPObject *grpObj = grp.GetGRPData();
  
  if( !grpObj ){
    Error("No GRP object found!!");
  }

  if(!AliGeomManager::GetGeometry()){
    AliGeomManager::LoadGeometry();
  }
 
  AliTPCcalibDB *calib=AliTPCcalibDB::Instance();  
  if(!calib){
    Error("AliTPCcalibDB does not exist");
    return -ENOENT;
  }
    
  Double_t bz = AliTracker::GetBz();

  cout<<"\n\nBz field is set to "<<bz<<", time stamp is set to "<<grpObj->GetTimeEnd()<<endl<<endl;

  calib->SetExBField(bz);  
  calib->SetRun( runMin );
  calib->UpdateRunInformations( runMin );  

  
  AliHLTTPCClusterTransformation hltTransform;  
  TStopwatch timer;
  timer.Start();
  int err = hltTransform.Init( bz, grpObj->GetTimeEnd() );
  timer.Stop();
  cout<<"\n\n Initialisation: "<<timer.CpuTime()<<" / "<<timer.RealTime()<<" sec.\n\n"<<endl;
  
  if( err!=0 ){
    cerr << Form("Cannot retrieve offline transform from AliTPCcalibDB, AliHLTTPCClusterTransformation returns %d",err).Data()<<endl;
    return -1;
  }
  
  const AliHLTTPCFastTransform &fastTransform = hltTransform.GetFastTransform();
  cout<<"create config object.."<<endl;
  AliHLTTPCFastTransformObject configObj;
  cout<<"call FastTransform::WriteToObject().."<<endl;
  Int_t err = fastTransform.WriteToObject( configObj );
  if( err ){
    Error("Can not create transformation object, error code %d",err);
  }
  cout<<"Write to object ok."<<endl;

  {
    const AliHLTTPCSpline2D3DObject &spl = configObj.GetSplineNonConst(0,10,1);
    float x=0, y=0, z=0;
    spl.GetGridValue(1,x,y,z);
    cout<<"object to write to OCDB: n points "<<spl.GetNPoints()<<" x "<<x<<" y "<<y<<" z "<<z<<endl; 
  }

  TString path("HLT/ConfigTPC/TPCFastTransform");
  AliCDBPath cdbPath(path);
  AliCDBId   cdbId(cdbPath, runMin, runMax);
  AliCDBMetaData cdbMetaData;

  printf("Adding %s type OCDB object to %s [%d,%d] in %s \n",
	 configObj.ClassName(), 
	 path.Data(),
	 runMin, runMax, storage.Data());
  
  storage="local://$ALICE_ROOT/OCDB";
  man->SetDefaultStorage(storage);
  
  man->Put( (TObject*)configObj, cdbId, &cdbMetaData);

  /*
  cout<<" try to write to a file.."<<endl;
  configObj.SaveAs("bla.root");
  */
  cout<<" now try to read the object from OCDB.."<<endl;

  AliCDBEntry *entr = man->Get(cdbPath);
  AliHLTTPCFastTransformObject *obj = dynamic_cast<AliHLTTPCFastTransformObject*> (entr->GetObject());
  if( !obj ) {
    Error("can not read object from OCDB");
  }

  {
    const AliHLTTPCSpline2D3DObject &spl1 = obj->GetSpline( 0,10,1);
    float x=0, y=0, z=0;
    spl1.GetGridValue(1,x,y,z);    
    cout<<"object read from OCDB: npoints "<<spl1.GetNPoints()<<" x "<<x<<" y "<<y<<" z "<<z<<endl; 
  } 

}

