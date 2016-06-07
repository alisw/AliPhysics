// $Id$
/**
 * @file makeGlobalHistoConfigObject.C
 * @brief Creation of Global Histo component configuration objects in OCDB
 *
 * <pre>
 * Usage:
 *
 * aliroot -b -q $ALICE_ROOT/HLT/TPCLib/macros/makeTPCFastTransformOCDBObject.C'("uri", runmin, runmax)'
 *
 * </pre>
 *
 * Create an OCDB entry

 * Parameters: <br>
 * - uri   (opt)    the OCDB URI, default $ALICE_ROOT/OCDB   
 * - runmin (opt)   default 0
 * - runmax (opt)   default 999999999

 * @author sergey gorbunov
 * @ingroup alihlt_tutorial
 */

/*
aliroot -b -q $ALICE_ROOT/HLT/TPCLib/macros/makeTPCFastTransformOCDBObject.C'("local:///cvmfs/alice.gsi.de/alice/data/2010/OCDB/", 127941)'
//130704 )'

aliroot -q -b "$ALICE_ROOT/HLT/global/physics/macros/testconfigCalib.C"'("GLOBAL-flat-esd-converter","WriteAnalysisToFile=0 TPCcalibConfigString=TPCCalib:CalibTimeDrift EnableDebug=1 -pushback-period=1000")' $ALICE_ROOT/HLT/exa/recraw-local.C'("/hera/alice/local/filtered/alice/data/2010/LHC10e/000130704/raw/merged_highpt_12.root","local:///cvmfs/alice.gsi.de/alice/data/2010/OCDB/",'1', '1', "HLT", "chains=TPC-compression ignore

 */

void makeTPCFastTransformOCDBObject( const Char_t* cdbUri=NULL,
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

  const AliMagF * field = (AliMagF*) TGeoGlobalMagField::Instance()->GetField();
  calib->SetExBField(field);
  calib->SetRun( runMin );
  calib->UpdateRunInformations( runMin );

  AliHLTTPCClusterTransformation hltTransform;
  
  TStopwatch timer;
  timer.Start();
  //int err = hltTransform.Init( GetBz(), GetTimeStamp() );
  int err = hltTransform.Init( bz, grpObj->GetTimeEnd() );
  timer.Stop();
  cout<<"\n\n Initialisation: "<<timer.CpuTime()<<" / "<<timer.RealTime()<<" sec.\n\n"<<endl;
  if( err!=0 ){
    cerr << Form("Cannot retrieve offline transform from AliTPCcalibDB, AliHLTTPCClusterTransformation returns %d",err).Data()<<endl;
    return -1;
  }
  

  TString path("HLT/ConfigTPC/TPCFastTransform");

  // --------------------------------------
  // -- Create Config Object
  // --------------------------------------

  // here is the actual content of the configuration object
  
  
  const AliHLTTPCFastTransform &fastTransform = hltTransform.GetFastTransform();

  AliHLTTPCFastTransformObject configObj;
  
  Int_t err = fastTransform.WriteToObject( configObj );
  
  if( err ){
    Error("Can not create transformation object, error code %d",err);
  }

  // --------------------------------------
  // -- Fill Object
  // --------------------------------------
  
  AliCDBPath cdbPath(path);
  AliCDBId   cdbId(cdbPath, runMin, runMax);
  AliCDBMetaData cdbMetaData;
  
  printf("Adding %s type OCDB object to %s [%d,%d] in %s \n",
	 configObj.ClassName(), 
	 path.Data(),
	 runMin, runMax, storage.Data());

  
  if(1){ // for test proposes only
    storage="local://$ALICE_ROOT/OCDB";
    man->SetDefaultStorage(storage);
  }

  man->Put( (TObject*)&configObj, cdbId, &cdbMetaData);
  
  
}

