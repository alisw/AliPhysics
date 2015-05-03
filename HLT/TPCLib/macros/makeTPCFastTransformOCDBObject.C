// $Id$
/**
 * @file makeGlobalHistoConfigObject.C
 * @brief Creation of Global Histo component configuration objects in OCDB
 *
 * <pre>
 * Usage:
 *  aliroot -b -q makeTPCFastTransformOCDBObject'("uri", runmin, runmax)'
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
 if(!AliGeomManager::GetGeometry()){
    AliGeomManager::LoadGeometry();
  }

  AliTPCcalibDB *calib=AliTPCcalibDB::Instance();  
  if(!calib){
    HLTError("AliTPCcalibDB does not exist");
    return -ENOENT;
  }
    
  pCalib->SetExBField(-5.);
 
  calib->SetRun( runMin );
  calib->UpdateRunInformations( runMin );

  AliHLTTPCClusterTransformation hltTransform;
  
  TStopwatch timer;
  timer.Start();
  //int err = hltTransform.Init( GetBz(), GetTimeStamp() );
  int err = hltTransform.Init( -5., 0 );
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
  
  const AliHLTTPCFastTransformObject &configObj = hltTransform.GetFastTransform();

  // --------------------------------------
  // -- Fill Object
  // --------------------------------------
  
  if ( 0 ) {
    cerr << "Error : No configuration object created" << endl;
    return;
  }
    
  AliCDBPath cdbPath(path);
  AliCDBId   cdbId(cdbPath, runMin, runMax);
  AliCDBMetaData cdbMetaData;
  man->Put( &configObj, cdbId, &cdbMetaData);

  printf("Adding %s type OCDB object to %s [%d,%d] in %s \n",
	 configObj.ClassName(), 
	 path.Data(),
	 runMin, runMax, storage.Data());
}

