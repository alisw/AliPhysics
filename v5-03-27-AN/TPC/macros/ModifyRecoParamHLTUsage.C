/* ModifyRecoParamHLTUsage
 *  Changes the Input for Reconstruction ( TPC RAW data or HLT TPC clusters ) 
 *    1 -> only TPC raw/sim data
 *    2 -> if present TPC raw/sim data, otherwise HLT clusters
 *    3 -> only HLT clusters
 *    4 -> if present HLT clusters, otherwise TPC raw/sim data
 *
 * Usage : aliroot -b -q ModifyRecoParamHLTUsage.C'("/lustre/alice/alien/alice/data/2011/OCDB/TPC/Calib/RecoParam/Run136844_999999999_v2_s0.root",3,"local:///tmp/ocdb/")'
 *
 */

void ModifyRecoParamHLTUsage( const Char_t* lastOCDBEntry, 
			      Int_t iHLTusage,
			      const Char_t* newStoragePath = "local:///tmp/ocdb", 
			      const Char_t* author = "Jochen Thaeder", 
			      const Char_t* alirootVersion = "05-01-Release" ) {
  
  // -- Get RecoParam
  // -------------------------------------------------------------------
  TFile inFile(lastOCDBEntry);
  if (!inFile) {
    printf("File %s could not be found!\n", lastOCDBEntry);
    return -1;
  }

  AliCDBEntry * entry = ( AliCDBEntry *) inFile.Get("AliCDBEntry");
  if (!entry) {
    printf("File %s does not contain AliCDBEntry object!\n", lastOCDBEntry);
    return -2;
  }

  TObjArray * arr = entry->GetObject();

  AliTPCRecoParam * parHighFlux = NULL;
  AliTPCRecoParam * parLowFlux  = NULL;
  AliTPCRecoParam * parLaserFlux  = NULL;
  AliTPCRecoParam * parCosmicFlux  = NULL;

  for (Int_t idx = 0 ; idx < arr->GetEntries() ; ++idx) {
    AliTPCRecoParam *param = (AliTPCRecoParam*) arr->At(idx);

    TString name(param->GetName());

    if ( name.Contains("Low") )
      parLowFlux = param;
    else if ( name.Contains("High") )
      parHighFlux = param;
    else if ( name.Contains("Laser") )
      parLaserFlux = param;
    else if ( name.Contains("Cosmic") )
      parCosmicFlux = param;
  }

  // -- Set TPC RecoParam : UseHLTClusters
  // -------------------------------------------------------------------
  parLaserFlux->SetUseHLTClusters(1); // LASER use always raw data

  parHighFlux->SetUseHLTClusters(iHLTusage);
  parLowFlux->SetUseHLTClusters(iHLTusage);
  parCosmicFlux->SetUseHLTClusters(iHLTusage);

  // -- Write out
  // -------------------------------------------------------------------
  const Char_t *comment = "Modified RecoParam with HLT clusters usage";
  // -------------------------------------------------------------------
  
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible(author);
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion(alirootVersion); 
  metaData->SetComment(comment);

  AliCDBId id("TPC/Calib/RecoParam", 0, AliCDBRunRange::Infinity());

  AliCDBStorage * gStorage = AliCDBManager::Instance()->GetStorage(newStoragePath);
  gStorage->Put(arr, id, metaData);    

  return;
}

