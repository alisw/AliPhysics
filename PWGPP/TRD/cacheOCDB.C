void cacheOCDB(Int_t run, Int_t year)
{    
  AliCDBManager* ocdb = AliCDBManager::Instance();
  ocdb->SetDefaultStorage(Form("alien://folder=/alice/data/%d/OCDB?cacheFolder=%s/local", year, gSystem->ExpandPathName("$HOME")));
  ocdb->SetRun(run);
  // create geo manager
  if(!(obj = ocdb->Get(AliCDBPath("GRP", "Geometry", "Data")))){
    AliError("GEOMETRY failed initialization.");
  } else {
    AliGeomManager::SetGeometry((TGeoManager*)obj->GetObject());
    AliGeomManager::GetNalignable("TRD");
    AliGeomManager::GetNalignable("TPC");
    AliGeomManager::GetNalignable("ITS");
    AliGeomManager::ApplyAlignObjsFromCDB("ITS TPC TRD");
  }
  // get GRP/GRP/Data
  if(!(obj = ocdb->Get(AliCDBPath("GRP", "GRP", "Data")))){
    AliError("GRP failed initialization.");
  }
  // get HMPID/Calib/Nmean
  if(!(obj = ocdb->Get(AliCDBPath("HMPID", "Calib", "Nmean")))){
    AliError("HMPID Nmean failed initialization.");
  }
  //init magnetic field
  if(!TGeoGlobalMagField::Instance()->GetField() && !TGeoGlobalMagField::Instance()->IsLocked()){
    AliGRPManager grpManager;
    if(!grpManager.ReadGRPEntry()) AliError("Cannot get GRP entry"); 
    if(!grpManager.SetMagField()) AliError("Problem with magnetic field setup"); 
  }

  const Int_t nobj(14);
  const Char_t *tobj[nobj] = {"ChamberExB", 
                              "ChamberT0", 
                              "ChamberVdrift", 
                              "ChamberStatus", 
                              "ChamberGainFactor", 
                              "LocalGainFactor", 
                              "LocalVdrift", 
                              "LocalT0", 
                              "PRFWidth", 
                              "TrkAttach", 
                              "PIDLQ1D", 
                              "PHQ", 
                              "DCS",
                              "RecoParam"};
  for(Int_t iobj(0); iobj<nobj; iobj++){
    // get TRD object
    if(!(obj = ocdb->Get(AliCDBPath("TRD", "Calib", tobj[iobj])))){
      AliError(Form("TRD \"%s\" failed initialization.", tobj[iobj]));
    }
  }
}
