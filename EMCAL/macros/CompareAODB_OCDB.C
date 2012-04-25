// Compare the contents of the calibration, bad channels, etc in OCDB and AODB
// You need connexion to the grid to run this macro
// The OADB file can be in a place different than the default in aliroot

// Author : Gustavo Conesa Balbastre (LPSC-CNRS)

void CompareAODB_OCDB(Int_t run = 146806, TString pathOADB = "$ALICE_ROOT/OADB/EMCAL")
{  

  gSystem->Load("libOADB");
  TGrid::Connect("alien://");
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  man->SetRun(run);
  AliCDBStorage *storage = man->GetDefaultStorage();
  
  // Instantiate EMCAL geometry for the first time

  
  if     (run < 140000) geom = AliEMCALGeometry::GetInstance("EMCAL_FIRSTYEARV1");
  else if(run < 171000) geom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1"); 
  else                  geom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1"); 
      
  CheckBadChannels(run, pathOADB, storage, geom);
  //CheckEnergyCalibration(run, pathOADB, storage, geom);
  //CheckTimeCalibration(run, pathOADB, storage, geom);
  
}  

//___________________________
void CheckEnergyCalibration(Int_t run, TString pathOADB, AliCDBStorage * storage, AliEMCALGeometry *geom)
{
    // Check energy recalibration, Uncomplete
  
    AliOADBContainer *contRF=new AliOADBContainer("");
    contRF->InitFromFile(Form("%s/EMCALRecalib.root",pathOADB.Data()),"AliEMCALRecalib");
  
    const Int_t nSM = geom->GetNumberOfSuperModules();

    TH2I *h[12];
  
    TObjArray *recal=(TObjArray*)contRF->GetObject(run); 
    if(recal)
    {
      TObjArray *recalpass=(TObjArray*)recal->FindObject(pass);
      
      if(recalpass)
      {
        TObjArray *recalib=(TObjArray*)recalpass->FindObject("Recalib");
        
        if(recalib)
        {
          for (Int_t i=0; i < nSM; ++i) 
          {
            h[i] = (TH2F*)recalib->FindObject(Form("EMCALRecalFactors_SM%d",i));
            if (!h[i]) 
            {
              AliError(Form("Could not load EMCALRecalFactors_SM%d",i));
              continue;
            }
          } 
        }else printf("Energy recalibration OADB not found 1\n");
      }else printf("Energy recalibration OADB not found 2\n");
    }else printf("Energy recalibration OADB not found  3\n");
    
    // TO COMPLETE
}  

//_____________________
void CheckBadChannels(Int_t run, TString pathOADB, AliCDBStorage * storage, AliEMCALGeometry *geom)
{
  // Get the OCDB bad channels and compare to the OADB ones
  
  //const Int_t nSM = static_const (geom->GetNumberOfSuperModules());
  const Int_t nSM = static_cast<const Int_t> (geom->GetNumberOfSuperModules());

  // Access OCDB histograms
  AliCaloCalibPedestal* caloped  = (AliCaloCalibPedestal*) (storage->Get("EMCAL/Calib/Pedestals", run)->GetObject());
    
  // Access directly the OCDB file and not the latest version
  //TFile * f = TFile::Open("alien:///alice/data/2011/OCDB/EMCAL/Calib/Pedestals/Run145954_146856_v3_s0.root","READ");
  //AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");
  //AliCaloCalibPedestal * caloped = (AliCaloCalibPedestal *) cdb->GetObject();  
  
  TObjArray map = caloped->GetDeadMap();

  // Access OADB histograms
  TH2I *hbm[12];
  
  AliOADBContainer *contBC=new AliOADBContainer("");
  
  contBC->InitFromFile(Form("%s/EMCALBadChannels.root",pathOADB.Data()),"AliEMCALBadChannels"); 
  
  TObjArray *arrayBC=(TObjArray*)contBC->GetObject(run);
  
  if(arrayBC)
  {
    for (Int_t i=0; i < nSM; ++i) 
    {
      hbm[i]=(TH2I*)arrayBC->FindObject(Form("EMCALBadChannelMap_Mod%d",i));
      
      if (!hbm) 
      {
        AliError(Form("Can not get EMCALBadChannelMap_Mod%d",i));
        continue;
      }
    }
  } 
  else 
  { 
    printf("--- Bad map not available \n");
    return;
  }

  Int_t badMapOCDB = -1;
  Int_t badMapAODB = -1;
  Int_t iCol = -1, iRow = -1, iSM =-1, iMod = -1,iIphi =-1,iIeta = -1;    
  for(Int_t i=0;i < nSM*24*48; i++)
  {
    //printf("AbsID %d\n",i);
    geom->GetCellIndex(i,iSM,iMod,iIphi,iIeta); 
    geom->GetCellPhiEtaIndexInSModule(iSM,iMod, iIphi, iIeta,iRow,iCol);
    if(map.At(iSM))badMapOCDB = ((TH2F*)map.At(iSM))->GetBinContent(iCol,iRow);
    else badMapOCDB  = -1;
    if(hbm[iSM])   badMapAODB = hbm[iSM]->GetBinContent(iCol,iRow);
     else badMapAODB = -1;

    //if(badMapOCDB>0)
      if(badMapOCDB!=badMapAODB)  
      printf("DIFFERENT STATUS: ID %d, col %d, row %d, sm %d  OCDB %d, OADB %d\n",
             i,iCol,iRow,iSM,badMapOCDB, badMapAODB);

  }
  
}


