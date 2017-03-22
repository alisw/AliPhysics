/// \file CompareOADB_OCDB.C
/// \brief Compare contents in OADB and OCDB
///
/// Compare the contents of the calibration, bad channels, etc in OCDB and OADB
/// You need connexion to the grid to run this macro
/// The OADB file can be in a place different than the default in aliroot
///
/// \author Gustavo Conesa Balbastre, <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-CNRS
///
/// Main method
/// Executes the for the selected run and parameter
///
/// \param run: run numer
/// \param pathOADB: path location of OADB file
/// \param checkObject: Type of comparison, 1-bad channels, 2-energy, 3-energy online, 4-Time
/// \param pass: for some objects, pass string is needed
/// \param printAll: option to print messages
///
void CompareOADB_OCDB(Int_t run = 196000, TString pathOADB = "$ALICE_PHYSICS/OADB/EMCAL",
                      Int_t checkObject = 3, TString pass = "pass3", Bool_t printAll = kFALSE)
{  
  gSystem->Load("libOADB");
  TGrid::Connect("alien://");
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  //man->SetDefaultStorage("local://2013/OCDB");
  man->SetRun(run);
  AliCDBStorage *storage = man->GetDefaultStorage();
  
  // Instantiate EMCAL geometry for the first time

  AliEMCALGeometry * geom = 0;
  if     (run < 140000) geom = AliEMCALGeometry::GetInstance("EMCAL_FIRSTYEARV1"); // 2010
  else if(run < 198000) geom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");  // 2011-2012-2013
  else                  geom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM"); // Run2
  
  if     (checkObject == 0) CheckBadChannels            (run, pathOADB, storage, geom,     printAll);
  else if(checkObject == 1) CheckEnergyCalibration      (run, pathOADB, storage, geom,pass,printAll);
  else if(checkObject == 2) CheckEnergyOnlineCalibration(run, pathOADB, storage, geom,     printAll);
  if     (checkObject == 3) CheckTimeCalibration        (run, pathOADB, storage, geom,pass,printAll);
  else printf("non existing object option\n");
  
  printf("*** Comparisons ended *** \n");
}  

///
/// Check energy recalibration.
//__________________________
void CheckEnergyCalibration(Int_t run, TString pathOADB, AliCDBStorage * storage, AliEMCALGeometry *geom,
                            TString pass,Bool_t printAll = kTRUE)
{  
  AliOADBContainer *contRF=new AliOADBContainer("");
  contRF->InitFromFile(Form("%s/EMCALRecalib.root",pathOADB.Data()),"AliEMCALRecalib");
  
  const Int_t nSM = geom->GetNumberOfSuperModules();
  
  // Get the OCDB object
  
  AliEMCALCalibData* cparam = (AliEMCALCalibData*) (storage->Get("EMCAL/Calib/Data", run)->GetObject());

//  // Access directly the OCDB file and not the latest version
//  TFile * f = TFile::Open("alien:///alice/data/2011/OCDB/EMCAL/Calib/Data/Run144484_999999999_v2_s0.root","READ");
//  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");
//  AliEMCALCalibData* cparam2 = (AliEMCALCalibData*)  cdb->GetObject();

  // Get the OADB object
  
  TH2F *h[12];
  
  TObjArray *recal=(TObjArray*)contRF->GetObject(run);
  
  if(!recal)
  {
    printf("Energy recalibration OADB not found  3\n");
    return;
  }
  
  TObjArray *recalpass=(TObjArray*)recal->FindObject(pass);
  
  if(!recalpass)
  {
    printf("Energy recalibration OADB not found 2\n");
    return;
  }
  
  TObjArray *recalib=(TObjArray*)recalpass->FindObject("Recalib");
  
  if(!recalib)
  {
    printf("Energy recalibration OADB not found 1\n");
    return;
  }
  
  for (Int_t i=0; i < nSM; ++i)
  {
    h[i] = (TH2F*)recalib->FindObject(Form("EMCALRecalFactors_SM%d",i));
    if (!h[i])
    {
      AliError(Form("Could not load EMCALRecalFactors_SM%d",i));
      continue;
    }
  }
  
  // Do the comparison
  Float_t paramOCDB = -1;
  Float_t paramOADB = -1;
  Int_t nDiff = 0;
//  Int_t nDiff2 = 0;
  Int_t iCol = -1, iRow = -1, iSM =-1, iMod = -1,iIphi =-1,iIeta = -1;
  for(Int_t i=0;i < nSM*24*48; i++)
  {
    //printf("AbsID %d\n",i);
    geom->GetCellIndex(i,iSM,iMod,iIphi,iIeta);
    geom->GetCellPhiEtaIndexInSModule(iSM,iMod, iIphi, iIeta,iRow,iCol);
    
    Float_t paramOCDB = -1;
    if(cparam) paramOCDB = cparam->GetADCchannel(iSM,iCol,iRow);

//    Float_t paramOCDB2 = -1;
//    if(cparam) paramOCDB2 = cparam2->GetADCchannel(iSM,iCol,iRow);
    
    Float_t paramOADB = -1;
    if(h[iSM]) paramOADB = h[iSM]->GetBinContent(iCol,iRow);
    paramOADB*=0.0162; // Transformation into OCDB parameter, it will not work for all channels
    
    if    (printAll)
    printf("Recalibration parameter: ID %d, col %d, row %d, sm %d  OCDB %1.4f, OADB %1.4f\n",
           i,iCol,iRow,iSM,paramOCDB, paramOADB);
    else if(paramOADB > 0)
    {
      if (paramOCDB/paramOADB > 1.01 || paramOCDB/paramOADB < 0.99) )
      {
        printf("DIFFERENT recalibration parameter: ID %d, col %d, row %d, sm %d  OCDB %1.4f, OADB %1.4f, ratio OCDB/OADB %1.4f\n",//, old OCDB param %1.4f\n",
               i,iCol,iRow,iSM,paramOCDB, paramOADB,paramOCDB/paramOADB);//,paramOCDB2);
        nDiff++;
      }
      else if(paramOADB <= 0)
      {
        printf("DIFFERENT recalibration parameter: ID %d, col %d, row %d, sm %d  OCDB %1.4f, OADB %1.4f\n",//, old OCDB param %1.4f\n",
               i,iCol,iRow,iSM,paramOCDB, paramOADB);//,paramOCDB2);
        nDiff++;
      }
      
//      if(TMath::Abs(paramOCDB2-0.0162)> 0.0001)
//      {
//        printf("\t Different initial calib! %1.4f\n",paramOCDB2);
//        nDiff2++;
//      }
    }
  }
  
  if(!printAll) printf("Total number of different channels %d / %d\n",nDiff,nSM*24*48);//, origin %d\n",nDiff,nSM*24*48,nDiff2);
}

///
/// Check energy online calibration
//___________________________
void CheckEnergyOnlineCalibration(Int_t run, TString pathOADB, AliCDBStorage * storage, AliEMCALGeometry *geom,
                            Bool_t printAll = kTRUE)
{  
  AliOADBContainer *contRF=new AliOADBContainer("");
  contRF->InitFromFile(Form("%s/EMCALCalibOnlineRef.root",pathOADB.Data()),"AliEMCALRecalib");
  
  const Int_t nSM = geom->GetNumberOfSuperModules();
  
  // Get the OCDB object
  // Access directly the OCDB file and not the latest version
  TFile * f = 0;
  if     (run < 140000)
    f = TFile::Open("alien:///alice/data/2010/OCDB/EMCAL/Calib/Data/Run113461_999999999_v3_s0.root","READ");
  else if(run < 171000)
    f = TFile::Open("alien:///alice/data/2011/OCDB/EMCAL/Calib/Data/Run144484_999999999_v3_s0.root","READ");
  else if(run < 198000)
    f = TFile::Open("alien:///alice/data/2012/OCDB/EMCAL/Calib/Data/Run177115_999999999_v2_s0.root","READ");
  else {
    printf("Run not available\n");
    return;
  }
  
  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");
  AliEMCALCalibData* cparam = (AliEMCALCalibData*)  cdb->GetObject();
  
  // Get the OADB object
  
  TH2F *h[12];
  
  TObjArray *cal=(TObjArray*)contRF->GetObject(run);
  
  if(!cal)
  {
    printf("Energy online calibration OADB not found  3\n");
    return;
  }

  for (Int_t i=0; i < nSM; ++i)
  {
    h[i] = (TH2F*)cal->FindObject(Form("EMCALRecalFactors_SM%d",i));
    if (!h[i])
    {
      AliError(Form("Could not load EMCALRecalFactors_SM%d",i));
      continue;
    }
  }
  
  // Do the comparison
  Float_t paramOCDB = -1;
  Float_t paramOADB = -1;
  Int_t nDiff = 0;
  //  Int_t nDiff2 = 0;
  Int_t iCol = -1, iRow = -1, iSM =-1, iMod = -1,iIphi =-1,iIeta = -1;
  for(Int_t i=0;i < nSM*24*48; i++)
  {
    //printf("AbsID %d\n",i);
    geom->GetCellIndex(i,iSM,iMod,iIphi,iIeta);
    geom->GetCellPhiEtaIndexInSModule(iSM,iMod, iIphi, iIeta,iRow,iCol);
    
    Float_t paramOCDB = -1;
    if(cparam) paramOCDB = cparam->GetADCchannel(iSM,iCol,iRow);
    
    Float_t paramOADB = -1;
    if(h[iSM]) paramOADB = h[iSM]->GetBinContent(iCol,iRow);
    
    if    (printAll)
    printf("Online Calib: ID %d, col %d, row %d, sm %d  OCDB %1.4f, OADB %1.4f\n",
           i,iCol,iRow,iSM,paramOCDB, paramOADB);
    
    if (TMath::Abs(paramOCDB-paramOADB) > 0.001)//  || (TMath::Abs(paramOCDB-0.0162)> 0.001))
    {
      printf("Online Calib DIFFERENT: ID %d, col %d, row %d, sm %d  OCDB %1.4f, OADB %1.4f, diff OCDB-OADB %1.4f\n",
             i,iCol,iRow,iSM,paramOCDB, paramOADB,paramOCDB-paramOADB);
      nDiff++;
    }
    
    
  }
  
  if(!printAll) printf("Total number of different channels %d / %d\n",nDiff,nSM*24*48);//, origin %d\n",nDiff,nSM*24*48,nDiff2);
}

///
/// Get the OCDB bad channels and compare to the OADB ones.
//____________________
void CheckBadChannels(Int_t run, TString pathOADB, AliCDBStorage * storage, AliEMCALGeometry *geom, Bool_t printAll = kTRUE)
{  
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
  
  if(!arrayBC)
  {
    printf("--- Bad map not available \n");
    return;
  }
  
  for (Int_t i=0; i < nSM; ++i)
  {
    hbm[i]=(TH2I*)arrayBC->FindObject(Form("EMCALBadChannelMap_Mod%d",i));
    
    if (!hbm)
    {
      AliError(Form("Can not get EMCALBadChannelMap_Mod%d",i));
      continue;
    }
  }
  
  Int_t badMapOCDB = -1;
  Int_t badMapOADB = -1;
  Int_t iCol = -1, iRow = -1, iSM =-1, iMod = -1,iIphi =-1,iIeta = -1;    
  for(Int_t i = 0; i < nSM*24*48; i++)
  {
    //printf("AbsID %d\n",i);
    geom->GetCellIndex(i,iSM,iMod,iIphi,iIeta); 
    geom->GetCellPhiEtaIndexInSModule(iSM,iMod, iIphi, iIeta,iRow,iCol);
    if(map.At(iSM))badMapOCDB = ((TH2F*)map.At(iSM))->GetBinContent(iCol,iRow);
    else badMapOCDB  = -1;
    if(hbm[iSM])   badMapOADB = hbm[iSM]->GetBinContent(iCol,iRow);
     else badMapOADB = -1;

    if    (printAll && badMapOCDB > 0)
      printf("STATUS: ID %d, col %d, row %d, sm %d  OCDB %d, OADB %d\n",
             i,iCol,iRow,iSM,badMapOCDB, badMapOADB);
    else if(badMapOCDB!=badMapOADB)
      printf("DIFFERENT STATUS: ID %d, col %d, row %d, sm %d  OCDB %d, OADB %d\n",
             i,iCol,iRow,iSM,badMapOCDB, badMapOADB);

  }
}

///
/// Get the OCDB time calibration shifts and compare to the OADB ones.
//________________________
void CheckTimeCalibration(Int_t run, TString pathOADB, AliCDBStorage * storage, AliEMCALGeometry *geom, TString pass, Bool_t printAll = kTRUE)
{
  // Get the OCDB object

  AliEMCALCalibTime* cparam = (AliEMCALCalibTime*) (storage->Get("EMCAL/Calib/Time", run)->GetObject());

  // Access directly the OCDB file and not the latest version
//  TFile * f = TFile::Open("2010/OCDB/EMCAL/Calib/Time/Run122195_126437_v0_s1.root","READ");
//
//  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");
//  AliEMCALCalibTime* cparam = (AliEMCALCalibTime*)  cdb->GetObject();

  if(!cparam)
  {
    printf("OCDB not available\n");
    return;
  }

  // Access OADB
  AliOADBContainer *contTRF=new AliOADBContainer("");
  contTRF->InitFromFile(Form("%s/EMCALTimeCalib.root",pathOADB.Data()),"AliEMCALTimeCalib");
  
  TObjArray *trecal=(TObjArray*)contTRF->GetObject(run); 
  
  if(!trecal)
  {
    return;
    printf("\t array not available\n");
  }
  
  TObjArray *trecalpass=(TObjArray*)trecal->FindObject(pass);
  
  if(!trecalpass)
  {
    return;
    printf("\t pass not available\n");
  }
  
  Int_t iTower  = -1;
  Int_t iIphi   = -1, iIeta   = -1;
  Int_t iSupMod = -1, iSupModMax = -1;
  Int_t iphi = -1, ieta =-1;
  
  for (Int_t ibc = 0; ibc < 4; ++ibc) 
  {
    TH1F *h = (TH1F*)trecalpass->FindObject(Form("hAllTimeAvBC%d",ibc));
    
    printf("\t BC%d -  %p\n",ibc,h);
    
    if (!h) 
    {
      AliError(Form("Could not load hAllTimeAvBC%d",ibc));
      continue;
    }
    
    // Fill params here
        
    const Int_t nSM = static_cast<const Int_t> (geom->GetNumberOfSuperModules());

    for(Int_t absId = 0; absId < nSM*24*48; absId++)
    {
      geom->GetCellIndex(absId,iSupMod,iTower,iIphi,iIeta); 
      geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);    
      
      Float_t paramOADB = h->GetBinContent(absId);
            
      Float_t paramOCDB = cparam->GetTimeChannel(iSupMod,ieta,iphi,ibc);
      
      
      if    (printAll && TMath::Abs(paramOCDB-paramOADB) < 0.001)
        printf("OK: ID %d, col %d, row %d, sm %d  OCDB %d, OADB %d\n",
               absId,ieta,iphi,iSupMod,paramOCDB,paramOADB);
      else if(TMath::Abs(paramOCDB-paramOADB) > 0.001)
        printf("DIFFERENT VALUES: ID %d, col %d, row %d, sm %d  OCDB %d, OADB %d\n",
               absId,ieta,iphi,iSupMod,paramOCDB,paramOADB);
    }
    
    // end
    
  } // bunch crossing loop
}

