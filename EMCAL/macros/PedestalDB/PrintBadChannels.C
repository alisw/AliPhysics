///
/// \file PrintBadChannels.C
/// \ingroup EMCAL_BadMapDB
/// \brief Print energy calibration parameters in OCDB
///
/// Macro to print the values stored in the OCDB with AliCaloCalibPedestal, channels bad map,
/// either local file or alien file
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

///
/// Main method
///
/// \param file: full file path
///
void PrintBadChannels
(TString file = "$ALICE_ROOT/OCDB/EMCAL/Calib/Pedestals/Run0_999999999_v0_s0.root"
//		          "alien:///alice/data/2014/OCDB/EMCAL/Calib/Pedestals/Run0_999999999_v1_s0.root"
) 
{ 
  if(file.Contains("alien:///"))
    TGrid::Connect("alien://");

  TFile * f = TFile::Open(file,"READ");
  
  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");
  
  AliCaloCalibPedestal * caloped = (AliCaloCalibPedestal *) cdb->GetObject();
  
  TObjArray map = caloped->GetDeadMap();
  Int_t nHisto = map.GetEntries();
  AliEMCALGeometry * geom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");
  
  printf ("N SM %d, array size %d\n",geom->GetNumberOfSuperModules(),nHisto);

  Int_t nbadTotal  = 0;
  Int_t nhotTotal  = 0;
  Int_t nwarmTotal = 0;
  Int_t ndeadTotal = 0;
  
  for(Int_t iSM = 0; iSM < geom->GetNumberOfSuperModules(); iSM ++)
  {     
    if(iSM >=nHisto || !(TH2D*)map[iSM]))
    {
      printf("No entry for SM %d, skip it!\n",iSM);
      continue;
    }

    printf(">>> SM %d <<< Entries %d \n",iSM,((TH2D*)map[iSM])->GetEntries());

    Int_t nbad  = 0;
    Int_t nhot  = 0;
    Int_t nwarm = 0;
    Int_t ndead = 0;

    for(Int_t i = 0; i < ((TH2D*)map[iSM])->GetNbinsX() ; i++)
    {
      for(Int_t j = 0; j < ((TH2D*)map[iSM])->GetNbinsY() ; j++)
      {
        if(((TH2D*)map[iSM])->GetBinContent(i, j)!=AliCaloCalibPedestal::kAlive)
        {
          Int_t id = geom->GetAbsCellIdFromCellIndexes(iSM,j,i);
          printf("\t Bin (%d-%d) Id %d Content: %d \n",i,j,id,((TH2D*)map[iSM])->GetBinContent(i, j));
          nbad++;
        }	
        
        if(((TH2D*)map[iSM])->GetBinContent(i, j)==AliCaloCalibPedestal::kDead   ) ndead++;
        if(((TH2D*)map[iSM])->GetBinContent(i, j)==AliCaloCalibPedestal::kHot    ) nhot ++;
        if(((TH2D*)map[iSM])->GetBinContent(i, j)==AliCaloCalibPedestal::kWarning) nwarm++;
      }	
    }
    
    printf("Summary : --- dead %d --- hot %d --- warm %d --- bad %d\n",ndead,nhot,nwarm,nbad);  
    
    nbadTotal  += nbad  ;
    nhotTotal  += nhot  ;
    ndeadTotal += ndead ;
    nwarmTotal += nwarm ;
  }
  
  
  printf("All SM summary : --- dead %d --- hot %d --- warm %d --- bad %d\n",ndeadTotal,nhotTotal,nwarmTotal,nbadTotal);  
  
  printf("Total BAD %d\n", caloped->GetDeadTowerCount());
}
