// Example to check the contents of a bad channels OCDB file locally
// Author: Gustavo Conesa Balbastre (LPSC-Grenoble)

void PrintBadChannels(char * file = "$ALICE_ROOT/OCDB/EMCAL/Calib/Pedestals/Run0_999999999_v0_s0.root") 
{
  // Read status map
  TFile * f = new TFile(file,"READ");
  
  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");
  
  AliCaloCalibPedestal * caloped = (AliCaloCalibPedestal *) cdb->GetObject();
  
  TObjArray map = caloped->GetDeadMap();
  
  printf("MAP entries %d\n",map.GetEntries());
  
  AliEMCALGeometry * geom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1");
  
  Int_t nbadTotal  = 0;
  Int_t nhotTotal  = 0;
  Int_t nwarmTotal = 0;
  Int_t ndeadTotal = 0;
  
  for(Int_t iSM = 0; iSM < geom->GetNumberOfSuperModules(); iSM ++)
  { 
    Int_t nbad  = 0;
    Int_t nhot  = 0;
    Int_t nwarm = 0;
    Int_t ndead = 0;
    
    printf(">>> SM %d <<< Entries %d \n",iSM,((TH2D*)map[iSM])->GetEntries());
    
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
