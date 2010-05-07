// some global var/constants
const Int_t kNRCU = 2; // per SM
AliCaloAltroMapping *fMapping[4]; // 1 for each side (A/C) and each RCU (0/1), i.e. 2*2 total
//
const Bool_t kDebug = kFALSE;
const char *branchStr[] = {"A", "B"};

// help methods
Int_t GetHWAddress(Int_t branch, Int_t FEC, Int_t chip, Int_t chan);
void GetMapping();

// main method
void 
FEEToTower(Int_t iSM=0, Int_t iRCU=0, Int_t branch=0, Int_t FEC=1, Int_t chip=0, Int_t chan=0)
{
  // Get mapping file info, and clear arrays
  GetMapping();

  Int_t isect = iSM / 2; //
  Int_t iside = iSM % 2; // A or C side

  Int_t iRCUSide = iside*2 + iRCU;
  Int_t iDDL = iRCU + iSM * kNRCU;

  Int_t hwAddress = GetHWAddress(branch, FEC, chip, chan);

  Int_t icol = fMapping[iRCUSide]->GetPad(hwAddress);
  Int_t irow = fMapping[iRCUSide]->GetPadRow(hwAddress); 
  Int_t caloflag = fMapping[iRCUSide]->GetSector(hwAddress);
  
  // report channel info
  printf(" iSM %d iRCU %d (iDDL %d EqId %d) : branch %d (%s) FEC %d chip %d chan %d \n corresponds to \n iSM %d icol %d irow %d caloflag (igain) %d \n", 
	 iSM, iRCU, iDDL, iDDL + 0x1200, 
	 branch, branchStr[branch], FEC, chip, chan,
	 iSM, icol, irow, caloflag);

  return;
}

Int_t
GetHWAddress(Int_t branch, Int_t FEC, Int_t chip, Int_t chan)
{
  Int_t hwAddr = (branch << 11) | (FEC << 7) | (chip << 4) | chan; 
  return hwAddr;
}

void 
GetMapping()
{  
  TString sides[]={"A","C"};
  // Read mapping files from $ALICE_ROOT/CALO/mapping/*.data
  TString path = gSystem->Getenv("ALICE_ROOT");
  path += "/EMCAL/mapping/RCU";
  TString path2;
  for(Int_t j = 0; j < 2; j++){ // sides
    for(Int_t i = 0; i < 2; i++) { // RCU
      path2 = path;
      path2 += i;
      path2 += sides[j];
      path2 += ".data";
      if (kDebug) { printf("Mapping file: %s\n",path2.Data()); }
      fMapping[j*2 + i] = new AliCaloAltroMapping(path2.Data());
    }
  }
  return;
}

