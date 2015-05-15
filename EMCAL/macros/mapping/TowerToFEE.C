// some global var/constants
const Int_t kNRCU = 2; // per SM
AliCaloAltroMapping *fMapping[4]; // 1 for each side (A/C) and each RCU (0/1), i.e. 2*2 total
//
const Bool_t kDebug = kFALSE;
const char *branchStr[] = {"A", "B"};

// help methods
Int_t GetHWAddress(Int_t iside, Int_t icol, Int_t irow, Int_t igain);
void DecodeHWAddress(Int_t hwAddr, Int_t & branch, Int_t & FEC, Int_t & chip, Int_t & chan);
void GetMapping();

// main method
void 
TowerToFEE(Int_t iSM=0, Int_t icol=0, Int_t irow=0, Int_t igain=0)
{
  // Get mapping file info, and clear arrays
  GetMapping();

  Int_t isect = iSM / 2; //
  Int_t iside = iSM % 2; // A or C side

  Int_t hwAddress = 0;
  Int_t iRCU = 0;
  Int_t branch = 0;
  Int_t FEC = 0;
  Int_t chip = 0;
  Int_t chan = 0;

  hwAddress = GetHWAddress(iside, icol, irow, igain, iRCU);
  DecodeHWAddress(hwAddress, branch, FEC, chip, chan);  

  Int_t iDDL = iRCU + iSM * kNRCU;

  // report channel info
  printf(" iSM %d icol %d irow %d caloflag (igain) %d \n corresponds to \n iSM %d iRCU %d (iDDL %d EqId %d) : branch %d (%s) FEC %d chip %d chan %d\n", 
	 iSM, icol, irow, igain, 
	 iSM, iRCU, iDDL, iDDL + 0x1200, 
	 branch, branchStr[branch], FEC, chip, chan);

  return;
}

void 
DecodeHWAddress(Int_t hwAddr, Int_t & branch, Int_t & FEC, Int_t & chip, Int_t & chan)
{
  chan = hwAddr & 0xf;
  chip = (hwAddr >> 4) & 0x7;
  FEC = (hwAddr >> 7) & 0xf;
  branch = (hwAddr >> 11) & 0x1;
  return;
}

Int_t 
GetHWAddress(Int_t iside, Int_t icol, Int_t irow, Int_t igain, Int_t & iRCU)
{
  iRCU = -111;

  //RCU0
  if (0<=irow&&irow<8) iRCU=0; // first cable row
  else if (8<=irow&&irow<16 && 0<=icol&&icol<24) iRCU=0; // first half; 
  //second cable row
  //RCU1
  else if(8<=irow&&irow<16 && 24<=icol&&icol<48) iRCU=1; // second half; 
  //second cable row
  else if(16<=irow&&irow<24) iRCU=1; // third cable row

  // swap for odd=C side, to allow us to cable both sides the same
  Int_t iRCUSide = iRCU; 
 if (iside == 1) {
    iRCU = 1 - iRCU;
    iRCUSide = iRCU + 2; // to make it map file index
  }
  Int_t hwAddress = fMapping[iRCUSide]->GetHWAddress(irow, icol, igain);

  return hwAddress;
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

