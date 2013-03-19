// some global var/constants
const Int_t kNSM = 12; // 
const Int_t kNRCU = 2;
const Int_t kNDTC = 40; // links for full SRU (corresponding to two readout crates or RCUs)
AliCaloAltroMapping *fMapping[4]; // 1 for each side (A/C) and each RCU (0/1), i.e. 2*2 total
const Int_t kNBranch = 2;
const Int_t kNFEC = 10; // 0..9, when including LED Ref
const Int_t kNChip = 5; // really 0,2..4, i.e. skip #1
const Int_t kNChan = 16;
Float_t fMeanPed[kNSM][kNRCU][kNBranch][kNFEC][kNChip][kNChan];
Float_t fRmsPed[kNSM][kNRCU][kNBranch][kNFEC][kNChip][kNChan];
//
const int kNStrips = 24; // per SM
Int_t fHWAddrLEDRef[kNStrips][2]; // [2] is for Low/High gain

const Bool_t kDebug = kFALSE;
const Float_t kBadRMS = 20;

// help methods
void GetPedVal(const Int_t iSM, const Int_t igain, const TProfile2D *h2);
void GetPedValLEDRef(const Int_t iSM, const Int_t igain, const TProfile *h);
void PrintScript(const Int_t runno);
void Clear();
Int_t GetHWAddress(Int_t iside, Int_t icol, Int_t irow, Int_t igain);
Int_t GetHWAddressLEDRef(Int_t istrip, Int_t igain);
void DecodeHWAddress(Int_t hwAddr, Int_t & branch, Int_t & FEC, Int_t & chip, Int_t & chan);
void GetMapping();
void CreateMappingLEDRef();

// main method
void 
GeneratePedestalScriptSRU(const char * filename = "EMCALPED.root") // from DA node
{
  // get the DA info/object
  TFile *file = TFile::Open(filename); 
  //TMP  AliCaloCalibPedestal *emcCalibPedestal = AliCDBEntry->GetObject();
  if (kDebug) {
    file->ls();
  }
  int runno = emcCalibPedestal->GetRunNumber();

  // Get mapping file info, and clear arrays
  Clear();
  GetMapping();
  CreateMappingLEDRef();

  // Store the pedestal info
  for (Int_t iSM=0; iSM<kNSM; iSM++) {
    GetPedVal( iSM, 0, emcCalibPedestal->GetPedProfileLowGain(iSM) );
    GetPedVal( iSM, 1, emcCalibPedestal->GetPedProfileHighGain(iSM) );
    GetPedValLEDRef( iSM, 0, emcCalibPedestal->GetPedLEDRefProfileLowGain(iSM) );
    GetPedValLEDRef( iSM, 1, emcCalibPedestal->GetPedLEDRefProfileHighGain(iSM) );
  }

  // Generate the needed scripts
  PrintScript(runno);

}

void 
GetPedVal(const Int_t iSM, const Int_t igain, const TProfile2D *h2)
{
  Int_t isect = iSM / 2; //
  Int_t iside = iSM % 2; // A or C side
  Int_t nCols = h2->GetNbinsX();
  Int_t nRows = h2->GetNbinsY();
  if (kDebug) {
    printf("GetPedVal: iSM %d isect %d iside %d igain %d nRows %d nCols %d\n",
	   iSM, isect, iside, igain, nRows, nCols);
  }
  Int_t hwAddress = 0;
  Int_t iRCU = 0;
  Int_t branch = 0;
  Int_t FEC = 0;
  Int_t chip = 0;
  Int_t chan = 0;

  Int_t icol = 0;
  Int_t irow = 0;
  Int_t bin = 0;

  for (icol=0; icol<nCols; icol++) {
    for (irow=0; irow<nRows; irow++) {

      hwAddress = GetHWAddress(iside, icol, irow, igain, iRCU);
      DecodeHWAddress(hwAddress, branch, FEC, chip, chan);  
      bin = h2->FindBin(icol, irow);
            
      // store the values
      fMeanPed[iSM][iRCU][branch][FEC][chip][chan] = h2->GetBinContent(bin);
      fRmsPed[iSM][iRCU][branch][FEC][chip][chan] = h2->GetBinError(bin);

      // report bad RMS channels:
      if (h2->GetBinError(bin) > kBadRMS) {
	printf(" bad pedestal RMS: iSM %d icol %d irow %d igain %d iRCU %d branch %d FEC %d chip %d chan %d - mean %4.1f rms %4.1f\n", 
	       iSM, icol, irow, igain, iRCU, branch, FEC, chip, chan,
	       h2->GetBinContent(bin), h2->GetBinError(bin));
      }
    }
  }

  return;
}


void 
GetPedValLEDRef(const Int_t iSM, const Int_t igain, const TProfile *h)
{
  Int_t isect = iSM / 2; //
  Int_t iside = iSM % 2; // A or C side
  Int_t nStrips = h->GetNbinsX();

  if (kDebug) {
    printf("GetPedValLEDRef: iSM %d isect %d iside %d igain %d nStrips %d\n",
	   iSM, isect, iside, igain, nStrips);
  }
  Int_t hwAddress = 0;
  Int_t iRCU = 0; // always true for LED Ref FEE
  Int_t branch = 0;
  Int_t FEC = 0;
  Int_t chip = 0;
  Int_t chan = 0;

  Int_t icol = 0;
  Int_t irow = 0;
  Int_t bin = 0;

  for (int istrip=0; istrip<nStrips; istrip++) {

    hwAddress = GetHWAddressLEDRef(istrip, igain);
    DecodeHWAddress(hwAddress, branch, FEC, chip, chan);  
    bin = h->FindBin(istrip);
            
    // store the values
    fMeanPed[iSM][iRCU][branch][FEC][chip][chan] = h->GetBinContent(bin);
    fRmsPed[iSM][iRCU][branch][FEC][chip][chan] = h->GetBinError(bin);
  }

  return;
}

void 
PrintScript(const int runno)
{
  const char * sideStr[] = {"A","C"};
  const char * branchStr[] = {"A","B"};
  int VFPED = 0x06;
  int SRUAltroWrite = 0x40000000; 
  int Port = 0x1001; // 4097

  char dirname[40];
  sprintf(dirname,"scriptsSRU_Run%d", runno);

  char cmd[500];
  sprintf(cmd,"mkdir -p %s/SM{A,C}{0,1,2,3,4,5}", dirname);
  gSystem->Exec(cmd);

  sprintf(cmd,"echo \"cp info.txt GeneratePedestalScriptSRU.C %s/.\" > scp.sh", dirname);
  gSystem->Exec(cmd);
  sprintf(cmd,"echo \"scp -r -p %s aldaqacr50:pedestals/.\" >> scp.sh", dirname);
  gSystem->Exec(cmd);
  sprintf(cmd,"echo \"ssh -t aldaqacr50 \'cd pedestals; scp -r -p %s emc@alidcscom702:srudcs/scripts/pedestals/.\'\" >> scp.sh", dirname);
  gSystem->Exec(cmd);

  char filename[100];
  char scriptLine[200];

  Int_t iSM = 0;
  Int_t iRCU = 0;
  Int_t ibranch = 0;
  Int_t iFEC = 0;
  Int_t ichip = 0;
  Int_t ichan = 0;
  Int_t Ped = 0;
  Int_t iDTC = 0;

  for (iSM=0; iSM<kNSM; iSM++) {
    int iside = iSM % 2;
    int isect = iSM / 2;

    char IP[100];
    if (iSM == 0) { sprintf(IP, "10.160.132.100"); } // SMA0 
    else if (iSM == 1) { sprintf(IP, "10.160.132.102"); } // SMC0
    else if (iSM == 2) { sprintf(IP, "10.160.132.104"); } // SMA1 
    else if (iSM == 3) { sprintf(IP, "10.160.132.106"); } // SMC1
    else if (iSM == 4) { sprintf(IP, "10.160.132.108"); } // SMA2 
    else if (iSM == 5) { sprintf(IP, "10.160.132.110"); } // SMC2
    else if (iSM == 6) { sprintf(IP, "10.160.132.112"); } // SMA3 
    else if (iSM == 7) { sprintf(IP, "10.160.132.114"); } // SMC3
    else if (iSM == 8) { sprintf(IP, "10.160.132.116"); } // SMA4 
    else if (iSM == 9) { sprintf(IP, "10.160.132.118"); } // SMC4
    else if (iSM == 10) { sprintf(IP, "10.160.36.155"); } // SMA5
    else if (iSM == 11) { sprintf(IP, "10.160.36.156"); } // SMC5

    // only do instrumented parts..
    int activeDTC[kNDTC] = {0};
    for (iDTC=0; iDTC<kNDTC; iDTC++) {
      if (iDTC==10 || iDTC==20 || iDTC==30) { // skip TRU 
	activeDTC[iDTC] = 0;
      } 
      else {
	if (iSM<10) { // not special third SM
	  activeDTC[iDTC] = 1;
	}
	else {
	  if (iSM==10) { // SMA5
	    if (iDTC<14) { activeDTC[iDTC] = 1; }
	    else { activeDTC[iDTC] = 0; }
	  }
	  else if (iSM==11) { // SMC5
	    if (iDTC==0 || iDTC>=27) { activeDTC[iDTC] = 1; }
	    else { activeDTC[iDTC] = 0; }
	  }

	}
      }
    }

    // OK, let's generate the files for all active FECs/DTCs
    for (iDTC=0; iDTC<kNDTC; iDTC++) {
      if (activeDTC[iDTC] == 0) { continue; }
      sprintf(filename, "%s/SM%1s%d/set_ped_DTC%02d.txt", 
	      dirname, sideStr[iside], isect, iDTC);
      ofstream fout(filename);

      iRCU = iDTC / 20;
      ibranch = (iDTC % 20) / 10;
      iFEC = iDTC % 10;
      int ipos = iFEC + 10*ibranch;

      // write DTC file header..
      sprintf(scriptLine, "%s # IP\n%d       #UDP port", IP, Port);
      fout << scriptLine << endl;

      int dtcselUpper = 0;
      int dtcselLower = 0;
      if (iRCU == 0) {
        dtcselLower = (1 << ipos);
      }
      else { // crate == 1
        dtcselUpper = (1 << ipos);
      }
      sprintf(scriptLine, "%08x # DTC SEL Upper\n%08x # DTC SEL Lower",
	      dtcselUpper, dtcselLower);
      fout << scriptLine << endl;
      sprintf(scriptLine, "%08x # FEE GTL Address\n%08x # Branch %s, Card %d",
	      0x3, iFEC + 16*ibranch, branchStr[ibranch], iFEC);
      fout << scriptLine << endl;
      // end header

      for (ichip=0; ichip<kNChip; ichip++) { // ALTRO 0,2,3,4
	if (ichip!=1) {
	  for (ichan=0; ichan<kNChan; ichan++) {

	    if (iFEC!=0 || (ichan<8 || ichan>11)) {

	      Ped = TMath::Nint(fMeanPed[iSM][iRCU][ibranch][iFEC][ichip][ichan]);
	      // raise Ped value to max for channels with exceptionally large RMS
	      if (fRmsPed[iSM][iRCU][ibranch][iFEC][ichip][ichan] > kBadRMS) {
		printf(" bad pedestal RMS: iSM %d iRCU %d ibranch %d iFEC %d ichip %d ichan %d - raising from %d to 0x3ff\n", 
		       iSM, iRCU, ibranch, iFEC, ichip, ichan, Ped);
		Ped = 0x3ff;
	      }
	      // 
	      int writeAddr = SRUAltroWrite | (ibranch << 16) | (iFEC << 12) | (ichip << 9) | (ichan << 5) | VFPED; 
	      sprintf(scriptLine, "%08x # Branch %s, Card %d, Altro %d, Chan %d", 
		      writeAddr, branchStr[ibranch], iFEC, ichip, ichan);
	      fout << scriptLine << endl;
	      
	      int writeVal = Ped; 
	      sprintf(scriptLine, "%08x # Pedestal 0x%x = %d", 
		       writeVal, Ped, Ped);
	      fout << scriptLine << endl;
	      
	    }
	  }
	}
      } // chip

      fout.close();
    } // iDTC
  }// iSM

  return;
}

void 
Clear()
{
  for (Int_t iSM=0; iSM<kNSM; iSM++) {
    for (Int_t iRCU=0; iRCU<kNRCU; iRCU++) {
      for (Int_t ibranch=0; ibranch<kNBranch; ibranch++) {
	for (Int_t iFEC=0; iFEC<kNFEC; iFEC++) {
	  for (Int_t ichip=0; ichip<kNChip; ichip++) {
	    for (Int_t ichan=0; ichan<kNChan; ichan++) {
	      fMeanPed[iSM][iRCU][ibranch][iFEC][ichip][ichan] = 0;
	      fRmsPed[iSM][iRCU][ibranch][iFEC][ichip][ichan] = 0;
	    }
	  }
	}
      }
    }
  }

  for (int istrip=0; istrip<kNStrips; istrip++) {
    fHWAddrLEDRef[istrip][0] = 0; 
    fHWAddrLEDRef[istrip][1] = 0; 
  }

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

Int_t 
GetHWAddressLEDRef(Int_t istrip, Int_t igain)
{
  Int_t iRCU = 0; // for both sides; LED ref info is the same for both sides
  Int_t caloflag = 3; // AliCaloRawStreamV3::kLEDMonData;
 
  Int_t hwAddress = fHWAddrLEDRef[istrip][igain];

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

void 
CreateMappingLEDRef()
{ 
  Int_t iRCU = 0; // for both sides; LED ref info is the same for both sides
  Int_t caloflag = 3; // AliCaloRawStreamV3::kLEDMonData;

  Int_t maxAddr = 1 << 7; // LED Ref FEE is in FEC pos 0, i.e. addr space 0..127

  int nLEDRefFEEChan = 0;

  Int_t branch = 0;
  Int_t FEC = 0;
  Int_t chip = 0;
  Int_t chan = 0;
  for (int hwaddr = 0; hwaddr<maxAddr; hwaddr++) {

    DecodeHWAddress(hwaddr, branch, FEC, chip, chan);  
    if ( (chip!=1 && chip<kNChip) &&  // ALTROs 0,2,3,4
	 (chan<8 || chan>11) ) { // actual installed LED Ref FEE channels

      int istrip = fMapping[iRCU]->GetPad(hwaddr);
      int igain = fMapping[iRCU]->GetPadRow(hwaddr);
      int iflag = fMapping[iRCU]->GetSector(hwaddr);
      if (iflag == caloflag) {
	fHWAddrLEDRef[istrip][igain] = hwaddr; 
	nLEDRefFEEChan++;
      }
    }
  }

  if (kDebug) { cout << " nLEDRefFEEChan " << nLEDRefFEEChan << endl; }
}
