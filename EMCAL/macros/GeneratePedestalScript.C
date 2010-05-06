// some global var/constants
const Int_t kNSM = 4; // for first LHC run
const Int_t kNRCU = 2;
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
void PrintScript();
void Clear();
Int_t GetHWAddress(Int_t iside, Int_t icol, Int_t irow, Int_t igain);
Int_t GetHWAddressLEDRef(Int_t istrip, Int_t igain);
void DecodeHWAddress(Int_t hwAddr, Int_t & branch, Int_t & FEC, Int_t & chip, Int_t & chan);
void GetMapping();
void CreateMappingLEDRef();

// main method
void 
//GeneratePedestalScript(const char * filename = "alien/Run113790_113790_v1_s0.root") // 1st set
GeneratePedestalScript(const char * filename = "Run117756_117756_v1_s0.root") // 2nd set
{
  // get the DA info/object
  TFile *file = TFile::Open(filename); 
  //TMP  AliCaloCalibPedestal *emcCalibPedestal = AliCDBEntry->GetObject();
  if (kDebug) {
    file->ls();
  }

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
  PrintScript();

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
PrintScript()
{
  const char * sideStr[] = {"A","C"};
  const char * branchStr[] = {"A","B"};
  int VFPED = 0x06;
  int RCUWrite = 0x200000;

  char filename[100];
  char scriptLine[200];

  Int_t iSM = 0;
  Int_t iRCU = 0;
  Int_t ibranch = 0;
  Int_t iFEC = 0;
  Int_t ichip = 0;
  Int_t ichan = 0;
  Int_t Ped = 0;

  for (iSM=0; iSM<kNSM; iSM++) {
    int iside = iSM % 2;
    int isect = iSM / 2;
    for (iRCU=0; iRCU<kNRCU; iRCU++) {
      sprintf(filename, "setSM%1s%dRCU%d.scr", 
	      sideStr[iside], isect, iRCU);
      ofstream fout(filename);
      int nscriptLines = 0;

      for (ibranch=0; ibranch<kNBranch; ibranch++) {
	int firstFEC = 1; 
	if (ibranch==0 && iRCU==0) { // extra LED Ref FEE in this location
	  firstFEC = 0;
	}

	for (iFEC=firstFEC; iFEC<kNFEC; iFEC++) { // FEC 1..9 (or 0..9)
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
		  int writeAddr = (ibranch << 16) | (iFEC << 12) | (ichip << 9) |
		    (ichan << 5) | VFPED | RCUWrite;  
		  sprintf(scriptLine, "w 0x%04x 0x%06x  # Branch %s, Card %d, Altro %d, Chan %d", 
			  nscriptLines, writeAddr, branchStr[ibranch],
			  iFEC, ichip, ichan);
		  fout << scriptLine << endl;
		  nscriptLines++;
		  
		  int writeVal = (Ped | RCUWrite);
		  sprintf(scriptLine, "w 0x%04x 0x%06x  # Pedestal 0x%x = %d", 
			  nscriptLines, writeVal, Ped, Ped);
		  fout << scriptLine << endl;
		  nscriptLines++;
		  
		}
	      }
	    }
	  }
	}
      }
      // ending, with execute and update step..	    
      sprintf(scriptLine, "w 0x%04x 0x380000  # End of the sequence", 
	      nscriptLines);
      fout << scriptLine << endl;
      nscriptLines++;

      sprintf(scriptLine, "w 0x%04x 0x3F0000  # End of the instruction memory", 
	      nscriptLines);
      fout << scriptLine << endl;
      nscriptLines++;

      fout << "wait 100 us" << endl;
      fout << "w 0x5304 0x0       \# execute and update registers" << endl; 
      fout.close();
    } // iRCU
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
