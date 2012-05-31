/*
*/
// from AliEMCALGeoParams.h
static const int fgkEMCALRows = 24; // number of rows per module for EMCAL
static const int fgkEMCALCols = 48; // number of columns per module for EMCAL

//____________________________________________________________________
void WriteOCDBCalibMapAPD(const char * inputDBName, const char * inputMapName,
			  const char * outputFileName, const int swapSides) 
{ 

  gSystem->Load("AliEMCALCalibAPD_cxx");
  AliEMCALCalibAPD *calibAPD = new AliEMCALCalibAPD();

  calibAPD->ReadCalibAPDInfo(10000, inputDBName);
  int fNCalibAPD = calibAPD->GetNCalibAPD();
  AliEMCALCalibAPD::AliEMCALCalibAPDData * fCalib = calibAPD->GetCalibAPDData();

  gSystem->Load("AliEMCALMapAPD_cxx");
  AliEMCALMapAPD *mapAPD = new AliEMCALMapAPD();

  // assume we do this for one SuperModule at a time; can merge the 
  // output files with 'cat' or so separately 
  int nSM = 1;
  mapAPD->ReadMapAPDInfo(nSM, inputMapName);
  AliEMCALMapAPD::AliEMCALSuperModuleMapAPD * fMap = mapAPD->GetSuperModuleData();

  // set up ouput file
  ofstream outputFile(outputFileName);

  // let's loop over the files..
  int nFound = 0;
  int nNotFound = 0;
  int iCol = 0; 
  int iRow = 0;
  for (int icol=0; icol<fgkEMCALCols; icol++) {
    for (int irow=0; irow<fgkEMCALRows; irow++) {
      iCol = icol;
      iRow = irow;
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = fgkEMCALCols-1 - iCol;
	iRow = fgkEMCALRows-1 - iRow;
      }

      int apdMap = fMap[0].fAPDNum[icol][irow]; // 0 = nSM - 1
      int i = 0;
      int apdCalib = -1;
      while (i<fNCalibAPD && apdMap!=apdCalib) {
	apdCalib = fCalib[i].fAPDNum;
	i++;
      }

      if (apdCalib == apdMap) { // found!
	i--; // go back to what we found

	// should also calculate the HardWare Id.. note that the functionality
	// of the Tower2FEEMap and GetHardWareId calls have not been tested here..
	
	int ircu, ibranch, card, icsp;
	Tower2FEEMap(icol, irow, &ircu, &ibranch, &card, &icsp);
	int iHW = GetCSPAddress(ibranch, card, icsp);
	iHW |= (ircu << 12); // RCU not part of normal HW addresses, but we add an extra bit for it here, just in case it might be useful

	// print some info..
	outputFile << iCol << " " << iRow << " " << iHW 
		   << " " << fCalib[i].fAPDNum << " " << fCalib[i].fV30 
		   << " " << fCalib[i].fPar[0] << " " << fCalib[i].fPar[1] << " " << fCalib[i].fPar[2]
		   << " " << fCalib[i].fParErr[0] << " " << fCalib[i].fParErr[1] << " " << fCalib[i].fParErr[2]
		   << " " << fCalib[i].fBreakDown << " " << fCalib[i].fDarkCurrent << endl;
	
	nFound++;
      }
      else {
	cout << " APD " << apdMap << " could not be found! " << endl;
	// print some dummy info
	nNotFound++;
      }

    }
  }

  cout << " found " << nFound << " matches " << endl;
  cout << " did not find " << nNotFound << " APDs " << endl;

  // close down
  outputFile.close();
}

// from AliEMCALGeoParams.h:
Int_t GetCSPAddress(Int_t iBranch, Int_t iFEC, Int_t iCSP) const
{ return ( (iBranch<<11) | (iFEC<<7) | iCSP ); }; // 

// from DCSGenerateAPD.C; the method assumes that the columns and rows
// are labelled as we have them on the A side (normal case for calibrations
// before the installation)
const int NRCU = 2; // per SM
const int NBranch = 2; // per RCU
const int NFEC = 9; // per branch, labelled 1..9
const int NCSP = 32; // per FEC

void Tower2FEEMap(const int icol, const int irow,
                  int *ircu, int *ibranch, int *card, int *icsp)
{ /*
    If you are interested in where these magic numbers come from -
    See mapping info on
    http://dsilverm.web.cern.ch/dsilverm/mapping/emcal_mapping.html
    http://dsilverm.web.cern.ch/dsilverm/mapping/ppt/Coordinates_and_Mapping.pdf   */
  
  // each FEC covers a 4x8 tower area
  int C = irow/8; // Cable bundle
  int FEC = C*12 + icol/4; // FEC in 0..35 range
  
  *ircu = FEC / 18; // 18 FEC per RCU
  *ibranch = (FEC%18) / 9;
  *card = FEC % 9;
  
  // columns and rows within an FEC area
  int tCol = icol%4;
  int tRow = irow%8;
  
  // The mapping to CSP is a bit complicated so I also define two more help variables here..
  // which T-card?
  int TCard = tCol/2; // 0=Top (even StripModules), 1=Bottom (odd StripModules)
  int locCol = tCol%2;  // local column inside T-card
  
  *icsp = (7 - tRow) + locCol*16 + TCard*8;
}
