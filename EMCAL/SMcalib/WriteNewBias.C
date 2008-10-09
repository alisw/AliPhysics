/*
*/

const int kMaxHV = 395;

static const int fgkEmCalRows = 24; // number of rows per module for EMCAL
static const int fgkEmCalCols = 48; // number of columns per module for EMCAL

int previousVoltage[fgkEmCalCols][fgkEmCalRows];
Float_t gainFactor[fgkEmCalCols][fgkEmCalRows];
int biasVoltage[fgkEmCalCols][fgkEmCalRows];

//____________________________________________________________________
void WriteNewBias(const char * inputDBName, const char * inputMapName,
		  const char * previousVoltageFileName, const char * gainFactorFileName, 
		  const char * outputFileName) 
{ 

  ofstream outputFile(outputFileName);
  ReadFiles(previousVoltageFileName, gainFactorFileName);

  SetBiasVoltage(inputDBName, inputMapName);

  for (int icol=0; icol<fgkEmCalCols; icol++) {
    for (int irow=0; irow<fgkEmCalRows; irow++) {
      outputFile << icol << " " << irow << " " << biasVoltage[icol][irow] << endl;
    }
  }

  outputFile.close();
}

//____________________________________________________________________
void ReadFiles(const char * previousVoltageFileName, 
	       const char * gainFactorFileName)
{
  ifstream previousVoltageFile(previousVoltageFileName);
  ifstream gainFactorFile(gainFactorFileName);

  int icolp, irowp;
  int icolg, irowg;
  Float_t gFactor;
  int pVoltage;

  for (int icol=0; icol<fgkEmCalCols; icol++) {
    for (int irow=0; irow<fgkEmCalRows; irow++) {
      previousVoltageFile >> icolp >> irowp >> pVoltage;
      gainFactorFile >> icolg >> irowg >> gFactor;

      previousVoltage[icolp][irowp] = pVoltage;
      gainFactor[icolg][irowg] = gFactor;
    }
  }

  previousVoltageFile.close();
  gainFactorFile.close();

  return;
}

//____________________________________________________________________
void SetBiasVoltage(const char * inputDBName, const char * inputMapName) 
{
  gSystem->Load("AliEMCALCalibAPD_cxx");
  AliEMCALCalibAPD *calibAPD = new AliEMCALCalibAPD();

  calibAPD->ReadCalibAPDInfo(10000, inputDBName);
  int fNCalibAPD = calibAPD->GetNCalibAPD();
  AliEMCALCalibAPD::AliEMCALCalibAPDData * fCalib = calibAPD->GetCalibAPDData();


  gSystem->Load("AliEMCALMapAPD_cxx");
  AliEMCALMapAPD *mapAPD = new AliEMCALMapAPD();

  int nSM = 1;
  mapAPD->ReadMapAPDInfo(nSM, inputMapName);
  AliEMCALMapAPD::AliEMCALSuperModuleMapAPD * fMap = mapAPD->GetSuperModuleData();

  int nFound = 0;
  for (int icol=0; icol<fgkEmCalCols; icol++) {
    for (int irow=0; irow<fgkEmCalRows; irow++) {

      int apdMap = fMap[0].fAPDNum[icol][irow]; // 0 = nSM - 1
      int i = 0;
      int apdCalib = -1;
      while (i<fNCalibAPD && apdMap!=apdCalib) {
	apdCalib = fCalib[i].fAPDNum;
	i++;
      }

      if (apdCalib == apdMap) { // found!
	i--; // go back to what we dound

	// estimate what the new/target HV should be
	biasVoltage[icol][irow] = CalculateTargetHV(previousVoltage[icol][irow], 
						    gainFactor[icol][irow], 
						    fCalib[i].fPar[0], 
						    fCalib[i].fPar[1], 
						    fCalib[i].fPar[2]);
	nFound++;
      }
      else { // no calib info, just use old settings
	biasVoltage[icol][irow] = previousVoltage[icol][irow];
      }

    }
  }

  cout << " found " << nFound << " matches " << endl;
  return;
}

//____________________________________________________________________
int CalculateTargetHV(int initHV, Float_t gainChange,		       
		      Float_t par0, Float_t par1, Float_t par2 )
{
  printf("parameters p0:%g p1:%g p2:%g\n", par0, par1, par2);

  // figure out what new HV should be, 
  // if we want to adjust the gain by some factor
  Float_t initialHV = (Float_t) initHV;
  Float_t initialGain = par0 + par1 * exp(par2*initialHV);
  Float_t newGain = initialGain * gainChange; // = par0 + par1 * exp(par2*newHV);

  printf("initialGain:%g newGain:%g\n", initialGain, newGain);

  Float_t fNewHV = -1;
  if ( par1>0 && par2>0 ) {
    fNewHV = log ( (newGain - par0)/par1 ) / par2;
  }
  // round-off
  int newHV = (int) (fNewHV + 0.5);
  
  // check results before returning..
  if (newHV == -1) {
   // conversion failed:  let's just keep the old custom value then
    newHV = initHV;
  }
  if (newHV>kMaxHV) {
    // we reached a too high voltage: let's keep the max then
   newHV = kMaxHV;
  }

  return newHV;
}
