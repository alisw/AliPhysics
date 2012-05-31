/*
*/

static const int fgkEmCalRows = 24; // number of rows per module for EMCAL
static const int fgkEmCalCols = 48; // number of columns per module for EMCAL

Float_t biasVoltage[fgkEmCalCols][fgkEmCalRows];

//____________________________________________________________________
void WriteBiasV30(const char * inputDBName, const char * inputMapName,
		  const int defaultVoltage, const char * outputFileName) 
{ 

  ofstream outputFile(outputFileName);
  SetBiasVoltage(inputDBName, inputMapName, defaultVoltage);

  for (int icol=0; icol<fgkEmCalCols; icol++) {
    for (int irow=0; irow<fgkEmCalRows; irow++) {
      outputFile << icol << " " << irow << " " << biasVoltage[icol][irow] << endl;
    }
  }

  outputFile.close();
}

//____________________________________________________________________
void SetBiasVoltage(const char * inputDBName, const char * inputMapName,
		    const int defaultVoltage) 
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
  int nNotFound = 0;
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
	i--; // go back to what we found
	biasVoltage[icol][irow] = fCalib[i].fV30;
	nFound++;
      }
      else {
	biasVoltage[icol][irow] = defaultVoltage;
	cout << " APD " << apdMap << " could not be found! " << endl;
	nNotFound++;
      }

    }
  }

  cout << " found " << nFound << " matches " << endl;
  cout << " did not find " << nNotFound << " APDs " << endl;

  return;
}

