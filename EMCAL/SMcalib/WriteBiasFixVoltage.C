/*
*/

static const int fgkEmCalRows = 24; // number of rows per module for EMCAL
static const int fgkEmCalCols = 48; // number of columns per module for EMCAL

//____________________________________________________________________
void WriteBiasFixVoltage(const int biasVoltage, const char * txtFileName) 
{ 

  ofstream outputFile(txtFileName);

  for (int icol=0; icol<fgkEmCalCols; icol++) {
    for (int irow=0; irow<fgkEmCalRows; irow++) {
      outputFile << icol << " " << irow << " " << biasVoltage << endl;
    }
  }

  outputFile.close();
}
