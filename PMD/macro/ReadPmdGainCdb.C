/**************************************************************************
 * PMD Gain OCDB Reader, this is a template to be used for the reading of
 * Gain ocdb entries. The custom ocdb to keep chain flag for ESD analysis
 * can use this macro. This is written for aditya. 
 *           
 * Auther: Satyajit Jena <sjena@cern.ch>
 * Wed Dec 17 16:08:38 CET 2014
 * In case it doesn't work, please let me know. 
 **************************************************************************/


const Int_t kDet = 2;
const Int_t kMod = 24;
const Int_t kRow = 48;
const Int_t kCol = 96;

void ReadPmdGainCdb(const Char_t* filename="GiveTheOcdbFileName.root") {
  
  TFile *file = new TFile(filename);
  
  if (!file) {
    AliError("TFile","No file presnt - please provide correct file");
    return;
  }
  
  AliCDBEntry *fCdbEntrObj = (AliCDBEntry*)file->Get("AliCDBEntry");
  if (!fCdbEntrObj) {
    AliError("AliCDBEntry","No CDB Entry");
    return;
  }
  
  AliPMDCalibData *fGainObj = (AliPMDCalibData*)fCdbEntrObj->GetObject();
  if (!fGainObj) {
    AliError("AliPMDPedestal","No Pedestal Object Present");
    return;
  }
    
  // Read Pdestal  Entries
  
  for (int i = 0; i < kDet; i++) {
    for (int j = 0; j < kMod; j++) {
      for (int k = 0; k < kRow; k++) {
        for (int l = 0; l < kCol; l++) {
	  Float_t    gain  = fGainObj->GetGainFact(i,j,k,l);
	  Printf("%4d %4d %4d %4d %10.5f", i,j,k,l,gain);
	}
      }
    }
  }

}
