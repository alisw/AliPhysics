/**************************************************************************
 * PMD Pdestal Reader, pedestal is kept in a folded number with two distint
 * values, mean and rms. This macro reads and unfold the mean and RMS. You 
 *          can further modify it to check and plot the values.  
 *            Auther:    Satyajit Jena |  sjena@cern.ch
 *                Mon Nov 22 19:54:27 CET 2010
 *                     
 * In case it doesn't work, please let me know. 
 **************************************************************************/


const Int_t kDet = 2;
const Int_t kMod = 24;
const Int_t kRow = 48;
const Int_t kCol = 96;

void ReadPmdPedestal(const Char_t* filename="GiveTheOcdbFileName.root") {
  
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
  
  AliPMDPedestal *fPedObj = (AliPMDPedestal*)fCdbEntrObj->GetObject();
  if (!fPedObj) {
    AliError("AliPMDPedestal","No Pedestal Object Present");
    return;
  }
    
  // Read Pdestal  Entries
  
  for (int i = 0; i < kDet; i++) {
    for (int j = 0; j < kMod; j++) {
      for (int k = 0; k < kRow; k++) {
        for (int l = 0; l < kCol; l++) {
	  
	  Int_t    folded  = fPedObj->GetPedMeanRms(i,j,k,l);
	  Int_t    temprms = Int_t(folded%100);
	  Double_t rms     = Double_t(temprms/10.);
	  Double_t mean    = Double_t((folded - temprms)/1000.0);
	  
	  Printf("%4d %4d %4d %4d %10.5f %10.5f", i,j,k,l,mean,rms);
	}
      }
    }
  }
  //


}
