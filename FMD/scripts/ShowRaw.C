//
// Script to read a raw data file, and dump it to std::cout 
//
#include <iomanip>

void 
ShowRaw(Int_t det=2,  bool verbose=false, Int_t event=0) 
{
  TString file(Form("raw%d/FMD_%d.ddl", event, AliFMD::kBaseDDL + det - 1));

  std::cout << "Reading raw data file " << file << std::endl;
  
  TH1* h = new TH1F("rawData", "Raw Data", 90, 0, 90);
  
  
  // This method creates a text file containing the same information
  // stored in an Altro file. The information in the text file is
  // organized pad by pad and and for each pad it consists in a
  // sequence of bunches (Bunch length +2, Time bin of the last
  // amplitude sample in the bunch, amplitude values) It is used
  // mainly //for debugging

  AliAltroBuffer buff(file.Data(),0);
  Int_t numWords,padNum,rowNum,secNum=0;
  Int_t value = 0;
  Int_t zero  = 0;
  // if (!buff.ReadDataHeader()) {
  // std::cout<< file << " isn't a valid data file!" << std::endl;
  // }
  
  while(buff.ReadTrailerBackward(numWords,padNum,rowNum,secNum)){
    if (verbose) 
      std::cout << "Ring: " << (secNum == 0 ? 'I' : 'O') 
		<< " Sector: " << std::setw(2) << rowNum 
		<< " Strip:  " << std::setw(3) << padNum 
		<< " Words:  " << std::setw(4) << numWords << std::endl;
    if (numWords == 0) zero++;
    if (numWords % 4){
      if (verbose) 
	std::cout << "Skipping trailer of " 
		  << (4 - numWords % 4) << " words" << std::endl;
      for(Int_t j = 0; j < (4 - numWords % 4); j++)
	value=buff.GetNextBackWord(); 
    }//end if
    for(Int_t i = 0; i <numWords; i++) {
      value=buff.GetNextBackWord();
      if (verbose) {
	std::cout << std::setw(5) <<  value << std::flush;
	if (i % 16 == 15) std::cout << std::endl;
      }
      h->Fill(value);
    }//end for
    if (verbose)
      std::cout << std::endl;
    if (zero > 1) {
      std::cout << "Error: Read zero channels - should not happen" 
		<< std::endl;
      break;
    }
  }//end while
  h->Draw();
  return;
}
