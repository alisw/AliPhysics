#include <iomanip>

void
MakeFakeDigits()
{
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(0);
  AliFMDParameters* param = AliFMDParameters::Instance();
  param->Init(kFALSE,(AliFMDParameters::kPulseGain|
		      AliFMDParameters::kPedestal|
		      AliFMDParameters::kDeadMap|
		      AliFMDParameters::kSampleRate|
		      AliFMDParameters::kAltroMap|
		      AliFMDParameters::kStripRange));
  // param->SetPedestal(100);
  
  TFile*      file = TFile::Open("FMD.Digits.root", "RECREATE");
  TDirectory* dir  = file->mkdir("Event0");
  
  TTree*        tree   = new TTree("TreeD", "Digits container");
  TClonesArray* digits = new TClonesArray("AliFMDDigit");
  tree->Branch("FMD", &digits);
  
  Int_t i = 0;
  
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nrng = (d == 1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nrng; ir++) { 
      Char_t   r    = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ?  20 :  40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      for (UShort_t s = 0; s < nsec; s++) { 
	for (UShort_t t = 0; t < nstr; t++) { 
	  std::cout << "FMD" << d << r << "[" << std::setfill('0') 
		    << std::setw(2) << s << "," << std::setw(3) << t 
		    << "] " <<  std::setfill(' ') << " " 
		    << std::setw(5) << i << "/51200 " 
		    << std::setw(3) << int(i*100./51200) << "%\r"
		    << std::flush;

	  UShort_t nsam = param->GetSampleRate(d,r,s,t);
	  AliFMDDigit* digit = new((*(digits))[i]) AliFMDDigit(d,r,s,t);
	  i++;
	  UShort_t c = param->GetPedestal(d,r,s,t);
	  for (UShort_t z = 0; z < nsam; z++) { 
	    digit->SetCount(z, c);
	    if (nsam > 1 && z < 1) continue;
	    if ((t > 32 && t < 64) || (t > 96 && t < 98)) 
	      digit->SetCount(z, c + 100);
	  }
	}
      }
    }
  }
  std::cout << "Done making fake digits" << std::endl;
  dir->cd();
  tree->Fill();
  tree->Write();
  file->Close();
}

	    
