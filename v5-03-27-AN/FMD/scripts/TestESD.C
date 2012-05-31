#include <iomanip>
/**
 * Script to test for compatibility to read/write ESD information 
 */

Float_t EtaValue(UShort_t d, Char_t r, UShort_t t)
{
  return (d * 1000 + (r == 'I' || r == 'i' ? 0 : 1) * 100 + 0.001 * t);
}
Float_t MultValue(UShort_t d, Char_t r, UShort_t s, UShort_t t)
{
  return (EtaValue(d, r, t) + s);
}

void 
FillESD(AliESDFMD* esd)
{
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nRng = (d == 1 ? 1 : 2);
    for (UShort_t i = 0; i < nRng; i++) {
      Char_t   r    = (i == 0 ? 'I' : 'O');
      UShort_t nSec = (i == 0 ?  20 :  40);
      UShort_t nStr = (i == 0 ? 512 : 256);
      for (UShort_t s = 0; s < nSec; s++) { 
	for (UShort_t t = 0; t < nStr; t++) {
	  if (s == 0) esd->SetEta(d, r, 0, t, EtaValue(d, r, t));
	  esd->SetMultiplicity(d, r, s, t, MultValue(d, r, s, t));
	}
      } // for s
    } // for i 
  } // for d
}
	  
void 
WriteESD(const char* fileName)
{
  TFile*     file = TFile::Open(fileName, "RECREATE");
  TTree*     tree = new TTree("T", "T");
  AliESDFMD* fmd  = new AliESDFMD();
  tree->Branch("FMD", "AliESDFMD", &fmd);
  
  for (UShort_t i = 0; i < 10; i++) { 
    FillESD(fmd);
    tree->Fill();
  }
  file->Write();
  file->Close();
}

Bool_t 
PrintOne(UShort_t d, Char_t r, UShort_t s, UShort_t t, 
	 Float_t  m, Float_t e)
{
  Float_t em = MultValue(d, r, s, t);
  Float_t ee = EtaValue(d, r, t);
  if (m != em || ee != e) {
    std::cerr << "FMD" << d << r << '[' 
	      << std::setw(2) << s  << ',' 
	      << std::setw(3) << t  << "]: " 
	      << std::setw(8) << m  << " (" 
	      << std::setw(8) << em << ' ' 
	      << (m == em ? "ok" : "bad") << ") @ "
	      << std::setw(8) << e << " ("
	      << std::setw(8) << ee << ' ' 
	      << (ee == e ? "ok" : "bad") << ')' << std::endl;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t
PrintESD(AliESDFMD* esd)
{

  Bool_t ret = kTRUE;
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nRng = (d == 1 ? 1 : 2);
    for (UShort_t i = 0; i < nRng; i++) {
      Char_t   r    = (i == 0 ? 'I' : 'O');
      UShort_t nSec = (i == 0 ?  20 :  40);
      UShort_t nStr = (i == 0 ? 512 : 256);
      for (UShort_t s = 0; s < nSec; s++) { 
	for (UShort_t t = 0; t < nStr; t++) {
	  Float_t m = esd->Multiplicity(d, r, s, t);
	  Float_t e = esd->Eta(d, r, s, t);
	  if (!PrintOne(d, r, s, t, m, e)) ret = kFALSE;
	}
      }
    }
  }
  return ret;
}
  
void
ReadESD(const char* fileName)
{
  TFile*     file = TFile::Open(fileName, "READ");
  TTree*     tree = static_cast<TTree*>(file->Get("T"));
  AliESDFMD* fmd  = new AliESDFMD();
  tree->SetBranchAddress("FMD", &fmd);

  Bool_t ret = kTRUE;
  for (UShort_t i = 0; i < 10; i++) {
    fmd->Clear();
    tree->GetEntry(i);
    if (!PrintESD(fmd)) ret = kFALSE;
  }
  file->Close();
  if (!ret) 
    std::cerr << "There have been errors!" << std::endl;
  else 
    std::cout << "All correct" << std::endl;

}


void
TestESD(const char* fileName=0)
{
  if (!fileName) { 
    WriteESD("esd_test.root");
    return;
  }
  ReadESD(fileName);
}
