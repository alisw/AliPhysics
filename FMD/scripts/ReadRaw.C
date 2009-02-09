/** @file    ReadRaw.C
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Tue Mar 28 12:39:08 2006
    @brief   Script to read raw data 
*/
/** @ingroup FMD_script
    @brief Read raw data into a TClonesArray - for testing 
 */
void
ReadRaw(const char* file=0, Int_t evno=0, bool old=false)
{
  TString        src(file ? file : "");
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetRun(0);
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliLog::SetModuleDebugLevel("FMD", 1);
  AliFMDParameters::Instance()->Init();
  AliFMDParameters::Instance()->UseRcuTrailer(!old);
  AliFMDParameters::Instance()->UseCompleteHeader(!old);
  AliRawReader*                   r = 0;
  if (src.IsNull())               {
    std::cout << "Reading via AliRawReaderFile" << std::endl;
    r = new AliRawReaderFile(0);
  }
  else if (src.EndsWith(".root")) { 
    std::cout << "Reading via AliRawReaderRoot" << std::endl;
    r = new AliRawReaderRoot(src.Data(), evno);
  }
  else if (src.EndsWith(".raw"))  {
    std::cout << "Reading via AliRawReaderDate" << std::endl;
    r = new AliRawReaderDate(src.Data());
  }
  else {
    std::cerr << "Unknown raw type for source " << src 
	      << " assuming simulated raw files in directory " << src 
	      << std::endl;
    r = new AliRawReaderFile(src);
  }
  AliFMDRawReader* fr = new AliFMDRawReader(r, 0);
  TClonesArray*    a  = new TClonesArray("AliFMDDigit", 0);
  fr->ReadAdcs(a);

  std::cout << "Read " << a->GetEntriesFast() << " digits" << std::endl;
  
  bool read[3][2][40][512];
  for (UShort_t det = 0; det < 3; det++) {
    for (UShort_t rng = 0; rng < 2; rng++) { 
      for (UShort_t sec = 0; sec < 40; sec++) {
	for (UShort_t str = 0; str < 512; str++) { 
	  read[det][rng][sec][str] = false;
	}
      }
    }
  }
  

  TIter next(a);
  AliFMDDigit* d = 0;
  while ((d = static_cast<AliFMDDigit*>(next()))) {
    UShort_t det = d->Detector() - 1;
    UShort_t rng = d->Ring() == 'I' ? 0 : 1;
    UShort_t sec = d->Sector();
    UShort_t str = d->Strip();
    read[det][rng][sec][str] = true;
  }
  const UShort_t lineLength = 64;
  for (UShort_t det = 0; det < 3; det++) {
    for (UShort_t rng = 0; rng < 2; rng++) { 
      if (det == 0 && rng == 1) continue;
      Char_t   rid  = rng == 0 ? 'I' : 'O';
      UShort_t nsec = rng == 0 ?  20 :  40;
      UShort_t nstr = rng == 0 ? 512 : 256;
      std::cout << "FMD" << det+1 << rid << ":" << std::endl;
      for (UShort_t sec = 0; sec < nsec; sec++) {
	std::cout << " Sector " << sec << "\n" << std::flush;
	for (UShort_t str = 0; str < nstr; str++) { 
	  bool on = read[det][rng][sec][str];
	  if (str % lineLength == 0) std::cout << "  ";
	  std::cout << (on ? '+' : '-');
	  if (str % lineLength == lineLength-1) std::cout << "\n";
	}
      }
    }
  }
  // a->ls();
}
//____________________________________________________________________
//
// EOF
//
