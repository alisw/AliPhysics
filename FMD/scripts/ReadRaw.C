/** @file    ReadRaw.C
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Tue Mar 28 12:39:08 2006
    @brief   Script to read raw data 
*/
/** @ingroup FMD_script
    @brief Read raw data into a TClonesArray - for testing 
 */
void
ReadRaw(const char* src=0, Int_t nEv=0, Int_t skip=0)
{
  AliLog::SetModuleDebugLevel("FMD", 10);

  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetRun(0);
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  AliRawReader*    reader    = AliRawReader::Create(src);
  AliFMDRawReader* fmdReader = new AliFMDRawReader(reader, 0);
  TClonesArray*    array     = new TClonesArray("AliFMDDigit", 0);

  Int_t evCnt = 0;
  while (reader->NextEvent()) {
    if (nEv > 0 && (evCnt-skip) > nEv) break;
    evCnt++;
    array->Clear();
    fmdReader->ReadAdcs(array);

    std::cout << "Event # " << evCnt << std::endl;

    std::cout << "Read " << array->GetEntriesFast() << " digits" << std::endl;
  
#if 0
    AliFMDBoolMap read(0);
    read.Reset(kFALSE);

    TIter next(array);
    AliFMDDigit* digit = 0;
    while ((digit = static_cast<AliFMDDigit*>(next()))) {
      UShort_t d = digit->Detector();
      Char_t   r = digit->Ring();
      UShort_t s = digit->Sector();
      UShort_t t = digit->Strip();
      read(d,r,s,t) = true;
    }
    const UShort_t lineLength = 64;
    for (UShort_t d = 1; d <= 3; d++) {
      UShort_t nr = (d == 1 ? 1 : 2);
      for (UShort_t q = 0; q < nr; q++) { 
	Char_t   r  = q == 0 ? 'I' : 'O';
	UShort_t ns = q == 0 ?  20 :  40;
	UShort_t nt = q == 0 ? 512 : 256;
	std::cout << "FMD" << d << r << ":" << std::endl;
	for (UShort_t s = 0; s < ns; s++) {
	  std::cout << " Sector " << s << "\n" << std::flush;
	  for (UShort_t t = 0; t < nt; t++) { 
	    bool on = read(d,r,s,t);
	    if (t % lineLength == 0) std::cout << "  ";
	    std::cout << (on ? '+' : '-');
	    if (t % lineLength == lineLength-1) std::cout << "\n";
	  }
	}
      }
    }
#endif
  }
  // a->ls();
}
//____________________________________________________________________
//
// EOF
//
