//
// Test of the "NextSignal" interface of the AliFMDRawReader
//
// w/ZS:	08000059807016.10.root
// w/o ZS:	08000053726000.10.root

Int_t ReadEvent(AliRawReader* reader)
{
  if (!reader)              return -1;
  if (!reader->NextEvent()) return -1;
  Int_t evno = reader->GetEventIndex();
  // if (evno % 10 == 0) std::cout << "." << std::flush;
  std::cout << "In event # " << evno << "/" 
	    << reader->GetNumberOfEvents() << std::endl;  
  return evno;
}

void
TestRawReader(const char* file=
	      //"/media/cholm_data/alice/data/08000053726000.10.root" // w/o ZS
	      "/media/cholm_data/alice/data/08000059807016.10.root" // w/ZS
	      )
{
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetRun(59807);
  AliFMDParameters::Instance()->Init();
  AliLog::SetModuleDebugLevel("FMD", 2);

  AliRawReader* rawReader = AliRawReader::Create(file);
  if (!rawReader) return;
  
  AliFMDRawReader* fmdReader = new AliFMDRawReader(rawReader, 0);
  
  AliFMDBoolMap seen;
  
  Int_t evno = -1;
  while ((evno = ReadEvent(rawReader)) >= 0) { 
    if (evno <  25) continue;
    if (evno > 120) break;
    
    UShort_t det, sec, str, fac;
    Short_t adc;
    Bool_t zs;
    Char_t rng;
    seen.Reset(kFALSE);
    
    while (fmdReader->NextSignal(det, rng, sec, str, adc, zs, fac)) { 
      Printf("FMD%d%c[%02d,%03d]: %4d", det,rng,sec,str,adc);
      seen(det,rng,sec,str) = kTRUE;
    }
#if 0
    for (det = 1; det <= 3; det++) { 
      UShort_t nrng = (det == 1 ? 1 : 2);
      for (UShort_t ir = 0; ir < nrng; ir++) { 
	Char_t   rng  = (ir == 0 ? 'I' : 'O');
	UShort_t nsec = (ir == 0 ?  20 :  40);
	UShort_t nstr = (ir == 0 ? 512 : 256);
	for (sec = 0; sec < nsec; sec++) { 
	  for (str = 0; str < nstr; str++) {
	    if (!seen(det,rng,sec,str)) 
	      Printf("FMD%d%c[%02d,%03d] not seen", det,rng,sec,str);
	  }
	}
      }
    }
#endif
  }
}

    
    
