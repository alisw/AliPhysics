//____________________________________________________________________
//
// $Id$
//
// Test of AliFMDAltro{Reader,Writer}
//
/** @ingroup FMD_simple_script
 */
UShort_t 
MakeADC(UShort_t d, UShort_t q, UShort_t s, UShort_t t)
{
  return t & 0x3FF;
}
void 
MakeData(TClonesArray& array, AliFMDUShortMap& map)
{
  size_t i = 0;
  for (size_t d = 1; d <= 3; d++) { 
    size_t nRng = (d == 1 ? 1 : 2);
    for (size_t q = 0; q < nRng; q++) { 
      Char_t r    = (q == 0 ? 'I' : 'O');
      size_t nSec = (q == 0 ?  20 :  40);
      size_t nStr = (q == 0 ? 512 : 256);
      for (size_t s = 0; s < nSec; s++) { 
	for (size_t t = 0; t < nStr; t++) { 
	  UShort_t adc = MakeADC(d, q, s, t);
	  new (array[i++]) AliFMDDigit(d, r, s, t, adc, adc, adc, adc);
	  map(d,r,s,t) = adc;
	}
      }
    }
  }
}
void 
CompareData(TClonesArray& in, const AliFMDUShortMap& map)
{
  size_t i = 0;
  for (size_t i = 0; i < in.GetEntries(); i++) { 
    AliFMDDigit* inD  = static_cast<AliFMDDigit*>(in.At(i));
    // AliFMDDigit* outD = static_cast<AliFMDDigit*>(out.At(i));

    UShort_t out = map(inD->Detector(),
		       inD->Ring(),
		       inD->Sector(),
		       inD->Strip());
    if (out != inD->Counts()) {
      std::cout << "Entries " << i << " does not match up - expected " 
		<< out << "\n  ";
      inD->Print();
    }
  }
}

void
TestRawIO()
{
  AliCDBManager*    cdb   = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);

  AliFMDParameters* param = AliFMDParameters::Instance();
  param->Init();

  TClonesArray    out("AliFMDDigit");
  AliFMDUShortMap map(0);
  MakeData(out, map);

  AliFMDRawWriter writer(0);
  writer.WriteDigits(&out);

  // gSystem->mkdir("raw0");
  // gSystem->Rename("FMD_3072.ddl", "raw0/FMD_3072.ddl");
  // gSystem->Rename("FMD_3073.ddl", "raw0/FMD_3073.ddl");
  // gSystem->Rename("FMD_3074.ddl", "raw0/FMD_3074.ddl");

  AliRawReader* raw = AliRawReader::Create("./");
  AliFMDRawReader reader(raw, 0);
  TClonesArray    in("AliFMDDigit");
  raw->NextEvent();
  reader.ReadAdcs(&in);

  std::cout << "Got " << in.GetEntries() << std::endl;

  CompareData(in, map);
}
//____________________________________________________________________
//
// EOF
//
