//
// $Id$
//
// Small script to test consistency of writing and reading raw data.
//
/** Check raw data I/O
    @ingroup FMD_simple_script
 */
void
RawTest() 
{
  gRandom->SetSeed(12345);
  Int_t   sampleRate   = 4;
  Int_t   channelWidth = 128;
  Float_t shapingTime  = 5;
  UInt_t  maxAdc       = (1 << 10);
  UInt_t  threshold    = 0; // (1 << 8);
  TArrayI outData(sampleRate * channelWidth); 
  
  Float_t lastTotalCharge = 0;
  Int_t   ok = 0;
  for (Int_t channel = 0; channel < channelWidth; channel++) {
    Float_t totalCharge = gRandom->Uniform(0, 1);

    for (Int_t sample = 0; sample < sampleRate; sample++) {
      Float_t time   = Float_t(sample) / sampleRate + 1./sampleRate;
      Float_t charge = (totalCharge + (lastTotalCharge - totalCharge)
			* TMath::Exp(-shapingTime * time));
      UInt_t  adc    = channel; // UInt_t(maxAdc * charge);      
      outData[channel * sampleRate + sample] = adc;
      if (adc > threshold) ok++;
    }
    lastTotalCharge = totalCharge;
  }
  std::cout << "Total of " << outData.fN << " samples of which " 
	    << ok << " of them are above threshold (" << threshold 
	    << ")" << std::endl;
  
  { 
    AliAltroBuffer buffer("FMD_4096.ddl", new AliFMDAltroMapping());
    buffer.WriteDataHeader(kTRUE, kFALSE);
    buffer.WriteChannel(0, 0, 0, outData.fN, outData.fArray, threshold);
    buffer.Flush();
    buffer.WriteDataHeader(kFALSE, kFALSE);
  }
  
  AliRawReader* reader = new AliRawReaderFile(-1);
  if (!reader) {
    std::cerr << "Failed to make AliRawReader" << endl;
    return 0;
  }
  AliFMDRawStream input(reader); // , sampleRate);
  reader->Select(12); // AliFMDParameters::kBaseDDL >> 8);
  
  Int_t    oldDDL      = -1;
  Int_t    count       = 0;
  UShort_t detector    = 1; // Must be one here
  UShort_t oldDetector = 0;
  // Loop over data in file 
  Bool_t   next       = kTRUE;

  // local Cache 
  TArrayI counts(10);
  counts.Reset(-1);
  Int_t offset = 0;
  UInt_t ddl    = 0;
  UInt_t rate   = 0;
  UInt_t last   = 0;
  UInt_t hwaddr = 0;
  UShort_t data[2048];

  AliFMDParameters* pars = AliFMDParameters::Instance();
  TArrayI inputData(sampleRate * channelWidth); 
  Bool_t isGood = true;
  while (next) {
    isGood = input.ReadChannel(ddl, hwaddr, last, data);
    ddl = 0;
    
    std::cout << Form("Read channel %p 0x%x of size %d", ddl, hwaddr, last)
	      << std::endl;
    UShort_t det, sec, str;
    Char_t   ring;
    if (!pars->Hardware2Detector(ddl, hwaddr, det, ring, sec, str)) {
      std::cerr << Form("Failed to get detector id from DDL %d "
			"and hardware address 0x%x", ddl, hwaddr) << std::endl;
      continue;
    }
    rate           = pars->GetSampleRate(det, ring, sec, str);
    Int_t stripMin = pars->GetMinStrip(det, ring, sec, str);
    Int_t stripMax = pars->GetMaxStrip(det, ring, sec, str);
    std::cout << Form("DDL 0x%04x, address 0x%03x maps to FMD%d%c[%2d,%3d]", 
		      ddl, hwaddr, det, ring, sec, str) << std::endl;
    
    // Loop over the `timebins', and make the digits
    for (size_t i = 0; i < last; i++) {
      Int_t in  = data[i];
      Int_t out = outData[channel * sampleRate + sample];
      std::cout << "[\t" << channel << ",\t" << sample  << "]\t" 
		<< out << "\t" << in << std::flush; 
      if (out >= threshold && in != out) std::cout << "\tBad" << std::flush;
    }
#if 0
    next = input.Next();

    if (!next) break;

    Int_t  channel = input.Strip();
    Int_t  sample  = input.Sample();
    inputData[channel * sampleRate + sample] = input.Count();
    count++;

    Int_t in  = inputData[channel * sampleRate + sample];
    Int_t out = outData[channel * sampleRate + sample];
    std::cout << "[\t" << channel << ",\t" << sample  << "]\t" 
	      << out << "\t" << in << std::flush; 
    if (out >= threshold && in != out) std::cout << "\tBad" << std::flush;
    std::cout << std::endl;    
#endif
  }

  std::cout << "Read " << count << " values" << std::endl;
#if 0
  for (Int_t channel = channelWidth - 1; channel > 0; channel--) {
    for (Int_t sample =  sampleRate - 1; sample > 0; sample--) {
      Int_t in  = inputData[channel * sampleRate + sample];
      Int_t out = outData[channel * sampleRate + sample];
      std::cout << "[\t" << channel << ",\t" << sample  << "]\t" 
		<< out << "\t" << in << std::flush; 
      if (out >= threshold && in != out) std::cout << "\tBad" << std::flush;
      std::cout << std::endl; 
    }
  }
#endif
}

//____________________________________________________________________
//
// EOF
//
  
    
	
    
      
  
    
  
