#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <TStopwatch.h>
  #include "AliRawReaderRoot.h"
  #include "AliRawReaderDate.h"
  #include "AliAltroRawStreamFast.h"
  #include "AliLog.h"
#endif


void altroFast(const char *fileName)
{
  //  AliLog::SetGlobalLogLevel(AliLog::kFatal);

  //  AliRawReader *reader = new AliRawReaderRoot(fileName);
  AliRawReader *reader = new AliRawReaderDate(fileName);
  reader->Reset();

  TStopwatch timer;
  timer.Start();

  AliAltroRawStreamFast* stream = new AliAltroRawStreamFast(reader);
  stream->SelectRawData("TPC");

  while (reader->NextEvent()) {

    while (stream->NextDDL()) {

      while (stream->NextChannel()) {

	while (stream->NextBunch()) {
	  const UInt_t *adc = stream->GetSignals();
	  for(UInt_t i = stream->GetStartTimeBin(); i <= stream->GetEndTimeBin(); i++) {
	    // cout i - timebin, *adc - ADC signal, ...
	    adc++;
	  }
	}
      }
    }
  }

  timer.Stop();
  timer.Print();

  delete stream;

  return;
}
