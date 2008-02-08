// #include "TSystem.h"
// #include "TMath.h"
// #include "AliRawReader.h"
// #include "AliRawReaderRoot.h"
// #include "AliTPCRawStreamFast.h"
// #include "AliTPCCalibPedestal.h"

void fillPedestal( ){
  AliTPCCalibPedestal *signal = new AliTPCCalibPedestal;
    AliRawReader *rawReader = new AliRawReaderRoot("data.root");
    if ( !rawReader ) return;

    rawReader->RewindEvents();
    AliTPCRawStreamFast *strFast = new AliTPCRawStreamFast(rawReader);

    while (rawReader->NextEvent())  signal->ProcessEventFast(strFast);
    signal->Analyse();
    signal->DumpToFile("calibPedestal.root");


    delete strFast;
    delete rawReader;
    delete signal;
}
