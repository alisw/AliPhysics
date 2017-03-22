AliTPCCalibXXX* fillCalibObject( const Char_t *filename ){
  AliTPCCalibXXX *calibObject = new AliTPCCalibXXX;
  AliRawReader *rawReader = new AliRawReaderRoot(filename);
  if ( !rawReader ) return 0x0;
  rawReader->RewindEvents();
  //loop over all events
  while (rawReader->NextEvent())
    events = calibObject->ProcessEvent(rawReader);

  //Analyse accumulated data
  calibObject->Analyse();
  //Write Calibration class to file
  calibObject->DumpToFile("CalibXXXDataFile.root");
  //return the Calibration Object
  return calibObject;
}
