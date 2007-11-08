// $Id$
/*
 * Configuration macro defining a chain running the HLT Conformal mapping
 * tracker embedded into AliRoot simulation.
 * The reconstruction is done from the TPC digits.
 *
 * The output is written to an ESD file.
 *
 * Matthias.Richter@ift.uib.no
 */
{
  int iMinSlice=0; 
  int iMaxSlice=35;
  int iMinPart=0;
  int iMaxPart=5;
  TString writerInput;
  for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
    TString trackerInput;
    for (int part=iMinPart; part<=iMaxPart; part++) {
      TString arg, publisher, cf;

      // digit publisher components
      arg.Form("-slice %d -partition %d", slice, part);
      publisher.Form("DP_%02d_%d", slice, part);
      AliHLTConfiguration pubconf(publisher.Data(), "TPCDigitPublisher", NULL , arg.Data());

      // cluster finder components
      cf.Form("CF_%02d_%d", slice, part);
      AliHLTConfiguration cfconf(cf.Data(), "TPCClusterFinderUnpacked", publisher.Data(), "pp-run timebins 446");
      if (trackerInput.Length()>0) trackerInput+=" ";
      trackerInput+=cf;
    }
    TString tracker;
    // tracker finder components
    tracker.Form("TR_%02d", slice);
    AliHLTConfiguration trackerconf(tracker.Data(), "TPCSliceTracker", trackerInput.Data(), "pp-run bfield 0.5");

    if (writerInput.Length()>0) writerInput+=" ";
    writerInput+=tracker;
  }

  // the esd writer configuration
  AliHLTConfiguration esdwconf("esd-writer", "TPCEsdWriter"   , writerInput.Data(), "-datafile AliHLTESDs.root");
}
