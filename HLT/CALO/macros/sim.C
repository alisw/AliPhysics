void sim(Int_t nev=1, bool dophos = true, bool doemcal = true, bool dotm = true) {
  
  if(!dophos && !doemcal && !dotm)
    {
      cout << "No detectors selected for simulation, exiting..." << endl;
      return;
    }

  AliSimulation simulator;
  simulator.SetMakeDigitsFromHits(":");
  simulator.SetMakeSDigits(":");

  TString sdigits = ""; 
  TString rawdata = "";
  if(dophos) 
    {
      sdigits += "PHOS";
      rawdata += "PHOS";
    }
  if(doemcal) 
    {
      if(sdigits != "")
	{
	  sdigits += " ";
	  rawdata += " ";
	}
      sdigits += "EMCAL";
      rawdata += "EMCAL";
    }
  if(dotm) 
    {
      if(rawdata != "")
	{
	  rawdata += " ";
	}

      simulator.SetMakeDigitsFromHits("TPC");
      rawdata += "TPC";
    }

  if(sdigits != "") 
    {
      cout << "SetMakeSDigits("<< sdigits << ")" << endl;
      simulator.SetMakeSDigits(sdigits);
    }
  
  cout << "sdigits: " << sdigits <<endl;
  cout << "rawdata: " << rawdata <<endl;
  simulator.SetWriteRawData(rawdata);
  //simulator.SetWriteRawData("ALL");
  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));
  
  simulator.SetRunQA(":") ; 

  
  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
