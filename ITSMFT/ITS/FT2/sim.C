void sim(Int_t nev=3,Int_t runnumber) {
	
	gSystem->Exec(" rm itsSegmentations.root ");
	
	
	AliSimulation simulator;
	simulator.SetMakeSDigits("");
	simulator.SetMakeDigitsFromHits("");
	simulator.SetRunNumber(runnumber);
	//
	simulator.SetDefaultStorage("raw://");
	simulator.SetSpecificStorage("GRP/*/*","alien://folder=/alice/data/2010/OCDB");
	
	
	simulator.SetSpecificStorage("ITS/Align/Data", "alien://folder=/alice/simulation/LS1_upgrade/Ideal");
	simulator.SetSpecificStorage("ITS/Calib/RecoParam", "alien://folder=/alice/simulation/LS1_upgrade/Ideal");
	simulator.SetSpecificStorage("ITS/Calib/SimuParam", "alien://folder=/alice/simulation/LS1_upgrade/Ideal");
	simulator.SetDefaultStorage("alien://folder=/alice/data/2010/OCDB");
	simulator.SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/simulation/2008/v4-15-Release/Ideal/");
	simulator.SetSpecificStorage("TPC/Align/Data", "alien://folder=/alice/simulation/2008/v4-15-Release/Ideal/");
	simulator.SetSpecificStorage("TPC/Calib/ClusterParam", "alien://folder=/alice/simulation/2008/v4-15-Release/Ideal/");
	simulator.SetSpecificStorage("TPC/Calib/RecoParam", "alien://folder=/alice/simulation/2008/v4-15-Release/Residual/");
	simulator.SetSpecificStorage("TPC/Calib/TimeGain", "alien://folder=/alice/simulation/2008/v4-15-Release/Ideal/");
	simulator.SetSpecificStorage("TPC/Calib/AltroConfig", "alien://folder=/alice/simulation/2008/v4-15-Release/Ideal/");
	simulator.SetSpecificStorage("TPC/Calib/TimeDrift", "alien://folder=/alice/simulation/2008/v4-15-Release/Ideal/");
	simulator.SetSpecificStorage("TPC/Calib/Correction", "alien://folder=/alice/simulation/2008/v4-15-Release/Ideal/");
	simulator.Print("");
	
	//
	// Vertex and Mag.field from OCDB
	//
	simulator.UseVertexFromCDB();
	simulator.UseMagFieldFromGRP();
	simulator.SetRunHLT("");
	simulator.SetRunQA(":");
	//
	
	// The rest
	
	//
	
	printf("Before simulator.Run(nev);\n");
	simulator.Run(nev);
	printf("After simulator.Run(nev);\n");
}
