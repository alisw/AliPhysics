
// macro to test AliITSSAPTracker in offline mode using already reconstructed data 
// The input should be either
// 1) directory containing the output of the reconstruction, e.g. ppbench/
// (we need AliESDs.root to obtain the input SPD vertex and ITS.RecPoints.root)
// 2) for multiple input files: text file with paths to AliESDs.root files, 
// e.g. /data1/LHC10h8/137366/003/AliESDs.root etc.
// Full reconstruction output should be in these directories (including galice.root)

void TestITSSAP(const char *datapath = "lst.txt", int maxEv=-1) {
  //
  TString dtPath = datapath; 
  //
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("alien://folder=/alice/data/2011/OCDB");
  man->SetSpecificStorage("ITS/Align/Data","alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
  man->SetSpecificStorage("ITS/Align/SPDSparseDead","alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
  man->SetSpecificStorage("MUON/Align/Data","alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
  man->SetSpecificStorage("TPC/Calib/TimeGain","alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
  man->SetSpecificStorage("TPC/Calib/ClusterParam","alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
  man->SetSpecificStorage("TPC/Calib/AltroConfig","alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
  man->SetSpecificStorage("TPC/Calib/Correction","alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
  man->SetSpecificStorage("TPC/Align/Data","alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
  man->SetSpecificStorage("TPC/Calib/TimeDrift","alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
  man->SetSpecificStorage("TPC/Calib/RecoParam","alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
  

  //  man->SetDefaultStorage("local:///home/shahoian/ALICE/Aliroot/OCDB");
  //  man->SetSpecificStorage("ITS/Align/Data","alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  //  man->SetSpecificStorage("ITS/Align/Data","alien://folder=/alice/simulation/2008/v4-15-Release/Ideal");
  //  man->SetSpecificStorage("ITS/Align/Data","local:///alice/simulation/2008/v4-15-Release/Ideal");
  //man->SetSpecificStorage("ITS/Align/Data","local:///alice/simulation/2008/v4-15-Release/Residual");
  //
  TString inpData;
  if (!(dtPath.EndsWith(".txt")||dtPath.EndsWith(".dat"))) {
    inpData = Form("%s/AliESDs.root",dtPath.Data());
    //    gSystem->Exec(Form("ln -s -f %s/geometry.root ./",dtPath.Data()));
  }
  else inpData = dtPath;
  printf("InputData : %s\n",inpData.Data());
  gSystem->Load("libAliHLTITS");
  gROOT->ProcessLine(".L Process.C+");
  Process(inpData.Data(), maxEv);
  //
}
