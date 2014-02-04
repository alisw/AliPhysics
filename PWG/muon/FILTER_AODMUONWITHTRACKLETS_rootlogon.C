{
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libProofPlayer");
  gSystem->Load("libPhysics");
  gSystem->Load("libMatrix");
  gSystem->Load("libMinuit");
  gSystem->Load("libXMLParser");
  gSystem->Load("libGui");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libCDB");
  gSystem->Load("libSTEER");
  gSystem->Load("libANALYSISalice");
  
  //    gSystem->Load("libCORRFW");
  //    gSystem->Load("libPWGmuon");
  
  gROOT->LoadMacro("AliAnalysisNonMuonTrackCuts.cxx+");
  gROOT->LoadMacro("AliAnalysisNonPrimaryVertices.cxx+");
  gROOT->LoadMacro("AliAODMuonReplicator.cxx+");
  gROOT->LoadMacro("AliAnalysisTaskESDMuonFilter.cxx+");
  gROOT->LoadMacro("AliAnalysisTaskAOD2MuonAOD.cxx+");
}
