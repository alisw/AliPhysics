{

// sample rootlogon to be used with AliMuonAccEffSubmitter. Feel free to adapt to your needs

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
gSystem->Load("libCORRFW");

gSystem->Load("libPWGmuon");

gSystem->Load("libMUONcore");
gSystem->Load("libMUONmapping");
gSystem->Load("libMUONcalib");
gSystem->Load("libMUONgeometry");
gSystem->Load("libMUONtrigger");  
gSystem->Load("libRAWDatabase");
gSystem->Load("libMUONraw");
gSystem->Load("libMUONbase");
gSystem->Load("libMUONshuttle");
gSystem->Load("libMUONrec");
gSystem->Load("libMUONgraphics");
  
gSystem->Load("libPWGmuondep");

gSystem->Load("libEVGEN");

gSystem->SetIncludePath("-I. -I$ALICE_INSTALL/include -I$ALICE_ROOT/PWG/muon -I$ALICE_ROOT/PWG/muondep -I$ALICE_ROOT/MUON");
  
gStyle->SetPalette(1);

}

