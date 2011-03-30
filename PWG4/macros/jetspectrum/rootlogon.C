{

Printf("%s \n>>> Setting style",(char*)__FILE__);

 gStyle->SetPalette(1);
 gStyle->SetCanvasColor(10);
 gStyle->SetHistFillColor(10);
 gStyle->SetHistFillStyle(0);
 gStyle->SetOptTitle(0);
 gStyle->SetOptStat(0);
 gStyle->SetPadLeftMargin(0.17);
 gStyle->SetPadBottomMargin(0.2);
 gStyle->SetPadTickX(1);
 gStyle->SetPadTickY(1);
 gStyle->SetAxisColor(1, "X");
 gStyle->SetAxisColor(1, "Y");
 gStyle->SetAxisColor(1, "Z");
 gStyle->SetLabelColor(1, "X");
 gStyle->SetLabelColor(1, "Y");
 gStyle->SetLabelColor(1, "Z");
 gStyle->SetTickLength(0.03, "X");
 gStyle->SetTickLength(0.03, "Y");
 gStyle->SetTickLength(0.03, "Z");
 gStyle->SetTitleXSize(0.07);
 gStyle->SetTitleYSize(0.07);
 gStyle->SetNdivisions(505, "X");
 gStyle->SetNdivisions(505, "Y");
 gStyle->SetNdivisions(505, "Z");
 gStyle->SetTitleXOffset(1.0);
 gStyle->SetTitleYOffset(1.0);
 gStyle->SetLabelOffset(0.02, "X");
 gStyle->SetLabelOffset(0.02, "Y");
 gStyle->SetLabelOffset(0.02, "Z");
 gStyle->SetLabelSize(0.05, "X");
 gStyle->SetLabelSize(0.05, "Y");
 gStyle->SetLabelSize(0.05, "Z");

 gROOT->ForceStyle();

 // Printf("Adding include path $ALICE_ROOT/PWG4/JetTasks");
 gInterpreter->AddIncludePath("$ALICE_ROOT/PWG4/JetTasks");
 gInterpreter->AddIncludePath("$ALICE_ROOT/include");
 gInterpreter->AddIncludePath("$HOME/root");

 // Printf("Loading macro Normalize2D");
 // gROOT->LoadMacro("Normalize2D.C");

 Printf("Loading train libs");
 gSystem->Load("libTree");
 gSystem->Load("libPhysics");
 gSystem->Load("libHist");
 gSystem->Load("libVMC");
 gSystem->Load("libCGAL");
 gSystem->Load("libfastjet");

 gSystem->Load("libsiscone");
 gSystem->Load("libSISConePlugin");
 gSystem->Load("libSTEERBase");
 gSystem->Load("libESD");
 gSystem->Load("libAOD");
 gSystem->Load("libOADB");
 gSystem->Load("libANALYSIS");
 gSystem->Load("libANALYSISalice");
 gSystem->Load("libCORRFW");
 gSystem->Load("libPWG3base");
 gSystem->Load("libPWG3muon");
 gSystem->Load("libJETAN");
 gSystem->Load("libFASTJETAN");
 gSystem->Load("libPWG4JetTasks");
 gSystem->Load("libPWG4JCORRAN");
 gSystem->Load("libEMCALUtils");
 gSystem->Load("libPHOSUtils");
 gSystem->Load("libPWG4PartCorrBase");
 gSystem->Load("libPWG4PartCorrDep");
 gSystem->Load("libPWG4GammaConv");
 gSystem->Load("libPWG4omega3pi");

}
