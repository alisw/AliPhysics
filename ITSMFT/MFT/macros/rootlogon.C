{

  TString includePath = "-I${ROOTSYS}/include ";
  
  printf("\nLoading ALICE settings...\n\n");
  includePath        += "-I${ALICE_ROOT}/STEER ";
  includePath        += "-I${ALICE_ROOT}/STEER/STEERBase ";
  includePath        += "-I${ALICE_BUILD}/include ";
  includePath        += "-I${ALICE_ROOT}/RAW ";
  includePath        += "-I${ALICE_ROOT}/FASTSIM ";
  includePath        += "-I${ALICE_ROOT}/EVGEN ";
  includePath        += "-I${ALICE_ROOT}/SHUTTLE/TestShuttle ";
  includePath        += "-I${ALICE_ROOT}/ITS ";
  includePath        += "-I${ALICE_ROOT}/MFT ";
  includePath        += "-I${ALICE_ROOT}/MUON ";
  includePath        += "-I${ALICE_ROOT}/MUON/mapping ";
  includePath        += "-I${ALICE_ROOT}/RAW ";
  includePath        += "-I${ALICE_ROOT}/CORRFW ";
  includePath        += "-I${ROOTSYS}/net/alien/inc ";

  gSystem->SetIncludePath(includePath.Data());

  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS") ;
  gSystem->Load("libANALYSISalice") ;
  gSystem->Load("libCORRFW") ;
  gSystem->Load("libMFTbase") ;
  gSystem->Load("libMFTsim") ;
  gSystem->Load("libMFTrec") ;
  
  //================= define the style =================== //
  
  // divisions: axis marks outside the figure (minus sign)
  gStyle->SetTickLength(-0.03,"X");  gStyle->SetTickLength(-0.03,"Y"); 
  gStyle->SetNdivisions(508,"X");    gStyle->SetNdivisions(506,"Y"); // set ndvx (y)
  
  // dimensions and position of fonts connected to axes  
  gStyle->SetLabelSize(0.05,"X");    gStyle->SetLabelSize(0.05,"Y");  
  gStyle->SetLabelOffset(0.02,"X");  gStyle->SetLabelOffset(0.04,"Y");  
  
  // dimension of the canvas and of the figure (pad) inside it
  gStyle->SetCanvasDefW(600);        gStyle->SetCanvasDefH(600);
  gStyle->SetPadBottomMargin(0.16);  gStyle->SetPadLeftMargin(0.18);    
  gStyle->SetPadTopMargin(0.08);     gStyle->SetPadRightMargin(0.06);    
  gStyle->SetTitleXSize(0.05);       gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleXOffset(1.5);      gStyle->SetTitleYOffset(1.9);
  gStyle->SetStatW(0.16);            gStyle->SetStatH(0.16);
  gStyle->SetFitFormat("9.6g"); 
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1);
  
  gStyle->SetPalette(1,0);  // rainbow colors
  gStyle->SetHistLineWidth(2); 
  
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  
}

