// $Id$

{
  printf("\nWELCOME to the magic world of HLT\n\n"); 

  gSystem->Load("$(ROOTSYS)/lib/libPhysics");
  gSystem->Load("$(ROOTSYS)/lib/libEG");

  Int_t saveErrIgLevel=gErrorIgnoreLevel;
  gErrorIgnoreLevel=kFatal; //dont report errors
  if(gSystem->Load("$(ROOTSYS)/lib/libMC")==-1){
    gSystem->Load("$(ROOTSYS)/lib/libGeom");
    gSystem->Load("$(ROOTSYS)/lib/libVMC");
  }

  if(1)
    {
      gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libSTEER");
      gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libCONTAINERS");
      if(gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTPC")!=0){
	gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTCPbase");
	gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTCPrec");
	gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTCPsim");
	gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTCPfast");
      }
      cout<<"TPC libraries loaded"<<endl;
    }

  if(1)
    {
      gSystem->Load("$(ALIHLT_LIBDIR)/libAliHLTSrc");
      gSystem->Load("$(ALIHLT_LIBDIR)/libAliHLTMisc");
      gSystem->Load("$(ALIHLT_LIBDIR)/libAliHLTHough");
      gSystem->Load("$(ALIHLT_LIBDIR)/libAliHLTComp");
      cout<<"HLT libraries loaded"<<endl;
    }

  gSystem->Load("MakePileup_C.so");
  gSystem->Load("Read_C.so");
  gSystem->Load("runhough_C.so");
  gSystem->Load("runrowhough_C.so");
  gSystem->Load("deconvclusters_C.so");
  gSystem->Load("runtracker_C.so");
  gErrorIgnoreLevel=saveErrIgLevel;

  gSystem->SetIncludePath(" -DROWHOUGHPARAMS -Duse_reconstruction -Duse_newio -D do_mc -I$ALIHLT_TOPDIR/hough -I$ALIHLT_TOPDIR/src -I$ALIHLT_TOPDIR/comp -I$ALIHLT_TOPDIR/misc -I$ALICE_ROOT/include/ -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER ");

  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetPadColor(10);
  gStyle->SetPadBorderSize(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleSize(0.036,"X");
  gStyle->SetTitleSize(0.036,"Y");
  gStyle->SetTitleSize(0.036,"Z");
  gStyle->SetLabelSize(0.036,"X");
  gStyle->SetLabelSize(0.036,"Y");
  gStyle->SetLabelSize(0.036,"Z");
  gStyle->SetTitleOffset(1.2,"X");
  gStyle->SetTitleOffset(1.3,"Y");
  gStyle->SetTitleOffset(1.3,"Z");
  gStyle->SetTitleColor(10);
}
