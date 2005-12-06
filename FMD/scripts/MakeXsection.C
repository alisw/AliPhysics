//
// $Id$
//
// Script to generate cross-section tables for some particles in some
// mediums.   The output is put in xsec.root. 
// 
// It uses the script GetXsection.C and Compile.C to compile the
// aforementioned script. 
//
// Note, that VMC _must_ be the TGeant3TGeo VMC. 
//
void
MakeXsection()
{
  gROOT->ProcessLine(".x Compile.C(\"$ALICE_ROOT/FMD/scripts/GetSection.C\"");
  gAlice->InitMC("$(ALICE_ROOT)/FMD/Config.C");
  TFile* file = TFile::Open("xsec.root", "RECREATE");
  GetXsection("FMD_Si$", "pi+");
  file->Close();
}
//____________________________________________________________________
//
// EOF
// 

