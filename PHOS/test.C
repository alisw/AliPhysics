void test (Int_t nevent=1) 
{
 gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I.");
 gROOT->ProcessLine(".x $ALICE_ROOT/PHOS/testsim.C(nevent)") ; 
 gROOT->ProcessLine(".x $ALICE_ROOT/PHOS/testsimglobal.C++") ; 
}
