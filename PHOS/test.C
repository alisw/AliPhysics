void test (Int_t nevent=1) 
{
 gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I.");
 gROOT->ProcessLine(".x testsim.C(nevent)") ; 
 gROOT->ProcessLine(".x testsimglobal.C++") ; 
}
