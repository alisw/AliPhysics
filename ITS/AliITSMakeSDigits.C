void AliITSMakeSDigits(Int_t n){
//gInterpreter->ProcessLine(".L ($ALICE_ROOT)/ITS/AliITSHits2SDigits.C");

gROOT->LoadMacro("$ALICE_ROOT/ITS/AliITSHits2SDigits.C");

AliITSHits2SDigits(0,n-1);
}
