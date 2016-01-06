void includeMacro(){

  // macro to set the necessary include directories to compile AliTPCcalibAlignInterpolationMacro.C

  gSystem->SetIncludePath("-I$ALICE_PHYSICS/include -I$ALICE_ROOT/include"); 

  return;

}
