// $Id$
//
// Configuration macro for primary event generation for gun test (in vmctest).
//
// By I. Hrivnacova, IPN Orsay

void Config()
{
  // AliRoot setup
  commonConfig(kFALSE);

  // AliRoot event generator
  genGunConfig();
}  
