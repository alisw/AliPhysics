R__LOAD_LIBRARY(liblhapdf)
R__LOAD_LIBRARY(libEGPythia6)
R__LOAD_LIBRARY(libpythia6)
R__LOAD_LIBRARY(libAliPythia6)
R__LOAD_LIBRARY(libHIJING)
R__LOAD_LIBRARY(libTHijing)
R__LOAD_LIBRARY(libTTherminator)

#include "AliGenerator.h"
#include "AliGenTherminator.h"

AliGenerator *AddMCGenTherminator()
{  
// User defined generator  
  AliGenTherminator* gener = new AliGenTherminator(-1);

  gener->SetModel("SingleFreezeOut"); // setting for Cracow model
  gener->SetEventNumberInFile(500);
  gener->SetTau(9.0);
  gener->SetRhoMax(11.4);
  gener->SetMiuB(0.0);
//  gener->SetPtHardMin(2.3); this function is not derived in AliGenTherminator or its base classes ... 
  
  return gener;
}
  
  
