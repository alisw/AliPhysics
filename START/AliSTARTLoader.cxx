//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// Loader for START digits and RecPoints inherit from AliBaseDataLoader     //
// STARTDigits is TObject  consists time flight signal from each PMTtube,   //
// mean time right and left array (start signal) and time differnce         //
// (vertex z position)                                                      // 
// START RecPoints is TObject with mean time (start signal)                 // 
// and evrtex position  (cm)                                                //
/////////////////////////////////////////////////////////////////////////////
#include "AliSTARTLoader.h"
#include "AliConfig.h"
#include "AliRunLoader.h"


ClassImp(AliSTARTLoader)

/*****************************************************************************/ 
void AliSTARTLoader::InitObjectLoaders()
{
// use an object instead of a tree for digits and rec points

  // D I G I T S  
  if (fDataLoaders->At(kDigits)) {
    delete fDataLoaders->Remove(fDataLoaders->At(kDigits));
  }
  AliDataLoader* dl = new AliDataLoader(fDetectorName + ".Digits.root","START_D", "Digits","O");//we want to have object data not tree
  AliTaskLoader* tl = new AliTaskLoader(fDetectorName + AliConfig::Instance()->GetDigitizerTaskName(),dl,AliRunLoader::GetRunDigitizer(),kTRUE);
  dl->SetBaseTaskLoader(tl);
  fDataLoaders->AddAt(dl,kDigits);

  // R E C O N S T R U C T E D   P O I N T S, here: V E R T E X
  if (fDataLoaders->At(kRecPoints)) {
    delete fDataLoaders->Remove(fDataLoaders->At(kRecPoints));
  }
  dl = new AliDataLoader(fDetectorName + ".RecPoints.root","START_V", "Reconstructed Points","O");//we want to have object data not tree
  tl = new AliTaskLoader(fDetectorName + AliConfig::Instance()->GetReconstructionerTaskName(),dl,AliRunLoader::GetRunReconstructioner(),kTRUE);
  dl->SetBaseTaskLoader(tl);
  fDataLoaders->AddAt(dl,kRecPoints);  
}
