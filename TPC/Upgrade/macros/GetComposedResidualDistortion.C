#include <TString.h>
#include <TVectorD.h>
#include <TObjArray.h>
#include <TFile.h>

#include <AliTPCComposedCorrection.h>
#include <AliTPCCorrectionLookupTable.h>

#include <AliToyMCEventGenerator.h>


AliTPCComposedCorrection* GetComposedResidualDistortion(TString fluctuationMap, TString averageMap, Bool_t rescale=kTRUE)
{
  //
  //
  //


  TFile fFluct(fluctuationMap);
  AliTPCCorrectionLookupTable *corrFluct = (AliTPCCorrectionLookupTable*)fFluct.Get("map");
  fFluct.Close();

  TFile fAverage(averageMap);
  AliTPCCorrectionLookupTable *corrAverage = (AliTPCCorrectionLookupTable*)fAverage.Get("map");
  fAverage.Close();

  TObjArray *arrMaps = new TObjArray(2);
  arrMaps->Add(corrFluct);   // distortion with the fluctuation Map
  arrMaps->Add(corrAverage); // correction with the average Map

  // create the composed correction
  // if the weight are set to +1 and -1, the first map will be responsible for the distortions
  // The second map for the corrections
  AliTPCComposedCorrection *residualDistortion = new AliTPCComposedCorrection(arrMaps, AliTPCComposedCorrection::kQueueResidual);
  TVectorD weights(2);
  weights(0)=+1.;
  weights(1)=-1.;
  if (rescale) {
    Float_t dummy=0;
    weights(1)=-AliToyMCEventGenerator::GetSCScalingFactor(corrFluct, corrAverage,dummy);
  }
  residualDistortion->SetWeights(&weights);

  return residualDistortion;
}
