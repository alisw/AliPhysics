#include "AliHLTEmcalElectronMonitorComponent.h"
#include "AliHLTEmcalElectronMonitor.h"
#include "AliHLTScalars.h"

#include "TFile.h"
#include "TString.h"

AliHLTEmcalElectronMonitorComponent::AliHLTEmcalElectronMonitorComponent() :
  fRootFileName("EmcalElectrontrigger_histos.root"),
  fPushFraction(10),
  fLocalEventCount(0),
  fVerbose(0),
  fHistoPtr(NULL)
{

  // default constructor

}
//____________________________________________________________________________________________________________

AliHLTEmcalElectronMonitorComponent::~AliHLTEmcalElectronMonitorComponent()
{

  // default destructor

}
//____________________________________________________________________________________________________________

int AliHLTEmcalElectronMonitorComponent::DoInit(int argc, const char **argv)
{
  // initialize

  fHistoPtr = new AliHLTEmcalElectronMonitor();
  for (int i = 0; i < argc; i++) {
    if (!strcmp("-roothistofilename", argv[i]))
      fRootFileName = argv[i+1];
    if (!strcmp("-pushfraction", argv[i]))
      fPushFraction = atoi(argv[i+1]);
    if (!strcmp("-verbose", argv[i]))
      fVerbose = atoi(argv[i+1]);
  }

  return 0;

}
//____________________________________________________________________________________________________________

int AliHLTEmcalElectronMonitorComponent::Deinit()
{
  // de-initialize

  if (fHistoPtr) {
    delete fHistoPtr;
    fHistoPtr = NULL;
  }
  return 0;

}
//____________________________________________________________________________________________________________
const char* AliHLTEmcalElectronMonitorComponent::GetComponentID()
{
  // component id
  
  return "EmcalElectronMonitor";

}
//____________________________________________________________________________________________________________

void AliHLTEmcalElectronMonitorComponent::GetInputDataTypes(vector<AliHLTComponentDataType> &list)
{
  // define input data types
  
  list.clear();
  list.push_back(kAliHLTDataTypeEventStatistics|kAliHLTDataOriginHLT);

}
//____________________________________________________________________________________________________________

AliHLTComponentDataType AliHLTEmcalElectronMonitorComponent::GetOutputDataType()
{

  // return output data types
  
  return kAliHLTDataTypeHistogram | kAliHLTDataOriginEMCAL;

}
//____________________________________________________________________________________________________________

void AliHLTEmcalElectronMonitorComponent::GetOutputDataSize(unsigned long &constBase, double &inputMultiplier)
{

  // calculate output data size
  
  constBase = 0;
  inputMultiplier = 100;

}
//____________________________________________________________________________________________________________

int AliHLTEmcalElectronMonitorComponent::DoEvent(const AliHLTComponentEventData &evtData, const AliHLTComponentBlockData *blocks,
					   AliHLTComponentTriggerData &/*trigData*/, AliHLTUInt8_t */*outputPtr*/, AliHLTUInt32_t &/*size*/,
					   std::vector<AliHLTComponentBlockData> &/*outputBlocks*/)
{

  // do event

  const AliHLTComponentBlockData *iter = NULL;
  UInt_t specification = 0;
  
  for (unsigned long ij = 0; ij < evtData.fBlockCnt; ij++) {
    AliHLTScalars *scalarPtr = NULL;
    iter = blocks + ij;
    if (fVerbose) PrintComponentDataTypeInfo(iter->fDataType);

    if (iter->fDataType == kAliHLTDataTypeEventStatistics)  { 
      scalarPtr = reinterpret_cast<AliHLTScalars*>(iter->fPtr);
    }
    else {
      if (fVerbose)  HLTWarning("Electron Monitor: Data block does not contain event stats - check if flag is set for histograming for AliHLTEmcalElectronTrigger \n");
    }
    
    specification |= iter->fSpecification;
    
    if (scalarPtr)
      fHistoPtr->MakeHisto(scalarPtr);
  }

  fLocalEventCount++;
    
  TFile rootHistFile(fRootFileName, "RECREATE");
  
  fHistoPtr->GetHistograms()->Write();

  if (fLocalEventCount%fPushFraction == 0) {
    if (fVerbose) cout << "Emcal Electron Monitor: pushback done at " << fLocalEventCount << " evens " << endl;
    PushBack(fHistoPtr->GetHistograms(), kAliHLTDataTypeTObjArray | kAliHLTDataOriginEMCAL, specification);
  }
  
  return 0;
}
//____________________________________________________________________________________________________________

AliHLTComponent* AliHLTEmcalElectronMonitorComponent::Spawn()
{
  // spawn

  return new AliHLTEmcalElectronMonitorComponent();

}
