/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//-----------------------------------------------------------------------------
/// \class AliHLTMUONFullTrackerComponent
/// 
/// Component class for full tracker, see the detail description 
/// of the full tracker in the full tracker header file
///  \author :Indranil Das, email : indra.das@saha.ac.in | indra.ehep@gmail.com , Saha Institute of Nuclear Physics
//-----------------------------------------------------------------------------


#if __GNUC__== 3
using namespace std;
#endif

#include "AliHLTMUONFullTrackerComponent.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"

#include "AliHLTDefinitions.h"

#include "AliHLTMUONConstants.h"

#include "AliHLTMUONTracksBlockStruct.h"

ClassImp(AliHLTMUONFullTrackerComponent)

AliHLTMUONFullTrackerComponent::AliHLTMUONFullTrackerComponent() :
  AliHLTMUONProcessor(),
  fTracker(NULL)
{
  // see header file for class documentation
  
}

AliHLTMUONFullTrackerComponent::~AliHLTMUONFullTrackerComponent()
{
  // see header file for class documentation
  
  if (fTracker != NULL) delete fTracker;
}

const char* AliHLTMUONFullTrackerComponent::GetComponentID()
{
  // see header file for class documentation
  return AliHLTMUONConstants::FullTrackerId();
}

void AliHLTMUONFullTrackerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  /* in order to be backward compatible we have to keep the old code, at
   * least for a while. Remember to use the new const kAliHLTVoidDataType
   * if you are using a more recent AliRoot version (from Jan 07)
   list.push_back(kAliHLTAnyDataType); // We do not have any requirements for our input data type(s).
  */

  assert( list.empty() );
  list.push_back( AliHLTMUONConstants::TriggerRecordsBlockDataType() );
  list.push_back( AliHLTMUONConstants::RecHitsBlockDataType() );

}

AliHLTComponentDataType AliHLTMUONFullTrackerComponent::GetOutputDataType()
{
  // see header file for class documentation
  /* in order to be backward compatible we have to keep the old code, at
   * least for a while. Remember to use the new const kAliHLTVoidDataType
   * if you are using a more recent AliRoot version (from Jan 07)
   return kAliHLTVoidDataType;
  */
  return kAliHLTMultipleDataType;

}

int AliHLTMUONFullTrackerComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& list)
{
  /// Inherited from AliHLTComponent. Returns the output data types.
	
  assert( list.empty() );
  list.push_back( AliHLTMUONConstants::TracksBlockDataType() );
  return list.size();
}

void AliHLTMUONFullTrackerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase = sizeof(AliHLTMUONTracksBlockStruct) + 1024*1024;
  inputMultiplier = 1;

}



// Spawn function, return new instance of this class
AliHLTComponent* AliHLTMUONFullTrackerComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTMUONFullTrackerComponent;
}

int AliHLTMUONFullTrackerComponent::DoInit( int argc, const char** argv )
{
  // perform initialization
  bool useFast = false;
  int result = AliHLTMUONProcessor::DoInit(argc, argv);
  if (result != 0) return result;
  
  
  for (int i = 0; i < argc; i++){
	
    // To keep the legacy behaviour we need to have the following check
    // for -cdbpath here, before ArgumentAlreadyHandled.
    if (ArgumentAlreadyHandled(i, argv[i])) continue;

    if (strcmp(argv[i], "-usefast") == 0)
      {
	useFast = true;
      }
    
  }

  
  if (fTracker == NULL)
  {
    try
    {
      fTracker = new AliHLTMUONFullTracker;
    }
    catch (const std::bad_alloc&)
    {
      HLTError("Could not allocate a new AliHLTMUONFullTracker object. Ran out of memory.");
      return -ENOMEM;
    }
  }
  
  result = FetchMappingStores();
  fTracker->Init();
  fTracker->FastTracking(useFast);
  
  return 0;
}

int AliHLTMUONFullTrackerComponent::DoDeinit()
{
  // see header file for class documentation
  
  if (fTracker != NULL)
  {
    delete fTracker;
    fTracker = NULL;
  }
  
  return 0;
}

int AliHLTMUONFullTrackerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  // Process an event
  unsigned long totalSize = 0;
  AliHLTUInt32_t specification = 0;  // Contains the output data block spec bits.
  bool resultOk = true;
  // Loop over all input blocks in the event

  HLTDebug("Processing event %llu with %u input data blocks.",
	 evtData.fEventID, evtData.fBlockCnt
	 );

  //if(evtData.fBlockCnt==3) return 0;

    AliHLTMUONTracksBlockWriter block(outputPtr, size);  
    if (not block.InitCommonHeader())
    {
      Logging(kHLTLogError,
	      "AliHLTMUONFullTrackerComponent::DoEvent",
	      "Buffer overflow",
	      "The buffer is only %d bytes in size. We need a minimum of %d bytes.",
	      size, sizeof(AliHLTMUONTracksBlockWriter::HeaderType)
	      );
      if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
      size = 0; // Important to tell framework that nothing was generated.
      resultOk = false;
      return -ENOBUFS;
    }

  
  for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++){

    //if(evtData.fBlockCnt==3) continue;

    HLTDebug("Handling block: %u, with fDataType = '%s', fPtr = %p and fSize = %u bytes.",
  	     n, DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fPtr, blocks[n].fSize
  	     );

    if (blocks[n].fDataType == AliHLTMUONConstants::RecHitsBlockDataType()){
      specification |= blocks[n].fSpecification;
			
      AliHLTMUONRecHitsBlockReader inblock(blocks[n].fPtr, blocks[n].fSize);
      if (not BlockStructureOk(inblock)){
  	if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
	resultOk = false;
  	continue;
      }

      fTracker->SetInput(AliHLTMUONUtils::SpecToDDLNumber(blocks[n].fSpecification), inblock.GetArray(), inblock.Nentries());	
      

    }else if (blocks[n].fDataType == AliHLTMUONConstants::TriggerRecordsBlockDataType()){
      specification |= blocks[n].fSpecification;
      AliHLTMUONTriggerRecordsBlockReader inblock(blocks[n].fPtr, blocks[n].fSize);
      if (not BlockStructureOk(inblock)){
  	if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
	resultOk = false;
  	continue;
      }
      
      fTracker->SetInput(AliHLTMUONUtils::SpecToDDLNumber(blocks[n].fSpecification), inblock.GetArray(),inblock.Nentries());	
    }//check if trigger block
  }//loop over blocks array of rechit and trigrecs

  AliHLTUInt32_t nofTracks = 0;
  if( resultOk ){
    
    nofTracks   = block.MaxNumberOfEntries();
    //  if (evtData.fBlockCnt!=3){
    resultOk  = fTracker->Run(int(evtData.fEventID),block.GetArray(), nofTracks);
  }

  if (resultOk){
    assert( nofTracks <= block.MaxNumberOfEntries() );
    block.SetNumberOfEntries(nofTracks);
		
    HLTDebug("Number of reconstructed tracks found is %d\n", nofTracks);
    HLTDebug("sizeof  %d\n", sizeof(AliHLTMUONMansoTrackStruct));
    HLTDebug("Bytes Used  is %d\n",block.BytesUsed());    
    HLTDebug("specification is %d\n", specification);

    AliHLTComponentBlockData bd;
    FillBlockData(bd);
    bd.fPtr = outputPtr;
    bd.fOffset = 0;
    bd.fSize = block.BytesUsed();
    bd.fDataType = AliHLTMUONConstants::TracksBlockDataType() ;
      
    bd.fSpecification = specification;
    outputBlocks.push_back(bd);
    totalSize = block.BytesUsed();
  }else{
    HLTDebug("Error while processing the full tracker algorithm.");
    if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
    size = totalSize; // Must tell the framework how much buffer space was used.
    //return -EIO;
    return 0;
  }

  // Finally we set the total size of output memory we consumed.
  size = totalSize;
  return 0;

}


int AliHLTMUONFullTrackerComponent::Configure(const char* arguments)
{
  // see header file for class documentation
  int iResult=0;
  if (!arguments) return iResult;
  HLTInfo("parsing configuration string \'%s\'", arguments);

  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;

      // -config1
      if (argument.CompareTo("-config1")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("got \'-config1\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());

	// -config2
      } else if (argument.CompareTo("-config2")==0) {
	HLTInfo("got \'-config2\'");
      } else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTMUONFullTrackerComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation
  int iResult=0;
  const char* path="HLT/ConfigSample/FullTrackerComponent";
  const char* defaultNotify="";
  if (cdbEntry) {
    path=cdbEntry;
    defaultNotify=" (default)";
  }
  if (path) {
    HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
	HLTInfo("received configuration object string: \'%s\'", pString->GetString().Data());
	iResult=Configure(pString->GetString().Data());
      } else {
	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("can not fetch object \"%s\" from CDB", path);
    }
  }
  return iResult;
}

int AliHLTMUONFullTrackerComponent::ReadPreprocessorValues(const char* modules)
{
  // see header file for class documentation
  int iResult=0;
  TString detectors(modules!=NULL?modules:"");
  HLTInfo("read preprocessor values for detector(s): %s", detectors.IsNull()?"none":detectors.Data());
  return iResult;
}
