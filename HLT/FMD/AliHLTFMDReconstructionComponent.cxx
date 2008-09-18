
#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTSystem.h"
#include "AliRawReaderMemory.h"
#include "AliFMDRawReader.h"
#include "AliESDFMD.h"
#include "AliHLTDataTypes.h"
#include "TClonesArray.h"
#include "AliHLTFMDReconstructionComponent.h"
#include "AliHLTDefinitions.h"
#include "AliCDBManager.h"
#include "AliFMDParameters.h"
#include "AliHLTDataTypes.h"
#include "TH1F.h"
#include "AliLog.h"
#include "AliGeomManager.h"
#include "AliFMDGeometry.h"
#include <cstdlib>
#include <cerrno>

ClassImp(AliHLTFMDReconstructionComponent)

//_____________________________________________________________________

AliHLTFMDReconstructionComponent::AliHLTFMDReconstructionComponent()
  : fRunNumber(0)
{
  
}

//_____________________________________________________________________

AliHLTFMDReconstructionComponent::~AliHLTFMDReconstructionComponent()
{
   // see header file for class documentation
}

//_____________________________________________________________________
AliHLTFMDReconstructionComponent::AliHLTFMDReconstructor::AliHLTFMDReconstructor() 
  : AliFMDReconstructor()
{
  fDiagnostics = kFALSE;
  
}

//_____________________________________________________________________
void AliHLTFMDReconstructionComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
      list.push_back(kAliHLTDataTypeDDLRaw);
      
}

//______________________________________________________________________
AliHLTComponentDataType AliHLTFMDReconstructionComponent::GetOutputDataType()
{
      return kAliHLTDataTypeESDTree;
}

//______________________________________________________________________
void AliHLTFMDReconstructionComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  
  constBase = 2000000;
  inputMultiplier = 1;
}

//______________________________________________________________________

AliHLTComponent* AliHLTFMDReconstructionComponent::Spawn()
{
  
  return new AliHLTFMDReconstructionComponent;
}

//______________________________________________________________________
int AliHLTFMDReconstructionComponent::DoInit( int argc, const char** argv )
{
  char* cpErr;
  Int_t i = 0;
  while ( i < argc )
    {
      if ( !strcmp( argv[i], "run_number" ) ||
	   !strcmp( argv[i], "-run_number" ))
	{
	  if ( i+1>=argc )
	    {
	      HLTError("Missing run_number parameter");
	      return -EINVAL;
	    }
	  fRunNumber = strtoul( argv[i+1], &cpErr, 0 );
	  if ( *cpErr )
	    {
	      HLTError("Cannot convert run_number parameter '%s'", argv[i+1] );
	      return -EINVAL;
	    }
	  HLTInfo("Run Number set to %lu %%", fRunNumber );
	  i += 2;
	  continue;
	}
      HLTError("Unknown option '%s'", argv[i] );
      return -EINVAL;
    }
  
  AliLog::EnableDebug(kFALSE);

  if (GetRunNo() == kAliHLTVoidRunNo) {
    AliCDBManager* cdb = AliCDBManager::Instance();
    cdb->SetRun(fRunNumber);
  }
  AliFMDParameters::Instance()->Init();
  
  if (AliGeomManager::GetGeometry() == NULL)
    AliGeomManager::LoadGeometry();
  
  AliFMDGeometry* geo = AliFMDGeometry::Instance();
  geo->Init();
  geo->InitTransformations();
  
  return 0;
}

//______________________________________________________________________
int AliHLTFMDReconstructionComponent::DoDeinit()
{
  
  
  return 0;
}

//______________________________________________________________________
int AliHLTFMDReconstructionComponent::DoEvent( const AliHLTComponentEventData& evtData, 
					       const AliHLTComponentBlockData* blocks, 
					       AliHLTComponentTriggerData& /*trigData*/, 
					       AliHLTUInt8_t* /*outputPtr*/, 
					       AliHLTUInt32_t& size, 
					       vector<AliHLTComponentBlockData>& /*outputBlocks*/ )
{
  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )){
    size=0;
    return 0;
  }
  // -- Iterator over Data Blocks --
  const AliHLTComponentBlockData* iter = NULL;
  
  AliESDEvent  esd;
  esd.CreateStdContent();
    
  // Process an event
  //unsigned long totalSize = 0;
  // Loop over all input blocks in the event
  
  AliHLTFMDReconstructor rec;
  rec.Init();
  for ( unsigned long n = 0; n < evtData.fBlockCnt; n++ ) {   
      
    iter = blocks + n;
    
    // -- Check for the correct data type
    if ( iter->fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginFMD) )
      continue;
    
    AliRawReaderMemory* reader = new AliRawReaderMemory();
    
    reader->SetMemory((UChar_t*) iter->fPtr, iter->fSize );
    
    AliHLTUInt32_t spec = iter->fSpecification;
    
    Int_t id = 3072;
    for ( Int_t ii = 0; ii < 24 ; ii++ ) {
      if ( spec & 0x00000001 ) {
	id += ii;
	break;
      }
      spec = spec >> 1 ;
    }
    
    reader->SetEquipmentID( id );
    
    AliFMDRawReader* fmdReader  = new AliFMDRawReader(reader,0);
    TClonesArray*    digitArray = new TClonesArray("AliFMDDigit",0);
    
    digitArray->Clear();
    fmdReader->ReadAdcs(digitArray);
    rec.ReconstructDigits(digitArray);
    
  }
  
  Int_t iResult = 0;
  
  esd.SetFMDData(rec.GetFMDData());
  
  iResult=PushBack(&esd, kAliHLTDataTypeESDObject|kAliHLTDataOriginFMD, 0);
  
  return 0;
}

