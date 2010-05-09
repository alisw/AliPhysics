/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Base class fro anlyzing EMCAL raww data
// Further documentation found in base class
// --------------
// --------------
// --------------
// --------------



#include "AliHLTEMCALRawAnalyzerComponent.h"
#include "AliHLTEMCALMapper.h"
#include "AliHLTEMCALDefinitions.h"
#include "AliHLTCaloChannelDataHeaderStruct.h"
//#include "unistd.h"


#include  "TStopwatch.h"
TStopwatch  fgWatch; //CRAP PTH


AliHLTEMCALRawAnalyzerComponent::AliHLTEMCALRawAnalyzerComponent() : AliHLTCaloRawAnalyzerComponentv3("EMCAL")
{
  
   cout << __FILE__ << __FUNCTION__ << __LINE__ <<  endl;
}


AliHLTEMCALRawAnalyzerComponent::~AliHLTEMCALRawAnalyzerComponent()
{

}



void 
AliHLTEMCALRawAnalyzerComponent::GetInputDataTypes( vector <AliHLTComponentDataType>& list)
{
  list.clear();
  list.push_back( AliHLTEMCALDefinitions::fgkDDLRawDataType   | kAliHLTDataOriginEMCAL );
}



AliHLTComponentDataType
AliHLTEMCALRawAnalyzerComponent::GetOutputDataType()
{
  //comment
  return AliHLTEMCALDefinitions::fgkChannelDataType;
}


void
AliHLTEMCALRawAnalyzerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )
{
  //comment
  constBase = sizeof(AliHLTCaloChannelDataHeaderStruct);
  inputMultiplier = 0.5;
}


void 
AliHLTEMCALRawAnalyzerComponent::DoInit() 
{
  cout << __FILE__ << __FUNCTION__ << __LINE__ <<  endl;

  //  fgWatch.Start();
 
}

/*
struct AliHLTComponentDataType
  {
    AliHLTUInt32_t fStructSize;                            /// Size of this structure in bytes.
    char fID[kAliHLTComponentDataTypefIDsize];             /// Data type identifier.
    char fOrigin[kAliHLTComponentDataTypefOriginSize];     /// Subsystem or detector origin of the data.
  };
*/

bool 
AliHLTEMCALRawAnalyzerComponent::CheckInputDataType(const AliHLTComponentDataType &datatype)
{
  // Cheking if datatype is the correct one before processing 
  //    ////cout << __FILE__ << __LINE__ << "  :  fID  = " << datatype.fID <<  " : fOrigin = " <<  datatype.fOrigin << endl;
  //  ////cout << __FILE__ << __LINE__ << "fgkDDLRawDataType->fID = " << AliHLTEMCALDefinitions::fgkDDLRawDataType.fID <<
  //  "fgkDDLRawDataType->fOrigin = " << AliHLTEMCALDefinitions::fgkDDLRawDataType.fOrigin << endl;


  if ( datatype  == AliHLTEMCALDefinitions::fgkDDLRawDataType  )
     {
       return true;
     }
   else
     {
       //    return true;
       return false;
     }
}


void 
AliHLTEMCALRawAnalyzerComponent::InitMapping( const int specification )
{
  //-------------
  if ( fMapperPtr == 0 )
    {
      fMapperPtr =  new   AliHLTEMCALMapper( specification );
    }

  if(fMapperPtr->GetIsInitializedMapping() == false )
    {
      //      HLTError("%d:%d, ERROR, mapping not initialized ", __FILE__, __LINE__ );
      exit(-2);
    }
}


int 
AliHLTEMCALRawAnalyzerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& /*trigData*/, 
					 AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  
  //  cout << __FILE__ << __FUNCTION__ << __LINE__ <<  endl;
  
   if(!IsDataEvent())
   {
      size = 0;
      return 0;
   }

  static int evntcnt = 0;
  static double wlast = -1;
  static double wcurrent = 0;
  evntcnt  ++;
  
  if( evntcnt %1000 == 0  )
    {
      ////cout << __FILE__ << __LINE__ << " : Processing event "  << evntcnt   << endl; 
      wlast =  wcurrent;
      wcurrent = fgWatch.RealTime();
      ////cout << wlast << ":" << wcurrent << endl;
      cout << __FILE__ << __LINE__ << "The event rate is " <<  1000/( wcurrent  -  wlast ) << "  Hz" << endl; 
      fgWatch.Start(kFALSE); 
      //     wlast =  fgWatch.RealTime(); 
    }
  
  /*
  if( evntcnt %100 == 0  )
    {
      
      ////cout << __FILE__ << __LINE__ << " : Processing event "  << evntcnt   << endl; 
      wlast =  wcurrent;
      wcurrent = fgWatch.RealTime();
      ////cout << wlast << ":" << wcurrent << endl;
      ////cout << __FILE__ << __LINE__ << "The event rate is " <<  100/( wcurrent  -  wlast ) << "  Hz" << endl; 
      fgWatch.Start(kFALSE); 
      //     wlast =  fgWatch.RealTime(); 
    }
  */
  

  Int_t blockSize          = -1;
  UInt_t totSize           = 0;
  const AliHLTComponentBlockData* iter = NULL; 
  unsigned long ndx;

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      //   //cout << __FILE__ << __LINE__ <<  "ndx =" << ndx  << endl;
      iter = blocks+ndx;
      if(  ! CheckInputDataType(iter->fDataType) )
	{
	  //	  //cout << __FILE__ << __LINE__ <<  "  continue" << endl; 
	  continue;
	}
      else
	{
	  //	  //cout << __FILE__ << __LINE__ <<  "  else" << endl;  

	  InitMapping( iter->fSpecification); 
	  blockSize = DoIt(iter, outputPtr, size, totSize); // Processing the block
	  
	  if(blockSize == -1) // If the processing returns -1 we are out of buffer and return an error msg.
	    {
	      //	     //cout << __FILE__ << __LINE__ <<  "  return -ENOBUFS " << endl;   
	      return -ENOBUFS;
	    }
	  
	  totSize += blockSize; //Keeping track of the used size
	  AliHLTComponentBlockData bdChannelData;
	  FillBlockData( bdChannelData );
	  bdChannelData.fOffset = 0; //FIXME
	  bdChannelData.fSize = blockSize;
	  
	  //	  bdChannelData.fDataType = AliHLTPHOSDefinitions::fgkChannelDataType;
	  bdChannelData.fDataType = AliHLTEMCALDefinitions::fgkChannelDataType;

	  bdChannelData.fSpecification = iter->fSpecification;
	  outputBlocks.push_back(bdChannelData);
	  outputPtr += blockSize; //Updating position of the output buffer
	}

      fCaloEventCount++; 
      size = totSize; //telling the framework how much buffer space we have used.
    }

  
return 0;
  
}//end DoEvent

