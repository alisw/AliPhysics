/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>                *
 *          for The ALICE Off-line Project.                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTCOMPHuffmanAltroCalibComponent.cxx
    @author Jochen Thaeder
    @date   
    @brief  A calibration component for the Huffman code creation.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTCOMPHuffmanAltroCalibComponent.h"
#include "AliHLTCOMPHuffmanAltro.h"
#include "AliHLTCompDefinitions.h"
#include "AliHLTStdIncludes.h"
#include "AliHLTReadoutList.h"
#include "TFile.h" // necessary for HuffmanData writing

ClassImp(AliHLTCOMPHuffmanAltroCalibComponent)

AliHLTCOMPHuffmanAltroCalibComponent::AliHLTCOMPHuffmanAltroCalibComponent()
  :
  fHuffmanCompressor(NULL),
  fHuffmanData(NULL),
  fOrigin(kAliHLTVoidDataOrigin),
  fRunNumber(0),
  fSpecification(0),
  fTablePath(),
  fNRCUTrailerWords(0) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTCOMPHuffmanAltroCalibComponent::~AliHLTCOMPHuffmanAltroCalibComponent() {
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTCOMPHuffmanAltroCalibComponent::GetComponentID() {
  // see header file for class documentation

  return "COMPHuffmanTrainer";
}

void AliHLTCOMPHuffmanAltroCalibComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.clear(); 
  list.push_back( kAliHLTDataTypeDDLRaw );
}

AliHLTComponentDataType AliHLTCOMPHuffmanAltroCalibComponent::GetOutputDataType() {
  // see header file for class documentation
  AliHLTComponentDataType dt=AliHLTCompDefinitions::fgkHuffmanAltroCalDataType;
  if (!fOrigin.IsNull()) dt=dt|fOrigin.Data();
  return dt;
 
}

void AliHLTCOMPHuffmanAltroCalibComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
  // see header file for class documentation
  constBase = sizeof(AliHLTCOMPHuffmanData);
  inputMultiplier = (0.0);
}

AliHLTComponent* AliHLTCOMPHuffmanAltroCalibComponent::Spawn() {
  // see header file for class documentation

  return new AliHLTCOMPHuffmanAltroCalibComponent();
}  


Int_t AliHLTCOMPHuffmanAltroCalibComponent::ScanArgument( Int_t argc, const char** argv ) {
  // see header file for class documentation

  Int_t iResult = 0;
  TString argument = " ";
  TString parameter = " ";

  if ( !argc ) 
    return -EINVAL;

  argument = argv[iResult];
  
  if ( argument.IsNull() ) 
    return -EINVAL;

  // get data origin
  if ( argument.CompareTo("-origin") == 0 ) 
    {

      if ( ++iResult >= argc  ) 
	{
	  iResult = -EPROTO;
	}
      else 
	{
	  fOrigin = argv[1];

	  while(fOrigin.Length() <  kAliHLTComponentDataTypefOriginSize)
	    {
	      fOrigin.Append(" ");
	    }
	
	  HLTInfo("Origin is set to %s.", fOrigin.Data());	  

	}
    }

  else
    {
      // get run number   
      if ( argument.CompareTo("-runnumber") == 0 ) 
	{
	  
	if ( ++iResult >= argc  ) 
	  {
	    iResult = -EPROTO;
	  }
	else 
	  {
	    parameter = argv[1];
	    
	    // get run number
	 
	    fRunNumber =  atoi(parameter.Data());
	       
	    HLTInfo( "Run number is set to %d (Dec) = %X (Hex).", fRunNumber, fRunNumber ); 
	  }
	}
      else
	{
	  // get data specification
	  if(argument.CompareTo("-dataspec") == 0 ) 
	    {
	      if ( ++iResult >= argc  ) 
		{
		  iResult = -EPROTO;
		}
	      else 
		{
		  // get data specification
		  fSpecification = strtoul( argv[1], NULL, 16);
		  
		  HLTInfo( "Specification is set to %d (Dec) = %08X (Hex).", fSpecification, fSpecification ); 
		}
	    }
	  // get number of trailer words
	  else 
	    {
	      
	      if ( argument.CompareTo("-tablepath") == 0)
		{
		  if ( ++iResult >= argc  ) 
		    {
		      iResult = -EPROTO;
		    }
		  else 
		    {
		      // get table path
		      fTablePath = argv[1];
		      HLTInfo( "Path for Huffman table output is set to %s.", fTablePath.Data() ); 
		      if (!fTablePath.IsNull() && !fTablePath.EndsWith("/"))
			fTablePath+="/";		      
		    }
		}

	      else 
		{
		  if ( argument.CompareTo("-trailerwords") == 0 ) 
		    { 
		      
		      if ( ++iResult >= argc  ) 
			{
			  iResult = -EPROTO;
			}
		      else 
			{
			  parameter = argv[1];
			  if ( parameter.CompareTo("1") == 0 ) 
			    {
			      fNRCUTrailerWords = 1;
			      HLTInfo( "Number of trailer words is set to 1." );
			    }
			  else if ( parameter.CompareTo("2") == 0 ) 
			    {
			      fNRCUTrailerWords = 2;
			  HLTInfo( "Number of trailer words is set to 2." );
			    }
			  else if ( parameter.CompareTo("3") == 0 ) 
			    {
			      fNRCUTrailerWords = 3;
			  HLTInfo( "Number of trailer words is set to 3." );
			    }
			  else 
			    {
			      HLTError( "Invalid number of trailerwords: '%s'.", argv[1] );
			      iResult = -EPROTO;
			    }
			} 
		    }
		  else 
		    {
		      iResult = -EINVAL;
		    }
		}
	    }
	}
    }
  
  return iResult;
}

Int_t AliHLTCOMPHuffmanAltroCalibComponent::InitCalibration() {
  // see header file for class documentation
    
  // ** Create a calibration instance to train the Huffman code table
  if ( fHuffmanCompressor )
    return EINPROGRESS;
  
  if ( fHuffmanData )
    return EINPROGRESS;

  // create a new instance of HuffmanData to write results from training in
  fHuffmanData = new AliHLTCOMPHuffmanData();

  fHuffmanCompressor = new AliHLTCOMPHuffmanAltro(kTRUE, kTRUE, NULL, fNRCUTrailerWords);

  // initialise new training table
  fHuffmanCompressor->InitNewTrainingTable(); 
 
  return 0;
}

Int_t AliHLTCOMPHuffmanAltroCalibComponent::DeinitCalibration() {
  // see header file for class documentation

  if ( fHuffmanCompressor )
    delete fHuffmanCompressor; 
  fHuffmanCompressor = NULL;

  if ( fHuffmanData )
    delete fHuffmanData; 
  fHuffmanData = NULL;

  return 0;
}

/** function to do the calibration */
Int_t AliHLTCOMPHuffmanAltroCalibComponent::ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& /*trigData*/ ) {
  // see header file for class documentation
 
  if (evtData.fEventID==0) {
    // this is only to avoid missing parameter warning when compiling for non
    // debug. The parameter is used in the HLTDebug message only.
  }

  const AliHLTComponentBlockData* iter = NULL;

  //AliHLTUInt8_t slice, patch;
 
  // ** Loop over all input blocks and specify which data format should be read - only select Raw Data
  iter = GetFirstInputBlock( kAliHLTDataTypeDDLRaw );

  if ( iter != NULL ) do {
    
    // ** Print Debug output which data format was received
    HLTDebug ( "Event received - Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s", 
	       evtData.fEventID, evtData.fEventID, DataType2Text(iter->fDataType).c_str(), DataType2Text(kAliHLTDataTypeDDLRaw).c_str());

    
    TString blockorigin("");
    blockorigin.Insert(0, iter->fDataType.fOrigin, kAliHLTComponentDataTypefOriginSize);
    
    if (fOrigin.IsNull())
      {
	// if origin is not explicitly set by command line, take origin from data block
	fOrigin=blockorigin;
	HLTDebug("Origin of current data block set by block itself is %s.", blockorigin.Data());
      }
    else
      {
	// if set origin is not equal to current block origin, printout warning! 
	if(fOrigin.CompareTo(blockorigin)!=0) {
	  HLTWarning("Origin %s of current data block does not match origin set by command line argument %s.", blockorigin.Data(), fOrigin.Data());
	  continue;
	}
      }

    // ** Get DDL ID in order to tell the memory reader which slice/patch to use
    //fSlice = AliHLTCompDefinitions::GetMinSliceNr( *iter );
    //fPatch = AliHLTCompDefinitions::GetMinPatchNr( *iter );

    //HLTDebug ( "Input Raw Data - Slice/Patch: %d/%d.", fSlice, fPatch);

    // FIXME: set ddl no
    fHuffmanCompressor->AddInputData(reinterpret_cast<UChar_t*>(iter->fPtr), iter->fSize, 768);

    // only necessary for output in binary file
    //fHuffmanCompressor->SetSlice(fSlice);
    //fHuffmanCompressor->SetPatch(fPatch);
   
    fHuffmanCompressor->ProcessData();
  
    // ** Get next input block, with the same specification as defined in GetFirstInputBlock()
  } while ( (iter = GetNextInputBlock()) != NULL );

  // ** Get output specification
  // commented out for the moment to read spec in from command line argument
  //fSpecification = AliHLTCompDefinitions::EncodeDataSpecification( fSlice, fSlice, fPatch, fPatch );
  //fSpecification = fSlice<<24 | fSlice<<16 | fPatch<<8 | fPatch;

  // ** PushBack data to shared memory ... 

  // DATA TYE to DEFINE !!! XXXX
  PushBack( (TObject*) fHuffmanData, AliHLTCompDefinitions::fgkHuffmanAltroCalDataType|fOrigin.Data(), fSpecification);
 
  return 0;
} // Int_t AliHLTCOMPHuffmanAltroCalibComponent::ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData ) {


Int_t AliHLTCOMPHuffmanAltroCalibComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/ ) {
  // see header file for class documentation

  // create code from training table
  fHuffmanCompressor->CreateCodeTable();

  // write code table and occurrence table to HuffmanData instance
  fHuffmanCompressor->SetHuffmanData(fHuffmanData);

  TString rootfilename;
  if(fTablePath.IsNull() )
    {
      // if there is no explicit table path, take current path
      rootfilename.Form("huffmanData_%s_%08X_%08X.root", fOrigin.Data(), (unsigned)fRunNumber, (unsigned)fSpecification);      
    }
  else
    {
      rootfilename.Form("%shuffmanData_%s_%08X_%08X.root", fTablePath.Data(), fOrigin.Data(), (unsigned)fRunNumber, (unsigned)fSpecification);
    }
 
  TFile* huffmanrootfile = new TFile(rootfilename, "RECREATE");
  huffmanrootfile->WriteObject(fHuffmanData,"HuffmanData");
  huffmanrootfile->Write();
  huffmanrootfile->Close();

  // ** PushBack data to FXS ...
  // currently specification has to be put in by command line argument!
  Int_t dataspec = (Int_t) fSpecification;

  fHuffmanData->SetOCDBSpecifications(fOrigin, dataspec);
  static AliHLTReadoutList rdList(AliHLTReadoutList::kTPC);
  PushToFXS( (TObject*) fHuffmanData, "TPC", "HuffmanData", &rdList ) ;
  
  return 0;
} // Int_t AliHLTCOMPHuffmanAltroCalibComponent::ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData ) {
