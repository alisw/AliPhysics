// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


 
#include "AliHLTPHOSTreeMakerComponent.h"
#include "AliHLTPHOSTreeMaker.h"
#include "AliHLTPHOSProcessor.h"
#include "AliHLTPHOSDigitDataStruct.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TObject.h"
#include <fstream>
#include "TFile.h"
#include <sys/stat.h>
#include <sys/types.h>

/**
 * Tree maker component for PHOS HLT
 *
 * @file   AliHLTPHOSTreeMakerComponent.cxx
 * @author Oystein Djuvsland
 * @date   
 * @brief  A tree maker component for PHOS HLT
*/


// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

const AliHLTComponentDataType AliHLTPHOSTreeMakerComponent::fgkInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}};

AliHLTPHOSTreeMakerComponent gAliHLTPHOSTreeMakerComponent;

AliHLTPHOSTreeMakerComponent::AliHLTPHOSTreeMakerComponent() :
  AliHLTPHOSProcessor(),
  fDigitTreePtr(0),
  fTreeMakerPtr(0),
  fDirectory(0),
  fWriteInterval(1000)
{
  //See header file for documentation
}

AliHLTPHOSTreeMakerComponent::~AliHLTPHOSTreeMakerComponent()
{
  //See header file for documentation
}

int 
AliHLTPHOSTreeMakerComponent::Deinit()
{
  //See header file for documentation
  //  cout << "Writing file...";
  char filename [50];

  sprintf(filename, "%s/run%d_digitTree_%d.root", fDirectory, fRunNumber,(fPhosEventCount/fWriteInterval));
 
  TFile *outfile = new TFile(filename,"recreate");
  fDigitTreePtr->Write();
  delete outfile;
  outfile = 0;
  //  cout << "Done!\n";
  if(fDigitTreePtr) 
    {
      delete fDigitTreePtr;
      fDigitTreePtr = 0;
    }
  return 0;
}



const char*
AliHLTPHOSTreeMakerComponent::GetComponentID()
{
  //See header file for documentation
  return "PhosTreeMaker";
}

void
AliHLTPHOSTreeMakerComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
{ 
  //See header file for documentation
 //Get datatypes for input
  const AliHLTComponentDataType* pType=fgkInputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType); 
    pType++;
  }
}

AliHLTComponentDataType 
AliHLTPHOSTreeMakerComponent::GetOutputDataType()
{
  //See header file for documentation
  return AliHLTPHOSDefinitions::fgkAliHLTRootTreeDataType;
}

void 
AliHLTPHOSTreeMakerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  //See header file for documentation
  constBase = 30;
  inputMultiplier = 1;
}

int 
AliHLTPHOSTreeMakerComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
					AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* /*outputPtr*/, AliHLTUInt32_t& /*size*/,  //TODO: I think size should be set to zero when returning from this method if not data was written to the output buffer.
					std::vector<AliHLTComponentBlockData>& /*outputBlocks*/)

s{
  //Do event
  //See header file for documentation
  Bool_t digitEvent;
  Int_t nDigits = 0;
  Int_t totalDigits = 0;

  const AliHLTComponentBlockData* iter = 0;
  unsigned long ndx;

  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks + ndx;

      if ( iter->fDataType == AliHLTPHOSDefinitions::fgkAliHLTDigitDataType )

        {
          digitEvent = true;
          nDigits  = fTreeMakerPtr->MakeDigitArray ( reinterpret_cast<AliHLTPHOSDigitContainerDataStruct*> ( iter->fPtr ), totalDigits );
          totalDigits += nDigits;
	  //cout << totalDigits << endl;
          continue;
        }
      if ( iter->fDataType == AliHLTPHOSDefinitions::fgkAliHLTClusterDataType )
        {
          //
        }
    }

  fPhosEventCount++;

  fTreeMakerPtr->FillDigitTree();
  
  if( (fPhosEventCount%fWriteInterval == 0 ) && (fPhosEventCount != 0))
    {
      Write();
      ResetTrees();
    }

return 0;

}

int
AliHLTPHOSTreeMakerComponent::DoInit ( int argc, const char** argv )
{
  //See header file for documentation
  fTreeMakerPtr = new AliHLTPHOSTreeMaker();
  fDigitTreePtr = new TTree ( "digitTree", "Digits tree" );
  fDirectory = new char[50];

  for ( int i = 0; i < argc; i++ )
    {
      if ( !strcmp ( "-path", argv[i] ) )
        {
          strcpy ( fDirectory, argv[i+1] );
        }
      if ( !strcmp ( "-writeinterval", argv[i] ) )
	{
	  fWriteInterval = atoi(argv[i+1]);
	}
    }

  fTreeMakerPtr->SetDigitTree(fDigitTreePtr);
  //  cout << endl << "Run number is: " <<  fRunNumber  << "  -- Check that this is correct!!!\n" << endl;

  // fRunNumber

  return 0;
  
}


AliHLTComponent*
AliHLTPHOSTreeMakerComponent::Spawn()
{
  //See header file for documentation
  return new AliHLTPHOSTreeMakerComponent();
}

void
AliHLTPHOSTreeMakerComponent::Write()
{
  //See header file for documentation
  //  cout << "Writing file...";
  char filename [256];

  sprintf(filename, "%s/run%d_digitTree_%d.root", fDirectory, fRunNumber,(fPhosEventCount/fWriteInterval - 1));

  TFile *outfile = new TFile(filename,"recreate");
  fDigitTreePtr->Write();
  delete outfile;
  outfile = 0;
  //  cout << "Done!\n";
}

void
AliHLTPHOSTreeMakerComponent::ResetTrees()
{
  //See header file for documentation
  delete fDigitTreePtr;
  fDigitTreePtr = new TTree("digitTree", "Digits tree");
  fTreeMakerPtr->SetDigitTree(fDigitTreePtr);
}
    
  
  
  
  
  



