/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Author: Per Thomas Hille for the ALICE HLT Project.                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSFileWriterComponent.h"


AliHLTPHOSFileWriterComponent gAliHLTPHOSFileWriterComponent;

AliHLTPHOSFileWriterComponent::AliHLTPHOSFileWriterComponent():AliHLTFileWriter()
{
  //  cout <<"AliHLTPHOSFileWriterComponent: Creating new FILEWRTER"  << endl;
  
} 


AliHLTPHOSFileWriterComponent::~AliHLTPHOSFileWriterComponent()
{

}

const char* 
AliHLTPHOSFileWriterComponent::GetComponentID()
{
  //  cout <<"AliHLTPHOSFileWriterComponent::GetComponentID(): Returning ID"  << endl; 
  return "PhosFileWriter";
}

AliHLTComponent* AliHLTPHOSFileWriterComponent::Spawn()
{
  cout <<" AliHLTPHOSFileWriterComponent::Spawn()" << endl;
  return new AliHLTPHOSFileWriterComponent;
}

int 
AliHLTPHOSFileWriterComponent::DoInit( int argc, const char** argv )
{
   cout <<" AliHLTPHOSFileWriterComponent::DoInit()" << endl;
   return 0;

}


