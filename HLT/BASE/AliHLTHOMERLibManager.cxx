// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTHOMERLibManager.cxx
    @author Matthias Richter
    @date   
    @brief  dynamic HLT HOMER reader/writer generation and destruction.   */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cerrno>
#include <cassert>
#include "AliHLTHOMERLibManager.h"
#include "AliHLTHOMERReader.h"
#include "AliHLTHOMERWriter.h"
#include "TString.h"
#include "TSystem.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTHOMERLibManager)

AliHLTHOMERLibManager::AliHLTHOMERLibManager()
  :
  fLibraryStatus(0),
  fFctCreateReaderFromTCPPort(NULL),
  fFctCreateReaderFromTCPPorts(NULL),
  fFctCreateReaderFromBuffer(NULL),
  fFctDeleteReader(NULL),
  fFctCreateWriter(NULL),
  fFctDeleteWriter(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTHOMERLibManager::~AliHLTHOMERLibManager()
{
  // see header file for class documentation
}

AliHLTHOMERReader* AliHLTHOMERLibManager::OpenReader(const char* hostname, unsigned short port )
{
  // see header file for class documentation
  if (fLibraryStatus<0) return NULL;

  if (fLibraryStatus==0) {
    fLibraryStatus=LoadHOMERLibrary();
  }
  
  AliHLTHOMERReader* pReader=NULL;
  if (fFctCreateReaderFromTCPPort!=NULL && (pReader=(((AliHLTHOMERReaderCreateFromTCPPort_t)fFctCreateReaderFromTCPPort)(hostname, port)))==NULL) {
    //HLTError("can not create instance of HOMER reader (function %p)", fFctCreateReaderFromTCPPort);
  }
  
  return pReader;
}

AliHLTHOMERReader* AliHLTHOMERLibManager::OpenReader(unsigned int tcpCnt, const char** hostnames, unsigned short* ports)
{
  // see header file for class documentation
  if (fLibraryStatus<0) return NULL;

  if (fLibraryStatus==0) {
    fLibraryStatus=LoadHOMERLibrary();
  }
  
  AliHLTHOMERReader* pReader=NULL;
  if (fFctCreateReaderFromTCPPorts!=NULL && (pReader=(((AliHLTHOMERReaderCreateFromTCPPorts_t)fFctCreateReaderFromTCPPorts)(tcpCnt, hostnames, ports)))==NULL) {
    //HLTError("can not create instance of HOMER reader (function %p)", fFctCreateReaderFromTCPPorts);
  }
  
  return pReader;
}

AliHLTHOMERReader* AliHLTHOMERLibManager::OpenReaderBuffer(const AliHLTUInt8_t* pBuffer, int size)
{
  // see header file for class documentation
  if (fLibraryStatus<0) return NULL;

  if (fLibraryStatus==0) {
    fLibraryStatus=LoadHOMERLibrary();
  }
  
  AliHLTHOMERReader* pReader=NULL;
  if (fFctCreateReaderFromBuffer!=NULL && (pReader=(((AliHLTHOMERReaderCreateFromBuffer_t)fFctCreateReaderFromBuffer)(pBuffer, size)))==NULL) {
    //HLTError("can not create instance of HOMER reader (function %p)", fFctCreateReaderFromBuffer);
  }
  
  return pReader;
}

int AliHLTHOMERLibManager::DeleteReader(AliHLTHOMERReader* pReader)
{
  // see header file for class documentation
  if (fLibraryStatus<0) return fLibraryStatus;

  if (fLibraryStatus==0) {
    fLibraryStatus=LoadHOMERLibrary();
  }
  
  if (fFctDeleteReader!=NULL) {
    ((AliHLTHOMERReaderDelete_t)fFctDeleteReader)(pReader);
  }
  
  return 0;
}

AliHLTHOMERWriter* AliHLTHOMERLibManager::OpenWriter()
{
  // see header file for class documentation
  if (fLibraryStatus<0) return NULL;

  if (fLibraryStatus==0) {
    fLibraryStatus=LoadHOMERLibrary();
  }
  
  AliHLTHOMERWriter* pWriter=NULL;
  if (fFctCreateWriter!=NULL && (pWriter=(((AliHLTHOMERWriterCreate_t)fFctCreateWriter)()))==NULL) {
//     HLTError("can not create instance of HOMER writer (function %p)", fFctCreateWriter);
  }
  
  return pWriter;
}

int AliHLTHOMERLibManager::DeleteWriter(AliHLTHOMERWriter* pWriter)
{
  // see header file for class documentation
  if (fLibraryStatus<0) return fLibraryStatus;

  if (fLibraryStatus==0) {
    fLibraryStatus=LoadHOMERLibrary();
  }
  
  if (fFctDeleteWriter!=NULL) {
    ((AliHLTHOMERWriterDelete_t)fFctDeleteWriter)(pWriter);
  }
  
  return 0;
}

int AliHLTHOMERLibManager::LoadHOMERLibrary()
{
  // see header file for class documentation
  int iResult=-EBADF;
  const char* libraries[]={"libAliHLTHOMER.so", "libHOMER.so", NULL};
  const char** library=&libraries[0];
  do {
    TString libs = gSystem->GetLibraries();
    if (libs.Contains(*library) ||
	(gSystem->Load(*library)) >= 0) {
      iResult=1;
      break;
    }
  } while (*(++library)!=NULL);

  if (iResult>0 && *library!=NULL) {
    // print compile info
    typedef void (*CompileInfo)( char*& date, char*& time);
    CompileInfo fctInfo=(CompileInfo)gSystem->DynFindSymbol(*library, "CompileInfo");
    if (fctInfo) {
      char* date="";
      char* time="";
      (*fctInfo)(date, time);
      if (!date) date="unknown";
      if (!time) time="unknown";
      //HLTInfo("%s build on %s (%s)", *library, date, time);
    } else {
      //HLTInfo("no build info available for %s", *library);
    }

    fFctCreateReaderFromTCPPort=(void (*)())gSystem->DynFindSymbol(*library, ALIHLTHOMERREADER_CREATE_FROM_TCPPORT);
    fFctCreateReaderFromTCPPorts=(void (*)())gSystem->DynFindSymbol(*library, ALIHLTHOMERREADER_CREATE_FROM_TCPPORTS);
    fFctCreateReaderFromBuffer=(void (*)())gSystem->DynFindSymbol(*library, ALIHLTHOMERREADER_CREATE_FROM_BUFFER);
    fFctDeleteReader=(void (*)())gSystem->DynFindSymbol(*library, ALIHLTHOMERREADER_DELETE);
    fFctCreateWriter=(void (*)())gSystem->DynFindSymbol(*library, ALIHLTHOMERWRITER_CREATE);
    fFctDeleteWriter=(void (*)())gSystem->DynFindSymbol(*library, ALIHLTHOMERWRITER_DELETE);
    if (fFctCreateReaderFromTCPPort==NULL ||
	fFctCreateReaderFromTCPPorts==NULL ||
	fFctCreateReaderFromBuffer==NULL || 
	fFctDeleteReader==NULL ||
	fFctCreateWriter==NULL ||
	fFctDeleteWriter==NULL) {
      iResult=-ENOSYS;
    } else {
    }
  }
  if (iResult<0 || *library==NULL) {
    fFctCreateReaderFromTCPPort=NULL;
    fFctCreateReaderFromTCPPorts=NULL;
    fFctCreateReaderFromBuffer=NULL;
    fFctDeleteReader=NULL;
    fFctCreateWriter=NULL;
    fFctDeleteWriter=NULL;
  }

  return iResult;
}
