// @(#) $Id$

/** @file   AliHLTAgentSample.cxx
    @author Matthias Richter
    @date   
    @brief  Agent of the libAliHLTSample library
*/

#include "AliHLTAgentSample.h"
#include "AliHLTConfiguration.h"
#include "TSystem.h"

/** global instance for agent registration */
AliHLTAgentSample gAliHLTAgentSample;

const char* gAliHLTAgentSampleData="/tmp/testdata";
const char* gAliHLTAgentSampleOut="/tmp/hltout";

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTAgentSample)

AliHLTAgentSample::AliHLTAgentSample()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTAgentSample::~AliHLTAgentSample()
{
  // see header file for class documentation

  // delete the test data
  ofstream dump(gAliHLTAgentSampleData, ios::in);
  if (dump.good()) {
    TString arg("rm -f ");
    arg+=gAliHLTAgentSampleData;
    gSystem->Exec(arg.Data());
  }
}

int AliHLTAgentSample::CreateConfigurations(AliHLTConfigurationHandler* handler,
					  AliRunLoader* runloader) const
{
  // see header file for class documentation

  // create some test data
  ofstream dump(gAliHLTAgentSampleData, (ios::openmode)0);
  dump << "This is just some test data for the ALICE HLT analysis example";
  dump << "---- not copied" << endl;
  dump.close();

  if (handler) {
    // the publisher configuration for the test data
    TString arg("-datafile ");
    arg+=gAliHLTAgentSampleData;
    HLTDebug(arg.Data());
    handler->CreateConfiguration("fp1"  , "FilePublisher", NULL , arg.Data());

    // the configuration for the dummy component
    handler->CreateConfiguration("cp"   , "Dummy"        , "fp1", "output_percentage 80");

    // the writer configuration
    arg="-datafile "; arg+=gAliHLTAgentSampleOut;
    handler->CreateConfiguration("sink1", "FileWriter"   , "cp" , arg.Data());
  }
  return 0;
}

const char* AliHLTAgentSample::GetTopConfigurations(AliRunLoader* runloader) const
{
  // see header file for class documentation
  return "sink1";
}

const char* AliHLTAgentSample::GetRequiredComponentLibraries() const
{
  // see header file for class documentation
  return "libAliHLTUtil.so libAliHLTSample.so";
}
