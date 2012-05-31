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

/** @file   testAliHLTRootSchemaEvolutionComponent.C
    @author Matthias Richter
    @date   
    @brief  Test macro/program for the AliHLTRootSchemaEvolutionComponent
 */
int testAliHLTRootSchemaEvolutionComponent(const char* rawfile, int firstevent=0, int lastevent=-1) 
{
  AliRawReader* rawreader=AliRawReader::Create(rawfile);
  if (!rawreader) {
    cerr << "error: can not open raw file " << rawfile << endl;
    return -1;
  }

  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();
  pHLT->SetGlobalLoggingLevel(kHLTLogAll);
  pHLT->LoadComponentLibraries("libAliHLTUtil.so");

  AliHLTConfiguration publisher("hltout-publisher", "AliHLTOUTPublisher" , NULL, "");
  AliHLTConfiguration esdfilter("esdfilter", "BlockFilter" , "hltout-publisher", "-datatype ALIESDV0 'HLT ' -prescalar 2 -skip-events 3");  
  AliHLTConfiguration objfilter("objfilter", "BlockFilter" , "hltout-publisher", "-datatype ROOTTOBJ 'HLT '");  
  AliHLTConfiguration streamerinfo("streamerinfo", "ROOTSchemaEvolutionComponent" , "esdfilter objfilter", "");
  AliHLTConfiguration writer("writer", "FileWriter" , "streamerinfo", "-write-all-events -directory streamerinfo");

  pHLT->BuildTaskList("writer");
  AliHLTOfflineInterface::SetParamsToComponents(NULL, rawreader);

  if (!rawreader->GotoEvent(firstevent)) {
    cerr << "error: can not access event " << firstevent << endl;
    return 0;
  }
  int event=firstevent;
  do {
    pHLT->Run(1,0);
  } while ((++event<=lastevent || lastevent<0) && rawreader->NextEvent());
  pHLT->Run(-1);
}
