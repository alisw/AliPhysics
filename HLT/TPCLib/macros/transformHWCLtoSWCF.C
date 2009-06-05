
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>           *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** 
 * @file   transformHWCLtoSWCL.C
 * @author Kalliopi.Kanaki@ift.uib.no
 * @date   
 * @brief  Test macro for the AliHLTTPCHWClusterTransformComponent.cxx
 * 
 *
 * This macro tests the component that does the tranformation from the HW 
 * cluster format to the software cluster structure. 
 * For the macro to work, the user needs to use the FilePublisher component
 * to make the input files available. The data specification has to be given
 * and the input files have to be present.
 *
 * In order for the reconstruction to work properly at the absence of raw/ folders,
 * a directory has to be specified as the input argument of the function
 * and an empty raw0 folder has to be created in it.
 *
 * The component runs without arguments.
 * 
 */

void transformHWCLtoSWCL(const char* input="./"){
 
  if (!input) {
    cerr << "please specify input or run without arguments" << endl;
    return;
  }
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // 
  // init the HLT system in order to define the analysis chain below
  //
  gSystem->Load("libHLTrec.so");
  AliHLTSystem* gHLT=AliHLTReconstructorBase::GetInstance();

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
    
  int iMinSlice = 0; 
  int iMaxSlice = 0;
  int iMinPart  = 0;
  int iMaxPart  = 0;

  TString argument;

  argument.Form("-datatype 'HWCLUST1' 'TPC '  -datafile /scratch/FCF/PulserData_Limit_0xFF.fcf -dataspec 0x%02x%02x%02x%02x", iMinSlice, iMaxSlice, iMinPart, iMaxPart);
  
  AliHLTConfiguration pubconf("FP", "FilePublisher", NULL, argument.Data());
  
  AliHLTConfiguration hwconf("FCF", "TPCHWClusterTransform", "FP", "");
  
  AliHLTConfiguration clusDumpconf("sink1", "TPCClusterDump", "FCF", "-directory ClusterDump");
  
   
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstruction is switched off
  //
  AliReconstruction rec;
  rec.SetInput(input);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunLocalReconstruction("HLT");
  rec.SetRunTracking("");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");
  rec.SetFillESD("");
  rec.SetFillTriggerESD(false);
  rec.SetOption("HLT", "libAliHLTTPC.so loglevel=0x7c chains=sink1");
  rec.Run();
}
