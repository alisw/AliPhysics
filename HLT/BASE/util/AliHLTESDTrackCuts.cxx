//-*- Mode: C++ -*-
// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Jochen Thaeder <jochen@thaeder.de>                    *
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

/// @file   AliHLTESDTrackCuts.cxx
/// @author Jochen Thaeder <jochen@thaeder.de>
/// @brief  ESD track cuts used in the analysis of HLT data
///

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTESDTrackCuts.h"
#include "AliESDtrack.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTESDTrackCuts)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
AliHLTESDTrackCuts:: AliHLTESDTrackCuts(const Char_t* name, const Char_t* title) 
  : AliESDtrackCuts(name,title) {  
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  
}

//##################################################################################
AliHLTESDTrackCuts::~AliHLTESDTrackCuts() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 *                                     Selection
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Bool_t AliHLTESDTrackCuts::IsSelected(TObject* obj) {
  // see header file for class documentation

  AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*>(obj);
  if (!esdTrack)
    return kFALSE;

  return AcceptTrack(esdTrack);
}

/*
 * ---------------------------------------------------------------------------------
 *                           Standard Track Cut Definitions
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
AliHLTESDTrackCuts* AliHLTESDTrackCuts::GetStandardTrackCuts2010pp() {
  // see header file for class documentation

  //
  // !!! Be aware - this is not the final yet
  //

  // -- HLT adopted track cuts
  AliHLTESDTrackCuts* esdTrackCuts = new AliHLTESDTrackCuts("ALiHLTESDTrackCuts","HLT standard track cuts 2010 for pp");

  // -- turn off criteria
  esdTrackCuts->SetDCAToVertex2D(kFALSE);
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kFALSE);
  esdTrackCuts->SetRequireITSRefit(kFALSE);

  // -- CleanSample
  esdTrackCuts->SetMaxDCAToVertexXY(3.0);
  esdTrackCuts->SetMaxDCAToVertexZ(3.0);
  esdTrackCuts->SetEtaRange(-0.9,0.9);
  esdTrackCuts->SetMinNClustersTPC(60);

  // -- CleanSample Pt
  esdTrackCuts->SetPtRange(0.3,200.);         
  
  return esdTrackCuts;
}
