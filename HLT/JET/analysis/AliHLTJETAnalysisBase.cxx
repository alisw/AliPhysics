//-*- Mode: C++ -*-
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
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

/** @file   AliHLTJETAnalysisBase.cxx
    @author Jochen Thaeder
    @date   
    @brief  Base functionality for HLT JET analysis package
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
   using namespace std;
#endif

#include "AliHLTJETAnalysisBase.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETAnalysisBase)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
  AliHLTJETAnalysisBase::AliHLTJETAnalysisBase() {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}


//##################################################################################
AliHLTJETAnalysisBase::~AliHLTJETAnalysisBase() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 *                              Initialize static const
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
const Char_t *AliHLTJETAnalysisBase::fgkDeltaType[] = { 
  "all", "leading", 
  "Matched - all", "Matched -Leading" 
};

//##################################################################################
const Char_t *AliHLTJETAnalysisBase::fgkSpectraType[] = { 
  "Pythia",        "Pythia - Matched",        "Pythia - UnMatched",
  "Reconstructed", "Reconstructed - Matched", "Reconstructed - UnMatched",
  "Reco Leading",  "Reco Leading - Matched",  "Reco Leading - UnMatched"
};

//##################################################################################
const Char_t *AliHLTJETAnalysisBase::fgkPlotType[] = { 
  "all", "leading", 
};
