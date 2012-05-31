//-*- Mode: C++ -*-
// $Id: AliHLTJETConeHeader.cxx  $
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
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

/** @file   AliHLTJETConeHeader.cxx
    @author Jochen Thaeder
    @date   
    @brief  Header of the cone finder
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTJETConeHeader.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETConeHeader)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */
  
// #################################################################################
AliHLTJETConeHeader::AliHLTJETConeHeader()
  : 
  AliJetHeader("AliHLTJETConeHeader"),
  fJetCuts(NULL),
  fUseLeading(kFALSE) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

// #################################################################################
AliHLTJETConeHeader::~AliHLTJETConeHeader() {
  // see header file for class documentation

}
 
/*
 * ---------------------------------------------------------------------------------
 *                                    Initialize
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTJETConeHeader::Initialize() {
  // see header file for class documentation

  HLTInfo(" -= AliHLTJETHeader =- " );

  return 0;
}
