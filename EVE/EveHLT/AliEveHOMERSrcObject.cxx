// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007
// Author: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>                *
//         for The ALICE HLT Project.                                    *

//-*- Mode: C++ -*-

/** @file   AliEveHOMERSrcObject.cxx
    @author Jochen Thaeder
    @date
    @brief  Src Object for Src Mapping
*/

#if __GNUC__>= 3
   using namespace std;
#endif

#define use_aliroot
#define use_root
#define ROWHOUGHPARAMS
#define use_reconstruction
#define use_newio
#define ROOTVERSION    "unchecked"
#define ALIROOTVERSION "unchecked"
#define __ROOT__
#define USE_ALILOG
#define LINUX

#include "AliEveHOMERSrcObject.h"

//______________________________________________________________________________
//
// Translate HLT data-sources.

ClassImp(AliEveHOMERSrcObject)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
AliEveHOMERSrcObject::AliEveHOMERSrcObject( TString dataType, TString className, ULong_t specification ) :
  fDataType(dataType), 
  fClassName(className), 
  fSpecification(specification) {
  // This Class is a mapping object, which should allow, to 
  // add new sources out of the HLT.
  // This should be only needed by AliEveHOMERSrcTranslator
  
  printf ("%s - %s - %lu\n", fDataType.Data(),fClassName.Data(), fSpecification);

}

//##################################################################################
AliEveHOMERSrcObject::~AliEveHOMERSrcObject() {
  // The destructor

}

