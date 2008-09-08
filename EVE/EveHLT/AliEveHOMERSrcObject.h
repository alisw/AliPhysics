// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

//-*- Mode: C++ -*-
#ifndef ALIEVEHOMERSRCOBJECT_H
#define ALIEVEHOMERSRCOBJECT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliEveHOMERSrcObject.h
    @author Jochen Thaeder
    @date
    @brief  Src Object for Src Mapping
*/

#include "TString.h"
#include "TObject.h"

class AliEveHOMERSrcObject : public TObject
{
 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** constructor */
  AliEveHOMERSrcObject( TString dataType, TString className, ULong_t specification );
    
  /** destructor */
  virtual ~AliEveHOMERSrcObject();

  /*
   * ---------------------------------------------------------------------------------
   *                                Getter - public
   * ---------------------------------------------------------------------------------
   */

  /** Returns the HLT dataType */
  TString GetDataType() { return fDataType; } // Returns the HLT dataType

  /** Returns the HLT className */
  TString GetClassName() { return fClassName; } // Returns the HLT className

  /** Returns the HLT specification in HLT */
  ULong_t GetSpecification() { return fSpecification; } // Returns the HLT specification in HLT

  ///////////////////////////////////////////////////////////////////////////////////

private:

  AliEveHOMERSrcObject(const AliEveHOMERSrcObject&);            // Not implemented.
  AliEveHOMERSrcObject& operator=(const AliEveHOMERSrcObject&); // Not implemented.

  /*
   * ---------------------------------------------------------------------------------
   *                            Members - private
   * ---------------------------------------------------------------------------------
   */

  TString fDataType;   // Contains the HLT DataType
  TString fClassName;  // Contains the Classname in HLT
  ULong_t fSpecification; // Contains the HLT Specification
  
  ClassDef(AliEveHOMERSrcObject, 0); // Object for mapping of HLT data-sources.
};

#endif


