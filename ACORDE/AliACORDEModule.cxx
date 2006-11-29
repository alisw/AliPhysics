/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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

/* $Id: */

////////////////////////////////////////////////////////////////////////////
//
// ALICE Cosmic Ray Trigger
//
//  This class will provide the basic utilities to create the geometry of
//  the scintillatio array. This array is basically only the array
//  in the upper face of the magnet. The remaining arrays will be copies
//  of this array.
//
//   Authors:
//
//   Arturo Fernandez <afernand@fcfm.buap.mx>
//   Enrique Gamez    <egamez@fcfm.buap.mx>
//
////////////////////////////////////////////////////////////////////////////

#include "AliACORDEModule.h"

ClassImp(AliACORDEModule)

//_____________________________________________________________________________
AliACORDEModule::AliACORDEModule()
  : TNamed(),
    fScintillatorThickness(0),
    fScintillatorWidth(0),
    fScintillatorLength(0),
    fFrameThickness(0),
    fFrameWidth(0),
    fFrameLength(0),
    fNColumns(0),
    fNRows(0),
    fZGap(0),
    fXGap(0)
{
  //
  // Default constructor
  //
}

//_____________________________________________________________________________
AliACORDEModule::AliACORDEModule(const char* name, const char* title)
  : TNamed(name, title),
    fScintillatorThickness(1),
    fScintillatorWidth(19.7),
    fScintillatorLength(186),
    fFrameThickness(10),
    fFrameWidth(26),
    fFrameLength(300),
    fNColumns(2),
    fNRows(10),
    fZGap(100),
    fXGap(0)
{
  //
  // Standard constructor
  //
}

//_____________________________________________________________________________
AliACORDEModule::AliACORDEModule(const AliACORDEModule& mod)
  : TNamed(mod),
    fScintillatorThickness(mod.fScintillatorThickness),
    fScintillatorWidth(mod.fScintillatorWidth),
    fScintillatorLength(mod.fScintillatorLength),
    fFrameThickness(mod.fFrameThickness),
    fFrameWidth(mod.fFrameWidth),
    fFrameLength(mod.fFrameLength),
    fNColumns(mod.fNColumns),
    fNRows(mod.fNRows),
    fZGap(mod.fZGap),
    fXGap(mod.fXGap)
{
  //
  // Copy constructor
  //
}

//_____________________________________________________________________________
AliACORDEModule::~AliACORDEModule()
{
  //
  // Default destructor
  //
}

//_____________________________________________________________________________
AliACORDEModule& AliACORDEModule::operator=(const AliACORDEModule& mod)
{
  //
  // Asingment operator
  //
  fScintillatorThickness = mod.fScintillatorThickness;
  fScintillatorWidth = mod.fScintillatorWidth;
  fScintillatorLength = mod.fScintillatorLength;
  fFrameThickness = mod.fFrameThickness;
  fFrameWidth = mod.fFrameWidth;
  fFrameLength = mod.fFrameLength;
  fNColumns = mod.fNColumns;
  fNRows = mod.fNRows;
  fZGap = mod.fZGap;
  fXGap = mod.fXGap;
  return *this;
}
