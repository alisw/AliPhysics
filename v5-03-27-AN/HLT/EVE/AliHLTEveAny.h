//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTEVEANY_H
#define ALIHLTEVEANY_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTEveAny.h
/// @author Svein Lindal <slindal@fys.uio.no>
/// @date   
/// @brief  
///

#include "AliHLTEveBase.h"
class AliHLTHOMERBlockDesc;

class AliHLTEveAny : public AliHLTEveBase {

public:
  
  /** Constructor  **/
  AliHLTEveAny();

  /** Destructor **/
 ~AliHLTEveAny();

  /** Inherited form AliHLTEveBase */
  void ProcessBlock(AliHLTHOMERBlockDesc * block);

  /** inherited from AliHLTEveBase */
  void UpdateElements();
  
  /** inherited from AliHLTEveBase */
  void ResetElements();

private:
  
  /** copy constructor prohibited */
  AliHLTEveAny(const AliHLTEveAny&);
  /** assignment operator prohibited */
  AliHLTEveAny& operator = (const AliHLTEveAny &);

  void ProcessHistogram(AliHLTHOMERBlockDesc * block );

  ClassDef(AliHLTEveAny, 0);
};

#endif
