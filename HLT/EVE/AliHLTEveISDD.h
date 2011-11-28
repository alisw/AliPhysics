//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTEVEISDD_H
#define ALIHLTEVEISDD_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTEveISDD.h
/// @author Svein Lindal
/// @brief  ISDD Instance of Eve display processor

#include "AliHLTEveITS.h"


class AliHLTEveISDD : public AliHLTEveITS {

public:
  
  /** Constructor  **/
  AliHLTEveISDD();

  /** Destructor **/
 ~AliHLTEveISDD();

private:

  /** copy constructor prohibited */
  AliHLTEveISDD(const AliHLTEveISDD&);
  /** assignment operator prohibited */
  AliHLTEveISDD& operator = (const AliHLTEveISDD& );

  /** Inherited from AliHLTEveITS */
  void SetUpPointSet(TEvePointSet* ps);
  
  ClassDef(AliHLTEveISDD, 0);
};

#endif
