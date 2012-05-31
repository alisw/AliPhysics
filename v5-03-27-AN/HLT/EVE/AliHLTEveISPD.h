//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTEVEISPD_H
#define ALIHLTEVEISPD_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTEveISPD.h
/// @author Svein Lindal
/// @brief  SPD Instance of Eve display processor

#include "AliHLTEveITS.h"
class TEvePointSet;

class AliHLTEveISPD : public AliHLTEveITS {

public:
  
  /** Constructor  **/
  AliHLTEveISPD();

  /** Destructor **/
 ~AliHLTEveISPD();

private:

  /** copy constructor prohibited */
  AliHLTEveISPD(const AliHLTEveISPD&);
  /** assignment operator prohibited */
  AliHLTEveISPD& operator = (const AliHLTEveISPD& );

  /** Inherited from AliHLTEveITS */
  void SetUpPointSet(TEvePointSet * ps);
  
  ClassDef(AliHLTEveISPD, 0);
};

#endif
