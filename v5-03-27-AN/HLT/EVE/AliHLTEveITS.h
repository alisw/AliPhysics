//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTEVEITS_H
#define ALIHLTEVEITS_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTEveITS.h
/// @author Svein Lindal
/// @brief  ITS base class for the Eve display processors

#include "AliHLTEveBase.h"
#include "TString.h"
class TEveElementList;
class TEvePointSet;

class AliHLTEveITS : public AliHLTEveBase {

public :
  
  /** Constructor  **/
  AliHLTEveITS(TString name = TString("ITS"));

  /** Destructor **/
 ~AliHLTEveITS();

  /** Inherited from AliHLTEveBase */
  virtual void ProcessBlock(AliHLTHOMERBlockDesc * block);

  /** Inherited from AliHLTEveBase */
  virtual void UpdateElements();

  /** Inherited from AliHLTEveBase */
  virtual void ResetElements();


protected :

  /** Create new point set */
  TEvePointSet * CreatePointSet(TString name);
  
  /** Process the clusters block */
  void ProcessClusters(AliHLTHOMERBlockDesc * block, TEvePointSet * cont );
  
  /** Set up the look of the pointset, to be overridden in child instances (of one wishes) */
  virtual void SetUpPointSet(TEvePointSet * ps );


  TString fName;  //Detector (ITS, ISSD, ISPD, ISDD)
  TEvePointSet * fPointSet; //The pointset for the display


private :
  
  /** default constructor forbidden */
  //AliHLTEveITS();

  /** copy constructor prohibited */
  AliHLTEveITS(const AliHLTEveITS&);
  /** assignment operator prohibited */
  AliHLTEveITS& operator = (const AliHLTEveITS& );

  ClassDef(AliHLTEveITS, 0);
};

#endif
