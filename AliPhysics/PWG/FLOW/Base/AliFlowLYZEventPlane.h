/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliFlowLYZEventPlane_H
#define AliFlowLYZEventPlane_H

// AliFlowLYZEventPlane:
// Class to calculate the event plane and event weight from the LYZ method
// author: N. van der Kolk (kolk@nikhef.nl)

#include "TString.h"
#include "AliFlowVector.h"

class AliFlowEventSimple;
class TProfile;
class TFile;
class TList;

class AliFlowLYZEventPlane {
 public:
  AliFlowLYZEventPlane();
  virtual ~AliFlowLYZEventPlane();

  void Init();
  void CalculateRPandW(AliFlowVector aQ);

  Double_t GetWR() const  {return this->fWR; }
  Double_t GetPsi() const {return this->fPsi; }
  
  //input
  void       SetSecondRunList(TList* list) { this->fSecondRunList = list; }
  TList*     GetSecondRunList()            { return this->fSecondRunList; }
  
 private:
  
  AliFlowLYZEventPlane(const AliFlowLYZEventPlane& aAnalysis);             // copy constructor
  AliFlowLYZEventPlane& operator=(const AliFlowLYZEventPlane& aAnalysis);  // assignment operator
  
  TList*   fSecondRunList;   // list from Second LYZ run output
  Double_t fWR;              // event weight
  Double_t fPsi;             // reaction plane

  TProfile* fSecondReDtheta; // holds Re of Dtheta
  TProfile* fSecondImDtheta; // holds Im of Dtheta
  TProfile* fFirstr0theta;   // holds r0(theta)

  ClassDef(AliFlowLYZEventPlane, 0);          
};

#endif

