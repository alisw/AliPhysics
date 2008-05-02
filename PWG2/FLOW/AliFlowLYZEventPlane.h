/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliFlowLYZEventPlane_H
#define AliFlowLYZEventPlane_H


#include "TString.h"
#include "TVector2.h"

class AliFlowEventSimple;
class TProfile;
class TFile;

// AliFlowLYZEventPlane:
// Class to calculate the event plane and event weight from the LYZ method
// author: N. van der Kolk (kolk@nikhef.nl)

class AliFlowLYZEventPlane {
 public:
  AliFlowLYZEventPlane();
  virtual ~AliFlowLYZEventPlane();

  void Init();
  void CalculateRPandW(TVector2 fQ);

  Double_t GetWR() const  {return this->fWR; }
  Double_t GetPsi() const {return this->fPsi; }
  
  // input files
  void	   SetFirstRunFileName(TString name) 	
    { this->fFirstRunFileName = name ; }      // Sets input file name
  TString  GetFirstRunFileName() const		
    { return this->fFirstRunFileName ; }      // Gets output file name
  void     SetFirstRunFile(TFile* file)         
    { this->fFirstRunFile = file ; }          // Sets first run file

  void	   SetSecondRunFileName(TString name) 	
    { this->fSecondRunFileName = name ; }     // Sets input file name
  TString  GetSecondRunFileName() const		
    { return this->fSecondRunFileName ; }     // Gets output file name
  void     SetSecondRunFile(TFile* file)         
    { this->fSecondRunFile = file ; }         // Sets first run file



 private:
  
  TVector2 GetQ(AliFlowEventSimple* fEvent);

  TFile*   fFirstRunFile ;          //! pointer to file from first run
  TFile*   fSecondRunFile ;         //! pointer to file from second run
  TString  fFirstRunFileName;       //!
  TString  fSecondRunFileName;      //!


  TVector2 fQ;             // flow vector
  Double_t fWR;            // event weight
  Double_t fPsi;           // reaction plane

  TProfile* fSecondReDtheta; //!
  TProfile* fSecondImDtheta; //!
  TProfile* fFirstr0theta;   //!

  ClassDef(AliFlowLYZEventPlane, 0);          
};

#endif
