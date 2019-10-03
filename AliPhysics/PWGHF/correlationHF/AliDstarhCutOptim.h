#ifndef AliDstarhCutOptim_h
#define AliDstarhCutOptim_h

//
//  Class to crate a TTree for variable optimization
//
//-----------------------------------------------------------------------
//  Author M. Mazzilli
//  INFN Bari
//  marianna.mazzilli@cern.ch
//-----------------------------------------------------------------------

#include <iostream>
#include "TObject.h"
#include "TMath.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"

class AliDstarhCutOptim : public TObject
{
  public:
    AliDstarhCutOptim(); // default constructor
    Float_t  deltaInvMass;
    Float_t  cent;
    Float_t  pt;
    Float_t  dca;
    Float_t  cosThSt;
    Float_t  pTk;
    Float_t  pTpi;
    Float_t  d0k; 
    Float_t  d0pi;
    Float_t  d0xd0;
    Float_t  cosThPt;
    Float_t  normLxy;
   // Float_t  topom;
    
    ClassDef(AliDstarhCutOptim,1);
}; 

#endif
