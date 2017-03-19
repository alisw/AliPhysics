#ifndef AliD0hCutOptim_h
#define AliD0hCutOptim_h

//
//  Class to crate a TTree for variable optimization
//
//-----------------------------------------------------------------------
//  Author F.Colamaria
//  INFN Bari
//  fabio.colamaria@cern.ch
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

class AliD0hCutOptim : public TObject
{
  public:
    AliD0hCutOptim(); // default constructor
    Float_t  invMass;
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
    Float_t  topom;
    
    ClassDef(AliD0hCutOptim,1);
}; 

#endif
