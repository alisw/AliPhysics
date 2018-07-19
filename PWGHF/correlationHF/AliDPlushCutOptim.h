#ifndef AliDPlushCutOptim_h
#define AliDPlushCutOptim_h

//
//  Class to crate a TTree for variable optimization
//
//-----------------------------------------------------------------------
//  Author Shyam Kumar : copied from D0 class
//  IIT Bombay
//  shyam.kumar@cern.ch
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

class AliDPlushCutOptim : public TObject
{
  public:
    AliDPlushCutOptim(); // default constructor
    Float_t  invMass;    // "inv. mass [GeV]"
    Float_t  cent;       // Centrality by estimator set in Addtask of Dplus correlation
    Float_t  pt;        // Pt of DPlus meson
    Float_t  pTk; 	//"pTK [GeV/c]"
    Float_t  pTpi;      // "pTPi [GeV/c]"
    Float_t  d0k; 	// "d0K [cm]"
    Float_t  d0pi; 	// "d0Pi [cm]"
    Float_t  cosPA; 	// "cosThetaPoint"
    Float_t  sumd0square; // "Sum d0^2 (cm^2)"
    Float_t  dca;  	// "dca cut (cm)" distance between the tracks 
    Float_t dlXY; 	// "Normalized dec len XY scaled by P/Pt (cm)
    Float_t cosPAXY; 	// "cosThetaPointXY"
    Float_t  topom; 	// nsigma cut for impact parameter
    ClassDef(AliDPlushCutOptim,2);
}; 

#endif
