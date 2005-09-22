#ifndef ALIFASTMUONTRIGGEREFF_H
#define ALIFASTMUONTRIGGEREFF_H
/*  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <AliFastResponse.h>
#include <TString.h>
#include <TObjString.h>
#include <TFile.h>
#include <TH3.h>
#include <TROOT.h>
#include <stdlib.h>

// Debugging flag
//#define MYTRIGDEBUG

enum CutTupe {kLow, kHigh, kAny};

class AliFastMuonTriggerEff : public AliFastResponse {
    
 public:
    AliFastMuonTriggerEff();
    AliFastMuonTriggerEff(char* Name, char* Title) {;}
    virtual ~AliFastMuonTriggerEff(){;}
    virtual void    Init();
    virtual void    Evaluate(Float_t charge, Float_t pt, Float_t theta, Float_t phi,
			     Float_t& effLow, Float_t& effHigh, Float_t& effAny);
    virtual void    Evaluate(Float_t charge, Float_t pt, Float_t theta, Float_t phi,
			     Float_t& effLow, Float_t& effHigh);
    virtual Float_t Evaluate(Float_t charge, Float_t pt, Float_t theta, Float_t phi);
    virtual void    SetCut(Int_t cut = kLow);
    virtual Float_t Cut() {return fCut;}
    virtual Int_t   SetBkgLevel(Float_t Bkg=0.);
    virtual Int_t   ForceBkgLevel(Float_t Bkg=0.);
    virtual Float_t GetBkgLevel() {return fBkg;}
    Int_t LoadTables(Char_t *namet);  // Load trigger response tables
    void SetInt() {fInt=1;}
    void UnsetInt() {fInt=0;}
    Int_t GetInt() {return fInt;}
  protected:
    Double_t fPtMin;
    Double_t fPtMax;
    Double_t fDpt;                 // Delta_pt
    Int_t    fnptb;
    Double_t fPhiMin;              // lower limit for phi 
    Double_t fPhiMax;              // upper limit for phi
    Double_t fDphi;                // Delta_phi
    Int_t    fnphib;
    Double_t fThetaMin;            // lower limit for theta
    Double_t fThetaMax;            // upper limit for theta
    Double_t fDtheta;              // Delta_theta
    Int_t    fnthetab;
    Int_t   fCut;                 // Cut type (low/high)
    Int_t   fZones;               // Total number of zones
    TH3F*   fhEffAPt;             // Trig. prob. for Any Pt
    TH3F*   fhEffLPt;             // Trig. prob. for Low Pt
    TH3F*   fhEffHPt;             // Trig. prob. for High Pt
    TAxis*  fhLX;
    TAxis*  fhLY;
    TAxis*  fhLZ;
    Float_t fBkg;                 // Background level
    TString fTableTitle;          // Title of the LUT
    TString fDescription;         // Description of the table content
    Int_t fInt;                   // Interpolation flag (1 to interpolate)
  private:
    Int_t fibx;
    Int_t fiby;
    Int_t fibz;
    ClassDef(AliFastMuonTriggerEff,1)    // Fast Muon Trigger response
};

#endif 



