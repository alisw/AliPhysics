/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

  
#ifndef AliFlowLeeYangZerosMaker_H
#define AliFlowLeeYangZerosMaker_H

//#include <iostream>
//using namespace std; //required for resolving the 'cout' symbol

//#include "Riostream.h"
#include "TVector2.h"          //called explicitly
#include "AliFlowSelection.h"  //only with AliFlowSelection* called ???
#include "AliFlowTrack.h"      //only with AliFlowTrack* called ???

//class AliFlowTrack;          //why doesn't compile?
//class AliFlowSelection;      //why doesn't compile?

class AliFlowEvent;
class AliFlowLYZHist1; 
class AliFlowLYZHist2;
//class AliFlowLYZSummary;     //class not yet finished

class TH1F;
class TH1D;
class TProfile;
class TProfile2D;
class TObjArray;
class TFile;
class TTree;
class TComplex;
class Riostream;


// Description: Maker to analyze Flow by the LeeYangZeros method
//              One needs to do two runs over the data; 
//              First to calculate the integrated flow 
//              and in the second to calculate the differential flow
 
class AliFlowLeeYangZerosMaker {

 public:
 
   AliFlowLeeYangZerosMaker();                                            //default constructor
 
   AliFlowLeeYangZerosMaker(const AliFlowSelection* fFlowSelect);          //custum constructor

   AliFlowLeeYangZerosMaker(const AliFlowLeeYangZerosMaker&);              //copy constuctor (dummy)

   AliFlowLeeYangZerosMaker& operator=(const AliFlowLeeYangZerosMaker&);   //assignment operator (dummy)
   
   virtual  ~AliFlowLeeYangZerosMaker();                                   //destructor
 
   Bool_t    Init();    //defines variables and histograms
   Bool_t    Make(AliFlowEvent* fFlowEvent);    //calculates variables and fills histograms
   Bool_t    Finish();  //saves histograms

   Double_t  GetQtheta(AliFlowSelection* fFlowSelect, TObjArray* fFlowTracks, Float_t fTheta);

   void      SetFirstRun(Bool_t kt)              { this->fFirstRun = kt ; }
   Bool_t    GetFirstRun() const                 { return this->fFirstRun ; }
   void      SetUseSum(Bool_t kt)                { this->fUseSum = kt ; }
   Bool_t    GetUseSum() const                   { return this->fUseSum ; }
   void      SetDebug(Bool_t kt)                 { this->fDebug = kt ; }
   Bool_t    GetDebug() const                    { return this->fDebug ; }


   // Output 
   void	    SetHistFileName(TString name) 	{ this->fHistFileName = name ; } // Sets output file name
   TString  GetHistFileName() const		{ return this->fHistFileName ; } // Gets output file name
   TFile*   GetHistFile() const                 { return this->fHistFile ; }     // Gets output file
  
   // input for second run
   void	    SetFirstRunFileName(TString name) 	{ this->firstRunFileName = name ; } // Sets input file name
   TString  GetFirstRunFileName() const		{ return this->firstRunFileName ; } // Gets output file name
   void     SetFirstRunFile(TFile* file)        { this->firstRunFile = file ; }        // Sets first run file

 private:
   Bool_t   MakeControlHistograms(AliFlowEvent* fFlowEvent); 
   Bool_t   FillFromFlowEvent(AliFlowEvent* fFlowEvent);
   Bool_t   SecondFillFromFlowEvent(AliFlowEvent* fFlowEvent);

   //Double_t GetQtheta(AliFlowSelection* fFlowSelect, TObjArray* fFlowTracks, Float_t fTheta) const;
   TComplex GetGrtheta(AliFlowSelection* fFlowSelect, Float_t fR, Float_t fTheta);
   TComplex GetDiffFlow(AliFlowSelection* fFlowSelect, Float_t fR, Int_t theta); 
   Float_t  GetR0(TH1D* fHistGtheta);

  
   
#ifndef __CINT__
   TVector2  fQ;            // flow vector
   TVector2  fQsum;         // flow vector sum
   Float_t   fQ2sum;        // flow vector sum squared
   Double_t  fQtheta;       // flow vector projected on ref. angle theta
   
#endif /*__CINT__*/

   Int_t     fEventNumber;  // event counter
   Int_t     fMult;         // multiplicity
   Int_t     fNbins;        // number of bins
   Float_t   fTheta;        // ref. angle theta
   
   AliFlowEvent*      fFlowEvent;    //! pointer to AliFlowEvent
   AliFlowSelection*  fFlowSelect;   //! pointer to AliFlowSelection
   TObjArray*         fFlowTracks ;  //! pointer to the TrackCollection
   AliFlowTrack*      fFlowTrack ;   //! pointer to the AliFlowTrack
   //AliFlowLYZSummary* fLYZSummary;   //! 
   TTree*             fFlowTree;     //!

   Bool_t       fFirstRun ;         //! flag for lyz analysis: true=first run over data, false=second run 
   Bool_t       fUseSum ;           //! flag for lyz analysis: true=use sum gen.function, false=use product gen.function
   Bool_t       fDebug ;            //! flag for lyz analysis: more print statements

   TString      fHistFileName;      //!
   TFile*       fHistFile;          //!
   TFile*       fSummaryFile;       //!
   TDirectory * fFistLoop;          //!
   TString      firstRunFileName;   //!
   TFile*       firstRunFile;       //!
  // for single histograms
  TH1F*        fHistOrigMult;      //!
  TH1F*        fHistMult;          //!
  TH1F*        fHistQ;             //!
  TH1F*        fHistPt;            //!
  TH1F*        fHistEta;           //!
  TH1F*        fHistPhi;           //!
  TH1F*        fHistQtheta;        //!
   

  TProfile*    fHistProVthetaHar1;     //!
  TProfile*    fHistProVthetaHar2;     //!
  TProfile*    fHistProVetaHar1;       //!
  TProfile*    fHistProVetaHar2;       //!
  TProfile*    fHistProVPtHar1;        //!
  TProfile*    fHistProVPtHar2;        //!
  TH1D*        fHistVPtHar2;           //!
  TProfile*    fHistProVR0;            //!
  TH1D*        fHistVR0;               //!
  TProfile*    fHistProV;              //!
  TProfile*    fHistProR0thetaHar1;    //!
  TProfile*    fHistProR0thetaHar2;    //!
  TProfile*    fHistProReDenomHar1;    //!
  TProfile*    fHistProReDenomHar2;    //!
  TProfile*    fHistProImDenomHar1;    //!
  TProfile*    fHistProImDenomHar2;    //!
  TProfile*    fHistProReDtheta;       //!
  TProfile*    fHistProImDtheta;       //!
   
   
  //class AliFlowLYZHist1 defines the histograms: fHistProGtheta, fHistProReGtheta, fHistProImGtheta
  AliFlowLYZHist1* fHist1[2][5];       //!

  //class AliFlowLYZHist1 defines the histograms: fHistProReNumer, fHistProImNumer, fHistProReNumerPt,
  //fHistProImNumerPt, fHistProReNumer2D, fHistProImNumer2D.
  AliFlowLYZHist2* fHist2[2][5];       //!
 
  ClassDef(AliFlowLeeYangZerosMaker,0)  // macro for rootcint
    };
 
     
#endif
