/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//////////////////////////////////////////////////////////////////////////
// AliFlowEventSimpleMaker:
// Class to fill the AliFlowEventSimple with AliFlowTrackSimple objects
// author: N. van der Kolk (kolk@nikhef.nl),  Ante Bilandzic (anteb@nikhef.nl),  
// Raimond Snellings (Raimond.Snellings@nikhef.nl) 
//////////////////////////////////////////////////////////////////////////


#ifndef ALIFLOWEVENTSIMPLEMAKER_H
#define ALIFLOWEVENTSIMPLEMAKER_H

class AliFlowEventSimple;
class AliFlowTrackSimpleCuts;
class TTree;
class AliCFManager;
class AliMCEvent;
class AliESDEvent;
class AliAODEvent;
          
class AliFlowEventSimpleMaker {

 public:

  AliFlowEventSimpleMaker();             //constructor
  virtual ~AliFlowEventSimpleMaker();    //destructor

  void SetMCReactionPlaneAngle(Double_t fPhiRP)  { this->fMCReactionPlaneAngle = fPhiRP; } 
  //TTree
  AliFlowEventSimple* FillTracks(TTree* anInput, const AliFlowTrackSimpleCuts* rpCuts, const AliFlowTrackSimpleCuts* poiCuts);   //use own cut class
  //AliMCEvent
  AliFlowEventSimple* FillTracks(AliMCEvent* anInput);   //use own cuts
  AliFlowEventSimple* FillTracks(AliMCEvent* anInput, const AliCFManager* rpCFManager, const AliCFManager* poiCFManager ); //use CF(2x)
  //AliESDEvent
  AliFlowEventSimple* FillTracks(AliESDEvent* anInput);   //use own cuts
  AliFlowEventSimple* FillTracks(AliESDEvent* anInput,  const AliCFManager* rpCFManager, const AliCFManager* poiCFManager); //use CF(2x)
  //AliESDEvent & AliMCEvent
  AliFlowEventSimple* FillTracks(AliESDEvent* anInput, const AliMCEvent* anInputMc, Int_t anOption);  //use own cuts
  AliFlowEventSimple* FillTracks(AliESDEvent* anInput, const AliMCEvent* anInputMc, const AliCFManager* rpCFManager, const AliCFManager* poiCFManager, Int_t anOption);  //use CF(2x)
  // anOption = 0 : kine from ESD
  // anOption = 1 : kine from MC
  //AliAODEvent
  AliFlowEventSimple* FillTracks(AliAODEvent* anInput); //use own cuts
  AliFlowEventSimple* FillTracks(AliAODEvent* anInput, const AliCFManager* rpCFManager, const AliCFManager* poiCFManager);  //use CF(2x)
  
  void  SetNoOfLoops(Int_t noofl) {this->fNoOfLoops = noofl;}
  Int_t GetNoOfLoops() const      {return this->fNoOfLoops;} 
  
  void     SetEllipticFlowValue(Double_t elfv) {this->fEllipticFlowValue = elfv;}
  Double_t GetEllipticFlowValue() const        {return this->fEllipticFlowValue;} 
  
  void  SetMultiplicityOfEvent(Int_t multevnt) {this->fMultiplicityOfEvent = multevnt;}
  Int_t GetMultiplicityOfEvent() const         {return this->fMultiplicityOfEvent;} 
  
  void  SetMinMult(Int_t multmin) {this->fMinMult = multmin; }
  Int_t GetMinMult() const        {return this->fMinMult; }
  void  SetMaxMult(Int_t multmax) {this->fMaxMult = multmax; }
  Int_t GetMaxMult() const        {return this->fMaxMult; }

  void SetSubeventEtaRange(Double_t minA,Double_t maxA,Double_t minB,Double_t maxB) 
    {this->fEtaMinA = minA; this->fEtaMaxA = maxA;this->fEtaMinB = minB; this->fEtaMaxB = maxB;};

 private:
  AliFlowEventSimpleMaker(const AliFlowEventSimpleMaker& anAnalysis);            //copy constructor
  AliFlowEventSimpleMaker& operator=(const AliFlowEventSimpleMaker& anAnalysis); //assignment operator

  Double_t  fMCReactionPlaneAngle; // the angle of the reaction plane from the MC truth
  Int_t     fCount;                // counter for the number of events processed

  Int_t     fNoOfLoops;            // number of times to use the same particle (nonflow) 
  Double_t  fEllipticFlowValue;    // Add Flow. Must be in range [0,1].
  Int_t     fMultiplicityOfEvent;  // Set maximal multiplicity (when manipulating the event).

  Int_t     fMinMult;              // Minimum multiplicity from tracks selected using CORRFW
  Int_t     fMaxMult;              // Maximum multiplicity from tracks selected using CORRFW

  Double_t  fEtaMinA;              // minimum eta of subevent A eta range
  Double_t  fEtaMaxA;              // maximum eta of subevent A eta range
  Double_t  fEtaMinB;              // minimum eta of subevent B eta range
  Double_t  fEtaMaxB;              // maximum eta of subevent B eta range

  ClassDef(AliFlowEventSimpleMaker,1)    // macro for rootcint
};
 
     
#endif

