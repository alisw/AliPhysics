#ifndef ALIQUARKONIAACCEPTANCE_H
#define ALIQUARKONIAACCEPTANCE_H

/* $Id$ */
//===================================================================
//  class AliQUARKONIAACCEPTANCE                               
//  This class will provide the quarkonia decay acceptance in ALICE
//  for different resonances :
//    kJpsi
//    kPsiP
//    kUpsilon
//    kUpsilonP
//    kUpsilonPP
//  and for some vector mesons :
//    kPhi
//    kOmega
//  and different channels
//    kDimuon
//    kDielectron
//
//  Acceptance for the Dimuon channel is defined with respect to 
//  a flat distribution of quarkonia emited in the rapidity range 
//  -4 < y < -2.5. Acceptance is defined as both muon from the 
//  decay to be in the theta range 171. < theta < 178.
//                                                              
//   Gines MARTINEZ, Subatech, May 06   
//===================================================================
#include "TNamed.h"
class TH2F;
class TString;


class AliQuarkoniaAcceptance : public TNamed
{
 public:
    
  enum quarkonia{kJpsi, kPsiP, kUpsilon, kUpsilonP, kUpsilonPP, kOmega, kPhi};
  enum channel{kDimuon, kDielectron};
  
  AliQuarkoniaAcceptance(Int_t quarkoniaResonance=kJpsi, Int_t decayChannel=kDimuon);
  virtual ~AliQuarkoniaAcceptance();  
  void   Init(); 
  TH2F*  GetAcceptanceHisto() const;
  void   GetAcceptance(Float_t rap, Float_t pT, Double_t & accep, Double_t & error); 
  void   SetAcceptanceFileName(char * acceptanceFileName) { fAcceptanceFileName = acceptanceFileName; }
  void   SetQuarkoniaResonance(Int_t quarkoniaResonance=kJpsi) { fQuarkoniaResonance= quarkoniaResonance;}
  void   SetDecayChannel(Int_t decayChannel=kDimuon) { fDecayChannel = decayChannel;}
  
 protected: 
  AliQuarkoniaAcceptance(const AliQuarkoniaAcceptance& rhs);
  AliQuarkoniaAcceptance& operator=(const AliQuarkoniaAcceptance& rhs);
  
  TString        fAcceptanceFileName;      // Name of the acceptance root file
  Int_t          fQuarkoniaResonance;      // Resonance Acceptance
  Int_t          fDecayChannel;            // Studied decay channel
  TH2F *         fAcceptance;              // Acceptance histogram
  
 private:  
  
  ClassDef(AliQuarkoniaAcceptance,1)
    };
#endif

