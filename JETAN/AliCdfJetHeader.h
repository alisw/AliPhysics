#ifndef ALICDFJETHEADER_H
#define ALICDFJETHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-----------------------------------------------
// PxCone CDF Algorithm Jet finder header class
// Stores the parameters of CDF Jet finder
//-----------------------------------------------

#include "AliJetHeader.h"

class AliCdfJetHeader : public AliJetHeader
  {
  public:

    AliCdfJetHeader();
    virtual ~AliCdfJetHeader() { }

    Double_t GetRadius ()  { return fRadius; }
    Double_t GetPtMin  ()  { return fPtMin ; }
    Double_t GetPtMax  ()  { return fPtMax ; }
    Double_t GetEtaMin ()  { return fEtaMin ; }
    Double_t GetEtaMax ()  { return fEtaMax ; }
    Double_t GetPhiMin ()  { return fPhiMin ; }
    Double_t GetPhiMax ()  { return fPhiMax ; }


    // Setters

    void SetRadius ( Double_t f ) {fRadius = f;}
    void SetPtMin  ( Double_t f ) {fPtMin = f;}
    void SetPtMax  ( Double_t f ) {fPtMax = f;}
    void SetEtaMin ( Double_t f ) {fEtaMin = f;}
    void SetEtaMax ( Double_t f ) {fEtaMax = f;}
    void SetPhiMin ( Double_t f ) {fPhiMin = f;}
    void SetPhiMax ( Double_t f ) {fPhiMax = f;}


    // others

//     void PrintParameters() const ;

  protected:

   AliCdfJetHeader(const AliCdfJetHeader &jh);
	 AliCdfJetHeader& operator=(const AliCdfJetHeader &jh);

   // parameters of algorithm
   Double_t fRadius;      //  Cone radius

   // ranges of Pt,Eta and Phi cut ranges ; by default == 0
   Double_t fPtMin;       // minimum pt
   Double_t fPtMax;       // maximum pt
   Double_t fEtaMin;      // minimum eta
   Double_t fEtaMax;      // maximum eta
   Double_t fPhiMin;      // minimum phi
   Double_t fPhiMax;      // maximum phi

   ClassDef ( AliCdfJetHeader, 1 )

  };
#endif
