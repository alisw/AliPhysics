#ifndef ALIRSNVALUEEVENT_H
#define ALIRSNVALUEEVENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
////////////////////////////////////////////////////////////////////////////////
//
//  Values which depend on 4-momentum of the pair.
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRsnValue.h"

class AliRsnValueEvent : public AliRsnValue {
public:

   enum EType {
      kLeadingPt,       // transverse momentum of the event leading particle
      kMult,            // multiplicity computed as the number of tracks
      kMultMC,          // multiplicity from MC
      kMultESDCuts,     // multiplicity of good quality tracks
      kMultSPD,         // multiplicity from SPD
      kVz,              // Z position of event primary vertex
      kCentralityV0,    // event centrality (V0 method)
      kCentralityTrack, // event centrality (tracks method)
      kCentralityCL1,   // event centrality (CL1 method)
      kTypes            
   };

   AliRsnValueEvent(const char *name = "valEvent", EType type = kTypes);
   AliRsnValueEvent(const AliRsnValueEvent& copy);
   AliRsnValueEvent& operator=(const AliRsnValueEvent& copy);
   virtual ~AliRsnValueEvent() { }

   void             SetType(EType type)  {fType = type;}
   EType            GetType()     const  {return fType;}
   const char*      GetTypeName() const;

   virtual Bool_t   Eval(TObject *object);

protected:

   EType           fType;         //  type from enumeration

   ClassDef(AliRsnValueEvent, 1)  // AliRsnValueEvent class
};

#endif
