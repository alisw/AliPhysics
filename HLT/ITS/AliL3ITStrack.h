#ifndef ALIL3ITSTRACK_H
#define ALIL3ITSTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                  High Level Trigger ITS Track Class
//
//        Origin: Cvetan Cheshkov, CERN, Cvetan.Cheshkov@cern.ch 
//-------------------------------------------------------------------------


/*****************************************************************************
 *                          October 11, 2004                                 *
 * The class inherits all the properties of the off-line AliITStrackV2 class *
 * and in addition it contains an interface to the HLT ESD track             *
 *****************************************************************************/

#include <AliITStrackV2.h>

class AliESDHLTtrack;

class AliL3ITStrack : public AliITStrackV2 {
public:
  AliL3ITStrack();
  AliL3ITStrack(const AliL3ITStrack& t);
  AliL3ITStrack(AliESDHLTtrack& t, Double_t zvertex) throw (const Char_t *);
  AliL3ITStrack(const AliESDHLTtrack& t, Double_t zvertex) throw (const Char_t *);

  Int_t Compare(const TObject *o) const;

  // Set and get the pointer to the HLT ESD track
  AliESDHLTtrack *GetESDHLTtrack() const {return fESDHLTtrack; }
  void SetESDHLTtrack(AliESDHLTtrack *esdhlttrack) { fESDHLTtrack = esdhlttrack; }
  Bool_t GetPxPyPzAt(Double_t x,Double_t *p) const;

protected:
  void Set(const AliESDHLTtrack& t, Double_t zvertex) throw (const Char_t *);

  AliESDHLTtrack *fESDHLTtrack;   //! pointer to the connected ESD HLT track

  ClassDef(AliL3ITStrack,1)   //HLT ITS reconstructed track
};

#endif
