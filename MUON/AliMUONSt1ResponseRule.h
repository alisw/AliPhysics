#ifndef ALI_MUON_ST1_RESPONSE_RULE_H
#define ALI_MUON_ST1_RESPONSE_RULE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay
//
// Class AliMUONSt1ResponseRule
// -----------------------------
// Describes a response rule.
// A "rule" is defined as being a set of electronic filters to be applied 
// (ie. a set of AliMUONSt1ResponseParameter) and a set of cathode pads to 
// which these filters should be applied (set of AliMUONSt1ElectronicElement)

#include <TObject.h>
#include <TList.h>

class AliMpPad;

class AliMUONSt1ElectronicElement;
class AliMUONSt1ResponseParameter;

class AliMUONSt1ResponseRule : public TObject 
{
  public:
    AliMUONSt1ResponseRule();
    virtual ~AliMUONSt1ResponseRule();
  
    void   AddElement(AliMUONSt1ElectronicElement* element);
    void   AddParameter(AliMUONSt1ResponseParameter* param);
    Bool_t Contains(const AliMpPad& pad) const;
    TList* GetParameters() {return &fParameters;}

  private:
    TList  fElementList;// list of electronic elements to which this rule is applied
    TList  fParameters; // parameters for this rule

  ClassDef(AliMUONSt1ResponseRule,1) // A set of electronic elements and the linked electronic parameters 
};

#endif //ALI_MUON_ST1_RESPONSE_RULE_H
