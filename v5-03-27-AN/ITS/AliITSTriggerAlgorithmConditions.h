#ifndef AliITSTriggerAlgorithmConditions_H
#define AliITSTriggerAlgorithmConditions_H

////////////////////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                                         //
//                                                                                //
// Class for storing conditions data from Pixel Trigger (PIT) algorithms.         //
// This holds a sub set of the conditions data needed.                            //
// It is used by AliITSTriggerConditions, which holds all the information.        //
// AliITSTriggerConditions contains a TObjArray of this type.                     //
//                                                                                //
////////////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TString.h>
#include <TObjArray.h>
#include <TArrayI.h>

class AliITSTriggerAlgorithmConditions : public TObject {

 public:
  AliITSTriggerAlgorithmConditions();
  AliITSTriggerAlgorithmConditions(UShort_t id, const Char_t* label, const Char_t* descr);
  AliITSTriggerAlgorithmConditions(const AliITSTriggerAlgorithmConditions& cond);
  virtual ~AliITSTriggerAlgorithmConditions();
  AliITSTriggerAlgorithmConditions& operator=(const AliITSTriggerAlgorithmConditions& cond);

  virtual UShort_t      GetID() const {return fId;}
  virtual const Char_t* GetLabel() const {return fLabel.Data();}
  virtual const Char_t* GetDescription() const {return fDescription.Data();}

  virtual void          SetID(UShort_t id) {fId=id;}
  virtual void          SetLabel(const Char_t* label) {fLabel=label;}
  virtual void          SetDescription(const Char_t* descr) {fDescription=descr;}

  virtual void          ClearParams();

  virtual void          AddParam(const Char_t* name, Int_t value);

  virtual UShort_t      GetNumParam() const {return fNumParam;}
  virtual const Char_t* GetParamNameI(UShort_t index) const;
  virtual Int_t         GetParamValueI(UShort_t index) const;
  virtual Int_t         GetParamValueN(const Char_t* name) const;


 protected:
  UShort_t  fId;          // ID of output (1-10 for real system)
  TString   fLabel;       // label of output (ex: OSH1)
  TString   fDescription; // description of output
  UShort_t  fNumParam;    // Number of parameters used
  TObjArray fParamNames;  // list of parameter names (strings)
  TArrayI   fParamValues; // list of parameter values (integers)
  
  ClassDef(AliITSTriggerAlgorithmConditions,1) // Trigger algorithm conditions class
};

#endif
