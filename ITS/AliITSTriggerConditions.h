#ifndef AliITSTriggerConditions_H
#define AliITSTriggerConditions_H

////////////////////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                                         //
//                                                                                //
// Implementation of conditions data from Pixel Trigger (PIT)                     //
//                                                                                //
// The information is propagated from pixel trigger system to DCS file exchange   //
// server (text file format). The ReadFromTextFile method will populate this      //
// object with the values from the text file. Via a Preprocessor, this object     //
// can be stored in OCDB.                                                         //
//                                                                                //
// One can also manually create conditions data that may be interesting for       //
// simulation.                                                                    //
//                                                                                //
////////////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TObjArray.h>
#include <TString.h>
#include <TBits.h>

class AliITSTriggerConditions : public TObject{
 public:
    AliITSTriggerConditions();
    AliITSTriggerConditions(const AliITSTriggerConditions& cond);
    virtual ~AliITSTriggerConditions();
    AliITSTriggerConditions& operator=(const AliITSTriggerConditions& cond);

    virtual Bool_t        IsEqualTo(AliITSTriggerConditions *cond) const;

    virtual void          DumpAll() const;
    virtual void          ResetAll();

    virtual void          SetRunNumber(UInt_t num) {fRunNumber=num;}
    virtual UInt_t        GetRunNumber() const {return fRunNumber;}
    virtual void          SetFirmWareVersion(UShort_t num) {fFirmWareVersion=num;}
    virtual UShort_t      GetFirmWareVersion() const {return fFirmWareVersion;}
    virtual void          SetGlobalDescription(const Char_t* descr) {fGlobalDescription=descr;}
    virtual const Char_t* GetGlobalDescription() const {return fGlobalDescription.Data();}
    virtual void          SetVersionRegister(UShort_t num) {fVersionRegister=num;}
    virtual UShort_t      GetVersionRegister() const {return fVersionRegister;}
    virtual void          SetInputConditionsVersion(UShort_t num) {fInputConditionsVersion=num;}
    virtual UShort_t      GetInputConditionsVersion() const {return fInputConditionsVersion;}
    virtual void          SetParametersVersion(UShort_t num) {fParametersVersion=num;}
    virtual UShort_t      GetParametersVersion() const {return fParametersVersion;}

    virtual void          SetInActiveChip(UInt_t eq, UInt_t hs, UInt_t chip) 
      {fInActiveChips.SetBitNumber(GetChipKey(eq,hs,chip));}
    virtual void          ResetInActiveChips() {fInActiveChips.ResetAllBits();}
    virtual void          SetActiveChip(UInt_t eq, UInt_t hs, UInt_t chip) 
      {fInActiveChips.SetBitNumber(GetChipKey(eq,hs,chip),kFALSE);}
    virtual void          DumpInActiveChips() const;

    virtual Bool_t        IsChipActive(UInt_t eq, UInt_t hs, UInt_t chip) const 
      {return !IsChipInActive(eq,hs,chip);}
    virtual Bool_t        IsChipInActive(UInt_t eq, UInt_t hs, UInt_t chip) const 
      {return fInActiveChips.TestBitNumber(GetChipKey(eq,hs,chip));}
    virtual Bool_t        GetNextInActiveChip(Int_t& eq, Int_t& hs, Int_t& chip) const;

    virtual void          ClearAlgorithms();
    virtual void          ClearAlgoParamsI(UShort_t aIndex);
    virtual void          ClearAlgoParamsL(const Char_t* aLabel);

    virtual void          AddAlgo(UShort_t id, const Char_t* aLabel, const Char_t* aDescr);
    virtual void          AddAlgoParam(UShort_t id, const Char_t* pName, Int_t pValue);

    virtual UShort_t      GetNumAlgo() const {return fNumAlgo;}
    virtual Short_t       GetAlgoIndexL(const Char_t* aLabel) const;
    virtual Short_t       GetAlgoIDI(UShort_t aIndex) const;
    virtual const Char_t* GetAlgoLabelI(UShort_t aIndex) const;
    virtual const Char_t* GetAlgoDescriptionI(UShort_t aIndex) const;

    virtual Short_t       GetNumAlgoParamI(UShort_t aIndex) const;
    virtual const Char_t* GetAlgoParamNameII(UShort_t aIndex, UShort_t pIndex) const;
    virtual Int_t         GetAlgoParamValueII(UShort_t aIndex, UShort_t pIndex) const;
    virtual Int_t         GetAlgoParamValueIN(UShort_t aIndex, const Char_t* pName) const;
    virtual Short_t       GetNumAlgoParamL(const Char_t* aLabel) const;
    virtual const Char_t* GetAlgoParamNameLI(const Char_t* aLabel, UShort_t pIndex) const;
    virtual Int_t         GetAlgoParamValueLI(const Char_t* aLabel, UShort_t pIndex) const;
    virtual Int_t         GetAlgoParamValueLN(const Char_t* aLabel, const Char_t* pName) const;

    virtual void          ReadFromTextFile(const Char_t* fileName);

 protected:
    UInt_t    fRunNumber;              // Run number
    UShort_t  fFirmWareVersion;        // PIT Processing firmware version
    TString   fGlobalDescription;      // PIT Global description
    UShort_t  fVersionRegister;        // PIT Version register value
    UShort_t  fInputConditionsVersion; // PIT Input configuration version
    UShort_t  fParametersVersion;      // PIT Parameters version
    TBits     fInActiveChips;          // Map of PIT de-activated chips
    UShort_t  fNumAlgo;                // Number of algorithms used
    TObjArray fAlgoList;               // List of conditions for each algorithm used

    UInt_t  GetChipKey(Int_t eq, Int_t hs, Int_t chip) const;
    void    GetChipFromKey(UInt_t key, Int_t& eq, Int_t& hs, Int_t& chip) const;
    Bool_t  SplitStringIn2(TString orig, TString& word1, TString& word2, Char_t sep);
    TString GetStringBetween(TString orig, Char_t sep1, Char_t sep2);

    ClassDef(AliITSTriggerConditions,1) // Trigger conditions class
};

#endif
