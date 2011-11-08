/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// HFE correction framework container
// Contains many single containers
// Extra fuctionality like appending added
//
#ifndef ALIHFECONTAINER_H
#define ALIHFECONTAINER_H

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

#ifndef ROOT_TObjArray
#include <TObjArray.h>
#endif

#ifndef ROOT_TArrayD
#include <TArrayD.h>
#endif

class TArrayF;
template <class X>
class THnSparseT;
typedef class THnSparseT<TArrayF> THnSparseF;
class TCollection;
class TList;
class THashList;
class TString;
class AliCFContainer;

class AliHFEcontainer : public TNamed{
  public:
    AliHFEcontainer();
    AliHFEcontainer(const Char_t *name);
    AliHFEcontainer(const Char_t *name, UInt_t nVar);
    AliHFEcontainer(const AliHFEcontainer &ref);
    AliHFEcontainer& operator=(const AliHFEcontainer &ref);
    ~AliHFEcontainer();

    virtual Long64_t Merge(TCollection *coll);

    void CreateContainer(const Char_t *name, const Char_t *title, UInt_t nStep);
    void CreateCorrelationMatrix(const Char_t *name, const Char_t *title);
    AliCFContainer *GetCFContainer(const Char_t *name) const;
    THnSparseF *GetCorrelationMatrix(const Char_t *name) const;
    THashList *GetListOfCorrelationMatrices() const { return fCorrelationMatrices; }
    void FillCFContainer(const Char_t *name, UInt_t step, const Double_t * const content, Double_t weight = 1.) const;
    void FillCFContainerStepname(const Char_t *name, const Char_t *step, const Double_t *const content, Double_t weight = 1.) const;
    AliCFContainer *MakeMergedCFContainer(const Char_t *name, const Char_t *title, const Char_t *contnames) const;

    Int_t GetNumberOfCFContainers() const;
    Int_t GetNumberOfEvents() const { return fNEvents; };
    void NewEvent() { fNEvents++; };
    void SetNumberOfVariables(UInt_t nVar);
    inline void SetBinning(UInt_t var, UInt_t nBins, Double_t *content);
    void SetVariableName(UInt_t var, const Char_t *varname);
    void SetStepTitle(const Char_t *contname, const Char_t *steptitle, UInt_t step);
    void MakeLinearBinning(UInt_t var, UInt_t nBins, Double_t begin, Double_t end);
    void MakeLogarithmicBinning(UInt_t var, UInt_t nBins, Double_t begin, Double_t end);
    void MakeUserDefinedBinning(UInt_t var, UInt_t nBins, const Double_t *binning);

    virtual void Print(const Option_t * opt = NULL) const;

    struct AliHFEvarInfo : public TObject{
        AliHFEvarInfo();
        AliHFEvarInfo(const Char_t *name);
        AliHFEvarInfo(const AliHFEvarInfo &ref);
        AliHFEvarInfo &operator=(const AliHFEvarInfo &ref);
        ~AliHFEvarInfo();

        UInt_t GetNumberOfBins() const { return fBinning->GetSize() ? fBinning->GetSize() - 1 : 0; };
        const Double_t *GetBinning() const { return fBinning->GetArray(); };
        TString *GetVarName() const { return fVarName; };

        void SetVarName(const Char_t *name);
        void SetBinning(UInt_t nBins, const Double_t *binning);
      private:
        TString *fVarName;  // Variable Name
        TArrayD *fBinning;  // Binning
        
        ClassDef(AliHFEcontainer::AliHFEvarInfo, 1)  // Variable information
    };

  private:
    THashList *fContainers;     // TObjArray for Containers
    THashList *fCorrelationMatrices; // Container for Correlation Matrices
    TObjArray *fVariables;      // Variable Information
    UInt_t fNVars;              // Number of Variables
    Int_t fNEvents;             // Number of Events

    ClassDef(AliHFEcontainer, 1)  // HFE Efficiency Container
};

//__________________________________________________________________
void AliHFEcontainer::SetBinning(UInt_t var, UInt_t nBins, Double_t *content){
  //
  // Set Binning for a given variable
  //
  if(var >= fNVars) return;
  AliHFEvarInfo *inf = dynamic_cast<AliHFEvarInfo *>(fVariables->UncheckedAt(var));
  if(!inf) return;
  inf->SetBinning(nBins, content); 
}
#endif
