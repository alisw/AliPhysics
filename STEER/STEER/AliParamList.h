#ifndef ALIPARAMLIST_H
#define ALIPARAMLIST_H

#include <TNamed.h>

class TString;

class AliParamList : public TNamed
{
  public:
  AliParamList(Int_t n=0, const Double_t *parVal=0);
  AliParamList(const AliParamList& src);
  AliParamList& operator=(const AliParamList& src);
  virtual ~AliParamList();
  //
  Int_t         GetID()                const {return fID;}
  Int_t         GetNParams()           const {return fNPar;}
  Double_t*     GetParams()            const {return (Double_t*)fParams;}
  TString*      GetNames()             const {return (TString*) fNames;}
  Double_t      GetParameter(Int_t i)  const {return fParams[i];}
  const Char_t* GetParName(Int_t i)    const;
  //
  void          SetID(Int_t id)              {fID = id;}
  void          SetNParams(Int_t n);
  void          SetParName(Int_t i, const char* nm);
  void          SetParameter(Int_t i, Double_t v, const char* nm=0);
  void          SetParameters(const Double_t* vals) {for (int i=0;i<fNPar;i++) SetParameter(i,vals[i]);}
  //
  virtual void  Print(Option_t *opt="") const;
  //
 protected:
  Int_t     fID;       // user defined id
  Int_t	    fNPar;     // number of parameters
  TString*  fNames;    //[fNPar] parameter names
  Double_t* fParams;   //[fNPar] parameter values
  //
  ClassDef(AliParamList,1)
};

#endif
