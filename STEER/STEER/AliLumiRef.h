#include <TNamed.h>

class AliLumiRef : public TNamed
{
 public:
  
  AliLumiRef(const char* trigger="", const char* comment="", UInt_t run=0, float sig=0.f, float eff=1.f);
  ~AliLumiRef() {}
  AliLumiRef(const AliLumiRef& src);
  AliLumiRef &operator=(const AliLumiRef& src);
  virtual void Print(const Option_t *opt) const;

  
  const char* GetRefTrigger()  const {return GetName();}
  const char* GetComment()     const {return GetTitle();}
  UInt_t      GetRunStart()    const {return GetUniqueID();}
  float GetRefSigma()          const {return fRefSigma;}
  float GetRefEff()            const {return fRefEff;}

  void  SetRefTrigger(const char* trig) {SetName(trig);}
  void  SetComment(const char* cmt)     {SetTitle(cmt);}
  void  SetRunStart(UInt_t run)         {SetUniqueID(run);}
  void  SetRefSigma(float sig)          {fRefSigma = sig;}
  void  SetRefEff(float eff)            {fRefEff = eff;}

  Bool_t IsSortable() const {return kTRUE;}
  virtual Int_t Compare(const TObject* obj) const;

 protected:
  Float_t fRefSigma;       // ref. x-section
  Float_t fRefEff;         // ref. eff
  //
  ClassDef(AliLumiRef,1)
};
