#include <TNamed.h>

class AliVersion : public TNamed
{
public: 
  AliVersion();
  virtual ~AliVersion() {if (fgInstance==this) fgInstance=0;};
  static AliVersion* Instance();
  const TString& GetHash() const { return fHash; }
  const TString& GetTag() const { return fTag; }
  UInt_t GetSerial() const { return GetUniqueID(); }
  Bool_t IsSortable() const { return kTRUE; }
  Int_t Compare(const TObject* other) const;
  virtual void Print(Option_t *opt) const;
protected:
  //
  TString fHash;                                  //last hash
  TString fTag;                                   //revision/tag name
  static AliVersion* fgInstance;
  //
  AliVersion &operator=(const AliVersion& c);     //dummy assignment operator
  AliVersion(const AliVersion& c);                //dummy copy constructor
  //
  ClassDef(AliVersion,1);
};
