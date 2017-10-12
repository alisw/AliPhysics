#ifndef _ALIANALYSISNANOAODCUTSANDSETTERSJET_H_
#define _ALIANALYSISNANOAODCUTSANDSETTERSJET_H_

#include "AliNanoAODCustomSetter.h"

class TClonesArray;

class AliNanoAODSimpleSetterJet : public AliNanoAODCustomSetter
{
public:
  AliNanoAODSimpleSetterJet();
  virtual ~AliNanoAODSimpleSetterJet(){;}

  virtual void SetNanoAODHeader(const AliAODEvent * event   , AliNanoAODHeader * head ,TString varListHeader  );
  virtual void SetNanoAODTrack (const AliAODTrack * aodTrack, AliNanoAODTrack * spTrack);

  void SetArrayPythiaName(TString name) { fArrayPythiaName = name; } 

private:
  TClonesArray *fArrayPythia;
  TString       fArrayPythiaName;  

  ClassDef(AliNanoAODSimpleSetterJet, 1);

};

#endif /* _ALIANALYSISNANOAODCUTSANDSETTERS_H_ */
