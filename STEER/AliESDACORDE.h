#ifndef AliESDACORDE_H
#define AliESDACORDE_H

#include <TObject.h>

class AliESDACORDE : public TObject
{

 public:
  AliESDACORDE();
  AliESDACORDE(const AliESDACORDE&);
  AliESDACORDE(Int_t *ACORDESingleMuon,Int_t *ACORDEMultiMuon); 
  virtual ~AliESDACORDE() {};
  virtual void Copy(TObject &) const;
  void SetACORDEMultiMuon(Bool_t ACORDEMultiMuon[60]){for(Int_t i=0;i<60;i++){fACORDEMultiMuon[i]=ACORDEMultiMuon[i];}}
  
  void SetACORDESingleMuon(Bool_t ACORDESingleMuon[60]){for(Int_t i=0;i<60;i++){fACORDESingleMuon[i]=ACORDESingleMuon[i];}} 
  
	
  AliESDACORDE &operator=(const AliESDACORDE& source);
  
 protected:
  
  Bool_t	fACORDESingleMuon[60];	// array with the Single Muon hits in the 60 Acorde's Modules 
  Bool_t	fACORDEMultiMuon[60];   // array with the Multi Muon hits in the 60 Acorde's Modules

  ClassDef(AliESDACORDE,2)

};

#endif
