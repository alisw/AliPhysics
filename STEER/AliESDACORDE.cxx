#include "AliESDACORDE.h"

ClassImp(AliESDACORDE)

AliESDACORDE::AliESDACORDE():TObject()
{
 //Default constructor
	for(Int_t i=0;i<60;i++)
	{
		fACORDESingleMuon[i] = fACORDEMultiMuon[i] = 0;
	}
}


AliESDACORDE::AliESDACORDE(const AliESDACORDE &o)
  :TObject(o)

{	
	//Default constructor
	for(Int_t i=0;i<60;i++)
	{
		fACORDESingleMuon[i] = o.fACORDESingleMuon[i];
		fACORDEMultiMuon[i] = o.fACORDEMultiMuon[i];
	}
}


AliESDACORDE::AliESDACORDE(Int_t* MACORDESingleMuon, Int_t* MACORDEMultiMuon):TObject()
{

	//Constructor

	for(Int_t i=0;i<60;i++)
	{
		fACORDESingleMuon[i] = MACORDESingleMuon[i];
		fACORDEMultiMuon[i] = MACORDEMultiMuon[i];
	}
}

AliESDACORDE& AliESDACORDE::operator=(const AliESDACORDE& o)
{
	if(this==&o)return *this;
	TObject::operator=(o);

	// Assignment operator
	for(Int_t i=0; i<60; i++)
	{
		fACORDESingleMuon[i] = o.fACORDESingleMuon[i];
		fACORDEMultiMuon[i] = o.fACORDEMultiMuon[i];
	}
	
	return *this;
}

void AliESDACORDE::Copy(TObject &obj) const {
  
  // this overwrites the virtual TOBject::Copy()
  // to allow run time copying without casting
  // in AliESDEvent

  if(this==&obj)return;
  AliESDACORDE *robj = dynamic_cast<AliESDACORDE*>(&obj);
  if(!robj)return; // not an AliESDACRDE
  *robj = *this;

}


