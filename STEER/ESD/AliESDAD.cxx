

// Last update: Nov. 5th 2013 

#include "AliESDAD.h"

ClassImp(AliESDAD)

AliESDAD::AliESDAD():TObject()
{
 //Default constructor
	for(Int_t i=0;i<16;i++)
	{
		fADCellID[i] = 0;
	}
}


AliESDAD::AliESDAD(const AliESDAD &o)
  :TObject(o)

{	
	//Default constructor
	for(Int_t i=0;i<16;i++)
	{
		fADCellID[i] = o.fADCellID[i];
	}
}


AliESDAD::AliESDAD(Bool_t* MADBitCell):TObject()
{

	//Constructor

	for(Int_t i=0;i<16;i++)
	{
		fADCellID[i] = MADBitCell[i];
	}
}

AliESDAD& AliESDAD::operator=(const AliESDAD& o)
{
// Copy Constructor
	if(this==&o)return *this;
	TObject::operator=(o);

	// Assignment operator
	for(Int_t i=0; i<16; i++)
	{
		fADCellID[i] = o.fADCellID[i];
	}
	
	return *this;
}


Bool_t AliESDAD::GetADCell(Int_t i) const
{
	return fADCellID[i];
}

void AliESDAD::Copy(TObject &obj) const {
  
  // this overwrites the virtual TOBject::Copy()
  // to allow run time copying without casting
  // in AliESDEvent

  if(this==&obj)return;
  AliESDAD *robj = dynamic_cast<AliESDAD*>(&obj);
  if(!robj)return; // not an AliESDAD
  *robj = *this;

}


