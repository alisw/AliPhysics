#ifndef ALIPHICUT_H
#define ALIPHICUT_H

#include <TObject.h>
#include <TArrayD.h>
#include <Rtypes.h>

class AliPhiCut : public TObject {
public:
	AliPhiCut( Int_t energy = -1 );
	AliPhiCut ( Double_t* e, Double_t* w, Int_t n );
	~AliPhiCut();
	virtual void Print ( Option_t* option ) const;
	virtual Bool_t CheckCut ( const Double_t phi );

	virtual void CreateDefaultCut ( Int_t energy );
private:
	TArrayD widths;
	TArrayD edges;
	
	ClassDef ( AliPhiCut, 1 );
};

#endif // ALIPHICUT_H
