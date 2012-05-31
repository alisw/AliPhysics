//
// AliGenEposIsajetToPdgConverter.h
//
//  Helper class used by TEpos to insert EPOS internal objects and higher
//  resonances to PDG database.
//  Quark clusters has has unaltered code, unknown objects have code of
//  form 6xxxxxxxx, and higher resonances have codes from Particle Physics
//  Review '96
//
//  Created on: Aug 03, 2009
//      Author: Piotr Ostrowski, postrow@if.pw.edu.pl
//


#ifndef ALI_GEN_EPOS_ISAJET_TO_PDG_CONVERTER_H
#define ALI_GEN_EPOS_ISAJET_TO_PDG_CONVERTER_H

#include <Rtypes.h>

class TDatabasePDG;

class AliGenEposIsajetToPdgConverter
{
public:
	AliGenEposIsajetToPdgConverter();
	
	Int_t ConvertIsajetToPdg(Int_t isajet) const;

	virtual ~AliGenEposIsajetToPdgConverter();
	static void AddHigherResonances();

private:
	Int_t ExtendedMapping(Int_t isajet) const;
	void AddQuarkCluster(Int_t clusterCode) const ;
	Int_t AddUnknownObject(Int_t code) const;

	static Bool_t fgParticlesAdded; //flag indicating that resonances
					// have been alredy inserted
	ClassDef(AliGenEposIsajetToPdgConverter,1)
};

#endif /* ALI_GEN_EPOS_ISAJET_TO_PDG_CONVERTER_H */
