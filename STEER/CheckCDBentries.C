#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliCDBId.h"
#include <TString.h>
#include <TList.h>
#endif

void CheckCDBentries(const char* dest)
{
	AliCDBManager* cdb = AliCDBManager::Instance();
	const char* ref="local://$ALICE_ROOT/OCDB";
	cdb->SetDefaultStorage(ref);
	cdb->SetRun(0);
	AliCDBStorage* newstor = cdb->GetStorage(dest);
	// Missing here a check that newstor is a valid storage
	// otherwise exit
	TList* allentries = cdb->GetAll("*/*/*");
	TList* allnewentries = newstor->GetAll("*/*/*",0);
	Int_t nall = allentries->GetEntries();
	Int_t nallnew = allnewentries->GetEntries();
	Printf("Number of entries in reference OCDB %d  and in checked OCDB %d",nall, nallnew);
	TString missing;
	Int_t nMissing=0;
	if(nall!=nallnew)
	{
		AliCDBEntry *entry, *newentry;
		TString cdbpath;
		for(Int_t i=0; i<nall; i++)
		{
			entry = dynamic_cast<AliCDBEntry*>(allentries->At(i));
			cdbpath = ((AliCDBId)entry->GetId()).GetPath();
			newentry = newstor->Get(cdbpath.Data(),0);
			if(!newentry)
			{
				missing += cdbpath;
				missing.Insert(missing.Length(),'\n');
				nMissing++;
			}
		}
	}
	if(nMissing)
	{
		Printf("\n\nEntries missing in destination OCDB folder %s w.r.t. reference folder %s:",dest,ref);
		Printf("%s",missing.Data());
	}else{
		Printf("\n\nNo entry is missing in destination OCDB folder %s w.r.t. reference folder %s:",dest,ref);
	}
}
