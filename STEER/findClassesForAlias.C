#include "ARVersion.h"
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TNamed.h>
#include <THashList.h>
#include <TObjArray.h>
#include <TString.h>
#include <TObjString.h>
#endif

TObjArray* findClassesForAlias(THashList &list, const char* aliasName)
{

    TObjArray* matchingTrClasses = new TObjArray(2);
    TIter iter(&list);
    TNamed *n = 0;
    iter.Reset();
    while((n = dynamic_cast<TNamed*>(iter.Next()))){
	TString aliasList(n->GetTitle());
	if(aliasList.Contains(aliasName)){
	    TObjArray* arrAliases = aliasList.Tokenize(',');
	    Int_t nAliases = arrAliases->GetEntries();
	    for(Int_t i=0; i<nAliases; i++){
		TObjString *alias = (TObjString*) arrAliases->At(i);
		if(alias->String()==TString(aliasName)){
		    TObjString *trClass = new TObjString(n->GetName());
		    matchingTrClasses->Add(trClass);
		}
	    }
	}
    }

    if (matchingTrClasses->GetEntries() == 0){
	Printf("No entries for the trigger alias \"%s\" were found. Returning null pointer.", aliasName);
	return 0;
    }

    return matchingTrClasses;
}
