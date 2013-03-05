#include "ARVersion.h"
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TNamed.h>
#include <THashList.h>
#include <TObjArray.h>
#include <TString.h>
#include <TObjString.h>
#endif

TObjArray* findAliasesForClass(THashList &list, const char* className)
{

    TObjArray* matchingTrAliases = new TObjArray(2);
    TNamed *n = dynamic_cast<TNamed*>(list.FindObject(className));
    if(!n){
	Printf("No entry for a trigger-class named \"%s\"",className);
	return;
    }
    TString aliasList = n->GetTitle();
    TObjArray* arrAliases = aliasList.Tokenize(',');
    Int_t nAliases = arrAliases->GetEntries();
    for(Int_t i=0; i<nAliases; i++){
	TObjString *alias = dynamic_cast<TObjString*>(arrAliases->At(i));
	matchingTrAliases->Add(alias);
    }

    return matchingTrAliases;
}
