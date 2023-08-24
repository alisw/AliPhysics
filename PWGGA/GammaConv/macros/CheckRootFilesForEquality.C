// author: Stephan Stiefelmaier (stephan.friedrich.stiefelmaier@cern.ch) please contact in case of questions, difficulties, found issues, improvement etc
/* purpose: use this macro to compare two .root files if they contain identical objects. The following object types can be checked by the script:
- TH1* TH2* TH3*
- TProfile, TProfile2D
- TTree (only certain branch types so far, you get warned if a not supported branch type is in the Tree)
- TList
*/

/* todos:
- clean up
- go through comments
- ensure proper naming of not identical objects:
"Files differ. 2 objects have same name and type but differing content. There are 1049 identical objects. The objects that differ are the following:
THashList/GammaConvV1_2500/BrokenFiles (TTree)
THashList/GammaConvV1_2500 (TList)"
- go again through compare branches functions of trees if all elements get compared */

#include <iostream>
#include <string.h>
#include <vector>


#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

// todo: understand why dynamic_cast works without std::
using std::cout;
using std::endl;

Bool_t sameName(TObject* obj1, TObject* obj2){
    return strncmp(obj1->GetName(), obj2->GetName(), 100)==0;
}

Bool_t sameClassName(TObject* obj1, TObject* obj2){
    return strncmp(obj1->ClassName(), obj2->ClassName(), 100)==0;
}

template <class T>
Bool_t compareNEntries(T* obj1, T* obj2, TString fullName=""){

    if (!(obj1&&obj2)){
        cout << "compareNEntries(): nullptr!. fullName = " << fullName.Data() << endl;
        return kFALSE;
    }

    if (obj1->GetEntries() != obj2->GetEntries()){
        cout << "difference in " << fullName.Data();
        cout <<  "\n\tentries: " << obj1->GetEntries() << " versus " << obj2->GetEntries() << endl;
        return kFALSE;
    }
    return kTRUE;
}

template <class T>
Bool_t compareHistos(T* h1, T* h2, TString fullName=""){

    if (!(h1&&h2)){
        cout << "compareHistos(): nullptr!. fullName = " << fullName.Data() << endl;
        return kFALSE;
    }

    // compare first number of entries (fast)
    Bool_t lDifferenceFound = !compareNEntries(h1, h2, fullName);

    // uncomment this to only check for number of entries but not bin counts
    //~ if (lDifferenceFound){
        //~ return kFALSE;
    //~ }

    // if previous check is good, compare counts in each bin
    Int_t nBins = h1->GetNbinsX();
    if (h1->InheritsFrom("TH2")){
        nBins *= h1->GetNbinsY();
    }
    if (h1->InheritsFrom("TH3")){
        nBins *= h1->GetNbinsZ();
    }


    for (Int_t iBin=0; iBin<nBins; ++ iBin){
        if (h1->GetBinContent(iBin) != h2->GetBinContent(iBin)){
            if (!lDifferenceFound) cout << "difference in " << fullName.Data();
            cout <<  "\n\tbin content bin " << iBin << ": " << h1->GetBinContent(iBin) << " vs. " << h2->GetBinContent(iBin) << endl;
            lDifferenceFound = kTRUE;
            //~ return kFALSE;
        }
    }
    return !lDifferenceFound;
    //~ return kTRUE;
}


// todo: why does lArray1.GetSetupStatus() return < 0? Should be treated as error according to
//       https://root.cern.ch/doc/v608/classTTreeReader.html
//       -> ignoring for now.
// todo: I think handing over both branch readers and branches might be superflous
template<typename T>
Bool_t compareBranches_fT(TTreeReader &tr1, TBranch &b1, TTreeReader &tr2, TBranch &b2){

    if ((b1.GetNleaves() != b2.GetNleaves()) || (b1.GetNleaves() !=1 )){
        cout << "INFO: Differing number of leaves or number !=1, returning kFALSE.\n"
                "Change this to kTRUE if you want the branches to appear as if identical.\n";
        return kFALSE;
    }

    TObjArray *lLeaves1 = b1.GetListOfLeaves();
    TObjArray *lLeaves2 = b2.GetListOfLeaves();
    if (!(lLeaves1 && lLeaves2)){
        return kFALSE;
    }

    TLeaf *lLeaf1 = dynamic_cast<TLeaf*>(lLeaves1->At(0));
    TLeaf *lLeaf2 = dynamic_cast<TLeaf*>(lLeaves2->At(0));

    if (!(lLeaf1 && lLeaf2)){
        cout << "One of the leaves was a nullptr.\n";
        return kFALSE;
    }

    TString lCompleteName(Form("%s.%s", b1.GetName(), lLeaf1->GetName()));

    TTreeReaderArray<T> lArray1(tr1, lCompleteName.Data());
    TTreeReaderArray<T> lArray2(tr2, lCompleteName.Data());
    //~ if ((lArray1->GetSetupStatus() >=0) && (lArray2->GetSetupStatus()>=0)) {}

    if (!tr1.GetEntries() && !tr2.GetEntries()){
        cout << "both trees seem to be empty, they will appear as identical\n";
        return kTRUE;
    }

    if (tr1.GetEntries() != tr2.GetEntries()){
        cout << "differing number of entries\n";
        return kFALSE;
    }

    Bool_t lOneIt = kFALSE;
    // todo: check this is safe - the next line
    while (tr1.Next()){
        lOneIt = kTRUE;
        if (!tr2.Next()){
            cout << "treereader for second tree could not advance when treereader for first one could.\n";
            return kFALSE;
        }

        if (lArray1.GetSize() != lArray2.GetSize()){
            return kFALSE;
        }
        if (lArray1.GetSize() > 1){
            return kFALSE;
        }
        if ((lArray1.GetSize() == 1) && (lArray1[0] != lArray2[0])){
            cout << "differing entry in branch\n";
            tr1.Restart();
            tr2.Restart();
            return kFALSE;
        }
    }
    tr1.Restart();
    tr2.Restart();
    return lOneIt;
}

Bool_t compareBranchesTO(TTreeReader &tr1, TBranch &b1, TTreeReader &tr2, TBranch &b2){

    if ((b1.GetNleaves() != b2.GetNleaves()) || (b1.GetNleaves() !=1 )){
        cout << "INFO: Differing number of leaves or number !=1, returning kFALSE.\n"
                "Change this to kTRUE if you want the branches to appear as if identical.\n";
        return kFALSE;
    }

    TObjArray *lLeaves1 = b1.GetListOfLeaves();
    TObjArray *lLeaves2 = b2.GetListOfLeaves();
    if (!(lLeaves1 && lLeaves2)){
        return kFALSE;
    }

    TLeaf *lLeaf1 = dynamic_cast<TLeaf*>(lLeaves1->At(0));
    TLeaf *lLeaf2 = dynamic_cast<TLeaf*>(lLeaves2->At(0));

    if (!(lLeaf1 && lLeaf2)){
        cout << "One of the leaves was a nullptr.\n";
        return kFALSE;
    }

    TTreeReaderValue<TObject> lValue1(tr1, lLeaf1->GetName());
    TTreeReaderValue<TObject> lValue2(tr2, lLeaf1->GetName());

    ROOT::Internal::TTreeReaderValueBase &cVal(lValue1);
    //~ cout << "setup status = " << cVal.GetSetupStatus() << endl;

    // todo: check this is safe - the next line
    while (tr1.Next() && tr2.Next()){

        if (kTRUE) { // some meaningfull comparison between lValue1 and lValue2 should go here
            cout << "not same\n";
            tr1.Restart();
            tr2.Restart();
            return kFALSE;
        }
    }
    tr1.Restart();
    tr2.Restart();
    return kTRUE;
}


Bool_t compareTreesT(TTree *t1, TTree *t2, TString fullName=""){

    if (!(t1&&t2)){
        cout << "compareTrees(): nullptr!. fullName = " << fullName.Data() << endl;
        return kFALSE;
    }

    TObjArray *lBranches1 = t1->GetListOfBranches();
    if (!lBranches1){
        cout << "INFO: No branches in tree t1 " << fullName.Data() << endl;
        return  kTRUE;
    }

    TObjArray *lBranches2 = t2->GetListOfBranches();
    if (!lBranches2){
        cout << "no branches in tree t2 but there are branches in t1. -> not identical" << fullName.Data() << endl;
        return  kFALSE;
    }

    TTreeReader tr1(t1);
    TTreeReader tr2(t2);

    Bool_t lCheckedGood = kTRUE;
    for (TObject *iObj1 : *lBranches1){
        TBranch *iBranch1 = dynamic_cast<TBranch*>(iObj1);
        if (iBranch1){
            //~ cout << iBranch1->GetName() << " nLeaves " << iBranch1->GetNleaves() << endl;

            Bool_t lFoundIdenticalObject = kFALSE;
            for (TObject *iObj2 : *lBranches2){
                TBranch *iBranch2 = dynamic_cast<TBranch*>(iObj2);
                if (iBranch2 && sameName(iBranch1, iBranch2)){
                    //~ cout << "found same branches " << iBranch1->GetName() << endl;

                    // find out branch type
                    TClass *lClass1 = nullptr;
                    TClass *lClass2 = nullptr;
                    EDataType lType1, lType2;

                    Int_t s1 = iBranch1->GetExpectedType(lClass1, lType1);
                    Int_t s2 = iBranch2->GetExpectedType(lClass2, lType2);

                    if (s1 || s2){
                        cout << "could not obtain type information for\n"; // todo add proper name
                        return kFALSE;
                    }

                    if (lType1==lType2){
                        //~ cout << "INFO: Found branches with same name and type and one leaf each. Name is " << iBranch1->GetName() << "and type is " << lType1 << endl;

                        if(lType1==kFloat_t){
                            //~ cout << "kFloat_t\n";
                            lFoundIdenticalObject = compareBranches_fT<Float_t>(tr1, *iBranch1, tr2, *iBranch2);
                        }
                        else if (lType1==kUChar_t){
                            //~ cout << "kUChar_t\n";
                            lFoundIdenticalObject = compareBranches_fT<UChar_t>(tr1, *iBranch1, tr2, *iBranch2);
                        }
                        //~ else if (lType1==-1){
                            //~ cout << "no comparison for this branch type implemented, they will appear as if not identical.\n";
                            //~ /*lFoundIdenticalObject = compareBranchesTO(tr1, *iBranch1, tr2, *iBranch2);*/
                        //~ }
                        else{
                            cout << "INFO: In compareTreesT(): Found branches with name \'" << iBranch1->GetName() << "\': No comparison for this branch type implemented, they will appear as if not identical.\n";
                        }
                    }
                    //~ else { cout << "types differ"; }

                    if (lFoundIdenticalObject){
                        break;
                    }
                }
            }
            lCheckedGood &= lFoundIdenticalObject;
        }
        // todo: maybe remove
        else {
            cout << "INFO: Found non branch type object in branches of tree1 with name " << fullName.Data() << endl;
        }
    }
    return lCheckedGood;
}

Bool_t compareOtherTypes(TObject* obj1, TObject* obj2, TString& fullName){

    if (!(obj1&&obj2)){
        cout << "compareOtherTypes(): nullptr!. fullName = " << fullName.Data() << endl;
        return kFALSE;
    }

    if (obj1->InheritsFrom("TH1")){
        return compareHistos(dynamic_cast<TH1*>(obj1), dynamic_cast<TH1*>(obj2), fullName);
    }
    else if (obj1->InheritsFrom("TH2")){
        return compareHistos(dynamic_cast<TH2*>(obj1), dynamic_cast<TH2*>(obj2), fullName);
    }
    else if (obj1->InheritsFrom("TTree")){
        //~ return compareNEntries(dynamic_cast<TTree*>(obj1), dynamic_cast<TTree*>(obj2), fullName);
        return compareTreesT(dynamic_cast<TTree*>(obj1), dynamic_cast<TTree*>(obj2), fullName);

    }
    else{
        cout << "no comparison implemented for these objects: " << fullName.Data() << endl;
        return kFALSE;
    }
}

// todo: think again if all the naming with theHierarchy and lFullName is working is supposed and makes sense

/* Stratgey: Iterate over every element on theList1. For each element (lObj1) found,
 * loop over theList2 until an element (lObj2) is found with the same name.
 * If lObj1 and lObj2 have the same type, check if the type is TList and if it is,
 * call compareLists() again with these two objects as arguments.
 * If the type is not TList, call compareOtherTypes()
 *
 * If theList1 and theList2 come from a GetListOfKeys() call,
 * theFile1 and theFile2 are required as well. See comment over theFile1->Get( to see why.
 */
Bool_t compareLists(TList* theList1, TList* theList2,
                    std::vector<TString>& theOnBothAndEqual,
                    std::vector<TString>& theOnBothDiffering,
                    std::vector<TString>& theOnFirstNotOnSecond,
                    TString theHierachy, TString theIndent,
                    TFile* theFile1=nullptr, TFile* theFile2=nullptr){
    if (!theList1) cout << theIndent.Data() << "theList1 is null\n";
    if (!theList2) cout << theIndent.Data() << "theList2 is null\n";

    // todo: check that theList1 always delivers correct name for lNewHierachy
    TString lNewHierarchy = theHierachy  + theList1->GetName() + "/";
    Bool_t lCheckedGood = kTRUE;

    TIter lIter1(theList1);
    while (TObject *lObj1 = lIter1.Next()){

        // lObj1 is a TKey if theList1 was obtained from TFile::GetListOfKeys() call. This will bring its natural type
        if (theFile1) {
            lObj1 = theFile1->Get(lObj1->GetName());
        }
        TString lFullName = TString() + lNewHierarchy.Data() + lObj1->GetName() + " (" + lObj1->ClassName() + ")";

        Bool_t lFoundObjectWithSameNameAndType = kFALSE;
        Bool_t lFoundIdenticalObject = kFALSE;
        TIter lIter2(theList2);
        while (TObject *lObj2 = lIter2.Next()) {

            if (sameName(lObj1, lObj2)){
                if (theFile2) {
                    lObj2 = theFile2->Get(lObj2->GetName()); // same reasoning as above
                }

                if ((lFoundObjectWithSameNameAndType = sameClassName(lObj1, lObj2))){

                    if ((lFoundIdenticalObject = lObj1->InheritsFrom("TList") ?
                        compareLists(dynamic_cast<TList*>(lObj1),
                                     dynamic_cast<TList*>(lObj2),
                                     theOnBothAndEqual,
                                     theOnBothDiffering,
                                     theOnFirstNotOnSecond,
                                     lNewHierarchy,
                                     theIndent + "    "):
                        compareOtherTypes(lObj1,
                                          lObj2,
                                          lFullName)))
                    {
                        theOnBothAndEqual.push_back(lFullName);
                        // break inner while loop, which tries to find a matching lObj2 for lObj1
                        break;
                    }
                    else {
                        // same name and type but differing content
                        theOnBothDiffering.push_back(lFullName);
                    }
                }
                else {
                    cout << theIndent.Data() << "same name but different types\n";
                }
            }
        }
        if (!lFoundObjectWithSameNameAndType) {
            theOnFirstNotOnSecond.push_back(lFullName);
        }
        lCheckedGood &= lFoundIdenticalObject;
    }
    return lCheckedGood;
}



void CheckRootFilesForEquality(TString fname1="", TString fname2=""){

    gROOT->Reset();

    TFile* lFile1 = new TFile(fname1);
    TFile* lFile2 = new TFile(fname2);

    std::vector<TString> lOnBothAndEqual12;
    std::vector<TString> lOnBothDiffering12;

    std::vector<TString> lOnBothAndEqual21;
    std::vector<TString> lOnBothDiffering21;

    std::vector<TString> lOnFirstNotOnSecond;
    std::vector<TString> lOnSecondNotOnFirst;

    Bool_t lSame12 = kTRUE;
    Bool_t lSame21 = kTRUE;

    cout << "===============================================================\n";
    cout << "Current pwd is: " << gSystem->pwd() << endl;
    cout << "Comparing the following files:\n" << lFile1->GetName() << endl << lFile2->GetName() << endl;
    cout << "===============================================================\n";


    lSame12 = compareLists(lFile1->GetListOfKeys(),
                           lFile2->GetListOfKeys(),
                           lOnBothAndEqual12,
                           lOnBothDiffering12,
                           lOnFirstNotOnSecond,
                           "", "",
                           lFile1, lFile2);
    cout << "===============================================================\n";
    cout << "swapping files\n";
    cout << "===============================================================\n";
    lSame21 = compareLists(lFile2->GetListOfKeys(),
                           lFile1->GetListOfKeys(),
                           lOnBothAndEqual21,
                           lOnBothDiffering21,
                           lOnSecondNotOnFirst,
                           "", "",
                           lFile2, lFile1);

    cout << "===============================================================\n";
    cout << "Comparing of files finished. Printing starts now. (Un)comment to steer verbosity.\n";
    cout << "===============================================================\n";


    // uncomment this to get a list of all compared objects (that compared equal)
    //~ cout << "The following objects were found in both files and compared equal:\n";
    //~ for (auto &lName : lOnBothAndEqual12){ cout << lName.Data() << endl; }
    //~ cout << "Printing of equal objects finished.\n";


    if (lSame12 && lSame21){
        cout << "file1 and file2 are identical. They contain " << lOnBothAndEqual12.size() << " identical objects.\n";
    }
    else{
        cout << "Files differ. " << lOnBothDiffering12.size() << " objects have same name and type but differing content. There are " << lOnBothAndEqual12.size() << " identical objects. The objects that differ are the following:\n";

        for (auto &lName : lOnBothDiffering12){
            cout << lName.Data() << endl;
        }
        if (lOnFirstNotOnSecond.size()){
            cout << "The following objects are in file1 but not in file2:\n";
            for (auto &lName : lOnFirstNotOnSecond){
                cout << lName.Data() << endl;
            }
        }
        if (lOnSecondNotOnFirst.size()){
            cout << "The following objects are in file2 but not in file1:\n";
            for (auto &lName : lOnSecondNotOnFirst){
                cout << lName.Data() << endl;
            }
        }
    }
}
