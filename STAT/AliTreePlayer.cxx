/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/




#include "TStatToolkit.h"
#include "Riostream.h"
#include <iostream>
#include "TSystem.h"
#include "TNamed.h"
#include "TFile.h"
#include "TTree.h"
#include "TPRegexp.h"
#include "TFriendElement.h"
#include "AliExternalInfo.h"
#include "TTreeFormula.h"
#include "TTreeFormulaManager.h"
#include "AliTreePlayer.h"
#include "TEntryList.h"
#include "THn.h"
#include "TLegend.h"
#include "TBufferJSON.h"
#include <sstream>
#include <bitset>
#include "TTimeStamp.h"
using namespace std;
ClassImp(AliTreePlayer)
ClassImp(AliTreeFormulaF)


/// Default constructor
AliTreeFormulaF::AliTreeFormulaF() : TTreeFormula(), fTextArray(NULL), fFormatArray(NULL), fFormulaArray(NULL),fValue(""),fDebug(0) {
}

///
/// \param name
/// \param formula
/// \param tree
AliTreeFormulaF::AliTreeFormulaF(const char *name, const char *formula, TTree *tree, Int_t debug):
        TTreeFormula(), fTextArray(NULL), fFormatArray(NULL), fFormulaArray(NULL), fValue(""), fDebug(debug) {
  SetName(name);
  SetTitle(formula);
  SetTree(tree);
  Compile(formula);
  SetBit(kIsCharacter);

}

AliTreeFormulaF::~AliTreeFormulaF() {
  /// TODO  - check if not memory leak
  delete fTextArray;
  delete fFormulaArray;
  delete fFormatArray;
}

/// Compile Formatted expression
/// \param expression
/// \return
/// TODO - Error handling - invalidate expression in case of the compilation error
Int_t AliTreeFormulaF::Compile(const char *expression) {
  TTree *tree = GetTree();
  TString fquery = expression;
  //
  // 1. GetList of %format{variables}  to query
  //    format:   TStringf
  //    variable: TTreeFormula
  // 3. GetFormatted string
  //
  //Int_t sLength = fquery.Length();     // original format string
  //Int_t iVar = 0;
  //Int_t varBegin = -1;
  fFormulaArray = new TObjArray;
  fFormatArray = new TObjArray;
  fTextArray = new TObjArray;
  Int_t nVars = 0;
  Int_t lastI = 0;
  Int_t iChar = 0;
  for (iChar = 0; iChar < fquery.Length(); iChar++) {
    if (fquery[iChar] != '{') continue;
    Int_t delta0 = 1, deltaf = 0;
    for (delta0 = 1; delta0 + iChar < fquery.Length(); delta0++) if (fquery[iChar + delta0] == '}') break;
    for (deltaf = -1; deltaf + iChar >= 0; deltaf--) if (fquery[iChar + deltaf] == '%') break;
    TString stext(fquery(lastI, iChar + deltaf - lastI));
    TString sformat(fquery(iChar + deltaf + 1, -(deltaf + 1)));
    TString sformula(fquery(iChar + 1, delta0 - 1));
    if (fDebug&4) { ///TODO - add AliDebug
      printf("%d\t%d\t%d\n", iChar, deltaf, delta0);
      printf("%d\t%s\t%s\t%s\n", nVars, stext.Data(), sformat.Data(), sformula.Data());
    }
    fFormatArray->AddAtAndExpand(new TObjString(sformat.Data()), nVars);
    fTextArray->AddAtAndExpand(new TObjString(stext.Data()), nVars);
    fFormulaArray->AddAtAndExpand(new TTreeFormula(sformula.Data(), sformula.Data(), tree), nVars);
    nVars++;
    lastI = iChar + delta0 + 1;
  }
  TString stext(fquery(lastI, fquery.Length() - lastI));
  fTextArray->AddAtAndExpand(new TObjString(stext.Data()), nVars);
  return 0;
}

///
/// \param mode
/// \return
char *AliTreeFormulaF::PrintValue(Int_t mode) const {
   PrintValue(mode, 0, "");
   return NULL;
}

/// Overwrite TTreeFormula PrintValue
/// \param mode
/// \param instance
/// \param decform
/// \return
/// TODO  - proper treatment of the 64bits/32 bits integer
/// TODO  - speed up evaluation caching information (e.g format type )
///       - get rid of the root formatting string and use sprintf  formatting
///       - e.g using t
char *AliTreeFormulaF::PrintValue(Int_t mode, Int_t instance, const char *decform) const {
  std::stringstream stream;
  Int_t nVars = fFormulaArray->GetEntries();
  for (Int_t iVar = 0; iVar <= nVars; iVar++) {
    stream << fTextArray->At(iVar)->GetName();
    if (fDebug&2) cout<<"T"<<iVar<<"\t\t"<<fTextArray->At(iVar)->GetName()<<endl;
    if (fDebug&4) cout<<"F"<<iVar<<"\t\t"<<stream.str().data()<<endl;
    if (iVar < nVars) {
      TTreeFormula *treeFormula = (TTreeFormula *) fFormulaArray->At(iVar);
      Bool_t isHex=TString(fFormatArray->At(iVar)->GetName()).Contains("x") || TString(fFormatArray->At(iVar)->GetName()).Contains("X");
      if (isHex) { // TODO check 64bits/32 bits
        char buffer[64];
        Int_t value=treeFormula->EvalInstance();
        //sprintf(format,"%s",fFormatArray->At(iVar)->GetName());
        sprintf(buffer,"%X",(Int_t)value);
        stream<<buffer;
        continue;
      }
      if ((fFormatArray->At(iVar)->GetName()[0]=='b')) { //TODO  parse number of bits
        Long64_t value=treeFormula->EvalInstance();
        bitset<32> b(value);
        stream<<b.to_string();
        continue;
      }
      if ((fFormatArray->At(iVar)->GetName()[0]=='t')) { //TODO use formats
        Long64_t value=treeFormula->EvalInstance64();
        TTimeStamp stamp(value);
        stream<<stamp.AsString("s");
        continue;
      }
      TString value=treeFormula->PrintValue(0, instance, fFormatArray->At(iVar)->GetName());
      stream << value.Data();
      if (fDebug&&1) cout<<"F"<<iVar<<"\t\t"<<value.Data()<<"\t"<<fFormatArray->At(iVar)->GetName()<<endl;
    }
  }
  //const_cast<AliTreeFormulaF*>(this)->fValue=TString(stream.str().data());      /// TODO check compliation of mutable in old c++
  fValue=stream.str().data();
  return (char *) fValue.Data();
}
///
void        AliTreeFormulaF::UpdateFormulaLeaves(){
  if (fFormulaArray==NULL) return;
  return;
  Int_t nVars = fFormulaArray->GetEntries();
  for (Int_t iVar = 0; iVar <= nVars; iVar++) {
    TTreeFormula *treeFormula = (TTreeFormula *) fFormulaArray->At(iVar);
    treeFormula->UpdateFormulaLeaves();
  }
}

/// Dummy AliTreePlayerConstructor
/// \param name
/// \param title
AliTreePlayer::AliTreePlayer(const char *name, const char *title)
        :TNamed(name, title)
{
  //
  //
}



/// Selected metadata fil filing query
/// retun TObjArray of selected metadata
/// \param tree
/// \param query
/// \param verbose
/// \return         TObjArray of selected metadata
/*!
    Example usage:
    \code
    query="[AxisTitle]";                     // only metadata with metadata axis existing
    query="[class]";                         // only metadata with metadata class existing
    query="[class==\"Base&&Log&&Stat\"]";    // metadata class contains
    //
    AliTreePlayer::selectMetadata(treeLogbook, "[class==\"Logbook&&(Stat||Base)\"]",0)->Print();
    AliTreePlayer::selectMetadata(testTree, "[Axis==\"counts\"]",0)->Print();
    \endcode
 */
TObjArray  *  AliTreePlayer::selectMetadata(TTree * tree, TString query, Int_t verbose, TString *idList){
  /*
    query -  case sensitive matching is done using Contains method (e.g N(Axis)<=N(xis) )
          -  WILL  be  case insesitive 
  */
  //
  // 1. Parse query and map query to TFormula selection
  //
  TObjArray *queryArray=query.Tokenize("[=]");
  TObjArray *varArray=0;
  TVectorD vecParam;
  TFormula *pFormula=0;
  if (queryArray->GetEntries()>1) {
    TString  stringFormula="";
    varArray = TString(queryArray->At(1)->GetName()).Tokenize("\"&|!()[]+-/");  //find variable list
    vecParam.ResizeTo(varArray->GetEntriesFast());
    stringFormula=queryArray->At(1)->GetName();
    stringFormula.ReplaceAll("\"","");
    for (Int_t ivar=0; ivar<varArray->GetEntriesFast(); ivar++){
      stringFormula.ReplaceAll(varArray->At(ivar)->GetName(),TString::Format("x[%d]",ivar).Data());   // Logic should be improved - to match only "full strings"
      vecParam[ivar]=1;
    }
    pFormula= new TFormula("printMetadataFromula",stringFormula.Data());
    if (verbose&0x2) pFormula->Print();
  }
  //
  if (queryArray->GetEntriesFast()<=0 ){
    return 0;
  }else{
    if (queryArray->GetEntriesFast()>2) {
      return 0;
    }
  }
  //
  // 2.) Find selected data and add the to the list
  //
  TObjArray * metaData = (TObjArray*)tree->GetUserInfo()->FindObject("metaTable");
  Int_t entries = metaData->GetEntries();
  TObjArray *selected = new TObjArray(entries);
  for (Int_t ientry=0;ientry<entries; ientry++){
    TString index=metaData->At(ientry)->GetName();
    if (strstr(metaData->At(ientry)->GetName(),queryArray->At(0)->GetName())==NULL) continue;
    Bool_t isSelected=kTRUE;
    if (queryArray->GetEntriesFast()>1){
      for (Int_t ipar=0; ipar< vecParam.GetNrows(); ipar++){
        vecParam[ipar]=strstr(metaData->At(ientry)->GetTitle(),varArray->At(ipar)->GetName())!=NULL;
      }
      isSelected=pFormula->EvalPar(vecParam.GetMatrixArray());
      if (verbose&0x2){
        vecParam.GetMatrixArray();
      }
    }
    if (isSelected){
      selected->AddLast(metaData->At(ientry));
      if (idList){
        TString id=metaData->At(ientry)->GetName();
        id.Remove(id.Last('.'),id.Length());
        if (id.Length()>0) {
          (*idList) += id;
          (*idList) +=":";
        }
      }
      if (verbose&0x1){
        printf("%s:%s:%d\n",metaData->At(ientry)->GetName(),metaData->At(ientry)->GetTitle(),isSelected);
      }
    }
  }
  if (idList!=NULL){
    delete selected;
    return NULL;
  }
  return selected;
}


/// Find all branches containing "Data string" and having class base
///    search follow only branches - does not not enter into classes
///   == - for exact match
///    :  - for contain
///  Implemented using TFormula interpreter:
///    - TODO - not case sensitive (should be part of unit test)
/// \param tree     - input tree
/// \param query    - query string
/// \param verbose  - verbosity
/// \return
/*!
 ### Implementation using TFormula
\code
  TString query="([.name==Data])&&([.class==base]||[.class==Data])";
    ==>
    TFormula   (x[0])&&(x[1]||x[2])
\endcode
### Example cases :
\code
      AliTreePlayer::selectTreeInfo(treeTPC, "([.name:mean] && [.name:MIP])",0)->Print();
      AliTreePlayer:: selectTreeInfo(treeLogbook, "([.name:ata] && [.AxisTitle:size])",0)->Print();
    In case of problems unit test in separate file
        -- AliTreePlayerTest.C::testselectTreeInfo()
\endcode
 */
TObjArray * AliTreePlayer::selectTreeInfo(TTree* tree, TString query,Int_t verbose){
  // 1.) Make a TFormula replacing logical operands with values
  TObjArray *queryArray=query.Tokenize("&|()!");
  Int_t formulaEntries=queryArray->GetEntries();
  for (Int_t i=0;i<formulaEntries; i++){
    if (TString(queryArray->At(i)->GetName()).ReplaceAll(" ","").Length()<=0) queryArray->RemoveAt(i);
  }
  queryArray->Compress();
  formulaEntries=queryArray->GetEntries();
  if (formulaEntries<=0){
    ::Error("selectTreeInfo","Empty or wrong selection %s", query.Data());
    return 0;
  }

  TString   stringFormula=query;
  TVectorD  formValue(formulaEntries);     // boolen or weight value
  TVectorF  formType(formulaEntries);      // type of formula - name or metadata
  TVectorF  formExact(formulaEntries);     // type of formula - exact or contain
  TObjArray formMatch(formulaEntries);     // string to match
  TObjArray formMatchType(formulaEntries); // type to match
  for (Int_t ivar=0; ivar<formulaEntries; ivar++){
    TString formName=queryArray->At(ivar)->GetName();
    stringFormula.ReplaceAll(queryArray->At(ivar)->GetName(),TString::Format("x[%d]",ivar).Data());   // Logic should be improved - to match only "full strings"
    formName.ReplaceAll(" ","");
    formName.ReplaceAll("\t","");
    TObjArray *tokenArray=formName.Tokenize("[=:]");
    if (tokenArray->GetEntries()<2){
      ::Error("selectTreeInfo","Wrong subformula %s",formName.Data());
      ::Error("selectTreeInfo","Full formula was %s",query.Data());
      tokenArray->Print();
      queryArray->Print();
      delete tokenArray;
      delete queryArray;
      return 0;
    }
    formMatch[ivar]=tokenArray->At(1);
    formMatchType[ivar]=tokenArray->At(0);
    TString queryType(tokenArray->At(0)->GetName());
    queryType.ToLower();
    formType[ivar]=1;
    if (queryType.Contains(".name",TString::kIgnoreCase)){
      formType[ivar]=0;
    }
    if (formName.Contains(":")){
      formExact[ivar]=0;
    }else{
      formExact[ivar]=1;
    }
  }
  //
  TFormula *pFormula= new TFormula("printMetadataFormula",stringFormula.Data());
  if (verbose&0x4){
    ::Info("selectTreeInfo","Formula:");
    pFormula->Print();
    ::Info("selectTreeInfo","Query array:");
    queryArray->Print();
    ::Info("selectTreeInfo","To match array:");
    formMatch.Print();
    ::Info("selectTreeInfo","Exact type:");
    formExact.Print();
  }
  //
  // 2.) array loop and evaluate filters
  //
  TObjArray *selected = new TObjArray(10000); // TODO: we should get upper limit estimate from
  Int_t nTrees=1;
  if (tree->GetListOfFriends()!=NULL) nTrees+=tree->GetListOfFriends()->GetEntries();
  for (Int_t iTree=0; iTree<nTrees; iTree++){
    TTree * cTree = tree;
    if (iTree>0) cTree=tree->GetFriend(tree->GetListOfFriends()->At(iTree-1)->GetName());
    if (cTree==NULL) continue;
    for (Int_t itype=0; itype<2; itype++){
      //   TTree::GetListOfAliases() is a TList
      //   TTree::GetListOfBranches  is a TObjArray
      TSeqCollection* elemList=0;
      if (itype==0) elemList=cTree->GetListOfBranches();
      if (itype==1) elemList=cTree->GetListOfAliases();
      if (elemList==NULL) continue;
      Int_t elemEntries=elemList->GetEntries();
      for (Int_t ientry=0; ientry<elemEntries; ientry++){ // to be rewritten to TSeqColection iterator
        TString elemName(elemList->At(ientry)->GetName());
        //
        for (Int_t icheck=0; icheck<formulaEntries; icheck++){ // check logical expression
          formValue[icheck]=0;
          if (formType[icheck]==0){  // check the name
            if (formExact[icheck]==1) formValue[icheck]=(elemName==formMatch.At(icheck)->GetName());
            if (formExact[icheck]==0) formValue[icheck]=(strstr(elemName.Data(),formMatch.UncheckedAt(icheck)->GetName())!=NULL);
          }
          if (formType[icheck]==1){  // check corresponding metadata
            TObject* metaObject = TStatToolkit::GetMetadata(cTree,TString::Format("%s%s",elemName.Data(), formMatchType.UncheckedAt(icheck)->GetName()).Data());
            if (metaObject){
              TString metaName(metaObject->GetTitle());
              if (formExact[icheck]==1) formValue[icheck]=(metaName==formMatch.At(icheck)->GetName());
              if (formExact[icheck]==0) formValue[icheck]=(strstr(metaName.Data(),formMatch.UncheckedAt(icheck)->GetName())!=NULL);
            }
          }
        }
        Bool_t isSelected=pFormula->EvalPar(formValue.GetMatrixArray());
        if (isSelected){
          if (iTree==0) selected->AddLast(new TObjString(elemList->At(ientry)->GetName()));
          if (iTree>0) selected->AddLast(new TObjString(TString::Format("%s.%s",tree->GetListOfFriends()->At(iTree-1)->GetName(),elemList->At(ientry)->GetName())));
          if (verbose&0x1) printf("%s\n",elemName.Data());
        }
      }
    }
  }
  return selected;
}


///
/// \param tree         - input master tree pointer
/// \param infoType     - array in which information is queried - "branch" "alias" "metaData" and logical or
/// \param regExpFriend - regular expression patter, where to seek (default empty - master tree)
/// \param regExpTag    - regular expression tag to be queried - use non case sensitive Contain ( use as ^regExpTag$ to get full match )
/// \param verbose      - verbosity flag
/// \return             return selected string
TString AliTreePlayer::printSelectedTreeInfo(TTree *tree, TString infoType, TString regExpFriend, TString regExpTag,
                                             Int_t verbose) {
  /*
    Example usage:
      AliExternalInfo info;
      TTree * treeTPC = info.GetTree("QA.TPC","LHC15o","cpass1_pass1","QA.TRD;QA.TPC;QA.TOC;QA.TOF;Logbook");

      AliTreePlayer::printSelectedTreeInfo(treeTPC,"branch alias","","MIP",1);
      --  print infoType ="branches or alias" containing  regExpTag == MIP in the master tree 
      AliTreePlayer::printSelectedTreeInfo(treeTPC,"branch","QA.TRD","Eff",1);
      --  print infoType ="branches" containing  regExpTag == Eff in the QA.TRD tree
      AliTreePlayer::printSelectedTreeInfo(treeTPC,"branch",".*","run",1)
      --  print infoType ="branches" containing  regExpTag == Eff in the all  friend trees
      AliTreePlayer::printSelectedTreeInfo(treeTPC,"branch",".*","^run$",1)
      --  print infoType ="branches" containing  regExpTag == run in the all  friend trees
      
      In case arrays or function should be selected - use aliases:
        treeTPC->SetAlias("meanMIPvsSectorArray","meanMIPvsSector.fElements")'
        treeTPC->SetAlias("runTypeName","runType.GetName()")
        AliTreePlayer::printSelectedTreeInfo(treeTPC,"branch alias",".*","meanMIPvsSectorArray",1); ==> meanMIPvsSectorArray
        AliTreePlayer::printSelectedTreeInfo(treeTPC,"branch alias ",".*","Name$",1) ==> runTypeName
   */
  TString result = "";
  TList *treeFriends = tree->GetListOfFriends();
  Int_t ntrees = 1 + ((treeFriends != NULL) ? treeFriends->GetEntries() : 0);
  TPRegexp pregExpFriend = regExpFriend;
  TPRegexp pregExpTag = regExpTag;
  const char *dataTypes[3] = {"branch", "alias", "metaData"};
  for (Int_t itree = 0; itree < ntrees; itree++) {
    TTree *currentTree = 0;
    if (itree == 0) {
      currentTree = tree;
      if (pregExpFriend.Match(currentTree->GetName()) == 0) continue;
    } else {
      if (pregExpFriend.Match(treeFriends->At(itree - 1)->GetName()) == 0) continue;
      currentTree = ((TFriendElement * )(treeFriends->At(itree - 1)))->GetTree();
    }

    if (verbose & 0x1) {
      ::Info("printSelectedTreeInfo", "tree %s selected", currentTree->GetName());
    }
    for (Int_t iDataType = 0; iDataType < 3; iDataType++) {
      if (infoType.Contains(dataTypes[iDataType], TString::kIgnoreCase) == 0) continue;
      TList *selList = 0;
      if (iDataType == 0) selList = (TList *) currentTree->GetListOfBranches();
      if (iDataType == 1) selList = (TList *) currentTree->GetListOfAliases();
      if (iDataType == 2 && tree->GetUserInfo())
        selList = (TList * )(currentTree->GetUserInfo()->FindObject("metaTable"));

      if (selList == NULL) continue;
      Int_t selListEntries = selList->GetEntries();
      for (Int_t iEntry = 0; iEntry < selListEntries; iEntry++) {
        if (pregExpTag.Match(selList->At(iEntry)->GetName()) <= 0) continue;
        if (verbose & 0x1) {
          ::Info(" printSelectedTreeInfo", "%s.%s", currentTree->GetName(), selList->At(iEntry)->GetName());
        }
        if (result.Length() > 0) result += ":";
        if (itree > 0) {
          result += treeFriends->At(itree - 1)->GetName();
          result += ".";
        }
        result += selList->At(iEntry)->GetName();
      }
    }
  }
  return result;
}

///
/// \param tree            - input tree
/// \param what            - variables to export (any valid TTree formula or class (in case ob JSON export))
/// \param where           - selection criteria as in the TTree::Draw
/// \param firstentry      - first entry in tree to export
/// \param nentries        - number of entries to
/// \param outputFormat    - output format (csv,json, elastic json - for bulk export,html table)
/// \param outputName
/// \return                - status
Int_t  AliTreePlayer::selectWhatWhereOrderBy(TTree * tree, TString what, TString where, TString /*orderBy*/,  Int_t firstentry, Int_t nentries, TString outputFormat, TString outputName){
  //
  // Select entry similar to the SQL select
  //    tree instead - SQL FROM statement used
  //         - it is supposed that all information is contained in the tree and related friend trees
  //         - build tree with friends done in separate function
  //    what -
  // code inspired by the TTreePlay::Scan but another treatment of arrays needed
  // but wih few changes realted to different output formatting
  // NOTICE:
  // for the moment not all parameters used
  // Example usage:
  /*
    AliTreePlayer::selectWhatWhereOrderBy(treeTPC,"run:Logbook.run:meanMIP:meanMIPele:meanMIPvsSector.fElements:fitMIP.fElements","meanMIP>0", "", 0,10,"html","qatpc.html");
  */

  if (tree==NULL || tree->GetPlayer()==NULL){
    ::Error("AliTreePlayer::selectWhatWhereOrderBy","Input tree not defiend");
    return -1;
  }
  // limit number of entries - shorter version of the TTreePlayer::GetEntriesToProcess - but not fully correct ()
  if (firstentry +nentries >tree->GetEntriesFriend()) nentries=tree->GetEntriesFriend()-firstentry;
  if (tree->GetEntryList()){ //
    if (tree->GetEntryList()->GetN()<nentries) nentries=tree->GetEntryList()->GetN();
  }
  //
  Int_t  tnumber = -1;
  Bool_t isHTML=outputFormat.Contains("html",TString::kIgnoreCase );
  Bool_t isCSV=outputFormat.Contains("csv",TString::kIgnoreCase);
  Bool_t isElastic=outputFormat.Contains("elastic",TString::kIgnoreCase);
  Bool_t isJSON=outputFormat.Contains("json",TString::kIgnoreCase)||isElastic;

  //
  FILE *default_fp = stdout;
  if (outputName.Length()>0){
    default_fp=fopen (outputName.Data(),"w");
  }
  TObjArray *fArray=what.Tokenize(":");
  Int_t      nCols=fArray->GetEntries();
  TObjArray *fFormulaList       = new TObjArray(nCols+1);
  TTreeFormula ** rFormulaList  = new TTreeFormula*[nCols+1];
  TObjString **printFormatList  = new TObjString*[nCols];     // how to format variables
  TObjString **columnNameList   = new TObjString*[nCols];
  TObjString **outputFormatList = new TObjString*[nCols];     // root header formatting
  Bool_t isIndex[nCols];
  Bool_t isParent[nCols];
  TClass*  isClass[nCols];       // is object
  TPRegexp indexPattern("^%I"); // pattern for index
  TPRegexp parentPattern("^%P"); // pattern for parent index
  for (Int_t iCol=0; iCol<nCols; iCol++){
    isClass[iCol]=NULL;
    rFormulaList[iCol]=NULL;
    TObjArray * arrayDesc = TString(fArray->At(iCol)->GetName()).Tokenize(";");
    if (arrayDesc->GetEntries()<=0) {
      ::Error("AliTreePlayer::selectWhatWhereOrderBy","Invalid descriptor %s", arrayDesc->At(iCol)->GetName());
      //return -1;
    }
    TString  fieldName=arrayDesc->At(0)->GetName();    // variable content
    if (tree->GetBranch(fieldName.Data())){
      if (TString(tree->GetBranch(fieldName.Data())->GetClassName()).Length()>0) {
        isClass[iCol]=TClass::GetClass(tree->GetBranch(fieldName.Data())->GetClassName());
      }
    }
    if (fieldName.Contains(indexPattern)){             // variable is index
      isIndex[iCol]=kTRUE;
      indexPattern.Substitute(fieldName,"");
    }else{
      isIndex[iCol]=kFALSE;
    }
    if (fieldName.Contains(parentPattern)){             // variable is parent
      isParent[iCol]=kTRUE;
      parentPattern.Substitute(fieldName,"");
    }else{
      isParent[iCol]=kFALSE;
    }
    TTreeFormula * formula = NULL;
    TNamed* htmlTag=TStatToolkit::GetMetadata(tree,(fieldName+".html").Data());
    Bool_t isFormulaF= (fieldName[0]=='#') || htmlTag!=NULL;

    if (!isFormulaF) {
      formula=new TTreeFormula(fieldName.Data(), fieldName.Data(), tree);
    }else {
      TString fstring = fieldName;
      if (htmlTag) {
        fstring=htmlTag->GetTitle();
      } else {
        fieldName.Replace(0, 1, "");
        if (tree->GetAlias(fieldName.Data())) {
          fstring = tree->GetAlias(fieldName.Data());
        }
      }
      formula=new AliTreeFormulaF(fieldName.Data(), fstring.Data(), tree);
    }

    if (formula->GetTree()==NULL){
      ::Error("AliTreePlayer::selectWhatWhereOrderBy","Invalid formula %s, parsed from the original string %s",fieldName.Data(),what.Data());
      if (isJSON==kFALSE) return -1;
    }
    TString  printFormat="";                       // printing format          - column 1 - use default if not specified
    TString  colName=arrayDesc->At(0)->GetName();  // variable name in ouptut  - column 2 - using defaut ("variable name") as in input
    TString  outputFormat="";                      // output column format specification for TLeaf (see also reading using TTree::ReadFile)
    if (arrayDesc->At(1)!=NULL){  // format
      printFormat=arrayDesc->At(1)->GetName();
    }else{
      if (formula->IsInteger()) {
        printFormat="1.20";
      }else{
        printFormat="1.20";
      }
    }
    if (arrayDesc->At(2)!=NULL)  {
      colName=arrayDesc->At(2)->GetName();
    }else{
      colName=arrayDesc->At(0)->GetName();
    }
    //colName = (arrayDesc->GetEntries()>1) ? arrayDesc->At(2)->GetName() : arrayDesc->At(0)->GetName();
    if (arrayDesc->At(3)!=NULL){  //outputFormat (for csv)
      outputFormat= arrayDesc->At(3)->GetName();
    }else{
      outputFormat="/D";
      if (formula->IsInteger()) outputFormat="/I";
      if (formula->IsString())  outputFormat="/C";
    }
    fFormulaList->AddAt(formula,iCol);
    rFormulaList[iCol]=formula;
    printFormatList[iCol]=new TObjString(printFormat);
    outputFormatList[iCol]=new TObjString(outputFormat);
    columnNameList[iCol]=new TObjString(colName);
  }
  TTreeFormula *select = new TTreeFormula("Selection",where.Data(),tree);
  fFormulaList->AddLast(select);
  rFormulaList[nCols]=select;


  Bool_t hasArray = kFALSE;
  Bool_t forceDim = kFALSE;
  for (Int_t iCol=0; iCol<nCols; iCol++){
    if (rFormulaList[iCol]!=NULL) rFormulaList[iCol]->UpdateFormulaLeaves();
    if (rFormulaList[iCol]->GetManager()==NULL) continue;   // TODO - check of needed for AliTreeFormulaF
    // if ->GetManager()->GetMultiplicity()>0 mean there is at minimum one array
    switch( rFormulaList[iCol]->GetManager()->GetMultiplicity() ) {
      case  1:
      case  2:
        hasArray = kTRUE;
        forceDim = kTRUE;
        break;
      case -1:
        forceDim = kTRUE;
        break;
      case  0:
        break;
    }
  }


  // print header
  if (isHTML){
    fprintf(default_fp,"%s", "<table class=\"display\" cellspacing=\"0\" width=\"100%\">\n"); // add metadata info
    fprintf(default_fp,"\t<thead class=\"header\">\n"); // add metadata info
    fprintf(default_fp,"\t<tr>\n"); // add metadata info
    for (Int_t iCol=0; iCol<nCols; iCol++){
      TNamed *tHeadName = TStatToolkit::GetMetadata(tree,TString(columnNameList[iCol]->GetName())+".thead");
      TNamed *tooltipName = TStatToolkit::GetMetadata(tree,TString(columnNameList[iCol]->GetName())+".tooltip");
      TString sthName=(tHeadName!=NULL) ? tHeadName->GetTitle() : columnNameList[iCol]->GetName();
      TNamed *descriptionName = TStatToolkit::GetMetadata(tree,TString(columnNameList[iCol]->GetName())+".headerTooltip");
      fprintf(default_fp, "\t\t<th");
      if (tooltipName!=0)  fprintf(default_fp, " class=\"tooltip\" data-tooltip=\"%s\"", tooltipName->GetTitle());
      if (descriptionName!=0)  fprintf(default_fp, " title=\"%s\"", descriptionName->GetTitle());
      fprintf(default_fp, ">%s</th>\n",  sthName.Data()); // add metadata info
//      if (tooltipName==NULL){
//        fprintf(default_fp, "\t\t<th>%s</th>\n", sthName.Data()); // add metadata info
//      }else {
//        fprintf(default_fp, "\t\t<th class=\"tooltip\" data-tooltip=\"%s\">%s</th>\n", tooltipName->GetTitle(), sthName.Data()); // add metadata info
//      }
    }
    fprintf(default_fp,"\t</tr>\n"); // add metadata info
    fprintf(default_fp,"\t</thead>\n"); // add metadata info
    //
    fprintf(default_fp,"\t<tfoot class=\"header\">\n"); // add metadata info
    fprintf(default_fp,"\t<tr>\n"); // add metadata info
    for (Int_t iCol=0; iCol<nCols; iCol++){
      fprintf(default_fp,"\t\t<th>%s</th>\n",columnNameList[iCol]->GetName()); // add metadata info
    }
    fprintf(default_fp,"\t</tr>\n"); // add metadata info
    fprintf(default_fp,"\t</tfoot>\n"); // add metadata info
    fprintf(default_fp,"\t<tbody>\n"); // add metadata info
  }
  if (isCSV && outputFormat.Contains("csvroot")){
    // add header info
    for (Int_t iCol=0; iCol<nCols; iCol++){
      fprintf(default_fp,"%s%s",columnNameList[iCol]->GetName(), outputFormatList[iCol]->GetName());
      if (iCol<nCols-1)  {
        fprintf(default_fp,":");
      }else{
        fprintf(default_fp,"\n"); // add metadata info
      }
    }
  }
  if (isJSON&&!isElastic){
    fprintf(default_fp,"{\n\t \"tree\": [\n "); // add metadata info
  }

  Int_t selected=0;
  for (Int_t ientry=firstentry; ientry<firstentry+nentries; ientry++){
    Int_t entryNumber = tree->GetEntryNumber(ientry);
    if (entryNumber < 0) break;
    Long64_t localEntry = tree->LoadTree(entryNumber);
    //
    if (tnumber != tree->GetTreeNumber()) {
      tnumber = tree->GetTreeNumber();
      for(Int_t iCol=0;iCol<nCols;iCol++) {
        rFormulaList[iCol]->UpdateFormulaLeaves();
      }
      select->UpdateFormulaLeaves();
    }
    if (select) {
      //      if (select->EvalInstance(inst) == 0) {
      if (select->EvalInstance(0) == 0) {  // for the moment simplified version of selection - not treating "array" selection
        continue;
      }
    }
    selected++;
    // if json out
    if (isJSON){
      if (selected>1){
        if (isElastic) {
          fprintf(default_fp,"}\n{\"index\":{\"_id\": \"");
        }
        else{
          fprintf(default_fp,"},\n{\n");
        }
      }else{
        if (isElastic){
          fprintf(default_fp,"{\"index\":{\"_id\": \"");
        }else{
          fprintf(default_fp,"{\n");
        }
      }
      for (Int_t icol=0; icol<nCols; icol++){
        TString obuffer="";
        if (isClass[icol]){
          const char*bname=rFormulaList[icol]->GetName();
          TBranch *br = tree->GetBranch(bname);
          br->GetEntry(entryNumber); // not sure if the entry is loaded
          void **ppobject=(void**)br->GetAddress();
          if (ppobject) {
            obuffer = TBufferJSON::ConvertToJSON(*ppobject, isClass[icol]);
            if (isElastic) obuffer.ReplaceAll("\n","");
            //fprintf(default_fp,"%s",clbuffer.Data());
          }
        }
        if (rFormulaList[icol]->GetTree()==NULL) continue;
        Int_t nData=rFormulaList[icol]->GetNdata();
        if (isJSON==kFALSE){
          fprintf(default_fp,"\t\"%s\":",rFormulaList[icol]->GetName());
        }else{
          if (isIndex[icol]==kFALSE && isParent[icol]==kFALSE){
            TString fieldName(rFormulaList[icol]->GetName());
            if (isElastic) {
              if (isClass[icol]) fieldName.Remove(fieldName.Length() - 1);
              fieldName.ReplaceAll(".", "%_");
            }
            if (icol>0 && isIndex[icol-1]==kFALSE && isParent[icol-1]==kFALSE ){
              fprintf(default_fp,"\t,\"%s\":",fieldName.Data());
            }else{
              fprintf(default_fp,"\t\"%s\":",fieldName.Data());
            }
          }
        }
        if (nData<=1){
          if ((isIndex[icol]==kFALSE)&&(isParent[icol]==kFALSE)){
            if (isClass[icol]!=NULL){
              fprintf(default_fp,"%s",obuffer.Data());
            }else {
              if ( (isElastic || isJSON)  && rFormulaList[icol]->IsString()) {
                fprintf(default_fp, "\t\"%s\"",
                        rFormulaList[icol]->PrintValue(0, 0, printFormatList[icol]->GetName()));
              } else {
                fprintf(default_fp, "\t%s",
                        rFormulaList[icol]->PrintValue(0, 0, printFormatList[icol]->GetName()));
              }
            }
          }
          if (isIndex[icol]){
            fprintf(default_fp,"%s",rFormulaList[icol]->PrintValue(0,0,printFormatList[icol]->GetName()));
            if (isIndex[icol+1]){
              fprintf(default_fp,".");
            }else
            if (isParent[icol+1]==kFALSE){
              fprintf(default_fp,"\"}}\n{");
            }
          }
          if (isParent[icol]==kTRUE){
            fprintf(default_fp,"\", \"parent\": \"%s\"}}\n{",rFormulaList[icol]->PrintValue(0,0,printFormatList[icol]->GetName()));
          }
        }else{
          fprintf(default_fp,"\t[");
          for (Int_t iData=0; iData<nData;iData++){
            fprintf(default_fp,"%f",rFormulaList[icol]->EvalInstance(iData));
            if (iData<nData-1) {
              fprintf(default_fp,",");
            }   else{
              fprintf(default_fp,"]");
            }
          }
        }
        //	if (icol<nCols-1 && (isIndex[icol]==kFALSE && isParent[icol]==kFALSE) ) fprintf(default_fp,",");
        //fprintf(default_fp,"\n");
      }
    }

    //
    if (isHTML){
      fprintf(default_fp,"<tr>\n");
      for (Int_t icol=0; icol<nCols; icol++){
        Int_t nData=rFormulaList[icol]->GetNdata();
        if (nData<=1){
          fprintf(default_fp,"\t<td>%s</td>",rFormulaList[icol]->PrintValue(0,0,printFormatList[icol]->GetName()));
        }else{
          fprintf(default_fp,"\t<td>");
          for (Int_t iData=0; iData<nData;iData++){
            fprintf(default_fp,"%f",rFormulaList[icol]->EvalInstance(iData));
            if (iData<nData-1) {
              fprintf(default_fp,",");
            }else{
              fprintf(default_fp,"</td>");
            }
          }
        }
        fprintf(default_fp,"\n");
      }
      fprintf(default_fp,"</tr>\n");
    }
    //
    if (isCSV){   ///dummy
      for (Int_t icol=0; icol<nCols; icol++){  // formula loop
        Int_t nData=rFormulaList[icol]->GetNdata();
        if (nData<=1){
          TString qValue=rFormulaList[icol]->PrintValue(0,0,printFormatList[icol]->GetName());
          qValue.ReplaceAll("\t","");
          fprintf(default_fp,"%s\t",qValue.Data());
        }else{
          for (Int_t iData=0; iData<nData;iData++){  // array loo
            fprintf(default_fp,"%f",rFormulaList[icol]->EvalInstance(iData));
            if (iData<nData-1) {
              fprintf(default_fp,",");
            }else{
              fprintf(default_fp,"\t");
            }
          }
        }
      }
      fprintf(default_fp,"\n");
    }
    if (isCSV){   ///dummy
      Int_t nDataMax=1;
      for (Int_t icol=0; icol<nCols; icol++) {       // formula loop
        Int_t nData = rFormulaList[icol]->GetNdata();
        nDataMax=TMath::Max(nDataMax,nData);
      }
      for (Int_t iData=0; iData<nDataMax;iData++){  // array loop
        for (Int_t icol=0; icol<nCols; icol++){  // formula loop
          Int_t nData=rFormulaList[icol]->GetNdata();
          if (nData<=1){
            fprintf(default_fp,"%s\t",rFormulaList[icol]->PrintValue(0,0,printFormatList[icol]->GetName()));
          }else{
              fprintf(default_fp,"%f",rFormulaList[icol]->EvalInstance(iData));
              fprintf(default_fp,"\t");

          }
        }
        fprintf(default_fp,"\n");
      }
    }



  }
  if (isJSON) fprintf(default_fp,"}\t]\n}\n");
  if (isHTML){
    fprintf(default_fp,"\t</tbody>"); // add metadata info
    fprintf(default_fp,"</table>"); // add metadata info
  }

  if (default_fp!=stdout) fclose (default_fp);
  return selected;
}


/// Get the enumeration type from a string
/// \param stat     - enumuration tyle as a string
/// \return         Get the enumeration type from a string
Int_t AliTreePlayer::GetStatType(const TString &stat){

  if(!stat.CompareTo("median",TString::kIgnoreCase)){
    return kMedian;
  } else if(!stat.CompareTo("medianLeft",TString::kIgnoreCase)){
    return kMedianLeft;
  } else if(!stat.CompareTo("medianRight",TString::kIgnoreCase)){
    return kMedianRight;
  } else if(!stat.CompareTo("RMS",TString::kIgnoreCase)){
    return kRMS;
  } else if(!stat.CompareTo("Mean",TString::kIgnoreCase)){
    return kMean;
  }else if(!stat.CompareTo("Max",TString::kIgnoreCase)){
    return kMax;
  }else if(!stat.CompareTo("Min",TString::kIgnoreCase)){
    return kMin;
  } else if(stat.BeginsWith("LTMRMS",TString::kIgnoreCase)){ // Least trimmed RMS, argument is LTMRMS0.95 or similar
    return kLTMRMS;
  } else if(stat.BeginsWith("LTM",TString::kIgnoreCase)){ // Least trimmed mean, argument is LTM0.95 or similar
    return kLTM;
  } else {
    ::Error("GetStatType()","Cannot decode string \"%s\"."
            " Use one of \"median\", \"medianLeft\", \"medianRight\", \"RMS\", or \"Mean\". "
            " Also supported is \"LTM\", or \"LTMRMS\", which should be succeeded by a float like"
            " \"LTM0.95\" or an integer (interpreted as percentage) like \"LTMRMS95\" to specify "
            " the fraction of data to be kept."
            " Use a colon separated list like \"median:medianLeft:medianRight:RMS\""
            " as the fifth argument to the AddStatInfo().",stat.Data());
    return kUndef;
  }
}
///
/// \param treeLeft
/// \param treeRight
/// \param refQuery
/// \param deltaT
/// \param statString
/// \param maxEntries
void  AliTreePlayer::AddStatInfo(TTree* treeLeft,  TTree * treeRight , const TString refQuery, Double_t deltaT,
                                 const TString statString,
                                 Int_t maxEntries){

  //
  // 1.) Get variables from the right tree, sort them according to the coordinate
  //
  Int_t entries  = treeRight->Draw(refQuery.Data(),"","goffpara",maxEntries);
  if(entries<1){
    ::Error("AddStatInfo","No matching entries for query");
    return;
  }else if (entries>treeRight->GetEstimate()){
    treeRight->SetEstimate(entries*2);
    entries  = treeRight->Draw(refQuery.Data(),"","goffpara",maxEntries);
  }
  Int_t * indexArr = new Int_t[entries];
  TMath::Sort(entries, treeRight->GetV1(), indexArr,kFALSE);
  Double_t * coordArray  = new Double_t[entries];
  for (Int_t icoord=0; icoord<entries; icoord++) coordArray[icoord]=treeRight->GetV1()[indexArr[icoord]];
  //
  // 2.) Attach the median,.. of the variables to the left tree
  //
  // The first token is the coordinate
  TString var;
  Ssiz_t from=0;
  if(!refQuery.Tokenize(var,from,":")){
    ::Error("AddStatInfo","Cannot tokenize query \'%s\'. Use colon separated list"
            ,refQuery.Data());
    delete[]indexArr;indexArr=0;
    delete[]coordArray;coordArray=0;
    return;
  }
  Int_t entriesCoord =  treeLeft->GetEntries();
  TBranch * br = treeLeft->GetBranch(var.Data());
  Int_t coordValue;
  br->SetAddress(&coordValue);
  //TLeaf *coordLeaf = (TLeaf*)(br->GetListOfLeaves()->At(0));

  while(refQuery.Tokenize(var,from,":")){
    //   var="valueAnodeRaw.fElements[0]";
    TString stat;
    Ssiz_t fromStat=0;
    while(statString.Tokenize(stat,fromStat,":")){
      // stat = "median"; stat="medianLeft"; stat="medianRight";
      //      Int_t statType=TStatToolkit::GetStatType(stat);
      Int_t statType=GetStatType(stat);
      if(statType==kUndef)continue;
      //printf("\n\n\n--- StatType %d ---\n\n\n\n",statType);

      // In case of LMS or LMR determine the fraction
      Float_t frac=0.;
      if(statType==kLTM || statType==kLTMRMS){ // stat="LTM0.95" or "LTMRMS0.95" similar
        TString tmp= stat;
        tmp.ReplaceAll("LTMRMS","");
        tmp.ReplaceAll("LTM","");
        frac = tmp.Atof(); // atof returns 0.0 on error
        if(frac>1)frac/=100.; // allows "LTM95" for 95%
      }

      // Determine the offset in coordinate to the left and right
      Double_t leftOffset=-1e99,rightOffset=-1e99;
      if(statType==kMedian||statType==kRMS || statType==kMean
         || statType==kLTM || statType==kLTMRMS || statType==kMin || statType==kMax ){ // symmetric
        leftOffset=-deltaT;rightOffset=deltaT;
      } else if(statType==kMedianLeft){ //left
        leftOffset=-2.*deltaT;rightOffset=0.;
      } else if(statType==kMedianRight){ //right
        leftOffset=0.;rightOffset=2.*deltaT;
      }

      TString brName=Form("%s_%s",var.Data(),stat.Data());
      brName.ReplaceAll("[","_");
      brName.ReplaceAll("]","_");
      brName.ReplaceAll(".","_");
      Double_t statValue=0;
      if (treeLeft->GetDirectory()) treeLeft->GetDirectory()->cd();
      TBranch *brToFill = treeLeft->Branch(brName.Data(),&statValue, (brName+"/D").Data());
      TVectorD dvalues(entries);

      for (Int_t icoord=0; icoord<entriesCoord; icoord++){
        // Int_t icoord=0
        br->GetEntry(icoord);
        Double_t startCoord=coordValue+leftOffset;
        Double_t endCoord       =coordValue+rightOffset;
        //Double_t startCoord=coordLeaf->GetValue()+leftOffset;
        //Double_t endCoord	=coordLeaf->GetValue()+rightOffset;

        // Binary search finds last element smaller or equal
        Int_t index0=BinarySearchSmaller(entries, coordArray, startCoord)+1;
        Int_t index1=TMath::BinarySearch(entries, coordArray, endCoord) +1;
        statValue=0;
        if (index1>=0 && index0>=0){
          //if (values.capacity()<index1-index0)  values.reserve(1.5*(index1-index0));     //resize of needed
          if (index1-index0<=0){
            index0--;
            index1++;
            if (index0<0) index0=0;
            if (index1>=entries) index1=entries-1;
          }
          for (Int_t i=0; i<index1-index0; i++){
            dvalues[i]=treeRight->GetV2()[indexArr[i+index0]];
          }
          // Calculate
          if(statType==kMedian||statType==kMedianLeft||statType==kMedianRight){
            statValue=TMath::Median(index1-index0, dvalues.GetMatrixArray());
          } else if(statType==kRMS){
            statValue=TMath::RMS(index1-index0, dvalues.GetMatrixArray());
          } else if(statType==kMean){
            statValue=TMath::Mean(index1-index0, dvalues.GetMatrixArray());
          }else if(statType==kMin){
            statValue=TMath::MinElement(index1-index0, dvalues.GetMatrixArray());
          }else if(statType==kMax){
            statValue=TMath::MaxElement(index1-index0, dvalues.GetMatrixArray());
          } else if(statType==kLTM){
            TVectorD params(7);
            try {
              TStatToolkit::LTMUnbinned(index1 - index0, dvalues.GetMatrixArray(), params, frac);
            }
            catch (int e)
            {
              cout << "TStatToolkit::LTMUnbinned. Catch Exception Nr. " << e << '\n';
              continue;
            }
            statValue = params[1];
          } else if(statType==kLTMRMS){
            TVectorD params(7);
            try {
              TStatToolkit::LTMUnbinned(index1 - index0, dvalues.GetMatrixArray(), params, frac);
            }
            catch (int e)
            {
              cout << "TStatToolkit::LTMUnbinned. Catch Exception Nr. " << e << '\n';
              continue;
            }
            statValue = params[2];
          } else{
            ::Error("AddStatInfo()","String %s StatType %d not implemented",stat.Data(),statType);
          }
        }
        brToFill->Fill();
      } // Loop over coordinates
      brToFill->FlushBaskets();
    } // Loop over stat types
  } // Loop over variables
  treeLeft->Write();
  delete[]indexArr;indexArr=0;
  delete[]coordArray;coordArray=0;

}


/// Create list of histograms specified by selection
/// Should be rough equivalent of the "ALICE train" TTree->Draw();
///  a.) Data are read only once
///  b.) values expression are reused (evaluated only once)
///  c.) Axis labelling and names of variables extracted from the tree metadata (.AxisTitle)
/// * default cut
///   * default selection applied common for all histograms (can be empty)
///
/// * hisString : - semicolomn separated string
///   * his0;his1; ...; hisN
/// * histogram syntax:
///    * var0:var1:...:<#weight>>>hisName(bins0,min0,max0,bins1,min0,min, minValue,maxValue)
///    * Syntax:
///      * vari are histogramming expression
///      * weight (or cut) entry is optional
///        * default cut is always applied, weight is applied on top
///    * ranges syntax:
///      *  nbins,max,min where max and min are double or format strings
///        * in case format string % specified using (Fraction, mean,meanFraction, rms, rmsFraction)
///          *  %fraction.sigma
///          *  #cumulant
///          *  range for bin content can be specified in the same format (by default is not set)
/*!
##### CPU time to process one histogram or set of histograms (in particular case of esdTrack queries) is the same - and it is determined (90 %) by tree->GetEntry
\code
  THn * his0= (THn*)hisArray->At(0);
  his0->Projection(0)->Draw("");
  tree->SetLineColor(2);
  TStopwatch timer; tree->Draw("esdTrack.Pt()","(esdTrack.fFlags&0x40)>0&&esdTrack.fTPCncls>70","same",60000); timer.Print();
\endcode
*/

/// \param tree         - input tree
/// \param hisString    - selection string
/// \param defaultCut   - default selection applied common for all histograms (can be empty)
/// \param firstEntry   - first entry to process
/// \param lastEntry    - last entry to process
/// \param chunkSize    - chunk size
/// \param verbose      - verbosity
/// \return             - TObjArray of N-dimensional histograms
/*!
#### Example usage:
\code
    chunkSize=10000;
    verbose=7;
    chinput=gSystem->ExpandPathName("$NOTES/JIRA/PWGPP-227/data/2016/LHC16t/000267161/pass1_CENT_wSDD/filteredLocal.list");
    TString defaultCut="esdTrack.fTPCncls>70";
    TTree *tree=(TTree*)AliXRDPROOFtoolkit::MakeChain(chinput, "highPt", 0, 1000000000,0);
    TString hisString="";
    hisString+="esdTrack.Pt():#esdTrack.fTPCncls>70>>hisPtAll(100,0,30);";
    hisString+="esdTrack.GetAlpha():#esdTrack.fTPCncls>70>>hisAlpha(90,-3.2,3.2);";
    hisString+="esdTrack.GetTgl():#esdTrack.fTPCncls>70>>hisTgl(20,-1.2,1.2);";
    hisString+="esdTrack.Pt():esdTrack.GetAlpha():esdTrack.GetTgl():#esdTrack.fTPCncls>70>>hisPtPhiThetaAll(100,0,30,90,-3.2,3.2,20,-1.2,1.2);";
    hisString+="esdTrack.Pt():#(esdTrack.fFlags&0x4)>0>>hisPtITS(100,1,10);";
    hisString+="esdTrack.fIp.Pt():#(esdTrack.fFlags&0x4)>0>>hisPtTPCOnly(100,1,10);";
    TStopwatch timer; hisArray = AliTreePlayer::MakeHistograms(tree, hisString, "(esdTrack.fFlags&0x40)>0&&esdTrack.fTPCncls>70",0,60000,100000); timer.Print();
\endcode
 */
TObjArray  * AliTreePlayer::MakeHistograms(TTree * tree, TString hisString, TString defaultCut, Int_t firstEntry, Int_t lastEntry, Int_t chunkSize, Int_t verbose){
  //  Algorithm:
  //  1.) Analyze formula, create minimal formula expression
  //  2.) Event loop
  //  3.) Filling of histograms
  //

  const Int_t kMaxDim=10;
  if (tree== nullptr){
    ::Error("AliTreePlayer::MakeHistograms","Tree=0");
    return 0;
  }
  Int_t entriesAll=tree->GetEntries();
  if (chunkSize<=0) chunkSize=entriesAll;
  if (lastEntry>entriesAll) lastEntry=entriesAll;
  if (verbose&0x7){
    ::Info("AliTreePlayer::MakeHistograms","firstEntry %d", firstEntry);
    ::Info("AliTreePlayer::MakeHistograms","lastEntry %d", lastEntry);
    ::Info("AliTreePlayer::MakeHistograms","entriesAll %d", entriesAll);
  }
  //
  TObjArray *hisDescriptionList=hisString.Tokenize(";");
  Int_t nHistograms = hisDescriptionList->GetEntries();
  TObjArray * hisArray = new TObjArray(nHistograms);
  TObjArray * hisDescriptionArray=new TObjArray(nHistograms); // OWNER
  TObjArray * hisFormulaArray=new TObjArray(nHistograms);     // array of TFomula arrays - Not OWNER
  TObjArray * hisWeightArray=new TObjArray(nHistograms);     // array of TFomula arrays - Not OWNER
  TArrayI     hisDims(nHistograms);
  //
  Int_t nExpressions=hisString.CountChar(':')+hisString.CountChar(';')+1;
  TObjArray * formulaArray   = new TObjArray(nExpressions);    // array of all expressions  - OWNER
  TString queryString = "";
  Int_t hisSizeFull=0;
  //
  //  1.) Analyze formula, book list of TObjString
  //
  Bool_t isOK=kTRUE;
  for (Int_t iHis=0; iHis<nHistograms; iHis++){
    TString hisDescription = hisDescriptionList->At(iHis)->GetName();
    Int_t hisIndex=hisDescription.Index(">>");
    if (hisIndex<=0) {
      isOK=kFALSE;
      ::Error("AliTreePlayer::MakeHistograms","Invalid expression %s",hisDescription.Data());
      break;
    }else{
      hisDescriptionArray->AddAtAndExpand(new TObjString(((hisDescriptionList->At(iHis)->GetName()))+(hisIndex+2)),iHis);
    }
    hisDescription.Remove(hisIndex);
    TObjArray *hisDimArray=hisDescription.Tokenize(":");
    Int_t nDims=hisDimArray->GetEntries();
    if (nDims<=0){
      isOK=kFALSE;
      ::Error("AliTreePlayer::MakeHistograms","Invalid description %s",hisDescription.Data());
      delete hisDimArray;
      break;
    }
    TObjArray * varArray = new TObjArray(nDims);
    if (hisDimArray->At(nDims-1)->GetName()[0]=='#'){
      TString formulaName=&((hisDimArray->At(nDims-1)->GetName())[1]);
      nDims-=1;
      TObjString *tFormula = new TObjString(formulaName.Data());
      hisWeightArray->AddAt(tFormula,iHis);
      if (formulaArray->FindObject(formulaName.Data())==NULL){
        formulaArray->AddLast(tFormula);
        varArray->AddAt(tFormula,nDims);
      }
    }
    for (Int_t iDim=0; iDim<nDims;iDim++){
      TObjString *tFormula =  (TObjString*) (formulaArray->FindObject(hisDimArray->At(iDim)->GetName()));
      if (tFormula==NULL){
        tFormula = new TObjString(hisDimArray->At(iDim)->GetName());
        formulaArray->AddLast(tFormula);
      }
      varArray->AddAt(tFormula,iDim);
    }
    hisFormulaArray->AddAt(varArray,iHis);
    hisDims[iHis]=nDims;
    delete hisDimArray;
  }
  queryString="";
  Int_t nFormulas=formulaArray->GetEntries();
  for (Int_t iFor=0; iFor<nFormulas; iFor++){
    queryString+=formulaArray->At(iFor)->GetName();
    queryString+=":";
  }
  queryString+=formulaArray->At(nFormulas-1)->GetName();
  if (verbose&0x2) hisDescriptionArray->Print();
  if (verbose&0x4) formulaArray->Print();
  //
  //  2.) Event loop
  //
  Int_t tNumber=-1;
  for (Int_t bEntry=firstEntry; bEntry<lastEntry; bEntry+=chunkSize){  // chunks loop
    Int_t toQuery=TMath::Min(chunkSize, lastEntry-bEntry);
    Int_t qLength = tree->Draw(queryString,defaultCut,"goffpara",toQuery, bEntry); // query varaibles
    if (qLength>tree->GetEstimate()){
      tree->SetEstimate(qLength*1.5);
      qLength = tree->Draw(queryString,defaultCut,"goffpara",chunkSize, bEntry);
    }

    //  2.2 fill histograms if  not yet done
    if (hisArray->GetEntriesFast()==0){  // book histograms if not yet done
      for (Int_t iHis=0; iHis<nHistograms; iHis++){
        if (hisDescriptionArray->At(iHis)==NULL){
          ::Error("AliTreePlayer::MakeHistograms", "Empty description %d",iHis);
          continue;
        }
        TString hisDescription= hisDescriptionArray->At(iHis)->GetName();
        TString  varDecription=hisDescriptionList->At(iHis)->GetName();
        TObjArray * descriptionArray=hisDescription.Tokenize("(,)");
        TObjArray * varArray= TString(hisDescriptionList->At(iHis)->GetName()).Tokenize(":");
        Int_t nLength=descriptionArray->GetEntries();
        if ((nLength-1)/3 < hisDims[iHis]){
          ::Error("AliTreePlayer::MakeHistograms", "Histogram dimension Mismatch %s", hisDescriptionArray->At(iHis)->GetName());
          return NULL;
        }
        if (varArray->GetEntries()<hisDims[iHis]){
          ::Error("AliTreePlayer::MakeHistograms", "Variable mismatch %s", hisDescriptionArray->At(iHis)->GetName());
          return NULL;
        }
        TString hName(descriptionArray->At(0)->GetName());
        THnBase * his = 0;
        Int_t nBins[kMaxDim];
        Double_t xMin[kMaxDim], xMax[kMaxDim];
        for (Int_t iDim=0; iDim<hisDims[iHis]; iDim++){
          nBins[iDim]= TString(descriptionArray->At(3*iDim+1)->GetName()).Atoi();
          if (descriptionArray->At(3*iDim+2)->GetName()[0]!='%'){
            xMin[iDim]= TString(descriptionArray->At(3*iDim+2)->GetName()).Atof();
          }else{
            if (descriptionArray->At(3*iDim+2)->GetName()[1]=='A'){ // %A - alias describing range
              TTreeFormula falias("falias",&(descriptionArray->At(3*iDim+2)->GetName()[2]),tree);
              xMin[iDim]=falias.EvalInstance();
            }
          }
          if (descriptionArray->At(3*iDim+3)->GetName()[0]!='%'){
            xMax[iDim]= TString(descriptionArray->At(3*iDim+3)->GetName()).Atof();
          }else{
            if (descriptionArray->At(3*iDim+3)->GetName()[1]=='A'){ // %A - alias describing range
              TTreeFormula falias("falias",&(descriptionArray->At(3*iDim+3)->GetName()[2]),tree);
              xMax[iDim]=falias.EvalInstance();
            }
          }
          if (xMax[iDim]<=xMin[iDim]){
            ::Error("AliTreePlayer::MakeHistograms","Invalid hstogram range specification for histogram  %s: %s\t%s",hisDescription.Data(), \
		    descriptionArray->At(3*iDim+2)->GetName(), descriptionArray->At(3*iDim+3)->GetName() );
          }
        }
        THnF * phis = new THnF(hName.Data(),hName.Data(), hisDims[iHis],nBins, xMin,xMax);
        hisArray->AddAt(phis,iHis);
        if (verbose&0x1) {
          ::Info("AliTreePlayer::MakeHistograms","%s: size=%lld",hisDescription.Data(), phis->GetNbins());
        }
        hisSizeFull+= phis->GetNbins();
        for (Int_t iDim=0;iDim<hisDims[iHis]; iDim++){
          phis->GetAxis(iDim)->SetName(varArray->At(iDim)->GetName());
          phis->GetAxis(iDim)->SetTitle(varArray->At(iDim)->GetName());
          TNamed *axisTitle=TStatToolkit::GetMetadata(tree,TString::Format("%s.AxisTitle",varArray->At(iDim)->GetName()).Data());
          if (axisTitle)  phis->GetAxis(iDim)->SetTitle(axisTitle->GetTitle());
        }
      }
      if (verbose&0x1) {
        ::Info("AliTreePlayer::MakeHistograms","Total size=%d",hisSizeFull);
      }
    }
    //    2.3 fill histograms
    if (tree->GetVal(0)==NULL){
      ::Error(" AliTreePlayer::MakeHistograms","Info  not available"); // TODO - fix the logic
      continue;
    }
    Double_t values[kMaxDim];
    for (Int_t iHis=0; iHis<nHistograms; iHis++){
      Int_t indeces[kMaxDim+1];
      TObjArray *formulaArrayHis = (TObjArray*) (hisFormulaArray->At(iHis));
      for (Int_t iVec=0; iVec<formulaArrayHis->GetEntriesFast(); iVec++){
        indeces[iVec]= formulaArray->IndexOf(formulaArray->FindObject(formulaArrayHis->At(iVec)->GetName()));
      }
      Int_t indexW=-1;
      if (hisWeightArray->GetEntriesFast()>=iHis){
        if (hisWeightArray->UncheckedAt(iHis)!=NULL){
          if (hisWeightArray->UncheckedAt(iHis)->GetName()){
            indexW= formulaArray->IndexOf(formulaArray->FindObject(hisWeightArray->UncheckedAt(iHis)->GetName()));
          }else{
            ::Error("xxx","Problem to find %s", hisWeightArray->UncheckedAt(iHis)->GetName());
          }
        }
      }
      THnBase * his = (THnBase*) hisArray->UncheckedAt(iHis);
      for (Int_t cEvent=0; cEvent<qLength; cEvent++){
        for (Int_t iDim=0; iDim<hisDims[iHis]; iDim++){
          if (tree->GetVal(indeces[iDim])==NULL){
            ::Error(" AliTreePlayer::MakeHistograms","Info %d not avaliable",iDim); // TODO - fix the logic
            continue;
          }
          values[iDim]=tree->GetVal(indeces[iDim])[cEvent];
        }
        Double_t weight=(indexW<0)? 1: tree->GetVal(indexW)[cEvent];
        if (weight>0) his->Fill(values,weight);
      }
    }
  }
  //
  delete hisDescriptionArray;
  delete formulaArray;
  return hisArray;

}




TPad *  AliTreePlayer::DrawHistograms(TPad  * pad, TObjArray * hisArray, TString drawExpression,  TObjArray *keepArray, Int_t verbose){
  //
  //
  // Example usage:
  /*
    TPad *pad= 0;    TString drawExpression="";
    drawExpression="[1,1,1]:" 
    drawExpression+="hisPtAll(0,10)(0)(errpl);hisPtITS(0,10)(0)(err);hisPtPhiThetaAll(0,10,-3.2,3.2,-1.2,1.2)(0)(err):";
    drawExpression+="hisAlpha(-3.2,3.2)(0)(errpl);hisPtPhiThetaAll(0,10,-3.2,3.2,-1.2,1.2)(1)(err):";
    drawExpression+="hisTgl(-1,1)(0)(errpl);hisPtPhiThetaAll(0,10,-3.2,3.2,-1.2,1.2)(2)(err):";
    pad = AliTreePlayer::DrawHistograms(0,hisArray,drawExpression);
  */
  TString projType[8]={"f-mean","f-rms","f-ltm","f-ltmsigma","f-gmean","f-grms","f-median","f-gmean"};

  TObjArray *drawList=drawExpression.Tokenize(":");
  // structure pad
  TString padDescription=drawList->At(0)->GetName();
  if (pad==NULL){
    static Int_t counter=0;
    pad = new TCanvas(TString::Format("canvasCounter%d",counter).Data(), drawExpression,1000,800);
    counter++;
  }
  // divide pads
  Int_t nPads=0, nRows=0;
  TObjArray *padRows=padDescription.Tokenize("[](),");
  nRows=padRows->GetEntries();
  for (Int_t iRow=0; iRow<nRows; iRow++){
    Int_t nCols=TString(padRows->At(iRow)->GetName()).Atoi();
    for (Int_t iCol=0; iCol<nCols; iCol++){
      pad->cd();
      TPad * newPad=new TPad(Form("pad%d",nPads),Form("pad%d",nPads),iCol/Double_t(nCols),(nRows-iRow-1)/Double_t(nRows),(iCol+1)/Double_t(nCols),(nRows-iRow)/Double_t(nRows));
      newPad->Draw();
      nPads++;
      newPad->SetNumber(nPads);
    }
  }
  delete padRows;
  //
  //
  TPRegexp isPadOption("^%O");
  Bool_t isLogY=kFALSE;
  for (Int_t iPad=0; iPad<nPads; iPad++){
    Int_t nGraphs=0, nHistos=0;
    if (drawList->At(iPad+1)==NULL) break;
    //TVirtualPad  *cPad = 
    TVirtualPad *cPad = pad->cd(iPad+1);
    TLegend * legend = new TLegend(cPad->GetLeftMargin()+0.02,0.7,1-cPad->GetRightMargin()-0.02 ,1-cPad->GetTopMargin()-0.02, TString::Format("Pad%d",iPad));
    legend->SetNColumns(2);
    legend->SetBorderSize(0);
    TString padSetup=drawList->At(iPad+1)->GetName();
    TString padOption="";
    Bool_t isTimeX=kFALSE, isTimeY=kFALSE;

    if (padSetup.Contains(isPadOption)){
      padOption=TString(&padSetup[2]);
      padOption.Remove(padOption.First(";"));
      padOption.ToLower();
      if (padOption.Contains("logy")) {
        pad->cd(iPad+1)->SetLogy();
        isLogY=kTRUE;
      }
      if (padOption.Contains("logx")) {
        pad->cd(iPad+1)->SetLogx();
      }
      if (padOption.Contains("gridy")) pad->cd(iPad+1)->SetGridy();
      if (padOption.Contains("gridx")) pad->cd(iPad+1)->SetGridx();
      if (padOption.Contains("timex")) isTimeX=kTRUE;
      if (padOption.Contains("timey")) isTimeY=kTRUE;
      padSetup.Remove(0,padSetup.First(';')+1);
    }

    TObjArray * padDrawList= padSetup.Tokenize(";");
    Double_t hisMin=0, hisMax=-1;
    TH1 * hisToDraw=0;
    TGraphErrors * grToDraw=0;
    //
    for (Int_t ihis=0; ihis<padDrawList->GetEntries(); ihis++){
      TObjArray *hisDescription= TString(padDrawList->At(ihis)->GetName()).Tokenize("()");
      THn * his = (THn*)hisArray->FindObject(hisDescription->At(0)->GetName());
      if (his==NULL) continue;
      if (verbose&0x4){
        ::Info("AliTreePlayer::DrawHistograms","Pad %d. Processing his %s",iPad, hisDescription->At(0)->GetName());
      }
      Int_t ndim=his->GetNdimensions();
      TString rangeDescription(hisDescription->At(1)->GetName());
      Int_t ndimRange= (rangeDescription.CountChar(','));
      if (ndimRange>0) ndimRange=ndimRange/2+1;
      if (ndimRange>0){
        TObjArray *rangeArray=rangeDescription.Tokenize(",");
        for (Int_t iDim=0; iDim<ndimRange; iDim++){
          if (rangeArray->At(iDim*2)->GetName()[0]=='U') {
            Double_t min=TString(&(rangeArray->At(iDim*2)->GetName()[1])).Atof();
            Double_t max=TString(&(rangeArray->At(iDim*2+1)->GetName()[1])).Atof();
            if (verbose&0x8){
              ::Info("AliTreePlayer::DrawHistograms","Pad %d. %s.GetAxis(%d).SetRangeUser(%f,%f).",iPad,hisDescription->At(0)->GetName(), iDim, min,max);
            }
            his->GetAxis(iDim)->SetRangeUser(min,max);
          }else{
            Int_t min=TString((rangeArray->At(iDim*2)->GetName())).Atoi();
            Int_t max=TString((rangeArray->At(iDim*2+1)->GetName())).Atoi();
            if (verbose&0x8){
              ::Info("AliTreePlayer::DrawHistograms","Pad %d. %s.GetAxis(%d).SetRange(%d,%d).",iPad,hisDescription->At(0)->GetName(), iDim, min,max);
            }
            his->GetAxis(iDim)->SetRange(min,max);
          }
        }
        delete rangeArray;
      }
      TString drawOption = hisDescription->At(3)->GetName();
      drawOption.ToLower();
      TString projString=hisDescription->At(2)->GetName();
      Int_t nDims = projString.CountChar(',')+1;
      TH1 * hProj =0;
      TGraphErrors*gr=0;
      //
      if (nDims==1) {
        hProj=his->Projection(projString.Atoi());
        //hProj->SetName(Form("Pad%d_His%d",iPad,nHistos));  //
        nHistos++;
      }
      if (nDims==2) {
        Int_t dim0 = projString.Atoi();
        Int_t dim1 = TString(&(projString[2])).Atoi();
        TH2* his2D =his->Projection(dim0,dim1);
        for (Int_t iProj=0; iProj<8; iProj++){
          if (drawOption.Contains(projType[iProj])){
            gr=TStatToolkit::MakeStat1D(his2D,0,1.0,iProj,21+ihis,ihis+1);
            //gr->SetName(Form("gr_Pad%d_Graph%d",iPad,nGraphs));  //
            nGraphs++;
            // gr->SetTitle(Form("gr_Pad%d_Graph%d",iPad,nGraphs));  //
            gr->SetTitle(padDrawList->At(ihis)->GetName());
            gr->GetXaxis()->SetTitle(his2D->GetXaxis()->GetTitle());
            gr->GetYaxis()->SetTitle(his2D->GetYaxis()->GetTitle());
            drawOption.ReplaceAll(projType[iProj].Data(),"");
          }
        }
        delete his2D;
      }
      if (gr){
        Double_t grMinI=TMath::MinElement(gr->GetN(),gr->GetY())-3.*TMath::Median(gr->GetN(),gr->GetEY());
        Double_t grMaxI=TMath::MaxElement(gr->GetN(),gr->GetY())+3.*TMath::Median(gr->GetN(),gr->GetEY());
        if (hisMax<hisMin) {hisMin=grMinI;  hisMax=grMaxI;}
        if (hisMax<grMaxI) hisMax=grMaxI;
        if (hisMin>grMinI) hisMin=grMinI;
        if (ihis==0)  {
          grToDraw = gr;
          gr->Draw((drawOption+"a").Data());
          legend->AddEntry(gr,"","p");
        }
        else{
          gr->Draw(drawOption.Data());
          legend->AddEntry(gr,"", "p");
        }
        if (keepArray) keepArray->AddLast(gr);
      }
      if (hProj){
        hProj->SetMarkerColor(ihis+1);
        hProj->SetLineColor(ihis+1);
        hProj->SetMarkerStyle(21+ihis);
        if (keepArray) keepArray->AddLast(hProj);
        //
        if (hisMax<hisMin) {
          hisMin=hProj->GetMinimum();
          hisMax=hProj->GetMaximum();
        }
        if (hisMax<hProj->GetMaximum()) hisMax=hProj->GetMaximum();
        if (hisMin>hProj->GetMinimum()) hisMin=hProj->GetMinimum();
        if (ihis==0)  {
          hisToDraw = hProj;
          hProj->Draw((TString(hisDescription->At(3)->GetName())+"").Data());
          legend->AddEntry(hProj);
        }
        else{
          hProj->Draw((TString(hisDescription->At(3)->GetName())+"same").Data());
          legend->AddEntry(hProj);
        }
      }
    }
    pad->cd(iPad+1);
    if (hisToDraw!=NULL){
      hisToDraw->SetMaximum(hisMax+(hisMax-hisMin)/2.);
      if (hisMin<=0) hisMin=TMath::Min(0.01, hisMax*0.01);
      hisToDraw->SetMinimum(hisMin);
      if (isLogY){
        hisToDraw->SetMaximum(hisMax*TMath::Max(10.,(hisMax/hisMin)/4.));
      }

      if ((verbose&0x8)>0){
        ::Info("AliTreePlayer::DrawHistograms:","Pad %d. %s  SetMinimum(%f). SetMaximum(%f)",iPad,padDrawList->At(0)->GetName(), hisMin,hisMax+(hisMax-hisMin)/2.);
      }
      if (isTimeX) hisToDraw->GetXaxis()->SetTimeDisplay(1);
      if (isTimeY) hisToDraw->GetYaxis()->SetTimeDisplay(1);
      pad->cd(iPad+1)->Modified();
      pad->cd(iPad+1)->Update();
      legend->Draw("same");
    }
    if (grToDraw!=NULL){
      grToDraw->SetMaximum(hisMax+(hisMax-hisMin)/2.);
      grToDraw->SetMinimum(hisMin-(hisMax-hisMin)/3.);
      if (isTimeX) grToDraw->GetXaxis()->SetTimeDisplay(1);
      if (isTimeY) grToDraw->GetYaxis()->SetTimeDisplay(1);

      if ((verbose&0x8)>0){
        ::Info("AliTreePlayer::DrawHistograms:","Pad %d. %s  SetMinimum(%f). SetMaximum(%f)",iPad,padDrawList->At(0)->GetName(), hisMax+(hisMax-hisMin)/2.,hisMin-(hisMax-hisMin)/3.);
      }

      pad->cd(iPad+1)->Modified();
      pad->cd(iPad+1)->Update();
      legend->Draw("same");
    }
    pad->cd(iPad+1);
  }
  return pad;

}


///\brief Fill tree with information specified in varList of TTreeFormulas
/// Used to cache CPU consuming formulas
/// To consider:
///    1.) In case input tree is "flat" - not array output tree can be used as a friend ....
///    2.) In case of not flat tree user should made appropiate tree->SetEstimate to allocate buffers otherwise only partial results available
///
/// \param tree         - TTree with input
/// \param varList      - list of TTreeFormulas to export
/// \param outFile      -  output file name
/// \param outTree      - output tree name
/// \param selection    - tree selection
/// \param firstEntry   - first entry to export
/// \param nEntries     - number of nEntries to export
void AliTreePlayer::MakeCacheTree(TTree * tree, TString varList, TString outFile, TString outTree, TCut selection, Int_t nEntries, Int_t firstEntry, const char *fileMode){
  TTreeSRedirector *pcstream = new TTreeSRedirector(outFile,fileMode);
  if (tree->GetEstimate()<tree->GetEntries()) tree->SetEstimate(tree->GetEntries());
  Int_t estimate=tree->GetEstimate();
  Int_t entries=0;
  if (firstEntry>=0 && nEntries>0) {
    entries = tree->Draw(varList.Data(),selection,"goffpara",nEntries,firstEntry);
  }else{
    entries = tree->Draw(varList.Data(),selection,"goffpara");
  }
  TObjArray * varName=varList.Tokenize(":");
  const Int_t nVars=varName->GetEntries();
  Double_t vars[nVars];
  TTree *treeOut=NULL;
  if (entries>estimate) entries=estimate;
  for (Int_t iPoint=0; iPoint <entries; iPoint++){
    for (Int_t iVar=0; iVar<nVars; iVar++){
      vars[iVar]=tree->GetVal(iVar)[iPoint];
      if (iPoint==0) (*pcstream)<<outTree.Data()<<TString::Format("%s=",varName->At(iVar)->GetName()).Data()<<vars[iVar];
    }
    if (iPoint==0) {
      (*pcstream)<<outTree.Data()<<"\n";
      treeOut=((*pcstream)<<outTree.Data()).GetTree();
    }else{
      treeOut->Fill();
    }
  }
  delete pcstream;
}

/// \brief LoadTrees and and append friend trees
/// \param inputDataList          - command returning input data list - line with #tag:value for the metadata definition
/// \param chRegExp               - regular expression  - trees to include
/// \param chNotReg               - regular expression  - trees to exclude
/// \param inputFileSelection     - : separated input file selection
/// \param axisAlias              - axis names
/// \param axisTitle              - axis titles
/// \return                       - resulting tree
TTree* AliTreePlayer::LoadTrees(const char *inputDataList, const char *chRegExp, const char *chNotReg, TString inputFileSelection, TString axisAlias, TString axisTitle, Int_t verbose) {
  //
  TPRegexp regExp(chRegExp);
  TPRegexp notReg(chNotReg);
  TObjArray *regExpArray = inputFileSelection.Tokenize(":");
  TObjArray *residualMapList = gSystem->GetFromPipe(inputDataList).Tokenize("\n");
  Int_t nFiles = residualMapList->GetEntries();
  TTree *treeBase = 0;
  TObjArray *arrayAxisAlias = axisAlias.Tokenize(":");
  TObjArray *arrayAxisTitle = axisTitle.Tokenize(":");
  std::map<string,string> tagValue;
  ///
  Int_t nSelected=0;
  for (Int_t iFile = 0; iFile < nFiles; iFile++) {
    tagValue.clear();
    /// read metadata #tag:value information - signed by # at the beginning of the line
    for (Int_t jFile=iFile;jFile<nFiles; jFile++){
      TString name = residualMapList->At(jFile)->GetName();
      if (name[0]=='#') {
        Int_t first=name.First(":");
        TString tag=name(1,name.First(":")-1);
        TString value=name(name.First(":")+1,name.Length());
        tagValue[tag.Data()]=value.Data();
      }else{
        iFile=jFile;
        break;
      }
    }
    TString fileName = residualMapList->At(iFile)->GetName(); ///
    Bool_t isSelected = kTRUE;
    // check if the file was selected
    if (regExpArray->GetEntries() > 0) {
      isSelected = kFALSE;
      for (Int_t entry = 0; entry < regExpArray->GetEntries(); entry++) {
        TPRegexp reg(regExpArray->At(entry)->GetName());
        if (reg.Match(fileName)) isSelected |= kTRUE;
      }
    }
    if (!isSelected) continue;
    ::Info("LoadTrees","Load file\t%s", fileName.Data());
    TString description = "";
    TString option="";
    if (fileName.Contains("http")) option="cacheread";
    TFile *finput = TFile::Open(fileName.Data(),option.Data());
    if (finput == NULL) {
      ::Error("MakeResidualDistortionReport", "Invalid file name %s", fileName.Data());
      continue;
    }
    TList *keys = finput->GetListOfKeys();
    Int_t isLegend = kFALSE;
    ///
    for (Int_t iKey = 0; iKey < keys->GetEntries(); iKey++) {
      if (regExp.Match(keys->At(iKey)->GetName()) == 0) continue;     // is selected
      if (notReg.Match(keys->At(iKey)->GetName()) != 0) continue;     // is rejected
      TTree *tree = (TTree *) finput->Get(keys->At(iKey)->GetName()); // better to use dynamic cast
      if (treeBase == NULL) {
        TFile *finput2 = TFile::Open(fileName.Data(),option.Data());
        treeBase = (TTree *) finput2->Get(keys->At(iKey)->GetName());
      }
      TString fileTitle=tagValue["Title"];
      if (fileTitle.Length()){
        treeBase->AddFriend(tree, TString::Format("%s.%s", fileTitle.Data(), keys->At(iKey)->GetName()).Data());
      }else{
        treeBase->AddFriend(tree, TString::Format("%s", keys->At(iKey)->GetName()).Data());
      }
      Int_t entriesF = tree->GetEntries();
      Int_t entriesB = treeBase->GetEntries();
      if (entriesB == entriesF) {
        if (verbose>0) ::Info("InitMapTree", "%s\t%s.%s:\t%d\t%d", treeBase->GetName(), fileName.Data(), keys->At(iKey)->GetName(), entriesB, entriesF);
      } else {
        ::Error("InitMapTree", "%s\t%s.%s:\t%d\t%d", treeBase->GetName(), fileName.Data(), keys->At(iKey)->GetName(), entriesB, entriesF);
      }
    }
    /// add metadata to all files ....
  }
  if (treeBase) {
    for (Int_t iAxis = 0; iAxis < arrayAxisAlias->GetEntries(); iAxis++) {
      treeBase->SetAlias(arrayAxisAlias->At(iAxis)->GetName(), TString::Format("%s%dCenter", arrayAxisAlias->At(iAxis)->GetName(), iAxis + 1).Data());
      TStatToolkit::AddMetadata(treeBase, TString::Format("%s.AxisTitle", arrayAxisAlias->At(iAxis)->GetName()).Data(), arrayAxisTitle->At(iAxis)->GetName());
    }
  }
  return treeBase;
}



/// AddMetadata to the input tree - see https://alice.its.cern.ch/jira/browse/ATO-290
/// \param tree         - input tree
/// \param varTagName   - tag to register
/// \param varTagValue  - value
/// \return             - return has list
/// Currently supported metadata
///* Drawing  :
///  * <varName>.AxisTitle
///  *<varName>.Legend
///  * <varname>.Color
///  * <varname>.MarkerStyle
///* This metadata than can be used by the TStatToolkit
///  * TStatToolkit::MakeGraphSparse
///  * TStatToolkit::MakeGraphErrors
///
THashList*  AliTreePlayer::AddMetadata(TTree* tree, const char *varTagName,const char *varTagValue){
  if (!tree) return NULL;
  THashList * metaData = (THashList*) tree->GetUserInfo()->FindObject("metaTable");
  if (metaData == NULL){
    metaData=new THashList;
    metaData->SetName("metaTable");
    tree->GetUserInfo()->AddLast(metaData);
  }
  if (varTagName!=NULL && varTagValue!=NULL){
    TNamed * named = TStatToolkit::GetMetadata(tree, varTagName,NULL,kTRUE);
    if (named==NULL){
      metaData->AddLast(new TNamed(varTagName,varTagValue));
    }else{
      named->SetTitle(varTagValue);
    }
  }
  return metaData;
}

/// Get metadata description
/// In case metadata contains friend part - friend path os returned as prefix
/// Metadata are supposed to be added into tree using TStatToolkit::AddMetadata() function
/// \param tree          - input tree
/// \param varTagName    - tag name
/// \param prefix        - friend prefix in case metadata are in friend tree
/// \param fullMatch     - request full match in varTagName (to enable different metadata for tree and friend tree)
/// \return              - metadata description
/// TODO: too many string operations - to be speed up using char array arithmetic
TNamed* AliTreePlayer::GetMetadata(TTree* tree, const char *varTagName, TString * prefix,Bool_t fullMatch){
  //

  if (!tree) return 0;
  TTree * treeMeta=tree;

  THashList * metaData = (THashList*) treeMeta->GetUserInfo()->FindObject("metaTable");
  if (metaData == NULL){
    metaData=new THashList;
    metaData->SetName("metaTable");
    tree->GetUserInfo()->AddLast(metaData);
    return 0;
  }
  TNamed * named = (TNamed*)metaData->FindObject(varTagName);
  if (named || fullMatch) return named;

  TString metaName(varTagName);
  Int_t nDots= metaName.CountChar('.');
  if (prefix!=NULL) *prefix="";
  if (nDots>1){ //check if friend name expansion needed
    while (nDots>1){
      TList *fList= treeMeta->GetListOfFriends();
      if (fList!=NULL){
        Int_t nFriends= fList->GetEntries();
        for (Int_t kf=0; kf<nFriends; kf++){
          TPRegexp regFriend(TString::Format("^%s.",fList->At(kf)->GetName()).Data());
          if (metaName.Contains(regFriend)){
            treeMeta=treeMeta->GetFriend(fList->At(kf)->GetName());
            regFriend.Substitute(metaName,"");
            if (prefix!=NULL){
              (*prefix)+=fList->At(kf)->GetName();
              // (*prefix)+=""; /// FIX use prefix as it is
            }
          }
        }
      }
      if (nDots == metaName.CountChar('.')) break;
      nDots=metaName.CountChar('.');
    }
  }

  named = (TNamed*)metaData->FindObject(metaName.Data());
  return named;

}
/// change to next pad in canvas - used in TTree::Draw()
/// Example:
///   draw histogram fulfilling criteria
///   tree->Draw("histCM.Draw(\"colz\"):histCM.GetXaxis()->SetRangeUser(130,170)","!T0_isOK&&abs(peakPosition-150)<20&&AliTreePlayer::nextPad()","goffpara",100,0);
Int_t AliTreePlayer::nextPad(){
  /// used for Tree draw queries
  static int counter=0;
  if (!gPad) return kFALSE;
  counter++;
  gPad->GetMother()->cd(counter);
  if (gPad->GetNumber()<counter) counter=0;
  gPad->Update();
  return counter;
}