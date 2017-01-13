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


/*
  Set of functions to extend functionality of the TTreePlayer
  data source (tree+friend trees)  configured independently in AliExternalInfo

  Functionality:
  * function to support metadata and collumns annotation
  * filtering function (branch, aliases metadata)
  * select function  - export to (json, csv,  html, JIRA) + metadata (not yet implemented)
  ** subset of metadata for columns included in select export outputName.metadata
  * select function  - for root trees
  ** to be done using TTree::CopyTree  functionality swithing ON/OFF selected branches - metadata should be exported automatically (to check) 
  
  See  example usage in the test macro AliTreePlayerTest.C
  * AliTreePlayerTest.C::testAll();
  * AliTreePlayerTest.C::testSelectMetadata()
  * AliTreePlayerTest.C::testselectTreeInfo()
  * AliTreePlayerTest.C:testselectWhatWhereOrderByForTRD()

*/

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

ClassImp(AliTreePlayer)

//_____________________________________________________________________________
AliTreePlayer::AliTreePlayer(const char *name, const char *title)
         :TNamed(name, title)
{
  //
  //
}

//_____________________________________________________________________________



///
/// Selected metadata fil filing query
/// retun ObjArray of selected metadatas
/// param

			  
TObjArray  *  AliTreePlayer::selectMetadata(TTree * tree, TString query, Int_t verbose){
  /*
    query -  case sensitive matching is done using Contains method (e.g N(Axis)<=N(xis) )
          -  WILL  be  case insesitive 
   
    Example usage:
    query="[AxisTitle]";                     // only metadata with metadata axis existing
    query="[class]";                         // only metadata with metadata class existing
    query="[class==\"Base&&Log&&Stat\"]";    // metadata class contains
    //
    AliTreePlayer::selectMetadata(treeLogbook, "[class==\"Logbook&&(Stat||Base)\"]",0)->Print();
    AliTreePlayer::selectMetadata(testTree, "[Axis==\"counts\"]",0)->Print();
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
      if (verbose&0x1){
	printf("%s:%s:%d\n",metaData->At(ientry)->GetName(),metaData->At(ientry)->GetTitle(),isSelected);
      }
    }
  }
  return selected;
}

TObjArray * AliTreePlayer::selectTreeInfo(TTree* tree, TString query,Int_t verbose){
  /*
    Find all branches containing "Data string" and having class base
    search follow only branches - doesn not enter into classes
    - WILL be not case sensitive (should be part of unit test) 
    == - for exact match
    :  - for contain
    Implemented using TFormula interpeter:
    TString query="([.name==Data])&&([.class==base]||[.class==Data])";
    ==>
    TFormula   (x[0])&&(x[1]||x[2])

    Example cases :
      AliTreePlayer::selectTreeInfo(treeTPC, "([.name:mean] && [.name:MIP])",0)->Print();
      AliTreePlayer:: selectTreeInfo(treeLogbook, "([.name:ata] && [.AxisTitle:size])",0)->Print();
    In case of problems unit test in separate file 
        -- AliTreePlayerTest.C::testselectTreeInfo()
    //


  */
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




TString  AliTreePlayer::printSelectedTreeInfo(TTree*tree, TString infoType,  TString regExpFriend, TString regExpTag, Int_t verbose){
  //
  //
  //
  //   tree         - input master tree pointer
  //   infoType     - array in which information is queried - "branch" "alias" "metaData" and logical or 
  //   regExpFriend - regeular expression patter, where to seek (default empty - master tree) 
  //   regExpTag    - regula expression tag to be queried - use non case sensitive Contain ( use as ^regExpTag$ to get full match )
  //
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
  TString result="";
  TList * treeFriends = tree->GetListOfFriends();
  Int_t ntrees = 1+( (treeFriends!=NULL)?treeFriends->GetEntries():0 );
  TPRegexp  pregExpFriend=regExpFriend;
  TPRegexp  pregExpTag=regExpTag;
  const char* dataTypes[3]={"branch","alias", "metaData"};
  for (Int_t itree=0; itree<ntrees; itree++){
    TTree * currentTree = 0;
    if (itree==0){
      currentTree=tree;
      if (pregExpFriend.Match(currentTree->GetName())==0) continue;
    }else{
      if (pregExpFriend.Match(treeFriends->At(itree-1)->GetName())==0) continue;
      currentTree = ((TFriendElement*)(treeFriends->At(itree-1)))->GetTree();
    }

    if (verbose&0x1){
      ::Info("printSelectedTreeInfo","tree %s selected", currentTree->GetName());
    }
    for (Int_t iDataType=0; iDataType<3; iDataType++){
      if (infoType.Contains(dataTypes[iDataType], TString::kIgnoreCase)==0) continue;
      TList * selList = 0;
      if (iDataType==0) selList=(TList*)currentTree->GetListOfBranches();
      if (iDataType==1) selList=(TList*)currentTree->GetListOfAliases();
      if (iDataType==2 && tree->GetUserInfo()) selList=(TList*)(currentTree->GetUserInfo()->FindObject("metaTable"));

      if (selList==NULL) continue;
      Int_t selListEntries=selList->GetEntries();
      for (Int_t iEntry=0; iEntry<selListEntries; iEntry++){
	if (pregExpTag.Match(selList->At(iEntry)->GetName())<=0) continue;
	if (verbose&0x1){
	  ::Info(" printSelectedTreeInfo","%s.%s", currentTree->GetName(),selList->At(iEntry)->GetName());
	}
	if (result.Length()>0) result+=":";
	if (itree>0) 	{
	  result+=treeFriends->At(itree-1)->GetName(); 
	  result+=".";
	}
	result+=selList->At(iEntry)->GetName();
      }
    }
  }
  return result;
}


void AliTreePlayer::selectWhatWhereOrderBy(TTree * tree, TString what, TString where, TString /*orderBy*/,  Int_t firstentry, Int_t nentries, TString outputFormat, TString outputName){
  //
  // Select entry similar to the SQL select
  //    tree instead - SQL FROM statement used
  //         - it is supposed that all infomtation is contained in the tree and related friend trees
  //         - build tree with friends done in separate function
  //    what -
  // code inspired by the TTreePlay::Scan but another tretment of arrays needed
  // but wih few changes realted to different output formating
  // NOTICE:
  // for the moment not all parameters used
  // Example usage:
  /*
    AliTreePlayer::selectWhatWhereOrderBy(treeTPC,"run:Logbook.run:meanMIP:meanMIPele:meanMIPvsSector.fElements:fitMIP.fElements","meanMIP>0", "", 0,10,"html","qatpc.html");
    
  */
  if (tree==NULL || tree->GetPlayer()==NULL){
    ::Error("AliTreePlayer::selectWhatWhereOrderBy","Input tree not defiend");
    return;
  }
  // limit number of entries - shorter version of the TTreePlayer::GetEntriesToProcess - but not fully correct ()
  if (firstentry +nentries >tree->GetEntriesFriend()) nentries=tree->GetEntriesFriend()-firstentry;
  if (tree->GetEntryList()){ //
    if (tree->GetEntryList()->GetN()<nentries) nentries=tree->GetEntryList()->GetN();
  }
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

  for (Int_t iCol=0; iCol<nCols; iCol++){
    TObjArray * arrayDesc = TString(fArray->At(iCol)->GetName()).Tokenize(";");
    if (arrayDesc->GetEntries()<=0) {
      ::Error("AliTreePlayer::selectWhatWhereOrderBy","Invalid descriptor %s", arrayDesc->At(iCol)->GetName());
      return;
    }
    TString  fname=arrayDesc->At(0)->GetName();    // variable content
    TTreeFormula * formula = new TTreeFormula(fname.Data(), fname.Data(), tree);
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
    fFormulaList->AddLast(formula);
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
    rFormulaList[iCol]->UpdateFormulaLeaves();
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

  //
  Int_t  tnumber = -1;
  Bool_t isHTML=outputFormat.Contains("html",TString::kIgnoreCase );
  Bool_t isCSV=outputFormat.Contains("csv",TString::kIgnoreCase);
  Bool_t isElastic=outputFormat.Contains("elastic",TString::kIgnoreCase);
  Bool_t isJSON=outputFormat.Contains("json",TString::kIgnoreCase)||isElastic;

  // print header
  if (isHTML){
    fprintf(default_fp,"<table>"); // add metadata info
    fprintf(default_fp,"<tr>"); // add metadata info
    for (Int_t iCol=0; iCol<nCols; iCol++){
      fprintf(default_fp,"<th>%s</th>",columnNameList[iCol]->GetName()); // add metadata info
    }
    fprintf(default_fp,"<tr>"); // add metadata info
  }
  if (isCSV){
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
    // if json out
    if (isJSON){
      fprintf(default_fp,"{\n");
      for (Int_t icol=0; icol<nCols; icol++){
	Int_t nData=rFormulaList[icol]->GetNdata();
	fprintf(default_fp,"\t\"%s\":",rFormulaList[icol]->GetName());
	if (nData<=1){
	  fprintf(default_fp,"\t%f",rFormulaList[icol]->EvalInstance(0));
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
	if (icol<nCols-1) fprintf(default_fp,",");
	fprintf(default_fp,"\n");
      }
      fprintf(default_fp,"}\n");
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
    if (isCSV){
      for (Int_t icol=0; icol<nCols; icol++){  // formula loop
	Int_t nData=rFormulaList[icol]->GetNdata();
	if (nData<=1){
	  fprintf(default_fp,"%s\t",rFormulaList[icol]->PrintValue(0,0,printFormatList[icol]->GetName()));
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
  }
  if (isHTML){
    fprintf(default_fp,"</table>"); // add metadata info
  }
  if (default_fp!=stdout) fclose (default_fp);

}












// Get the enumeration type from a string
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
//__________________________________________________________
void  AliTreePlayer::AddStatInfo(TTree* treeLeft,  TTree * treeRight , const TString refQuery, Double_t deltaT,
		 const TString statString,
		 Int_t maxEntries){  
 
  
  //
  // 1.) Get variables from the right tree, sort them according to the coordinate
  //
  treeRight->SetEstimate(treeRight->GetEntries()*540);
  Int_t entries  = treeRight->Draw(refQuery.Data(),"","goffpara",maxEntries);
  if(entries<1){
    ::Error("AddStatInfo","No matching entries for query");
    return;
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
	 || statType==kLTM || statType==kLTMRMS){ // symmetric
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
      TBranch *brToFill = treeLeft->Branch(brName.Data(),&statValue, (brName+"/D").Data());
      TVectorD dvalues(entries/10+1);
				 
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
	  } else if(statType==kLTM){
	    TVectorD params(7);
	    TStatToolkit::LTMUnbinned(index1-index0, dvalues.GetMatrixArray()
				      ,params,frac);
	    statValue = params[1];
	  } else if(statType==kLTMRMS){
	    TVectorD params(7);
	    TStatToolkit::LTMUnbinned(index1-index0, dvalues.GetMatrixArray()
				      ,params,frac);
	    statValue = params[2];
	  } else{
	    ::Error("AddStatInfo()","String %s StatType %d not implemented",stat.Data(),statType);
	  }

	  // // debug
	  // printf("startCoord=%.3e, endCoord=%.3e, median=%.3e, index0=%d, index1=%d "
	  // 	 ,startCoord,endCoord,statValue,index0,index1);
	  // for(Int_t i=0;i<index1-index0;i++){
	  //   printf("dvalues[%d]=%.3e "
	  // 	   ,i,dvalues[i]);
	  // }
	  // printf("\n");
	    

	}
	brToFill->Fill();
      } // Loop over cordinates
      brToFill->FlushBaskets();
    } // Loop over stat types
  } // Loop over variables
  treeLeft->Write();
  delete[]indexArr;indexArr=0;
  delete[]coordArray;coordArray=0;

}
