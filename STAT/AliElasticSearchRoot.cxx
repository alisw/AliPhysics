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


#include <iostream>
#include "TSystem.h"
#include "TString.h"
#include "TError.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "AliElasticSearchRoot.h"

/// jg filters
std::map<std::string,std::string> AliElasticSearchRoot::fJQFilters;

/// default constructor
AliElasticSearchRoot::AliElasticSearchRoot():
        TObject(),
        fNode(""),
        fReadOnly(0),
        fVerbose(15),
        fCertificate("")
{
}

/// named cosntructor
/// \param node           - node description
/// \param readOnly       - not used yet
/// \param verbose        - verbosity
/// \param certification  - certification argument - will be prepended to  in curl queries
AliElasticSearchRoot::AliElasticSearchRoot(const std::string node, bool readOnly, Int_t verbose, std::string certification){
  //
  //
  fNode=node;
  fReadOnly=readOnly;
  fVerbose=verbose;
  fCertificate=certification;
    AliElasticSearchRoot::RegisterDefaultJQFilters();
}

/// print JQ filters selected by regular expression
void AliElasticSearchRoot::PrintJQFilters(TPRegexp reg) {
    typedef std::map<std::string,std::string>::const_iterator it_type;
    for(it_type iterator = fJQFilters.begin(); iterator != fJQFilters.end(); ++iterator) {
      if (reg.Match(iterator->first)==0) continue;
      std::cout << iterator->first << " " << iterator->second << "\n";
    }
    return;
}

/// Register default JQ filters - jq is like sed for JSON data
/// See jq filters:  https://stedolan.github.io/jq/
/// Example usage
/// \code
/// // Example usage of filters
/// index="alice_mc", type="passguess"
/// query=TString::Format("curl -s -XGET 'localhost:9200/%s/%s/_mapping' | %s",index,type,  TString::Format(AliElasticSearchRoot::fJQFilters["jqELayoutType"].data(),index,type,"text").Data());
/// result=gSystem->GetFromPipe(query);
/// \endcode
void AliElasticSearchRoot::RegisterDefaultJQFilters(){
    fJQFilters["jqELayoutJSON"] ="jq '.%s.mappings.%s.properties'"; // colored  filter to process elastic layout
    fJQFilters["jqELayoutType"] ="jq '.%s.mappings.%s.properties |. as $in| keys[]  | select($in[.].type==\"%s\")'"; // colored  filter to process elastic layout
    fJQFilters["jqELayoutTable"]="jq '.%s.mappings.%s.properties | keys[] as $k | \"\\\($k), \\\(.[$k] |.type)\"'";
    // Query filters
    fJQFilters["jqEQueryStat"]="jq '. | {hits_total:.hits.total, hits_max_score:.hits.max_score, _shards, timed_out, took}'";  // query stats  - extract stat
    fJQFilters["jqEQueryHeader"]="jq  -r '.hits|.hits|[.[]._source]|(map(keys) | add | unique) as $cols| $cols'";              // query layout - return array of variables
    fJQFilters["jqEQueryArrayNamed"]="jq  -r '.hits|.hits|[.[]._source][]'";            // query layout - return array of variables
    fJQFilters["jqEQueryArrayFlat"]="jq  -r '.hits|.hits|[.[]._source][][]'";           // query layout - return array of variables flat
    //
    fJQFilters["jqEQueryArrayIndexValue"]="jq  -r '.hits.hits|[.[]._source.%s]|(. | keys[]) as $i | $i, .[$i]"; // query return  element %s as index,value
    fJQFilters["jqEQueryArrayValue"]="jq  -r '.hits.hits|.[]._source.%s]"; // query return  element %s as value
}

void AliElasticSearchRoot::SetCertificate(std::string certification){
  fCertificate=certification;
}

///
/// \param option  - not used
void AliElasticSearchRoot::Print(Option_t* /*option*/) const{
  ::Info("AliElasticSearch::Print","Host:%s\tCertificate:%s\tReadOnly:%d",fNode.data(),fCertificate.data(),fReadOnly);
}

/// Bulk export to Elastic server
/// \param indexName    - index to export
/// \param jsonFile     - input bulk json file (e.g create by AliTreePlayer::select ...)
/// \return
Int_t AliElasticSearchRoot::ExportBinary(std::string indexName, std::string jsonFile, std::string /*option*/){
  //
  //
  //
  TString query=TString::Format("curl  -s  %s -XPOST '%s%s_bulk?pretty' --data-binary \"@%s\"",fCertificate.data(), fNode.data(),indexName.data(),jsonFile.data());
  if (fVerbose>0){
    ::Info("AliElasticSearchRoot:ExportBinary","%s%s %s",fNode.data(),indexName.data(),jsonFile.data());
  }
  ::Info("AliElasticSearchRoot:ExportBinary","%s",query.Data());
  return (gSystem->GetFromPipe(query.Data())).Length();
}

/// List indices in the Elastic server
void AliElasticSearchRoot::ls(Option_t *option) const {
  //
  // make indices list on server
  TString query=TString::Format("curl -s %s %s/_cat/indices?v | grep \"%s\"",fCertificate.data(), fNode.data(),option);
  if (fVerbose>0){
    ::Info("AliElasticSearchRoot::ls","%s",query.Data());
  }
  TString returnValue = gSystem->Exec(query.Data());
  ::Info("AliElasticSearchRoot::ls","%s",fNode.data());
  printf("%s\n",returnValue.Data());
}

/// GetIndexLayout projections for given data types
///  \param index     - elastic index name
///  \param type      - elastic index type
///  \param dataType  - data types to print - in case of NULL pointer print all
///  \param toScreen  - in case specified use nice printing to stdout (default only export to TString)
TString  AliElasticSearchRoot::GetIndexLayout(const char *index, const char *type, const char *dataType, Bool_t toScreen){
    TString query="";
    if (dataType){
        query=TString::Format("curl -s %s -XGET '%s/%s/%s/_mapping' | %s", fCertificate.data(), fNode.data(), index,type,  TString::Format(AliElasticSearchRoot::fJQFilters["jqELayoutType"].data(),index,type,dataType).Data());
    }else{
        query=TString::Format("curl -s %s -XGET '%s/%s/%s/_mapping' | %s", fCertificate.data(), fNode.data(), index,type,  TString::Format(AliElasticSearchRoot::fJQFilters["jqELayoutTable"].data(),index,type).Data());
    }
    if (fVerbose>0) {
        ::Info("AliElasticSearchRoot::GetIndexLayout()","%s",query.Data());
    }
    if (toScreen){
        gSystem->Exec(query);
    }
    return gSystem->GetFromPipe(query);
}
/// GetIndexLayeout - assuming flat structure (to be replace by jq implementation soon)
/// \param indexName  - old version
/// \return           - array of parameters
TObjArray *  AliElasticSearchRoot::GetIndexLayout(const char *indexName){
  //
  // simplified version of parser - working only for flat structure
  //
  TString query=TString::Format("curl -s -XGET '%s%s_mapping'",fNode.data(),indexName);
  if (fVerbose>0){
    ::Info("AliElasticSearchRoot::GetIndexLayout","%s",query.Data());
  }
  TString mapping  = gSystem->GetFromPipe(query.Data());
  //
  TObjArray * mappingArray = mapping.Tokenize("{},:\"");  
  Int_t nDesc=mappingArray->GetEntries();
  TObjArray *layoutArray = new TObjArray(nDesc/3);
  for (Int_t iEntry=0; iEntry<nDesc-1; iEntry++){
    if (TString(mappingArray->At(iEntry)->GetName()).EqualTo("type")){
      layoutArray->AddLast(new TNamed(mappingArray->At(iEntry-1)->GetName(), mappingArray->At(iEntry+1)->GetName()));
      iEntry+=1;
    }
  }
  return layoutArray;
}


/// select query  with similar interace to the TTree::Draw query
///   more options to be added soon -> Elastic json-> TTree,   GetVal())
/// \param indexName
/// \param fields          - fields to export  "*" select all fields
/// \param query           - elastic query in format as in TTree - internaly converted to painless script
/// \param first           - first entry to select
/// \param size            -
/// \param filterArray     - in verbose mode - provide list of jq filters
/// \return                - result of query is written to temporary query.json file
/// Example queries : see $AliRoot_SRC/STAT/Macros/AliElasticSearchRootTest.c:testSelect()

///\code
///   select=  pelastic->select("/alice/qatpc_test0/","*","run>240000&&run<246844&&meanMIP>50",0,1000,"jqEQueryStat,jqEQueryHeader");
///    select=  pelastic->select("/alice/qatpc_test0/","*","run>240000&&run<246844&&meanMIP>50",0,1000,"jqEQueryArrayNamed,jqEQueryStat,jqEQueryHeader");
///\endcode

 TString AliElasticSearchRoot::select(const char *indexName, TString fields, TString query, Int_t first, Int_t size, TString filterArray){
  //
  /*
    const char *indexName="/alice/logbook_test0/"
    TString fields="run:totalSubEvents";
    TString query="run>240000&&run<246844"
    verbose=7;
    first=0; size=1000;
  */
  
  AliElasticSearchRoot * el=this;                     //AliElasticSearchRoot * el=&elastic
  TObjArray * layout = el->GetIndexLayout(indexName);
  TString where=query;
  TObjArray *queryFields=where.Tokenize("><&|+=()! ");
  if (fVerbose&0x2) queryFields->Print();
  //
  {
    Int_t currentIndex=-1;
    for (Int_t ifield=0; ifield<queryFields->GetEntriesFast(); ifield++){
      currentIndex= where.Index(queryFields->At(ifield)->GetName(),currentIndex);
      //printf("%d\t%s\n",currentIndex, queryFields->At(ifield)->GetName());
      TObjString * field=(TObjString*)layout->FindObject(queryFields->At(ifield)->GetName());
      if (field!=NULL){
        if (fVerbose&0x1){
            printf("%d\t%s\t%s\n",ifield,field->GetName(), where.Data());
        }
        TString insertString=TString::Format("doc[\\u0027%s\\u0027].value",field->GetName()).Data();
        where.Replace(currentIndex,field->String().Length(),"");
        where.Insert(currentIndex,insertString);
      }
    } 
  }
  //
  where=TString::Format("\"query\":\n{\"bool\":{\n\"must\":{\n\"script\":{\n\"script\":{\n\"inline\": \"%s\",\n\"lang\":\"painless\"\n}\n}\n}\n}\n}",where.Data());
  //where=TString::Format("\"query\":\n{\"bool\":{\n\"must\":{\n\"script\":{\n\"script\":{\n\"inline\": \"%s\",\n\"lang\":\"groovy\"\n}\n}\n}\n}\n}",where.Data());
  if (fVerbose&0x4) printf("%s\n",where.Data());
  //
  TObjArray *fieldsArray = fields.Tokenize(":|");
  TString what =TString::Format("\"_source\" :[");
  for (Int_t iField=0; iField<fieldsArray->GetEntriesFast(); iField++){
    what+="\"";
    what+=fieldsArray->At(iField)->GetName();
    if (iField<fieldsArray->GetEntriesFast()-1)  {
      what+="\",";
    }else{
      what+="\"]";
    }
  }
  TString elQuery="";
  {
    elQuery=TString::Format("curl -s -XPOST \"%s%s_search?pretty\" -d\'\n{",el->fNode.data(), indexName);
    elQuery+=TString::Format("\"from\": %d,\n",first);
    elQuery+=TString::Format("\"size\": %d,\n",size);
    elQuery+=what;
    elQuery+=",\n";
    elQuery+=where;
    elQuery+="\n}\n'>query.json\n";
  }
  if ((fVerbose&0x1)>0)  ::Info("selectFromWhere\n","%s",elQuery.Data());
  TString result= gSystem->GetFromPipe(elQuery.Data());
    if (filterArray.Length()>0) {
        TObjArray * fArray = filterArray.Tokenize(",");
        Int_t fentries = fArray->GetEntries();
        for (Int_t ientry=0; ientry<fentries; ientry++){
            std::string filter=fJQFilters[fArray->At(ientry)->GetName()];
            if (filter.size()>0) {
                gSystem->Exec(TString::Format("cat query.json | %s", filter.data()).Data());
            }else{

            }
        }
    }
    if (1/*isFlat*/) {
        //TString result = gSystem->GetFromPipe(TString::Format("cat query.json | %s", filter.data()).Data());
    }
  return result;
  
}


/*
void jsontocsv(){
  //
  https://csvkit.readthedocs.io/en/0.9.1/install.html
  sudo apt-get install jq
}
*/
