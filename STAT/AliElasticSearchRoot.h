#ifndef ALIELASTICSEARCHROOT_H
#define ALIELASTICSEARCHROOT__H

/// \ingroup STAT
/// \class AliElasticSearchRoot
/// \brief AliElasticSearchRoot - C++ interace to the Elastic --> root tree like-   using gSystem->GetFromPipe, Exec, Elastic, painlsees scripting,  curl,  jq
///      - Originally planned usage of  higher level C++ inteface not possible for a moment
///      -- C++11 or better C++14 and ROOT6 needed (e.g https://github.com/QHedgeTech/cpp-elasticsearch/tree/master/src )
///      - Using older  C++  and Root 5 to represent heirarchiacal structure root trees to be used
///      - Query language similar to TTree::Draw  - suing elastic painless scriptin langaug (see query expamples)
/// ### Example usage (expecting acces to Elastic server with alice indexes)
///  - see also $AliRoot_SRC/STAT/Macros/AliElasticSearchRootTest.c
///  - to create indeces needed $AliRoot_SRC/STAT/Macros/AliElasticSearchRootTest.c:testExportBinary() can be used
/// \code
///  hostName="localhost:9200";            // local elastic
///   // or remote ElastichostName="https://localhost:9206";    //  connect to remote via tunnel e.g ssh mivanov@lxplus.cern.ch -f -N -L 9206:es-alice.cern.ch:9203
///  //
///  AliElasticSearchRoot elastic(hostName,kFALSE,1);
///  elastic.SetCertificate();           // default certification="--insecure --cacert CERNRootCertificationAuthority2.crt --netrc-file $HOME/.globus/.elastic/.netrc"
///  elastic.ls();                       // list all indeces
///  elastic.ls("alice");                // list incdeces containing alice
///  elastic.ls(".alice_mc");             // list indeces containing alice_mc
///  elastic.Print();
///  // GetIndexLayeout
///  elastic.GetIndexLayout("/alice/logbook_test0/")->Print();
///  elastic.GetIndexLayout("/alice_mc/passguess/")->Print();
///  // make a Tree like query
///  TString select="";
///  select=  elastic.select("/alice/qatpc_test0/","run:meanMIP:meanTPCncl","run>240000&&run<246844&&meanMIP>50",0,1000,"jqEQueryStat,jqEQueryHeader")
///  select=  elastic.select("/alice/logbook_test0/","run:totalSubEvents:ctpDuration","run>140000",0,10000,"jqEQueryStat,jqEQueryHeader")
///  select=  elastic.select("/alice_mc/passguess/","*","1",0,1000,"jqEQueryStat,jqEQueryHeader,jqEQueryArrayNamed");  // get list of all runs output in query.json
///  //
///  select=  elastic.select("/alice/qatpc_test0/","*","run>240000&&run<246844&&meanMIP>50",0,1000,0)
///
/// To do  - json- formatted cvs conversion to enable root like queries
/// \endcode
/// \author marian  Ivanov marian.ivanov@cern.ch

#include "TObject.h"
#include <map>
#include <string>
#include "TPRegexp.h"


class AliElasticSearchRoot : public TObject{

public:
    AliElasticSearchRoot();
    AliElasticSearchRoot(const std::string node, bool readOnly, Int_t fVerbose, std::string certification="");
    void SetCertificate(std::string certification="--cacert $HOME/.globus/CERNGridCertificationAuthority.crt --netrc-file $HOME/.globus/.elastic/.netrc");
    void Print(Option_t* option = "") const;
    virtual void ls(Option_t *option="") const;
    // queries
    Int_t        ExportBinary(std::string indexName, std::string jsonFile, std::string option);
    TObjArray *  GetIndexLayout(const char *indexName);
    TString      GetIndexLayout(const char *index, const char *type, const char *dataType, Bool_t display=kFALSE);
    TString select(const char *indexName, TString fields, TString query, Int_t first, Int_t size, TString filterArray="");
    // JQ filters
    static void RegisterDefaultJQFilters();
    static void PrintJQFilters(TPRegexp reg);
public:
    std::string fNode;
    Bool_t  fReadOnly;
    Int_t   fVerbose;
    std::string  fCertificate;  // certificates description - not tested yet
    static std::map<std::string,std::string> fJQFilters;
private:
    ClassDef(AliElasticSearchRoot,1);
};

#endif