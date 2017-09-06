#ifndef ALIELASTICSEARCHROOT_H
#define ALIELASTICSEARCHROOT_H

/// \ingroup STAT
/// \class AliElasticSearchRoot
/*!
 \brief AliElasticSearchRoot - C++ interface to the Elastic --> root tree like-   using gSystem->GetFromPipe, Exec, Elastic, painless scripting,  curl,  jq
      - Originally planned usage of  higher level C++ interface not possible for a moment
      -- C++11 or better C++14 and ROOT6 needed (e.g https://github.com/QHedgeTech/cpp-elasticsearch/tree/master/src )
      - Using older  C++  and Root 5 to represent hierarchical structure root trees to be used
      - Query language similar to TTree::Draw  - using elastic painless scripting language (see query examples)
      - jq -  jq is like sed for JSON data - you can use it to slice and filter and map and transform structured data
      --  see:       https://stedolan.github.io/jq/
      --- jq player  to check/test jq syntax https://jqplay.org/

### Example usage  (assuming access to Elastic server with alice indexes)

  - see also $AliRoot_SRC/STAT/Macros/AliElasticSearchRootTest.c
  - to create indices needed $AliRoot_SRC/STAT/Macros/AliElasticSearchRootTest.c:testExportBinary() can be used

\code
   hostName="localhost:9200";            // local elastic
   // hostName="https://188.184.85.222/krb"
   AliElasticSearchRoot elastic(hostName,kFALSE,1);
\endcode

#### Set certification

\code
   elastic.SetCertificate();           // default certification="--insecure --cacert CERNRootCertificationAuthority2.crt --netrc-file $HOME/.globus/.elastic/.netrc"
   elastic.SetCertificate("-u elastic:changeme"); // using username(elastic)":pwd(changeme) as in Elastic tutorials
   // to add Kerberos and negotiate
   elastic.SetCertificate("-k -u : --negotiate")

\endcode

#### Browse content

\code
   elastic.ls();                       // list all indices
   elastic.ls("alice");                // list indices containing alice
   elastic.ls(".alice_mc");             // list indices containing alice_mc
   elastic.Print();
   // GetIndexLayout
   elastic.GetIndexLayout("/alice/logbook_test0/")->Print();
   elastic.GetIndexLayout("/alice_mc/passguess/")->Print();
\endcode

#### TTree like queries

\code
   TString select="";
   select=  elastic.select("/alice/qatpc_test0/","run:meanMIP:meanTPCncl","run>240000&&run<246844&&meanMIP>50",0,1000,"jqEQueryStat,jqEQueryHeader")
   select=  elastic.select("/alice/logbook_test0/","run:totalSubEvents:ctpDuration","run>140000",0,10000,"jqEQueryStat,jqEQueryHeader")
   select=  elastic.select("/alice_mc/passguess/","*","1",0,1000,"jqEQueryStat,jqEQueryHeader,jqEQueryArrayNamed");  // get list of all runs output in query.json
   select=  elastic.select("/alice/qatpc_test0/","*","run>240000&&run<246844&&meanMIP>50",0,1000,0)
\endcode

  ### To do

  \author marian  Ivanov marian.ivanov@cern.ch
*/

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
