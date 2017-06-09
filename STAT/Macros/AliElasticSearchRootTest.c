///  ### test of AliEleasticSearchRoot - C++ interace to the Elastic using gSystem->GetFromPipe,Exec, curl and jq
///      Origignally planned usage of  higher level C++ inteface not possible for a moment
///           - C++11 or better C++14 and ROOT6 needed (e.g https://github.com/QHedgeTech/cpp-elasticsearch/tree/master/src )
///  \author marian  Ivanov marian.ivanov@cern.ch

///          In order to run the test Elastic server has to be running = default location hostName="localhost:9200"  - but it can pbe specified as test argument
/////  to run tests one by one
/////  \code
//      .L $AliRoot_SRC/STAT/Macros/AliElasticSearchRootTest.c
//    InitElastic();
//    testExportBinary();
//    testGetIndexLayout();
//    testSelect();
///// \endcode
//
///// To run full test
/////  \code
// aliroot -b -q $AliRoot_SRC/STAT/AliElasticSearchRoot.cxx+  $AliRoot_SRC/STAT/Macros/AliElasticSearchRootTest.c >AliElasticSearchRootTest.log
//    Parse output
//    cat AliElasticSearchRootTest.log | grep "I-AliElasticSearchRoot.*Key"
//I-AliElasticSearchRoot::AliElasticSearchRootTest: KeyValue-Begin
//I-AliElasticSearchRoot::ExportBinary: KeyValue-Begin
//I-AliElasticSearchRoot::ExportBinary: KeyValue-End
//I-AliElasticSearchRoot::GetIndexLayout: KeyValue-Begin
//I-AliElasticSearchRoot::GetIndexLayout: KeyValue-End
//I-AliElasticSearchRoot::select: KeyValue-Begin
//I-AliElasticSearchRoot::select: KeyValue-End
//I-AliElasticSearchRoot::AliElasticSearchRootTest: KeyValue-End
///// \endcode




char *hostName="localhost:9200";
AliElasticSearchRoot *pelastic=NULL;

void InitElastic();
void testExportBinary();
void testGetIndexLayout();
void testSelect();


void AliElasticSearchRootTest(char *phostName=0){
    if (phostName) hostName=phostName;
    ::Info("AliElasticSearchRoot::AliElasticSearchRootTest","KeyValue-Begin");
    InitElastic();
    testExportBinary();
    testGetIndexLayout();
    testSelect();
    ::Info("AliElasticSearchRoot::AliElasticSearchRootTest","KeyValue-EndOK");
}

void InitElastic(){
    pelastic=new AliElasticSearchRoot(hostName,kFALSE,1);
    pelastic->SetCertificate();           // default certification="--insecure --cacert CERNRootCertificationAuthority2.crt --netrc-file $HOME/.globus/.elastic/.netrc"
    pelastic->ls();
    pelastic->Print();
}

///  testExportBinary() AliElasticSearchRoot::ExportBinary
///  Get production tree and export al branches to the elastic key /alice_mc/passguess/
void testExportBinary(){
  // Example export MC production table
  ::Info("AliElasticSearchRoot::ExportBinary","KeyValue-Begin");
  AliExternalInfo info;
  TTree * tree = info.GetTreeMCPassGuess();
  TString exportInfo="%IprodName:";
  Int_t nBranches = tree->GetListOfBranches()->GetEntries();
  for (Int_t iBranch=0; iBranch<nBranches; iBranch++){
    exportInfo+= tree->GetListOfBranches()->At(iBranch)->GetName();
    if (iBranch<nBranches) exportInfo+=":";
  }
  Int_t entries = AliTreePlayer::selectWhatWhereOrderBy(tree,exportInfo.Data(),"1", "", 0,100000, "elastic","mctableAnchor.json");
  pelastic->ExportBinary("/alice_mc/passguess/","mctableAnchor.json","xxx");
  ::Info("AliElasticSearchRoot::ExportBinary","KeyValue-EndOK");
}


void testGetIndexLayout(){
    //
    // elastic.GetInde xLayout("/alice/logbook_test0/")->Print();
    ::Info("AliElasticSearchRoot::GetIndexLayout","KeyValue-Begin");
    pelastic->GetIndexLayout("alice_mc","passguess","text").Tokenize("\n")->Print();
    pelastic->GetIndexLayout("alice_mc","passguess",0,kTRUE).Tokenize("\n")->Print();
    pelastic->GetIndexLayout("alice","qatpc_test0","float").Tokenize("\n")->Print();
    pelastic->GetIndexLayout("alice","qatpc_test0",0,kTRUE).Tokenize("\n")->Print();
    ::Info("AliElasticSearchRoot::GetIndexLayout","KeyValue-EndOK");
}

void testSelect(){
    ::Info("AliElasticSearchRoot::select","KeyValue-Begin");
    TString select="";
    select=  pelastic->select("/alice/qatpc_test0/","*","run>240000&&run<246844&&meanMIP>50",0,1000,"jqEQueryStat,jqEQueryHeader");
    select=  pelastic->select("/alice/qatpc_test0/","*","run>240000&&run<246844&&meanMIP>50",0,1000,"jqEQueryArrayNamed,jqEQueryStat,jqEQueryHeader");
    ::Info("AliElasticSearchRoot::select","KeyValue-EndOK");
}
