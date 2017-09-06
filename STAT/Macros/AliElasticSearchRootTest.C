///  ### test of AliEleasticSearchRoot - C++ interace to the Elastic using gSystem->GetFromPipe,Exec, curl and jq
///      Origignally planned usage of  higher level C++ inteface not possible for a moment
///           - C++11 or better C++14 and ROOT6 needed (e.g https://github.com/QHedgeTech/cpp-elasticsearch/tree/master/src )
///  \author marian  Ivanov marian.ivanov@cern.ch

/*!
    In order to run the test Elastic server has to be running = default location hostName="localhost:9200"  - but it can be specified as test argument
    to run tests one by one
    \code
      .L $AliRoot_SRC/STAT/Macros/AliElasticSearchRootTest.C
       InitElastic();
       testGetIndexLayout();
       testSelect();
    \endcode
     To run full test
    \code
    aliroot -b -q $AliRoot_SRC/STAT/Macros/AliElasticSearchRootTest.C |tee AliElasticSearchRootTest.log
    Parse output
        cat AliElasticSearchRootTest.log | grep "I-AliElasticSearchRoot.*Key"
        I-AliElasticSearchRoot::AliElasticSearchRootTest: KeyValue-Begin
        I-AliElasticSearchRoot::ExportBinary: KeyValue-Begin
        I-AliElasticSearchRoot::ExportBinary: KeyValue-End
        I-AliElasticSearchRoot::GetIndexLayout: KeyValue-Begin
        I-AliElasticSearchRoot::GetIndexLayout: KeyValue-End
        I-AliElasticSearchRoot::select: KeyValue-Begin
        I-AliElasticSearchRoot::select: KeyValue-End
        I-AliElasticSearchRoot::AliElasticSearchRootTest: KeyValue-End
    \endcode
*/



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
    testExportClass();
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
///  - 1.) Get production tree and export al branches to the elastic key /alice_mc/passguess/
///  - 2.) Get QA/TPC and export all branches to the elastic key /alice/tpc_bulk
///  - 3.) Get QA.TPC and export selected branches to the elastic key /alice/tpx_test0
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
    //
    tree = info.GetTree("QA.TPC","LHC15n","pass2","Logbook");
    //TObjArray * branches = AliTreePlayer::selectMetadata(tree, "[class==\"dEdx||DCA\"]",0);
    // 2.) Example all branches of the TPC QA
    exportInfo="%Irun:";
    nBranches = tree->GetListOfBranches()->GetEntries();
    for (Int_t iBranch=0; iBranch<nBranches; iBranch++){
        exportInfo+= tree->GetListOfBranches()->At(iBranch)->GetName();
        exportInfo+=":";
    }

    entries = AliTreePlayer::selectWhatWhereOrderBy(tree,exportInfo.Data(),"1", "", 0,100000, "elastic","qatpcbulk.json");
    pelastic->ExportBinary("/alice/qatpcbulk/","qatpcbulk.json","xxx");
    // 3. export selected branches of the TPC QA
    TString exportInfo="%ILogbook.LHCperiod:%Ipass.GetName():%Irun:%Prun:";
    exportInfo+="period.GetName():pass.GetName():run:Logbook.run:";
    exportInfo+="meanTPCncl:meanMIP;1.5:meanMIPele;1.5:";
    exportInfo+="grdcar_neg_ASidePhi.fY:grdcar_neg_CSidePhi.fY:grdcar_pos_ASidePhi.fY:grdcar_pos_CSidePhi.fY";

    entries = AliTreePlayer::selectWhatWhereOrderBy(tree,exportInfo.Data(),"1", "", 0,100000, "elastic","qatpc_test0.json");
    pelastic->ExportBinary("/alice/qatpc_test0/","qatpc_test0.json","xxx");

}

void testExportClass() {
    ::Info("AliElasticSearchRoot::testExportClass","KeyValue-Begin");
    AliExternalInfo info;
    TTree *tree = info.GetTree("QA.TPC","LHC15n","pass2","Logbook");
    TString exportInfo="%Irun:";
    exportInfo+="period.GetName():";
    exportInfo+="period.:";
    exportInfo+="run:";
    exportInfo+="grdcar_neg_ASidePhi.";
    Int_t entries = AliTreePlayer::selectWhatWhereOrderBy(tree,exportInfo.Data(),"1", "", 0,100000, "elastic","qatpc_testClass.json");
    pelastic->ExportBinary("/alice/qatpc_testClass1/","qatpc_testClass.json","xxx");
    ::Info("AliElasticSearchRoot::testExportClass","KeyValue-End");
}

void testExportClass() {
    ::Info("AliElasticSearchRoot::testExportClass","KeyValue-Begin");
    AliExternalInfo info;
    TTree *tree = info.GetTree("QA.TPC","LHC15n","pass2","Logbook:QA.ITS:QA.TRD");
    //

}

void testGetIndexLayout(){
    //
    // elastic.GetInde xLayout("/alice/logbook_test0/")->Print();
    ::Info("AliElasticSearchRoot::GetIndexLayout","KeyValue-Begin");
    pelastic->GetIndexLayout("alice_mc","passguess","text").Tokenize("\n")->Print();
    pelastic->GetIndexLayout("alice_mc","passguess",0,kTRUE).Tokenize("\n")->Print();
    pelastic->GetIndexLayout("alice","qatpc_test0","float").Tokenize("\n")->Print();
    pelastic->GetIndexLayout("alice","qatpc_test0",0,kTRUE).Tokenize("\n")->Print();
    //
    ::Info("AliElasticSearchRoot::GetIndexLayout","KeyValue-EndOK");
}

void testSelect(){
    ::Info("AliElasticSearchRoot::select","KeyValue-Begin");
    TString select="";
    select=  pelastic->select("/alice/qatpc_test0/","*","run>240000&&run<246844&&meanMIP>50",0,1000,"jqEQueryStat,jqEQueryHeader");
    select=  pelastic->select("/alice/qatpc_test0/","*","run>240000&&run<246844&&meanMIP>50",0,1000,"jqEQueryArrayNamed,jqEQueryStat,jqEQueryHeader");

    select=  pelastic->select("/alice/qatpcbulk/","*","run>0&&meanMIP>0",0,1000,",jqEQueryHeader,jqEQueryStat");

    ::Info("AliElasticSearchRoot::select","KeyValue-EndOK");
}
