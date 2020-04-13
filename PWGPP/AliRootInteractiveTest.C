/*
 .L $AliPhysics_SRC/PWGPP/AliRootInteractive.h
 .x $AliPhysics_SRC/PWGPP/AliRootInteractiveTest.C

 */
TTree * tree=0;
TMatrixD *testMatrix=0;
TStopwatch timer;



TTree *  makeABCtree(Int_t nPoints){
    TTreeSRedirector *pcstream = new TTreeSRedirector("treeABCD.root","recreate");
    Double_t abcd[4];
    for (Int_t i=0; i<nPoints; i++){
        for (Int_t j=0; j<4; j++) abcd[j]=gRandom->Rndm();
        (*pcstream)<<"tree"<<
                   "A="<<abcd[0]<<
                   "B="<<abcd[1]<<
                   "C="<<abcd[2]<<
                   "D="<<abcd[3]<<
                   "\n";
    }
    delete pcstream;
    TFile *f = TFile::Open("treeABCD.root");
    return (TTree*)f->Get("tree");
}

void testBokehDrawArray(){
    tree= makeABCtree(1000);
    TString query = "A>0";
    TString figureArray= "["
                         "[['A'], ['D+A','C-A'], {\"size\": 1 , 'colorZvar':'D'}],"
                         "[['A'], ['C+A', 'C-A']],"
                         "[['A'], ['C']], ['table'] ]";

    TString widgets="query.xx(),slider.A(0,1,0.01,0.1,0.9),slider.B(0,1,0.01,0.1,0.9),"
                    "slider.C(0,1,0.01,0.1,0.9),slider.D(0,1,0.01,0.1,0.9)";
    TString options = "tooltips=[('VarA', '(@A)'), ('VarB', '(@B)'), ('VarC', '(@C)'), ('VarD', '(@D)')],"
                      "layout= '((0,1,2 ),(3, x_visible=1),commonX=1, x_visible=1,y_visible=0,plot_height=450,plot_width=1000)'";

    AliRootInteractive::treeBokehDrawArray("tree", query, figureArray, widgets, options);
}

void AliRootInteractiveTest(){
  testBokehDrawArray();
}
