/*
 Macro to make all multicanvases without space between them
 
 use the functions 
 -Convert3x3Matrix(TCanvas* c,TCanvas* &c_output, Bool_t transpose)
 -Convert3x2Matrix(TCanvas* c,TCanvas* &c_output, Bool_t transpose)
 give in input a canvas that is divided in 9 (6) sub canvases, and the functions Convert will remove the empty spaces
 
 */

//set limits to the margin of the total canvas (this numbers are already optimised, change them as you like, the macro will compute all the margins automatically)
/*Double_t LimitLeft =0.25; //limit from left margin (15% of size)
Double_t LimitBottom =0.15; //limit from bottom margin (15% of size)
Double_t LimitRight =0.05; //limit from right margin (5% of size)
Double_t LimitTop =0.05; //limit from top margin (5% of size)*/
TString inputdirectory = "";
Color_t ppcolor=kBlack;
Color_t pPbcolor=kRed;
Int_t markerstyle[2]={20,21};
//TString inputdirectory = "~/Analysis/HFCorrelations/testScript/ReflectedPlots/StdRebin/AllPlots/NiceStylePlots/Output_Plots/WeightedAverageDzeroDstarDplus";
//________________________________________________
void SetInputDirectory(TString strdir){
   inputdirectory=strdir;
}

//________________________________________________
void DoPlotInSingleCanvasNoSpaces(){
    
    //SetInputDirectory("~/Analysis/HFCorrelations/testScript");
    
    PlotppAverages("CanvasOfWeightedAverageDzeroDstarDplus_pp.root",kTRUE);
    PlotpPbAverages("CanvasOfWeightedAverageDzeroDstarDplus_pPb.root",kTRUE);
}


//________________________________________________
void PlotppAverages(TString filename, Bool_t transpose){
    
     inputdirectory.Form("Output_Plots/WeightedAverageDzeroDstarDplus");
    
    TFile * file = TFile::Open(Form("%s/%s",inputdirectory.Data(),filename.Data()));
    if(!file->IsOpen()){
        cout << "cannot open file" << endl;
        return;
    }
    
    TCanvas * c = (TCanvas*)file->Get("coutput");
    
     TCanvas * c_ouput = new TCanvas("c_ouput","c_ouput",0,0,1030,1000);
     Convert3x3Matrix(c,c_ouput,transpose);
    
    
     SaveCanvas(c_ouput,inputdirectory,"CanvasNoSpaces_WeightedAverageDzeroDstarDplus_pp");
}
//________________________________________________
void PlotpPbAverages(TString filename, Bool_t transpose){
    
    
      inputdirectory.Form("Output_Plots/WeightedAverageDzeroDstarDplus");
    
    TFile * file = TFile::Open(Form("%s/%s",inputdirectory.Data(),filename.Data()));
    if(!file->IsOpen()){
        cout << "cannot open file" << endl;
        return;
    }
    
    TCanvas * c = (TCanvas*)file->Get("coutput");
    
    TCanvas * c_ouput = new TCanvas("c_ouput","c_ouput",0,0,800,1000);
    
    Convert3x2Matrix(c,c_ouput,transpose);
    
    SaveCanvas(c_ouput,inputdirectory,"CanvasNoSpaces_WeightedAverageDzeroDstarDplus_pPb");
    
}
//________________________________________________
void PlotComparisonspp_pPb(){
    
    /* to be finished
    my idea was the following,
     
     1) define a dummy canvas, that you divide in 2,3
     2) in each pad of the canvas, draw the single canvas from the files plotComparison_WeightedAverage_pp_pPb_58_0.3to99.0.root
     3) pass the dummy canvas to the Convert3x2Matrix that will remove all the spaces properly
    */
    
    TCanvas * cdummy = new TCanvas("coutput","coutput",0,0,1200,800); // keep the name coutput
    cdummy->Divide(2,3);
    cdummy->cd(1);

    // problem is with this function, see description below
    ImportAndDraw("~/Analysis/HFCorrelations/testScript/ReflectedPlots/StdRebin/AllPlots/Averages/ComparisonPPtoPPB/plotComparison_WeightedAverage_pp_pPb_58_0.3to99.0.root","ptmidsub",0);
    
    /*
     do the same for the other figures
     */
    
    
    // pass the canvas cdummy to the function Convert3x2Matrix, that will do the same job as before
    TCanvas * c_ouput = new TCanvas("c_ouput","c_ouput",0,0,800,1200);
    Convert3x2Matrix(cdummy,c_ouput,kFALSE);
    
    SaveCanvas(c_ouput,inputdirectory,"CanvasNoSpaces_ppvspPb");
}

//________________________________________________
void ImportAndDraw(TString filename, TString canvasname, Int_t k){
    
      TFile * file = TFile::Open(filename.Data());
    if(!file->IsOpen()){
        cout << "cannot open file " << filename << endl;
    }
    
    TCanvas * c = (TCanvas *)file->Get(canvasname.Data());
    if(!c){
        cout << "cannot get canvas" << canvasname << endl;
        return;
    }
    
    // get list of primitives
    TList * list = c->GetListOfPrimitives();
    list->ls();
    
    TH1D * hpp = (TH1D*)list->At(4); // acces histogram

   
    hpp->Draw(); /* >>>>>>>>>     I GOT STUCK HERE - for some reason it is not drawing it!   <<<<<<<<<<<<<<<<<<<<<<<<<<*/
    cout<<"drawn" << endl;
  
    
    
}



//________________________________________________
void Convert3x3Matrix(TCanvas* c, TCanvas* &c_output,Bool_t transpose){
    
    if(!c){
        cout << "Cannot read canvas!" << endl; return;
    }

  TPad * pd1 = (TPad*)c->GetPad(1);
  TPad * pd2 =(TPad*) c->GetPad(4);
  TPad * pd3 =(TPad*) c->GetPad(7);
  TPad * pd4 =(TPad*) c->GetPad(2);
  TPad * pd5 =(TPad*) c->GetPad(5);
  TPad * pd6 =(TPad*) c->GetPad(8);
  TPad * pd7 =(TPad*) c->GetPad(3);
  TPad * pd8 =(TPad*) c->GetPad(6);
  TPad * pd9 =(TPad*) c->GetPad(9);

  Double_t xl,xu,yl,yu;
  Double_t marginLeft=0.09;
  Double_t marginRight=0.02;
  Double_t marginTop=0.02;
  Double_t marginBottom=0.06;
  Double_t marginLeftForXAxis=0.02;
  Double_t marginBottomForYAxis=0.02;
  innerPadWidth=(1-marginLeft-marginRight)/3.;// this is the width w/o margin, not the real pad width!!
  innerPadHeight=(1-marginTop-marginBottom)/3.;// this is the height w/o margin, not the real pad height, which differs between inner pads and pads at the "boarders"!!


    // Bottom row
    pd7->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd7->SetPad(0.,0.,innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd7->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd7->SetRightMargin(0.);
    pd7->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd7->SetTopMargin(0.);

    pd7->Modified();
    pd7->Update();

    pd8->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd8->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,0,2.*innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd8->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd8->SetRightMargin(0.);
    pd8->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd8->SetTopMargin(0.);

    pd8->Modified();
    pd8->Update();
    
    pd9->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd9->SetPad(2*innerPadWidth+marginLeft-marginLeftForXAxis,0,1.,innerPadHeight+marginBottom);
    pd9->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd9->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd9->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd9->SetTopMargin(0.);

    pd9->Modified();
    pd9->Update();

// Middle row
    pd4->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd4->SetPad(0.,innerPadHeight+marginBottom-marginBottomForYAxis,innerPadWidth+marginLeft,2*innerPadHeight+marginBottom);
    pd4->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd4->SetRightMargin(0.);
    pd4->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd4->SetTopMargin(0.);

    pd4->Modified();
    pd4->Update();

    pd5->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd5->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,2.*innerPadWidth+marginLeft,2*innerPadHeight+marginBottom);
    pd5->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd5->SetRightMargin(0.);
    pd5->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd5->SetTopMargin(0.);

    pd5->Modified();
    pd5->Update();
    
    pd6->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd6->SetPad(2*innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,1.,2.*innerPadHeight+marginBottom);
    pd6->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd6->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd6->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd6->SetTopMargin(0.);

    pd6->Modified();
    pd6->Update();

    // Top Row
    pd1->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd1->SetPad(0,2.*innerPadHeight+marginBottom-marginBottomForYAxis,innerPadWidth+marginLeft,1.);
    pd1->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd1->SetRightMargin(0.);
    pd1->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd1->SetTopMargin(marginTop/(innerPadHeight+marginTop));

    pd1->Modified();
    pd1->Update();

    pd2->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd2->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,2.*innerPadHeight+marginBottom-marginBottomForYAxis,2.*innerPadWidth+marginLeft,1);
    pd2->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    //    pd2->SetLeftMargin(0.);
    pd2->SetRightMargin(0);
    pd2->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd2->SetTopMargin(marginTop/(innerPadHeight+marginTop));

    pd3->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd3->SetPad(2.*innerPadWidth+marginLeft-marginLeftForXAxis,2.*innerPadHeight+marginBottom-marginBottomForYAxis,1,1);
    //    pd3->SetLeftMargin(0.);
    pd3->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd3->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd3->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd3->SetTopMargin(marginTop/(innerPadHeight+marginTop));

    scaleHeightPads=pd1->GetHNDC();
    scaleWidthPads=pd1->GetWNDC();
   
    c_ouput->cd();
    pd1->Draw();
    pd2->Draw();
    pd3->Draw();
    pd4->Draw();
    pd5->Draw();
    pd6->Draw();
    pd7->Draw();
    pd8->Draw();
    pd9->Draw();

    
}

//________________________________________________
void Convert3x2Matrix(TCanvas* c,TCanvas* &c_output, Bool_t transpose){
    
    if(!c){
        cout << "Cannot read canvas!" << endl; return;
    }
    
    //default is 3x2, transposed is 2x3
    
    Double_t horizontalwidth =0;
    Double_t verticalwidth = 0;

//    TList * list = c->GetListOfPrimitives();
    // list->ls();
    
    TPad * pd1 = (TPad*)c->GetPad(1);
    TPad * pd2 =(TPad*) c->GetPad(4);
    TPad * pd3 =(TPad*) c->GetPad(2);
    TPad * pd4 =(TPad*) c->GetPad(5);
    TPad * pd5 =(TPad*) c->GetPad(3);
    TPad * pd6 =(TPad*) c->GetPad(6);
    
   // SetPadStyle(pd1);

    Double_t xl,xu,yl,yu;
    Double_t marginLeft=0.15;
    Double_t marginRight=0.03;
    Double_t marginTop=0.02;
    Double_t marginBottom=0.05;
    innerPadWidth=(1-marginLeft-marginRight)/2.;// this is the width w/o margin, not the real pad width!!
    innerPadHeight=(1-marginTop-marginBottom)/3.;// this is the height w/o margin, not the real pad height, which differs between inner pads and pads at the "boarders"!!
    Printf("innerPadHeight: %f",innerPadHeight);
    Printf("innerPadWidth: %f",innerPadWidth);
    Double_t marginLeftForXAxis=0.02;
    Double_t marginBottomForYAxis=0.02;

    // Bottom row
    pd5->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd5->SetPad(0.,0.,innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd5->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd5->SetRightMargin(0.);
    pd5->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd5->SetTopMargin(0.);

        pd5->Modified();
        pd5->Update();

    pd6->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd6->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,0,1.,innerPadHeight+marginBottom);
    pd6->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));//0.02/(1.-innerPadWidth-marginLeft));
    pd6->SetRightMargin(marginRight/(innerPadWidth+marginRight+marginLeftForXAxis));
    pd6->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd6->SetTopMargin(0.);

        pd6->Modified();
        pd6->Update();


    // Middle Row
    pd3->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd3->SetPad(0.,innerPadHeight+marginBottom-marginBottomForYAxis,innerPadWidth+marginLeft,2.*innerPadHeight+marginBottom);
    pd3->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd3->SetRightMargin(0.);
    pd3->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd3->SetTopMargin(0.);

        pd3->Modified();
        pd3->Update();

    pd4->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd4->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,1.,2.*innerPadHeight+marginBottom);
    pd4->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));//0.02/(1.-innerPadWidth-marginLeft));
    pd4->SetRightMargin(marginRight/(innerPadWidth+marginRight+marginLeftForXAxis));
    pd4->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd4->SetTopMargin(0.);

        pd4->Modified();
        pd4->Update();


    // Top Row
    pd1->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd1->SetPad(0,2.*innerPadHeight+marginBottom-marginBottomForYAxis,innerPadWidth+marginLeft,1.);
    pd1->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd1->SetRightMargin(0.);
    pd1->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd1->SetTopMargin(marginTop/(innerPadHeight+marginTop+marginBottomForYAxis));

    pd1->Modified();
    pd1->Update();
    
   
    scaleHeightPads=pd1->GetHNDC();
    scaleWidthPads=pd1->GetWNDC();

    pd2->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd2->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,2.*innerPadHeight+marginBottom-marginBottomForYAxis,1,1);
    pd2->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));//0.02/(1.-innerPadWidth-marginLeft));
    pd2->SetRightMargin(marginRight/(innerPadWidth+marginRight+marginLeftForXAxis));
    pd2->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis+marginTop));
    pd2->SetTopMargin(marginTop/(innerPadHeight+marginTop+marginBottomForYAxis));
   
    c_ouput->cd();
    pd1->Draw();
    pd2->Draw();
    pd3->Draw();
    pd4->Draw();
    pd5->Draw();
    pd6->Draw();
  

}



//________________________________________________*********
void SaveCanvas(TCanvas * c, TString directory, TString name){
    
    TString outputDir = "";//Plots/15_May/";
    //    gSystem->Exec(Form("mkdir %s",outputDir.Data()));
    
    if(directory != ""){outputDir += directory;
        TString exec = "mkdir -p ";
        exec += outputDir;
        cout << exec << endl;
        gSystem->Exec(exec.Data());
    }
    
    
    TString plotsout = "";//"Canvas_pT_05_";
    plotsout += name;
    c->SaveAs(Form("%s/%s.root",outputDir.Data(),plotsout.Data()));
    c->SaveAs(Form("%s/%s.png",outputDir.Data(),plotsout.Data()));
    c->SaveAs(Form("%s/%s.pdf",outputDir.Data(),plotsout.Data()));
    c->SaveAs(Form("%s/%s.eps",outputDir.Data(),plotsout.Data()));
    
}


/*
//________________________________________________
void PlotComparisonspp_pPb1(){
    
    
    TFile * file1 = TFile::Open("~/Analysis/HFCorrelations/testScript/ReflectedPlots/StdRebin/AllPlots/Averages/ComparisonPPtoPPB/plotComparison_WeightedAverage_pp_pPb_58_0.3to99.0.root");
    if(!file1->IsOpen()){
        cout << "cannot open file 1" << endl; return;
    }
    TPad*c1 = (TPad*)file1->Get("ptmidsub");
    if(!c1){
        cout << "cannot open canvas 1 " << endl; return;
    }
    
    TFile * file2 = TFile::Open("~/Analysis/HFCorrelations/testScript/ReflectedPlots/StdRebin/AllPlots/Averages/ComparisonPPtoPPB/plotComparison_WeightedAverage_pp_pPb_816_0.3to99.0.root");
    if(!file2->IsOpen()){
        cout << "cannot open file 2" << endl; return;
    }
    TPad*c2 = (TPad*)file2->Get("pthighsub");
    if(!c2){
        cout << "cannot open canvas 2 " << endl; return;
    }
    //=======================
    
    TFile * file3 = TFile::Open("~/Analysis/HFCorrelations/testScript/ReflectedPlots/StdRebin/AllPlots/Averages/ComparisonPPtoPPB/plotComparison_WeightedAverage_pp_pPb_58_0.3to1.0.root");
    if(!file3->IsOpen()){
        cout << "cannot open file 3" << endl; return;
    }
    TPad*c3 = (TPad*)file3->Get("ptmidsub");
    if(!c3){
        cout << "cannot open canvas 3 " << endl; return;
    }
    
    TFile * file4 = TFile::Open("~/Analysis/HFCorrelations/testScript/ReflectedPlots/StdRebin/AllPlots/Averages/ComparisonPPtoPPB/plotComparison_WeightedAverage_pp_pPb_816_0.3to1.0.root");
    if(!file4->IsOpen()){
        cout << "cannot open file 4" << endl; return;
    }
    TPad*c4 = (TPad*)file4->Get("pthighsub");
    if(!c4){
        cout << "cannot open canvas 4 " << endl; return;
    }
    //=======================
    
    
    TFile * file5 = TFile::Open("~/Analysis/HFCorrelations/testScript/ReflectedPlots/StdRebin/AllPlots/Averages/ComparisonPPtoPPB/plotComparison_WeightedAverage_pp_pPb_58_1.0to99.0.root");
    if(!file5->IsOpen()){
        cout << "cannot open file 5" << endl; return;
    }
    TPad*c5 = (TPad*)file5->Get("ptmidsub");
    if(!c5){
        cout << "cannot open canvas 5 " << endl; return;
    }
    
    TFile * file6 = TFile::Open("~/Analysis/HFCorrelations/testScript/ReflectedPlots/StdRebin/AllPlots/Averages/ComparisonPPtoPPB/plotComparison_WeightedAverage_pp_pPb_816_1.0to99.0.root");
    if(!file6->IsOpen()){
        cout << "cannot open file 6" << endl; return;
    }
    TPad*c6 = (TPad*)file6->Get("pthighsub");
    if(!c6){
        cout << "cannot open canvas 6 " << endl; return;
    }
    //=======================
    
    
    TCanvas * c = new TCanvas("coutput","coutput",0,0,1200,800);
    c->Divide(2,3);
    c->cd(1); c1->Draw();
    c->cd(2); c2->Draw();
    c->cd(3); c3->Draw();
    c->cd(4); c4->Draw();
    c->cd(5); c5->Draw();
    c->cd(6); c6->Draw();
    
    
    Convert3x2Matrix(c,kTRUE);
}*/


void MergePPandPPbInSingleCanvas(TString strFilePP="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015Ago26Draft2/ReflectedPlots/StdRebin/AllPlots/NiceStylePlots/Output_Plots/WeightedAverageDzeroDstarDplus/CanvasNoSpaces_WeightedAverageDzeroDstarDplus_pp.root",TString strFilePPb="/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015Ago26Draft2/ReflectedPlots/StdRebin/AllPlots/NiceStylePlots/Output_Plots/WeightedAverageDzeroDstarDplus/CanvasNoSpaces_WeightedAverageDzeroDstarDplus_pPb.root"){

  gStyle->SetOptStat(0000);
  TFile *fpp=(TFile*)TFile::Open(strFilePP.Data(),"READ");
  TCanvas *cPP=(TCanvas*)fpp->Get("c_ouput");
  //  cPP->Draw();
  TFile *fpPb=(TFile*)TFile::Open(strFilePPb.Data(),"READ");
  TCanvas *cpPb=(TCanvas*)fpPb->Get("c_ouput");
  //  cpPb->Draw();
  Double_t rangesY[3]={14.,11.,6.8};
  TPad *pd;
  TH1D *h;

  /*TPad *pd=(TPad*)cpPb->GetPad(1);
  TH1D *h=(TH1D*)pd->FindObject("hDraw");
  rangesY[0]=h->GetMaximum();

  pd=(TPad*)cpPb->GetPad(3);
  h=(TH1D*)pd->FindObject("hDraw");
  rangesY[1]=h->GetMaximum();

  pd=(TPad*)cpPb->GetPad(5);
  h=(TH1D*)pd->FindObject("hDraw");
  rangesY[2]=h->GetMaximum();
  */
  
  TH1D **hpp=new TH1D*[9];
  TH1D **hpPb=new TH1D*[6];
  TH1D **hDraw=new TH1D*[9];

  TGraphAsymmErrors **grpp=new TGraphAsymmErrors*[9];
  TGraphAsymmErrors **grpPb=new TGraphAsymmErrors*[6];
  TLatex **tlscalepPb=new TLatex*[6];
  TLatex **tlscalepp=new TLatex*[9];

  TLatex **tlptD=new TLatex*[9];
  TLatex **tlptAssoc=new TLatex*[9];

  Int_t padorderingPP[9];

  for(Int_t j=1;j<=9;j++){
    TPad *pd=(TPad*)cPP->GetPad(j);
    TList *lpd=(TList*)pd->GetListOfPrimitives();
    Int_t binD,binass;
    for(Int_t il=0;il<lpd->GetEntries();il++){
      TObject* obj=lpd->At(il);
      TString strObjName=obj->ClassName();
      if(strObjName.Contains("TH1D")){	
	TString strhName=obj->GetName();
	if(strhName.EqualTo("hDraw")){
	  hDraw[j-1]=(TH1D*)((TH1D*)obj)->Clone();;
	  hDraw[j-1]->GetYaxis()->CenterTitle();
	  hDraw[j-1]->GetXaxis()->CenterTitle();
	  hDraw[j-1]->GetYaxis()->SetTitle("#frac{1}{#it{N}_{D}} #frac{d#it{N}^{assoc}}{d#Delta#varphi} (rad^{-1})");
	}
	else{
	  hpp[j-1]=(TH1D*)((TH1D*)obj)->Clone();
	  hpp[j-1]->SetName("fhDaveragepp");
	  hpp[j-1]->SetLineColor(ppcolor);
	  hpp[j-1]->SetMarkerColor(ppcolor);
	  hpp[j-1]->SetMarkerStyle(markerstyle[0]);
	}
      }
      else if(strObjName.Contains("TGraphAsymmErrors")){
	grpp[j-1]=((TGraphAsymmErrors*))((TGraphAsymmErrors*)obj)->Clone();;
	grpp[j-1]->SetName("grDaveragepp");
	grpp[j-1]->SetLineColor(ppcolor);
	grpp[j-1]->SetMarkerColor(ppcolor);
	grpp[j-1]->SetMarkerStyle(markerstyle[0]);
      }
      else if(strObjName.Contains("TLatex")){
	//	Printf("found latex :D ");
	TLatex *tl=(TLatex*)obj;// HEREEEEE 
	TString str=tl->GetTitle();
	if(str.Contains("<5")){
	  binD=0;
	  tlptD[j-1]=(TLatex*)((TLatex*)tl)->Clone();;
	  TString strlatTitle=tlptD[j-1]->GetTitle();
	  strlatTitle.ReplaceAll("3.0","3");
	  strlatTitle.ReplaceAll("5.0","5");
	  tlptD[j-1]->SetTitle(strlatTitle.Data());
	}
	else if(str.Contains("<8")){
	  tlptD[j-1]=(TLatex*)((TLatex*)tl)->Clone();;
	  TString strlatTitle=tlptD[j-1]->GetTitle();
	  strlatTitle.ReplaceAll("5.0","5");
	  strlatTitle.ReplaceAll("8.0","8");
	  tlptD[j-1]->SetTitle(strlatTitle.Data());
	  binD=1;
	}
	else if(str.Contains("<16")){
	  tlptD[j-1]=(TLatex*)((TLatex*)tl)->Clone();;
	  TString strlatTitle=tlptD[j-1]->GetTitle();
	  strlatTitle.ReplaceAll("8.0","8");
	  strlatTitle.ReplaceAll("16.0","16");
	  tlptD[j-1]->SetTitle(strlatTitle.Data());
	  binD=2;
	}
	if(str.Contains(">0.3")){
	  tlptAssoc[j-1]=(TLatex*)((TLatex*)tl)->Clone();;
	  binass=0;
	}
	else if(str.Contains("<1")){
	  tlptAssoc[j-1]=(TLatex*)((TLatex*)tl)->Clone();;
	  binass=1;
	}
	else if(str.Contains(">1")){
	  tlptAssoc[j-1]=(TLatex*)((TLatex*)tl)->Clone();;
	  binass=2;
	}
      }
    }
    for(Int_t il=0;il<lpd->GetEntries();il++){
      TObject* obj=lpd->At(il);
      TString strObjName=obj->ClassName();
      if(strObjName.Contains("TLatex")){
	TLatex *tl=(TLatex*)obj;// HEREEEEE 
	TString str=tl->GetTitle();
	if(str.Contains("scale")){
	  str.ReplaceAll("uncertainty","uncertainty (pp)");
	  tl->SetTitle(str.Data());
	  tl->SetY(0.68);//y+0.05);
	  tl->SetX(0.23);
	  tlscalepp[binass+3*(binD)]=(TLatex*)((TLatex*)tl)->Clone();
	}
	else if(binD!=0){
	  
	  Double_t y=tl->GetY();
	  tl->SetY(y+0.05);	  
	}
	
      }
    }
    
    hDraw[j-1]->GetYaxis()->SetRangeUser(0,rangesY[binass]);
    Printf("Setting pp index: %d,%d, to %d",binass,binD,j);
    padorderingPP[binass+3*binD]=j;
    
  }



  for(Int_t j=1;j<=6;j++){
    TPad *pd=(TPad*)cpPb->GetPad(j);
    TList *lpd=(TList*)pd->GetListOfPrimitives();
    Int_t binass=-1;
    Int_t binD=-1;
    TH1D *htemp;
    TGraphAsymmErrors *grtemp;
    TLatex *tltemp;
    for(Int_t il=0;il<lpd->GetEntries();il++){
      TObject* obj=lpd->At(il);
      TString strObjName=obj->ClassName();

      if(strObjName.Contains("TH1D")){	
	TString strhName=obj->GetName();
	if(strhName.EqualTo("hDraw")){
	  continue;
	  //	  TH1D *h=(TH1D*)obj;
	  //	  h->GetYaxis()->SetRangeUser(0,rangesY[(j-1)/3]);
	}
	else{
	  htemp=(TH1D*)obj;
	  htemp->SetName("fhDaveragepPb");
	  htemp->SetLineColor(pPbcolor);
	  htemp->SetMarkerColor(pPbcolor);
	  htemp->SetMarkerStyle(markerstyle[1]);
	}
      }
      else if(strObjName.Contains("TGraphAsymmErrors")){
	grtemp=(TGraphAsymmErrors*)obj;
	grtemp->SetName("grDaveragepPb");
	grtemp->SetLineColor(pPbcolor);
	grtemp->SetMarkerColor(pPbcolor);
	grtemp->SetMarkerStyle(markerstyle[1]);
      }
      else if(strObjName.Contains("TLatex")){
	//	Printf("found latex :D ");
	TLatex *tl=(TLatex*)obj;// HEREEEEE 
	TString str=tl->GetTitle();
	if(str.Contains("<5")){
	  Printf("SHOULD NOT HAPPEN;");
	  return;
	}
	else if(str.Contains("<8")){
	  binD=1;
	}
	else if(str.Contains("<16")){
	  binD=2;
	}
	if(str.Contains(">0.3")){
	  binass=0;
	}
	else if(str.Contains("<1")){
	  binass=1;
	}
	else if(str.Contains(">1")){
	  binass=2;
	}
	else if(str.Contains("scale")){
	  Double_t y=tl->GetY();
	  str.ReplaceAll("uncertainty","uncertainty (p-Pb)");
	  tl->SetTitle(str.Data());
	  tl->SetY(0.58);
	  tl->SetX(0.23);
	  tltemp=tl;
	}
      }
    }
    hpPb[3*(binD-1)+binass]=htemp;
    grpPb[3*(binD-1)+binass]=grtemp;    
    tlscalepPb[3*(binD-1)+binass]=tltemp;
  }

  Int_t nskip=0;  
  cPP->SetCanvasSize(900,900);
  cPP->Draw();

  //  TCanvas *cPPcp=(TCanvas*)cPP->Clone("cPPandPPb");
  for(Int_t jD=0;jD<=2;jD++){
    for(Int_t jass=0;jass<=2;jass++){
          
      cPP->cd(padorderingPP[jD*3+jass]);
      Printf("Adding histo ptD=%d,ptAss=%d of p-Pb to canvas %d",jD,jass,padorderingPP[jass+3*jD]);
      hDraw[padorderingPP[jD*3+jass]-1]->Draw();
      Printf("hdraw added");
      hpp[padorderingPP[jD*3+jass]-1]->Draw("same");
      Printf("hpp added");
      grpp[padorderingPP[jD*3+jass]-1]->Draw("p2");
      Printf("grpp added");
      if(jD==0){
	tlscalepp[padorderingPP[jD*3+jass]-1]->SetX(0.32);
	if(jass==0){
	  tlscalepp[padorderingPP[jD*3+jass]-1]->SetY(0.25);
	  TLatex *tlALICE=new TLatex(0.8,0.84,"ALICE");
	  tlALICE->SetNDC();
	  tlALICE->SetTextFont(43);
	  tlALICE->SetTextSize(24);
	  tlALICE->Draw();
	  TLatex *tlAvD=new TLatex(0.27,0.84,"Average D^{0}, D^{+}, D^{*+}");
	  tlAvD->SetNDC();
	  tlAvD->SetTextFont(43);
	  tlAvD->SetTextSize(22);
	  tlAvD->Draw();
	  TLegend *leg=new TLegend(0.255,0.56,0.4,0.8);
	  leg->AddEntry(hpp[0],"pp, #sqrt{#it{s}}=7 TeV, |#it{y}^{D}|<0.5","p");
	  leg->AddEntry(hpPb[0],"p-Pb, #sqrt{#it{s}_{NN}}=5.02 TeV,","p");//"p-Pb, #sqrt{s_{NN}}=5.02 TeV,-0.96<y^{D}_{cms}<0.04","p");
	  leg->AddEntry((TObject*)0,"-0.96<#it{y}^{D}_{cms}<0.04","");
	  leg->SetFillStyle(0);
	  leg->SetTextFont(43);
	  leg->SetTextSize(20);
	  leg->Draw();
	  TLatex *tlDEta=new TLatex(0.27,0.51,"|#Delta#eta|<1");//(0.5*(1-gPad->GetLeftMargin())+gPad->GetLeftMargin(),0.34,"|#Delta#eta|<1");
	  tlDEta->SetNDC();
	  tlDEta->SetTextAlign(12);
	  tlDEta->SetTextFont(43);
	  tlDEta->SetTextSize(19);
	  tlDEta->Draw();
	}
      }
      //      tlptD[padorderingPP[jD*3+jass]-1]->SetTextFont(63);
      tlptD[padorderingPP[jD*3+jass]-1]->SetTextFont(63);// do not why with 43 it appears in bold
      tlptD[padorderingPP[jD*3+jass]-1]->SetTextSize(19);
      tlptD[padorderingPP[jD*3+jass]-1]->SetY(0.86);      
      tlptD[padorderingPP[jD*3+jass]-1]->SetX(0.5*(1-gPad->GetLeftMargin())+gPad->GetLeftMargin());   
      //->SetX(0.5+gPad->GetLeftMargin());      
      tlptD[padorderingPP[jD*3+jass]-1]->SetTextAlign(22);
      tlptD[padorderingPP[jD*3+jass]-1]->Draw();

      tlptAssoc[padorderingPP[jD*3+jass]-1]->SetTextFont(63);
      tlptAssoc[padorderingPP[jD*3+jass]-1]->SetTextSize(19);
      tlptAssoc[padorderingPP[jD*3+jass]-1]->SetNDC();
      tlptAssoc[padorderingPP[jD*3+jass]-1]->SetY(0.78);
      tlptAssoc[padorderingPP[jD*3+jass]-1]->SetX(0.5*(1-gPad->GetLeftMargin())+gPad->GetLeftMargin());   
      tlptAssoc[padorderingPP[jD*3+jass]-1]->SetTextAlign(22);   
      tlptAssoc[padorderingPP[jD*3+jass]-1]->Draw();
      if(jD==0&&jass==0){
	//	tlptD[padorderingPP[jD*3+jass]-1]->SetY(0.51);
	//	tlptAssoc[padorderingPP[jD*3+jass]-1]->SetY(0.43);
	tlptD[padorderingPP[jD*3+jass]-1]->SetY(0.43);
	tlptAssoc[padorderingPP[jD*3+jass]-1]->SetY(0.34);
	tlscalepp[padorderingPP[jD*3+jass]-1]->SetY(0.12);
      }
      tlscalepp[padorderingPP[jD*3+jass]-1]->Draw();
      tlscalepp[padorderingPP[jD*3+jass]-1]->SetTextFont(63);
      tlscalepp[padorderingPP[jD*3+jass]-1]->SetTextSize(19);
      if(jD>=1){
	Printf("Adding histo %p",hpPb[jass+(jD-1)*3]);
	hpPb[jass+(jD-1)*3]->Draw("same");
	grpPb[jass+(jD-1)*3]->Draw("p2");
	tlscalepPb[jass+(jD-1)*3]->SetTextFont(63);
	tlscalepPb[jass+(jD-1)*3]->SetTextSize(19);
	tlscalepPb[jass+(jD-1)*3]->Draw();

      }
    }
  }
  cPP->Update();
  //  cPPcp->Draw();
  cPP->SaveAs("CanvasNoSpaces_WeightedAverageDzeroDstarDplus_ppAndpPb.root");
  cPP->SaveAs("CanvasNoSpaces_WeightedAverageDzeroDstarDplus_ppAndpPb.eps");
  cPP->SaveAs("CanvasNoSpaces_WeightedAverageDzeroDstarDplus_ppAndpPb.png");
  cPP->SaveAs("CanvasNoSpaces_WeightedAverageDzeroDstarDplus_ppAndpPb.pdf");
}
