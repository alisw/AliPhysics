/*

Macro to perform fits of the Laser Central electrode data
Several fit methods implemented

1. RebuildData() - transform arbitrary layeut of the Input data to the internal format
   StoreData();  - The data tree expected in file inname (see variable bellow)
   StoreTree();  - Modify inname and xxside and tcor in order to transform data

2. MakeFit();    - Make a fit of the data - already in internal format    
   StoreData();  - Store
   StoreTree();

3. LoadViewer(); - Browse the fit parameters

4. 
 

.x ~/rootlogon.C
gSystem->AddIncludePath("-I$ALICE_ROOT/TPC -I$ALICE_ROOT/STAT");
gSystem->Load("libSTAT.so");
.L $ALICE_ROOT/TPC/CalibMacros/AnalyzeLaser.C+




*/
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TStatToolkit.h"
#include "AliTPCCalibViewer.h"
#include "AliTPCCalibViewerGUI.h"
#include "AliTPCPreprocessorOnline.h"
//
//Define interesting variables - file names
//
char * inname = "treeCE08_05-07.root";  // input file with tree
//
// variable name definition in input tree - change it according the input
//
TString qaside("CE_A00_Q_05");
TString taside("CE_A00_T_05");
TString raside("CE_A00_rms_05");
TString qcside("CE_C00_Q_05");
TString tcside("CE_C00_T_05");
TString rcside("CE_C00_rms_05");
//
// correction variable - usually Pulser time
//
TString tcor("(sector%36>30)*2");

//
char * fname  = "treefitCE.root";       // output file with tree
char * oname  = "fitCE.root";           // output file with CalPads fit 

//
//
// Input CalPads
//
AliTPCCalPad *calPadIn  = 0;            // original time pad
AliTPCCalPad *calPadF1  = 0;            // original time pad - fit plane
AliTPCCalPad *calPadF2  = 0;            // original time pad - fit parabola
AliTPCCalPad *calPadQIn = 0;            // original Q pad
AliTPCCalPad *calPadQF1 = 0;            // original Q pad
AliTPCCalPad *calPadQF2 = 0;            // original Q pad

AliTPCCalPad *calPadCor = 0;            // base correction CalPad
AliTPCCalPad *calPadOut = 0;            // outlyer CalPad
//
// cuts
//
const Float_t tThr=0.5;                // max diff in the sector
const Float_t qThr0=0.5;               // max Q diff in the sector
const Float_t qThr1=2;                 // max Q diff in the sector

//
//
// fit Cal Pads
AliTPCCalPad *calPad0   = 0;            // global fit 0 - base
AliTPCCalPad *calPad1   = 0;            // global fit 1 - common behavior rotation -A-C
AliTPCCalPad *calPad2   = 0;            // gloabl fit 2 - CE missalign rotation     A-C
//
AliTPCCalPad *calPadInOut = 0;          // misaalign in-out
AliTPCCalPad *calPadLX    = 0;          // local x missalign
AliTPCCalPad *calPadTL   = 0;           // tan 
AliTPCCalPad *calPadQ     = 0;          // time (q)  correction
AliTPCCalPad *calPadGXY   = 0;          // global XY missalign (drift velocity grad) 
AliTPCCalPad *calPadOff   = 0;          // normalization offset fit

//
// working variables
//
AliTPCCalibViewerGUI * viewer=0;   //viewerGUI
AliTPCCalibViewer    *makePad=0;   //viewer
TTree * tree=0;                    // working tree

void LoadViewer();
void RebuildData();   // transform the input data to the fit format 
void MakeFit();       // make fits
//
//   internal functions
//
void MakeAliases0();  // Make Aliases 0  - for user tree
void MakeAliases1();  // Make Aliases 1  - for default tree   
void LoadData();      // Load data
void StoreData();     // store current data
void StoreTree();     // store fit data in the output tree


void AnalyzeLaser(){
  //
  //
  //
  LoadViewer();
  MakeAliases1();
}


void MakeFit(){
  //
  //
  TStatToolkit stat;
  Int_t npoints;
  Double_t chi2;
  TVectorD vec0,vec1,vec2;
  TMatrixD mat;
  TString fstring="";
  //
  //Basic  correction
  //
  fstring+="side++";        // offset on 2 different sides              //1
  fstring+="(1/qp)++";      // Q -threshold effect correction           //2
  fstring+="(qp)++";        // Q -threshold effect correction           //3
  fstring+="(inn)++";       //  inner outer misalign   - common         //4 
  fstring+="(side*inn)++";  //                         - opposite       //5
  //
  fstring+="(gyr)++";       // drift velocity gradient - common         //6
  fstring+="(side*gyr)++";  //                         - opposite       //7
  fstring+="(gxr)++";       //  global x tilt          - common         //8
  fstring+="(side*gxr)++";  //                         - opposite       //9
  //
  fstring+="tl^2++";        // local phi correction                     //10
  //
  fstring+="(lxr)++";       // zr            angle      - common        //11
  fstring+="(side*lxr)++";  //                          - opposite      //12
  fstring+="(inn*lxr)++";   // inner outer angle        - common        //13             
  fstring+="(side*inn*lxr)++";//                        - opposite      //14
  fstring+="(lxr^2)++";       // zr          second     - common        //15
  fstring+="(side*lxr^2)++";  //                        - opposite      //16
  //
  TString *fit0 =stat.FitPlane(tree,"dt",fstring.Data(),"cutF&&cutCE",chi2,npoints,vec0,mat,0.90);
  tree->SetAlias("f0",fit0->Data());
  //
  // Common "deformation" tendencies
  //
  fstring+="(sin(atan2(gy.fElements,gx.fElements)))++";
  fstring+="(cos(atan2(gy.fElements,gx.fElements)))++";
  //
  fstring+="(sin(atan2(gy.fElements,gx.fElements)*2))++";
  fstring+="(cos(atan2(gy.fElements,gx.fElements)*2))++";
  fstring+="(sin(atan2(gy.fElements,gx.fElements)*3))++";
  fstring+="(cos(atan2(gy.fElements,gx.fElements)*3))++";
  //
  fstring+="(sin(atan2(gy.fElements,gx.fElements)*2))*lxr++";
  fstring+="(cos(atan2(gy.fElements,gx.fElements)*2))*lxr++";
  fstring+="(sin(atan2(gy.fElements,gx.fElements)*3))*lxr++";
  fstring+="(cos(atan2(gy.fElements,gx.fElements)*3))*lxr++";
  //

  TString *fit1 =stat.FitPlane(tree,"dt",fstring.Data(),"cutF&&cutCE",chi2,npoints,vec1,mat,0.95);
  tree->SetAlias("f1",fit1->Data());
  //
  // Central electrode "deformation"
  //
  fstring+="(side*sin(atan2(gy.fElements,gx.fElements)))++";
  fstring+="(side*cos(atan2(gy.fElements,gx.fElements)))++";
  //
  fstring+="(side*sin(atan2(gy.fElements,gx.fElements)*2))++";
  fstring+="(side*cos(atan2(gy.fElements,gx.fElements)*2))++";
  fstring+="(side*sin(atan2(gy.fElements,gx.fElements)*3))++";
  fstring+="(side*cos(atan2(gy.fElements,gx.fElements)*3))++";
  //
  fstring+="(side*sin(atan2(gy.fElements,gx.fElements)*2))*lxr++";
  fstring+="(side*cos(atan2(gy.fElements,gx.fElements)*2))*lxr++";
  fstring+="(side*sin(atan2(gy.fElements,gx.fElements)*3))*lxr++";
  fstring+="(side*cos(atan2(gy.fElements,gx.fElements)*3))*lxr++";
  
  TString *fit2 =stat.FitPlane(tree,"dt",fstring.Data(),"cutF&&abs(dt-f0)<0.7&&cutCE",chi2,npoints,vec2,mat,0.90);
  tree->SetAlias("f2",fit2->Data());
  //
  // Extract variables
  //
  TString tmpstr = fstring;
  TObjArray *arr = tmpstr.Tokenize("++");
  TString fitQ("0");       // q correction 
  TString fitLX("0");      // lx correction 
  TString fitInOut("0");   // inner-outer - match
  TString fitGXY("0");      // global xy fit
  TString fitOff("0");  // side offsets
  TString fitTL("0");  // side offsets
  //
  fitOff+="+";
  fitOff+=vec2[0];
  fitOff+="+side*";
  fitOff+=vec2[1];
  {
  for(Int_t i=0;i<arr->GetEntriesFast();i++){
    if (!arr->At(i)) continue;
    TString *fitstr = new TString(arr->At(i)->GetName());
    //
    Bool_t isQ      = fitstr->Contains("qp)");
    Bool_t isRot    = fitstr->Contains("sin(")+fitstr->Contains("cos(");
    Bool_t isLX     = fitstr->Contains("lxr");
    Bool_t isIn     = fitstr->Contains("inn");
    Bool_t isGXY    = fitstr->Contains("gxr")+fitstr->Contains("gyr");
    if (fitstr->Contains("tl^2")){
      fitTL+="+";
      fitTL+=(*fitstr)+"*";
      fitTL+=vec2[i+1];
    }
    if (isGXY){
      fitGXY+="+";
      fitGXY+=(*fitstr)+"*";
      fitGXY+=vec2[i+1];
    }
    if (isQ){
      //
      fitQ+="+";
      fitQ+=(*fitstr)+"*";
      fitQ+=vec2[i+1];
    }
    //
    if (isLX&&!isRot&&!isIn){
      fitLX+="+";
      fitLX+=(*fitstr)+"*";
      fitLX+=vec2[i+1];
    }
    //
    if (!isRot&&isIn){
      fitInOut+="+";
      fitInOut+=(*fitstr)+"*";
      fitInOut+=vec2[i+1];
    }
  }
  }
  //
  tree->SetAlias("fInOut",fitInOut.Data());
  tree->SetAlias("fLX",fitLX.Data());
  tree->SetAlias("fGXY",fitGXY.Data());
  tree->SetAlias("fOff",fitOff.Data());
  tree->SetAlias("fQ",fitQ.Data());
  tree->SetAlias("fTL",fitTL.Data());
  //
  //
  // fits
  // 
  calPad0 = makePad->GetCalPad("f0","cutCE", "ffit0");
  calPad1 = makePad->GetCalPad("f1","cutCE", "ffit1");
  calPad2 = makePad->GetCalPad("f2","cutCE", "ffit2");
  calPadInOut = makePad->GetCalPad("fInOut","cutCE", "fInOut");
  calPadLX    = makePad->GetCalPad("fLX","cutCE", "fLX");
  calPadTL    = makePad->GetCalPad("fTL","cutCE", "fTL");
  calPadQ      = makePad->GetCalPad("fQ","cutCE", "fQ");
  calPadGXY    = makePad->GetCalPad("fGXY","cutCE", "fGXY");
  calPadOff    = makePad->GetCalPad("fOff","cutCE", "fOff");  
}

void LoadViewer(){
  //
  // Load calib Viewer
  //
  TObjArray * array = AliTPCCalibViewerGUI::ShowGUI(fname);
  viewer = (AliTPCCalibViewerGUI*)array->At(0);
  makePad = viewer->GetViewer();
  tree = viewer->GetViewer()->GetTree();
  MakeAliases1();  
}






void RebuildData(){
  //
  // transform the input data to the fit format 
  //
  makePad = new AliTPCCalibViewer(inname);
  tree = makePad->GetTree();
  MakeAliases0(); //
  //
  calPadCor = makePad->GetCalPad("tcor","(cutA||cutC)", "tcor");
  calPadOut = makePad->GetCalPad("1","!((cutA||cutC)&&abs(ly.fElements/lx.fElements)<0.155)", "out");
  calPadIn  = makePad->GetCalPad("dt-tcor","(cutA||cutC)&&abs(ly.fElements/lx.fElements)<0.155","timeIn");
  calPadF1  = calPadIn->GlobalFit("timeF1", calPadOut,kTRUE,0,0.9);
  calPadQIn = makePad->GetCalPad("qa*(sector%36<18)+qc*(sector%36>17)","1","qIn");
  //
  // Update outlyer maps
  //
  for (Int_t isector=0;isector<72; isector++){
    for (UInt_t ich=0;ich<calPadIn->GetCalROC(isector)->GetNchannels();ich++){
      Float_t val0= calPadIn->GetCalROC(isector)->GetValue(ich);
      Float_t val1= calPadF1->GetCalROC(isector)->GetValue(ich);
      if (TMath::Abs(val0-val1)>tThr) calPadOut->GetCalROC(isector)->SetValue(ich,1);
    }
  }
  calPadF1  = calPadIn->GlobalFit("timeF1", calPadOut,kTRUE,0,0.9);
  calPadF2  = calPadIn->GlobalFit("timeF2", calPadOut,kTRUE,1,0.9);
  calPadQF1 = calPadQIn->GlobalFit("qF1", calPadOut,kTRUE,0,0.9);
  calPadQF2 = calPadQIn->GlobalFit("qF2", calPadOut,kFALSE,1);
  //
  // Update outlyer maps
  //
  for (Int_t isector=0;isector<72; isector++){
    for (UInt_t ich=0;ich<calPadIn->GetCalROC(isector)->GetNchannels();ich++){
      Float_t val0= calPadQIn->GetCalROC(isector)->GetValue(ich);
      Float_t val1= calPadQF2->GetCalROC(isector)->GetValue(ich);
      if (val1<0.1)  {
	calPadOut->GetCalROC(isector)->SetValue(ich,1);
	continue;
      }
      if (TMath::Abs(val0/val1)<qThr0) calPadOut->GetCalROC(isector)->SetValue(ich,1);
      if (TMath::Abs(val0/val1)>qThr1) calPadOut->GetCalROC(isector)->SetValue(ich,1);
    }
  }
  calPadF1  = calPadIn->GlobalFit("timeF1", calPadOut,kTRUE,0,0.9);
  calPadF2  = calPadIn->GlobalFit("timeF2", calPadOut,kTRUE,1,0.9);
  calPadQF1 = calPadQIn->GlobalFit("qF1", calPadOut,kTRUE,0,0.9);
  calPadQF2 = calPadQIn->GlobalFit("qF2", calPadOut,kFALSE,1);
}

void LoadData(){
  //
  // Get Data
  //
  TFile f(oname);
  calPadIn  = (AliTPCCalPad*)f.Get("timeIn");  // original time pad
  calPadF1  = (AliTPCCalPad*)f.Get("timeF1");  // original time pad - fit plane
  calPadF2  = (AliTPCCalPad*)f.Get("timeF2");  // original time pad - fit parabola
  //
  calPadQIn  = (AliTPCCalPad*)f.Get("qIn");  // original time pad
  calPadQF1  = (AliTPCCalPad*)f.Get("qF1");  // original time pad - fit plane
  calPadQF2  = (AliTPCCalPad*)f.Get("qF2");  // original time pad - fit parabola
  //
  calPadCor = (AliTPCCalPad*)f.Get("tcor");    // base correction CalPad
  calPadOut = (AliTPCCalPad*)f.Get("out");     // outlyer CalPad  
  //
  calPad0   = (AliTPCCalPad*)f.Get("ffit0");   // global fit 0 - base
  calPad1   = (AliTPCCalPad*)f.Get("ffit1");   // global fit 1 - common behavior rotation -A-C
  calPad2   = (AliTPCCalPad*)f.Get("ffit2");   // gloabl fit 2 - CE missalign rotation     A-C
  calPadInOut = (AliTPCCalPad*)f.Get("fInOut");// misaalign in-out
  calPadLX    = (AliTPCCalPad*)f.Get("fLX");   // local x missalign
  calPadTL    = (AliTPCCalPad*)f.Get("fTL");   // local y/x missalign
  calPadQ     = (AliTPCCalPad*)f.Get("fQ");    // time (q)  correction
  calPadGXY   = (AliTPCCalPad*)f.Get("fGXY");  // global XY missalign (drift velocity grad) 
  calPadOff   = (AliTPCCalPad*)f.Get("fOff");  // normalization offset fit
}

void StoreData(){
  //
  // Store data
  // 
  TFile * fstore = new TFile(oname,"recreate");
  if (calPadIn) calPadIn->Write("timeIn");   // original time pad
  if (calPadF1) calPadF1->Write("timeF1");   // original time pad - fit plane
  if (calPadF2) calPadF2->Write("timeF2");   // original time pad - fit parabola
  //
  if (calPadQIn) calPadQIn->Write("qIn");   // original time pad
  if (calPadQF1) calPadQF1->Write("qF1");   // original time pad - fit plane
  if (calPadQF2) calPadQF2->Write("qF2");   // original time pad - fit parabola
  //
  if (calPadCor) calPadCor->Write("tcor");   // base correction CalPad
  if (calPadOut) calPadOut->Write("out");    // outlyer CalPad  
  //
  if (calPad0)   calPad0->Write("ffit0");    // global fit 0 - base
  if (calPad1)   calPad1->Write("ffit1");    // global fit 1 - common behavior rotation -A-C
  if (calPad2)   calPad2->Write("ffit2");    // gloabl fit 2 - CE missalign rotation     A-C
  if (calPadInOut)calPadInOut->Write("fInOut");   // misaalign in-out
  if (calPadLX)  calPadLX->Write("fLX");     // local x missalign
  if (calPadTL)  calPadTL->Write("fTL");     // local y/x missalign
  if (calPadQ)   calPadQ->Write("fQ");       // time (q)  correction
  if (calPadGXY) calPadGXY->Write("fGXY");   // global XY missalign (drift velocity grad) 
  if (calPadOff) calPadOff->Write("fOff");   // normalization offset fit
  fstore->Close();
  delete fstore;
}

void StoreTree(){
  //
  //
  //
  AliTPCPreprocessorOnline * preprocesor = new AliTPCPreprocessorOnline;
  //
  if (calPadIn) preprocesor->AddComponent(calPadIn->Clone());
  if (calPadF1) preprocesor->AddComponent(calPadF1->Clone());   
  if (calPadF2) preprocesor->AddComponent(calPadF2->Clone());   
  //
  if (calPadQIn) preprocesor->AddComponent(calPadQIn->Clone());
  if (calPadQF1) preprocesor->AddComponent(calPadQF1->Clone());   
  if (calPadQF2) preprocesor->AddComponent(calPadQF2->Clone());   
  //
  if (calPadCor) preprocesor->AddComponent(calPadCor->Clone());   
  if (calPadOut) preprocesor->AddComponent(calPadOut->Clone());  
  //
  if (calPad0)   preprocesor->AddComponent(calPad0->Clone());
  if (calPad1)   preprocesor->AddComponent(calPad1->Clone());
  if (calPad2)   preprocesor->AddComponent(calPad2->Clone());
  if (calPadInOut)preprocesor->AddComponent(calPadInOut->Clone());
  if (calPadLX)  preprocesor->AddComponent(calPadLX->Clone());
  if (calPadTL)  preprocesor->AddComponent(calPadTL->Clone());
  if (calPadQ)   preprocesor->AddComponent(calPadQ->Clone());
  if (calPadGXY) preprocesor->AddComponent(calPadGXY->Clone());
  if (calPadOff) preprocesor->AddComponent(calPadOff->Clone());
  preprocesor->DumpToFile(fname);
  delete preprocesor;
}


void MakeAliases0(){
  //
  // Define variables and selection of outliers - for user defined tree
  //
  tree->SetAlias("tcor",tcor.Data());          // correction variable
  tree->SetAlias("ta",taside+".fElements");
  tree->SetAlias("tc",tcside+".fElements");
  tree->SetAlias("qa",qaside+".fElements");
  tree->SetAlias("qc",qcside+".fElements");
  tree->SetAlias("ra",raside+".fElements");
  tree->SetAlias("rc",rcside+".fElements");
  tree->SetAlias("side","1-(sector%36>17)*2");
  tree->SetAlias("dt","(ta)*(sector%36<18)+(tc)*(sector%36>17)+tcor");
  tree->SetAlias("cutA","qa>30&&qa<400&&abs(ta)<2&&ra>0.5&&ra<2");
  tree->SetAlias("cutC","qc>30&&qc<400&&abs(tc)<2&&rc>0.5&&rc<2");
  tree->SetAlias("cutF","(pad.fElements%4==0)&&(row.fElements%3==0)");
  tree->SetAlias("cutCE","V.out.fElements");
  //
  // fit param aliases
  //
  tree->SetAlias("inn","sector<36");
  tree->SetAlias("gxr","(gx.fElements/250.)"); //
  tree->SetAlias("gyr","(gy.fElements/250.)"); //
  tree->SetAlias("lxr","(lx.fElements-133.41)/250.");
  tree->SetAlias("qp","((sector%36<18)*sqrt(qa)/10.+(sector%36>17)*sqrt(qc)/10.)"); //
  tree->SetAlias("tl","(ly.fElements/lx.fElements)/0.17");  
}


void MakeAliases1(){
  //
  // Define variables and selection of outliers -for default usage
  //
  tree->SetAlias("tcor","tcor.fElements");          // correction variable  
  tree->SetAlias("side","1-(sector%36>17)*2");
  tree->SetAlias("dt","timeIn.fElements+tcor");
  //
  tree->SetAlias("cutA","out.fElements==1");
  tree->SetAlias("cutC","out.fElements==1");
  tree->SetAlias("cutF","(pad.fElements%5==0)&&(row.fElements%4==0)");
  tree->SetAlias("cutCE","out.fElements<0.5");
  //
  // fit param aliases
  //
  tree->SetAlias("inn","sector<36");
  tree->SetAlias("gxr","(gx.fElements/250.)"); //
  tree->SetAlias("gyr","(gy.fElements/250.)"); //
  tree->SetAlias("lxr","(lx.fElements-133.41)/250.");
  tree->SetAlias("qp","(sqrt(qIn.fElements)/10.)"); //
  tree->SetAlias("tl","(ly.fElements/lx.fElements)/0.17");  
}


