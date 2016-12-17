/*
  NDimensional formula browser. 
  Create browser  (set of sliders) for the formula visualization.
  
  Current status:
    Working prototype
        Used e.g for browsing of the TPC space charge distortion maps
    To do ? 
        Interface THn and TTree ?

  Example usage:   
  TFormula * formula = new TFormula("formula","[0]+[1]*(sin([2]*exp([3])))")
  formula->SetParName(0,"Offset");
  formula->SetParName(1,"Slope");
  formula->SetParName(2,"Freq.");
  formula->SetParName(3,"Exp. slope");
  formula->SetParameter(0,0);
  formula->SetParameter(1,1.);  formula->SetParameter(2,1.);
  formula->SetParameter(3,1.);
  AliNDFormulaBrowser::SetDefaultStyle()
  AliNDFormulaBrowser *browser = new AliNDFormulaBrowser(0, 0, formula, 1200, 1000); 


*/

#include <THn.h>
#include <TH2.h>
#include <TRandom.h>
#include <TMath.h>

#include <TPad.h>
#include <TCanvas.h>
#include <TRootEmbeddedCanvas.h>
#include <TGWindow.h>
#include <TGDoubleSlider.h> 
#include <TGSlider.h>
#include <TROOT.h> 
#include "TGTextBuffer.h"
#include "TGTextEntry.h"
#include "TGLabel.h"
#include "TFormula.h"
#include "TVectorD.h"
#include "TVectorF.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include <cstdio>
#include <map>
#include <string>
#include "TClass.h"
#include <RQ_OBJECT.h>
#include "AliNDFormulaBrowser.h"


Int_t nBins=200;

ClassImp(AliNDFormulaBrowser)

TLatex *AliNDFormulaBrowser::fgkLatex = 0;
TLegend *AliNDFormulaBrowser::fgkLegend = 0;
std::map<std::string,TFormula*>   AliNDFormulaBrowser::fgkFormulaMap;   // map of registered formulas 


void AliNDFormulaBrowser::SetDefaultStyle(){
  //  fgkFormulaMap= new map<std::string,TFormula*>;
  fgkLatex=new TLatex;
  fgkLatex->SetX(0.15);
  fgkLatex->SetY(0.85);
  fgkLatex->SetTextSize(0.03);
  //
  fgkLegend= new TLegend(0.75,0.75,0.89,0.89,"N-dim formula browser");
  fgkLegend->SetBorderSize(0);

}

AliNDFormulaBrowser::AliNDFormulaBrowser(const TGWindow *p, const TGWindow *main,TFormula *formula,
					 UInt_t width, UInt_t height):
  fFormula(formula),
  fFormulaParams(0)
{
   // Dialog used to test the different supported sliders.
  if (formula==NULL) return;
  fMain = new TGTransientFrame(p, main, width, height);
  fMain->Connect("CloseWindow()", "AliNDFormulaBrowser", this, "CloseWindow()");
  fMain->DontCallClose(); // to avoid double deletions.
  // use hierarchical cleaning
  fMain->SetCleanup(kDeepCleanup);
  fMain->ChangeOptions((fMain->GetOptions() & ~kVerticalFrame) | kVerticalFrame);
  fCanvas = new TRootEmbeddedCanvas("Ecanvas",fMain,200,200);
  fMain->AddFrame(fCanvas, new TGLayoutHints(kLHintsLeft | kLHintsExpandX| kLHintsExpandY, 0, 0, 3, 0));
  //
  
  fVframeFormula = new TGVerticalFrame(fMain, 0, 0, 0);
  fMain->AddFrame(fVframeFormula,  new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 3, 0));
  //
  //
  fFormula=formula;
  Int_t nPars = formula->GetNpar();
  fFormulaParams = new TVectorD(nPars);
  fParamFrame   = new TObjArray(nPars);
  fSliders=       new TObjArray(nPars);
  fCurrentEntry=  new TObjArray(nPars);
  fMinEntry=      new TObjArray(nPars);
  fMaxEntry=      new TObjArray(nPars);
  fParamLabels=   new TObjArray(nPars);
  //
  fFormulaFrame  = new TGGroupFrame(fVframeFormula,"Formula",kHorizontalFrame);
  fDrawFormula   = new TGTextEntry(fFormulaFrame,new TGTextBuffer(15) , TString(1+formula->GetName()).Hash());
  fDrawFormula->Connect("TextChanged(char*)", "AliNDFormulaBrowser", this, "DoText(char*)"); 
  fFormulaFrame->AddFrame(fDrawFormula,   new TGLayoutHints(kLHintsNormal | kLHintsExpandX,2,2,0,0));
  fDrawOption   = new TGTextEntry(fFormulaFrame,new TGTextBuffer(5) , TString(2+formula->GetName()).Hash());
  fDrawOption->Connect("TextChanged(char*)", "AliNDFormulaBrowser", this, "DoText(char*)"); 
  fFormulaFrame->AddFrame(fDrawOption,   new TGLayoutHints(kLHintsNormal,2,2,0,0));
  fFormulaEval   = new TGTextEntry(fFormulaFrame,new TGTextBuffer(5) , TString(3+formula->GetName()).Hash());
  fFormulaEval->Connect("TextChanged(char*)", "AliNDFormulaBrowser", this, "DoText(char*)"); 
  fFormulaFrame->AddFrame(fFormulaEval,   new TGLayoutHints(kLHintsNormal,2,2,0,0));
  //
  fFormulaNPoints   = new TGTextEntry(fFormulaFrame,new TGTextBuffer(5) , TString(4+formula->GetName()).Hash());
  fFormulaEval->Connect("TextChanged(char*)", "AliNDFormulaBrowser", this, "DoText(char*)"); 
  fFormulaFrame->AddFrame(fFormulaNPoints,   new TGLayoutHints(kLHintsNormal,2,2,0,0));
  fFormulaNLines   = new TGTextEntry(fFormulaFrame,new TGTextBuffer(5) , TString(5+formula->GetName()).Hash());
  fFormulaEval->Connect("TextChanged(char*)", "AliNDFormulaBrowser", this, "DoText(char*)"); 
  fFormulaFrame->AddFrame(fFormulaNLines,   new TGLayoutHints(kLHintsNormal,2,2,0,0));
  
  
  //
  for (Int_t iPar=0; iPar<nPars; iPar++){
    TGTextBuffer * buffer = 0;
     TString parName=formula->GetName();
     parName+=formula->GetParName(iPar);
     Int_t parID=parName.Hash();
     TGHorizontalFrame *hframe=new TGHorizontalFrame(fVframeFormula,20,20);
     TGLabel *label           =new TGLabel(hframe, TString::Format("p[%d]: %s", iPar, formula->GetParName(iPar)).Data());
     label->Connect("TextChanged(char*)", "AliNDFormulaBrowser", this, "DoText(char*)");
     hframe->AddFrame(label,   new TGLayoutHints(kLHintsNormal,2,2,0,0));
     TGTextEntry * textMin    =new TGTextEntry(hframe,new TGTextBuffer(5) , parID+2);
     textMin->Connect("TextChanged(char*)", "AliNDFormulaBrowser", this, "DoText(char*)");
     hframe->AddFrame(textMin, new TGLayoutHints(kLHintsNormal,2,2,0,0));
     TGTextEntry * textMax    =new TGTextEntry(hframe,new TGTextBuffer(5) , parID+3);
     textMax->Connect("TextChanged(char*)", "AliNDFormulaBrowser", this, "DoText(char*)");
     hframe->AddFrame(textMax, new TGLayoutHints(kLHintsNormal,2,2,0,0));
     TGTextEntry * text       =new TGTextEntry(hframe,new TGTextBuffer(5) , parID+1);
     text->Connect("TextChanged(char*)", "AliNDFormulaBrowser", this, "DoText(char*)");
     hframe->AddFrame(text,    new TGLayoutHints(kLHintsNormal,2,2,0,0));
     //
     //
     TGHSlider * slider = new TGHSlider(fVframeFormula, 100, kSlider1 | kScaleBoth, parID);
     slider->Connect("PositionChanged(Int_t)", "AliNDFormulaBrowser", this, "DoSlider(Int_t)");
     Double_t value=formula->GetParameter(iPar);
     slider->SetRange(0,nBins);
     slider->SetPosition(nBins/2);
     Double_t varMin=0, varMax=1;
     if (value>0){
       varMin=0;
       varMax=2*value;
     }
     buffer=textMin->GetBuffer();
     buffer->AddText(0, TString::Format("%.3f",varMin).Data());
     buffer=textMax->GetBuffer();
     buffer->AddText(0, TString::Format("%.3f",varMax).Data());
     buffer=text->GetBuffer();
     buffer->AddText(0, TString::Format("%.3f",(varMin+varMax)/2.).Data());

     //
     // store variables in array for later usage
     fParamFrame->AddAt(hframe,iPar);
     fParamLabels->AddAt(label, iPar);
     fMinEntry->AddAt(textMin,iPar);
     fMaxEntry->AddAt(textMax,iPar);
     fCurrentEntry->AddAt(text,iPar);
     fSliders->AddAt(slider,iPar); 
     //
   } 
   fVframeFormula->Resize(100, 100);
   //--- layout for buttons: top align, equally expand horizontally
   fBly = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 3, 0);
   //--- layout for the frame: place at bottom, right aligned
   fBfly1 = new TGLayoutHints(kLHintsTop | kLHintsRight, 20, 10, 15, 0);
   fVframeFormula->AddFrame(fFormulaFrame,fBly);
   for (Int_t iPar=0; iPar<nPars; iPar++){
     TGHorizontalFrame *hframe= (TGHorizontalFrame *)fParamFrame->At(iPar);
     TGHSlider * slider = (TGHSlider *)fSliders->At(iPar);
     fVframeFormula->AddFrame(hframe,fBly);
     fVframeFormula->AddFrame(slider,fBly);
   }
   //
   //
   fMain->SetWindowName("AliNDFromulaBrower");
   TGDimension size = fMain->GetDefaultSize();
   fMain->Resize(size);

   fMain->SetWMSize(size.fWidth, size.fHeight);
   fMain->SetWMSizeHints(size.fWidth, size.fHeight, size.fWidth, size.fHeight, 0, 0);
   fMain->SetMWMHints(kMWMDecorAll | kMWMDecorResizeH  | kMWMDecorMaximize |
                                     kMWMDecorMinimize | kMWMDecorMenu,
                      kMWMFuncAll |  kMWMFuncResize    | kMWMFuncMaximize |
                                     kMWMFuncMinimize,
                      kMWMInputModeless);

   // position relative to the parent's window
   fMain->CenterOnParent();
   fMain->MapSubwindows();
   fMain->MapWindow();
   gClient->WaitFor(fMain);
}

AliNDFormulaBrowser::~AliNDFormulaBrowser()
{
   // Delete dialog.
   fMain->DeleteWindow();  // deletes fMain
}

void AliNDFormulaBrowser::CloseWindow()
{
   // Called when window is closed via the window manager.
   delete this;
}

void AliNDFormulaBrowser::DoText(const char * /*text*/)
{
   // Handle text entry widgets.

   TGTextEntry *te = (TGTextEntry *) gTQSender;
   Int_t id = te->WidgetId();
   /*
   switch (id) {
      case HId1:
         fHslider1->SetPosition(atoi(fTbh1->GetString()));
         break;
      case VId1:
	fVslider1->SetPosition(atoi(fTbv1->GetString()));
         break;
      case HId2:
         fHslider2->SetPosition(atoi(fTbh2->GetString()));
         break;
      case VId2:
         fVslider2->SetPosition(atoi(fTbv2->GetString()),
                                     atoi(fTbv2->GetString())+2);
         break;
      default:
         break;
   }
   */
}

void AliNDFormulaBrowser::DoSlider(Int_t pos)
{
   // Handle slider widgets.

   Int_t id;
   TGFrame *frm = (TGFrame *) gTQSender;
   if (frm->IsA()->InheritsFrom(TGSlider::Class())) {
      TGSlider *sl = (TGSlider*) frm;
      id = sl->WidgetId();
   } else {
      TGDoubleSlider *sd = (TGDoubleSlider *) frm;
      id = sd->WidgetId();
   }

   char buf[32];
   sprintf(buf, "%d", pos);
   Int_t nPars = fFormula->GetNpar();
   for (Int_t iPar=0; iPar<nPars; iPar++){
     TString parName=fFormula->GetName();
     parName+=fFormula->GetParName(iPar);
     Int_t parID=parName.Hash();
     if (id==parID){
       TGHSlider * slider = (TGHSlider *)fSliders->At(iPar);
       TGTextEntry * text = (TGTextEntry *)fCurrentEntry->At(iPar);
       TGTextEntry * textMin = (TGTextEntry *)fMinEntry->At(iPar);
       TGTextEntry * textMax = (TGTextEntry *)fMaxEntry->At(iPar);
       //
       TGTextBuffer * buffer = text->GetBuffer();
       Double_t minPosS=slider->GetMinPosition();
       Double_t maxPosS=slider->GetMaxPosition();
       Double_t minPosT=atof(textMin->GetBuffer()->GetString());
       Double_t maxPosT=atof(textMax->GetBuffer()->GetString());
       Double_t valueS=(slider->GetPosition()-minPosS)/(maxPosS-minPosS);
       Double_t valueR=minPosT+valueS*(maxPosT-minPosT);       
       sprintf(buf, "%.4f", valueR);
       buffer->Clear();
       buffer->AddText(0, buf);
       text->SetCursorPosition(text->GetCursorPosition());
       text->Deselect();
       gClient->NeedRedraw(text);
       //printf("id\t%d\tpos\t%d\n",id, slider->GetPosition() );
     }
   }
   UpdateCanvas();
   UpdateFormula();
}

void AliNDFormulaBrowser::UpdateFormula(){
  //
  // 
  //
  for (Int_t iPar=0; iPar<fFormula->GetNpar(); iPar++){
    TGTextEntry * text = (TGTextEntry *)fCurrentEntry->At(iPar);
    (*fFormulaParams)[iPar]=atof(text->GetBuffer()->GetString());
  }
  Double_t value= fFormula->EvalPar(0,fFormulaParams->GetMatrixArray());
  TGTextBuffer * buffer = fFormulaEval->GetBuffer();
  buffer->Clear();
  buffer->AddText(0, TString::Format("%f",value).Data()); 
  fFormulaEval->Deselect();
  gClient->NeedRedraw(fFormulaEval);
}

void  AliNDFormulaBrowser::UpdateCanvas(){
  //
  // 
  //
  Int_t nPoints=200;
  Int_t nLines=3;
  const Int_t kMinPoints=100;
  const Int_t kMinGraphs=2;
  TLegend *legend = (TLegend*)fgkLegend->Clone(); // 
  //
  if (fFormulaNPoints->GetBuffer()->GetString()){
    Int_t nPoints2=atoi(fFormulaNPoints->GetBuffer()->GetString());
    nPoints=TMath::Max(nPoints2,kMinPoints);
  }
  if (fFormulaNLines->GetBuffer()->GetString()){
    Int_t nLines2=atoi(fFormulaNLines->GetBuffer()->GetString());
    nLines=TMath::Max(nLines2,kMinGraphs);
  }
  //
  Int_t nParams=fParamFrame->GetEntriesFast();
  TVectorF vecY(nPoints);
  TVectorF vecX(nPoints);
  TObjArray *formulaArray = TString(fDrawFormula->GetBuffer()->GetString()).Tokenize(":");
  Int_t entries=formulaArray->GetEntries();
  if (entries<2) return;
  TObjArray *formulaStack=TString(formulaArray->At(0)->GetName()).Tokenize(";");  // semicolomn separated list of functions
  Int_t nFormulas=formulaStack->GetEntries();
  TGraph *gr0=0;
  Float_t minF=0;
  Float_t maxF=-1;
  Int_t iPar = atoi(&((formulaArray->At(1)->GetName())[1]));
  Int_t jPar =-1;
  if (entries>2) jPar=atoi(&((formulaArray->At(2)->GetName())[1]));
  Int_t nGraphs=1;
  Double_t minPosT2=0, maxPosT2=0;

  for ( Int_t iForm=0; iForm<nFormulas; iForm++){
    TFormula *drawFormula=fFormula;
    if (fgkFormulaMap[formulaStack->At(iForm)->GetName()]) drawFormula=fgkFormulaMap[formulaStack->At(iForm)->GetName()];
    if (entries>1){
      if (jPar>0) {
	nGraphs=nLines;
	TGTextEntry * textMin = (TGTextEntry *)fMinEntry->At(jPar);
	TGTextEntry * textMax = (TGTextEntry *)fMaxEntry->At(jPar);
	minPosT2=atof(textMin->GetBuffer()->GetString());
	maxPosT2=atof(textMax->GetBuffer()->GetString());
      }
      for (Int_t igr=0; igr<nGraphs; igr++){
	Double_t parZ=(maxPosT2-minPosT2)*((igr+0.5)/Double_t((nGraphs)));
	if (iPar>=0){
	  fCanvas->Clear();
	  TGTextEntry * textMin = (TGTextEntry *)fMinEntry->At(iPar);
	  TGTextEntry * textMax = (TGTextEntry *)fMaxEntry->At(iPar);
	  Double_t minPosT=atof(textMin->GetBuffer()->GetString());
	  Double_t maxPosT=atof(textMax->GetBuffer()->GetString());
	  for (Int_t i=0; i<nPoints;i++){
	    vecX[i]=minPosT+i*(maxPosT-minPosT)/nPoints;
	    (fFormulaParams->GetMatrixArray())[iPar]= vecX[i];
	    if (jPar>0) (fFormulaParams->GetMatrixArray())[jPar]=parZ;
	    vecY[i]=drawFormula->EvalPar(0,fFormulaParams->GetMatrixArray());
	  }
	  TGraph *gr = new TGraph(nPoints, vecX.GetMatrixArray(), vecY.GetMatrixArray()); 
	  fCanvas->GetCanvas()->cd();
	  gr->SetLineWidth(2.);
	  gr->SetLineColor(igr+1);
	  gr->SetLineStyle(iForm+1);
	  if (minF>maxF){
	    minF=TMath::MinElement(gr->GetN(),  vecY.GetMatrixArray());
	    maxF=TMath::MaxElement(gr->GetN(),  vecY.GetMatrixArray());
	  }else{
	    minF=TMath::Min(TMath::MinElement(gr->GetN(),  vecY.GetMatrixArray()),minF);
	    maxF=TMath::Max(TMath::MaxElement(gr->GetN(),  vecY.GetMatrixArray()),maxF);	  
	  }
	  if (igr==0&&iForm==0) {
	    gr->Draw("al");
	    gr0=gr;
	  }
	  gr->Draw("l");
	  if (igr==0){
	    legend->AddEntry(gr, formulaStack->At(iForm)->GetName(),"lp");
	  }
	}
      }
      gr0->SetMaximum(maxF+(maxF-minF)*0.3);
      gr0->SetMinimum(minF-(maxF-minF)*0.3);
    }
  }

  Int_t iRow=0;
  for (Int_t iParam=0; iParam<nParams; iParam++){
    if (iParam!=iPar &&iParam!=jPar) {
      fgkLatex->DrawLatexNDC(fgkLatex->GetX(), fgkLatex->GetY()-fgkLatex->GetTextSize()*1.2*iRow, 
			     TString::Format("p_{%d} %s=%.5f", iParam, fFormula->GetParName(iParam),(*fFormulaParams)[iParam]).Data());      
      iRow++;
    }	
    TGTextEntry * textMin = (TGTextEntry *)fMinEntry->At(iParam);
    TGTextEntry * textMax = (TGTextEntry *)fMaxEntry->At(iParam);
    
    if (iParam==iPar) {
      fgkLatex->DrawLatexNDC(fgkLatex->GetX(), fgkLatex->GetY()-fgkLatex->GetTextSize()*1.2*iRow, 
			     TString::Format("x=p_{%d}: %s <%s,%s>", iParam, fFormula->GetParName(iParam),textMin->GetBuffer()->GetString() ,textMax->GetBuffer()->GetString()).Data());  
      gr0->GetXaxis()->SetTitle(fFormula->GetParName(iParam));
      iRow++;
    }
    if (iParam==jPar) {
      fgkLatex->DrawLatexNDC(fgkLatex->GetX(), fgkLatex->GetY()-fgkLatex->GetTextSize()*1.2*iRow, 
			     TString::Format("z=p_{%d}: %s <%s,%s>", iParam, fFormula->GetParName(iParam),textMin->GetBuffer()->GetString(),textMax->GetBuffer()->GetString()).Data());      
      iRow++;      
    }
  }
  if (jPar>0){
    for (Int_t igr=0; igr<nGraphs; igr++){
      fgkLatex->SetTextColor(igr+1);
      Double_t value= minPosT2+(maxPosT2-minPosT2)*igr/Double_t(nGraphs);
      fgkLatex->DrawLatexNDC(fgkLatex->GetX()+0.05, fgkLatex->GetY()-fgkLatex->GetTextSize()*1.2*iRow, 
			     TString::Format("gr[%d]=%.4f", igr, value).Data());      
      iRow++;
    }
    fgkLatex->SetTextColor(1);
  }
  legend->Draw();
  
  gPad->Update();
  gClient->NeedRedraw(fCanvas);
  delete formulaArray; 
}





