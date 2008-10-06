/***************************************************************************
 *   Copyright (C) 2007 by Filimon Roukoutakis                             *
 *   Filimon.Roukoutakis@cern.ch                                           *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "UIQA.h"
#include <AmoreDA.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TList.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TCollection.h>
#include <TGButton.h>
#include "AliTPCCalPad.h"
#include "AliTPCdataQA.h"
#include "AliTPCCalibCE.h"
#include "AliTPCCalibPulser.h"
#include "AliTPCCalibPedestal.h"
#include "AliTPCCalibViewer.h"
#include "AliTPCCalibViewerGUI.h"
#include "AliTPCPreprocessorOnline.h"
#include "AliCDBEntry.h"
ClassImp(amore::TPC::ui::UIQA)

#include <iostream>
#include <sstream>
#include <TCanvas.h>

namespace amore {

namespace TPC {

namespace ui {

using amore::subscriber::Subscribe;




UIQA::UIQA()  :
  fMapCalibObjects(0x0),
  fListCalibObjInfo(0x0), 
  fListGuiObjects(0x0),
  fAmoreDA(new amore::da::AmoreDA(amore::da::AmoreDA::kReceiver))
  {

 Construct(); // Temporary but important!!! Do not forget to put this call in the constructor for the time being!
 fCycle=0; 
}


UIQA::~UIQA()
{
  if ( fAmoreDA )       delete fAmoreDA;
  if ( fMapCalibObjects )delete fMapCalibObjects;
  if ( fListCalibObjInfo ) delete fListCalibObjInfo;
  if ( fListGuiObjects ) delete fListGuiObjects;
  printf("DA destructor called\n");
}

void UIQA::Construct() { // The custom GUI is constructed here. gRootFrame is the container of the custom widgets.
  
 fTab=new TGTab(amore::ui::gRootFrame);
 amore::ui::gRootFrame->AddFrame(fTab);
 //
 //
 TGCompositeFrame* tabCont1 =fTab->AddTab("Expert");
 fExpert = tabCont1;
 
 fViewerGUI = new AliTPCCalibViewerGUI(tabCont1, 1000, 600, 0);
 tabCont1->AddFrame(fViewerGUI , new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
 //
  TGCompositeFrame* tempFrame=fTab->AddTab("OverThreshold");
 fEC[1]=new TRootEmbeddedCanvas("fEC0", tempFrame, 1000, 650);
 tempFrame->AddFrame(fEC[1]);
 fEC[1]->GetCanvas()->Divide(3,3);
 //
 //
 tempFrame=fTab->AddTab("Charge");
 fEC[2]=new TRootEmbeddedCanvas("fEC1", tempFrame, 1000, 650);
 tempFrame->AddFrame(fEC[2]);
 fEC[2]->GetCanvas()->Divide(2,3);
 //
 //
 TGCompositeFrame* tabCalib =fTab->AddTab("DA Calib.");
 SetupTabDACalib(tabCalib);
 //
 //
 amore::ui::gRootFrame->MapSubwindows();
 amore::ui::gRootFrame->Resize();
 amore::ui::gRootFrame->MapWindow();

 gROOT->SetStyle("Plain");
 gStyle->SetFillColor(10);
 gStyle->SetPadColor(10);
 gStyle->SetCanvasColor(10);
 gStyle->SetStatColor(10);
 
 gStyle->SetPalette(1,0);
 gStyle->SetNumberContours(30);
 gStyle->SetOptFit(111);
 
 gStyle->SetCanvasBorderMode(-1);
 gStyle->SetCanvasBorderSize(1);
 gStyle->SetCanvasColor(10);
 
 gStyle->SetFrameFillColor(10);
 gStyle->SetFrameBorderSize(1);
 gStyle->SetFrameBorderMode(-1);
 gStyle->SetFrameLineWidth(1);
 SubscribeMonitorObjects();
 //Perhaps move later to a button action
 RetrieveFromAmoreDB();
}

void UIQA::SetupTabDACalib(TGCompositeFrame *frame){
  if ( !fListGuiObjects ){
    fListGuiObjects = new TList;
    fListGuiObjects->SetOwner();
  }
  // Update Button
  TGTextButton *btnUpdate = new TGTextButton(frame, "&Update");
  frame->AddFrame(btnUpdate, new TGLayoutHints(kLHintsExpandX, 10, 10, 0, 0));
   //fBtnDraw->Connect("Clicked()", "AliTPCCalibViewerGUI", this, "DoTest(=\"fBtnDraw clicked\")");
  btnUpdate->Connect("Clicked()", "UIQA", this, "UpdateAmoreDBValues()");
  btnUpdate->SetToolTipText("Update values from the .");
  fListGuiObjects->Add(btnUpdate);
  //Update info boxes 
/*  TGGroupFrame *ldcPedestalGroup = new TGGroupFrame(ftabLeft0, "Pedestal/Noise Info", kVerticalFrame | kFitWidth | kFitHeight);
  frame->AddFrame(ldcPedestalGroup, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
  for (Int_t ildc=0;ildc<36;++ildc){
    char side='A';
    if (isec>17) side='C';
    TString amoreDAname(Form("ldc-TPC-%c%02d-%s",side,isec,"pedestals"));
    
    
  }*/
}

void UIQA::UpdateAmoreDBValues(){
  RetrieveFromAmoreDB();
}

void UIQA::SubscribeMonitorObjects() { // Before using any MonitorObject, a subscription should be made.

  //std::ostringstream stringStream;
 //amore::core::String_t sourceName="CEDA/", subscription; // The agent name acting as a source could be concatenated with all the objects it contains
 //subscription=sourceName+"CE";
 //Subscribe(subscription.c_str()); // Here you put a series of subscriptions where the string corresponds to the object name as published in the Publisher Module. As these names are internal to the QA framework, the recommended way of having consistency between AMORE and QA is to factor-out of QA the function that represents the histogram naming convention as a separate AliRoot class/function and use it from inside QA and AMORE.
 //...
  printf("UIQA::SubscribeMonitorObjects\n");
  amore::core::String_t sourceName="TPCQA/", subscription; 
  subscription=sourceName+"TPCRAW";
  Subscribe(subscription.c_str());
  printf("%s\n",subscription.c_str());
  //subscription=sourceName+"hist";
  //Subscribe(subscription.c_str());
  
}

void UIQA::Update() { 
  // This is executed after getting the updated contents of the subscribed MonitorObjects. Notice that the output of moInt[i] and moString[i] varies with time for a specific i because on the dqmAgent the "quality" check fails or succeeds. This is the essence of automatic data quality checks in AMORE. Try to use the moString[i] on a text widget to alert the shifter, or -depending of the value of moInt[i], 0 or 1- make part of the screen change color...
 std::ostringstream stringStream;
 // Example of accessing a normal TObject. The name is the name of the object in the QA framework


 static Int_t counter0 = 0, counter1=0;
 printf("UIQA::Update0\t%d\t%d\n",counter0,counter1); 

 amore::core::MonitorObjectTObject* ptr=0;
 AliTPCdataQA *tpcqa=0; 
 ptr=gSubscriber->At<amore::core::MOTObj>("TPCQA/TPCRAW");
 tpcqa=0;
 printf("Pointer - %p\n",ptr);
 if(ptr) {
   tpcqa=(AliTPCdataQA*)ptr->Object();
   printf("Pointertpcqa - %p\n",tpcqa);
   tpcqa->Print();
 }
 
 printf("Update1\t%d\t%d\n",counter0,counter1); 
 counter0++;
 if (!tpcqa) {
   MakeTree(0x0);
   return;
 }
 printf("Update2\t%d",counter1);
 counter1++;
 //fExpert->SetTitle(name);
 //
 // Over threshold
 //

 TCanvas *canvas  = fEC[1]->GetCanvas();
 if (tpcqa->GetNoThreshold()){
   canvas->cd(1);
   tpcqa->GetNoThreshold()->MakeHisto1D()->Draw();
   canvas->cd(2);
   tpcqa->GetNoThreshold()->MakeHisto2D(0)->Draw("colz");
   canvas->cd(3);
   tpcqa->GetNoThreshold()->MakeHisto2D(1)->Draw("colz");
 }
 //
 if (tpcqa->GetNTimeBins()){
   canvas->cd(4);
   tpcqa->GetNTimeBins()->MakeHisto1D()->Draw();
   canvas->cd(5);
   tpcqa->GetNTimeBins()->MakeHisto2D(0)->Draw("colz");
   canvas->cd(6);
   tpcqa->GetNTimeBins()->MakeHisto2D(1)->Draw("colz");
 }

 if (tpcqa->GetNPads()){
   canvas->cd(7);
   tpcqa->GetNPads()->MakeHisto1D()->Draw();
   canvas->cd(8);
   tpcqa->GetNPads()->MakeHisto2D(0)->Draw("colz");
   canvas->cd(9);
   tpcqa->GetNPads()->MakeHisto2D(1)->Draw("colz");
 }

 //
 // Mean charge
 //
 canvas  = fEC[2]->GetCanvas();
 if (tpcqa->GetMeanCharge()){
   canvas->cd(1);
   tpcqa->GetMeanCharge()->MakeHisto1D()->Draw();
   canvas->cd(3);
   tpcqa->GetMeanCharge()->MakeHisto2D(0)->Draw("colz");
   canvas->cd(5);
   tpcqa->GetMeanCharge()->MakeHisto2D(1)->Draw("colz");
   //
   canvas->cd(2);
   tpcqa->GetMaxCharge()->MakeHisto1D()->Draw();
   canvas->cd(4);
   tpcqa->GetMaxCharge()->MakeHisto2D(0)->Draw("colz");
   canvas->cd(6);
   tpcqa->GetMaxCharge()->MakeHisto2D(1)->Draw("colz");
 }
 MakeTree(tpcqa);

 char name[1000];
 sprintf(name,"Summary Event -  %d",tpcqa->GetEventCounter());
 fExpert->SetName(name);
 for (Int_t i=0;i<3;i++) printf("*********************************\n");
 printf("\n\nSummary Event -  %d\n\n",tpcqa->GetEventCounter());
 for (Int_t i=0;i<3;i++) printf("*********************************\n");
 

}

void UIQA::Process() {

}

void UIQA::StartOfCycle() {
}

void UIQA::EndOfCycle() {

}


void UIQA::MakeTree(AliTPCdataQA* ped){
  //
  //
  AliTPCPreprocessorOnline * preprocesor = new AliTPCPreprocessorOnline;

  //==============================
  //    QA Information
  //==============================
  if (ped) {
    if (ped->GetNLocalMaxima()) 
      preprocesor->AddComponent(new AliTPCCalPad(*(ped->GetNLocalMaxima())));
    if (ped->GetMaxCharge()) 
      preprocesor->AddComponent(new AliTPCCalPad(*(ped->GetMaxCharge())));  
    if (ped->GetMeanCharge()) 
      preprocesor->AddComponent(new AliTPCCalPad(*(ped->GetMeanCharge())));  
    if (ped->GetNoThreshold()) 
      preprocesor->AddComponent(new AliTPCCalPad(*(ped->GetNoThreshold())));
    if (ped->GetNTimeBins()) 
      preprocesor->AddComponent(new AliTPCCalPad(*(ped->GetNTimeBins())));
    if (ped->GetNPads()) 
      preprocesor->AddComponent(new AliTPCCalPad(*(ped->GetNPads())));
    if (ped->GetTimePosition()) 
      preprocesor->AddComponent(new AliTPCCalPad(*(ped->GetTimePosition())));   
  }
  //==============================
  //    Calibration Information
  //==============================
  TIter next(fMapCalibObjects);
  while ( TObject *o=next() ) preprocesor->AddComponent(o->Clone());

  char fname[10000];
  sprintf(fname,"QAtree%d.root",fCycle);
  preprocesor->DumpToFile(fname);
  fCycle++;
    //init viewer
  AliTPCCalibViewer *viewer = fViewerGUI->GetViewer();
  AliTPCCalibViewer *nviewer = new  AliTPCCalibViewer(fname, "calPads");
  fViewerGUI->Initialize(nviewer);
  delete preprocesor;
}

void UIQA::RetrieveFromAmoreDB(){
  if ( fMapCalibObjects ) {
    fMapCalibObjects->Delete();
    delete fMapCalibObjects;
  }
  if ( !fListCalibObjInfo ){
    fListCalibObjInfo = new TList;
    fListCalibObjInfo->SetOwner();
  }
  fMapCalibObjects = new TList;
  fMapCalibObjects->SetOwner();
  //define calibration objects, type and role
  TString calNames("Pedestals;Noise;PulserT0;PulserQ;PulserRMS;CET0;CEQ;CERMS");
  TString daTypes("pedestals;pedestals;pulser;pulser;pulser;CE;CE;CE");
  TString daRoles("ldc;ldc;ldc;ldc;ldc;mon-DA-TPC-0;mon-DA-TPC-0;mon-DA-TPC-0");
  TObjArray *calNameArr=calNames.Tokenize(";");
  TObjArray *daTypeArr=daTypes.Tokenize(";");
  TObjArray *daRoleArr=daRoles.Tokenize(";");
  Int_t nobj=calNameArr->GetEntries();
  
  for (Int_t ient=0; ient<nobj; ++ient){
    const TString &calName = ((TObjString*)calNameArr->At(ient))->GetString();  
    const TString &daType  = ((TObjString*)daTypeArr->At(ient))->GetString();
    const TString &daRole  = ((TObjString*)daRoleArr->At(ient))->GetString();
    AliTPCCalPad *calPad  = new AliTPCCalPad(calName.Data(),calName.Data());
    fMapCalibObjects->Add(calPad);
    if ( daRole.Contains("ldc") ) CollectFromLDCs(calPad,calName,daType);
    else if ( daRole.Contains("mon") ) CollectFromMon(calPad,calName,daType,daRole);
  }
  delete calNameArr;
  delete daTypeArr;
  delete daRoleArr;
}

void UIQA::CollectFromLDCs(AliTPCCalPad *calPad, const TString &calName, const TString &daType){
  for (int isec=0; isec<36; isec++){
    TObject *obj=0x0;
    Int_t res=0;
    char side='A';
    if (isec>17) side='C';
    TString amoreDAname(Form("ldc-TPC-%c%02d-%s",side,isec,daType.Data()));
    TString objName(Form("%s/%s",amoreDAname.Data(),calName.Data()));   
    res=fAmoreDA->Receive(objName.Data(),obj);
    if ( !res&&obj ){
      printf ("Found Pedestal Object in '%s'\n",objName.Data());
      TObjArray *arr=(TObjArray*)obj;
       //IROC
      AliTPCCalROC *rocMerge=(AliTPCCalROC*)arr->At(isec);
      if (rocMerge) calPad->SetCalROC(rocMerge,isec);
        //OROC
      rocMerge=(AliTPCCalROC*)arr->At(isec+36);
      if (rocMerge) calPad->SetCalROC(rocMerge,isec+36);
    }else{
      printf("Could not receive '%s'\n",objName.Data());
    }      
    // Get info message
    TString infoName(Form("%s/%s",amoreDAname.Data(),"Info"));   
    res=fAmoreDA->Receive(infoName.Data(),obj);
    if ( !res&&obj ){
      TObjString *ostr=(TObjString*)obj;
      printf("Calib. object DA info from '%s': '%s'\n",objName.Data(),ostr->GetName());
      //save to info list
      TObjString *infoStr = (TObjString*)fListCalibObjInfo->FindObject(amoreDAname);
      if ( !infoStr ) {
        infoStr = new TObjString;
        fListCalibObjInfo->Add(infoStr);
      }
      infoStr->GetString()=ostr->GetString();
    }
    
  }    
}

void UIQA::CollectFromMon(AliTPCCalPad *calPad, const TString &calName, const TString &daType, const TString &mon){
  TObject *obj=0x0;
  Int_t res=0;
  TString amoreDAname(Form("%s-%s",mon.Data(),daType.Data()));
  TString objName(Form("%s/%s",amoreDAname.Data(),calName.Data()));   
  res=fAmoreDA->Receive(objName.Data(),obj);
  if ( !res&&obj ){
    printf ("Found Pedestal Object in '%s'\n",objName.Data());
    TObjArray *arr=(TObjArray*)obj;
    for (Int_t isec=0; isec<72; isec++){
      AliTPCCalROC *rocMerge=(AliTPCCalROC*)arr->At(isec);
      if (rocMerge) calPad->SetCalROC(rocMerge,isec);     
    }
  }else{
    printf("Could not receive '%s'\n",objName.Data());
  }  
    // Get info message
  TString infoName(Form("%s/%s",amoreDAname.Data(),"Info"));   
  res=fAmoreDA->Receive(infoName.Data(),obj);
  if ( !res&&obj ){
    TObjString *ostr=(TObjString*)obj;
    printf("Calib. object DA info from '%s': '%s'\n",objName.Data(),ostr->GetName());
    TObjString *infoStr = (TObjString*)fListCalibObjInfo->FindObject(amoreDAname);
    if ( !infoStr ) {
      infoStr = new TObjString;
      fListCalibObjInfo->Add(infoStr);
    }
    infoStr->GetString()=ostr->GetString();
  }
}

};

};

};
