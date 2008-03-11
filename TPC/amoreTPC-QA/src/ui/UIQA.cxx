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
#include "AliTPCCalPad.h"
#include "AliTPCdataQA.h"
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




UIQA::UIQA() {

 Construct(); // Temporary but important!!! Do not forget to put this call in the constructor for the time being!
 fCycle=0; 
}


UIQA::~UIQA()
{
}

void UIQA::Construct() { // The custom GUI is constructed here. gRootFrame is the container of the custom widgets.
  
 fTab=new TGTab(amore::ui::gRootFrame);
 amore::ui::gRootFrame->AddFrame(fTab);
 //
 //
 TGCompositeFrame* tabCont1 =fTab->AddTab("Expert");
 
 fViewerGUI = new AliTPCCalibViewerGUI(tabCont1, 1000, 600, 0);
 tabCont1->AddFrame(fViewerGUI , new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
 //
  TGCompositeFrame* tempFrame=fTab->AddTab("OverThreshold");
 fEC[0]=new TRootEmbeddedCanvas("fEC0", tempFrame, 1000, 650);
 tempFrame->AddFrame(fEC[0]);
 fEC[0]->GetCanvas()->Divide(3,3);
 //
 //
 tempFrame=fTab->AddTab("Charge");
 fEC[1]=new TRootEmbeddedCanvas("fEC1", tempFrame, 1000, 650);
 tempFrame->AddFrame(fEC[1]);
 fEC[1]->GetCanvas()->Divide(2,3);
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
}

void UIQA::SubscribeMonitorObjects() { // Before using any MonitorObject, a subscription should be made.

  //std::ostringstream stringStream;
 //amore::core::String_t sourceName="CEDA/", subscription; // The agent name acting as a source could be concatenated with all the objects it contains
 //subscription=sourceName+"CE";
 //Subscribe(subscription.c_str()); // Here you put a series of subscriptions where the string corresponds to the object name as published in the Publisher Module. As these names are internal to the QA framework, the recommended way of having consistency between AMORE and QA is to factor-out of QA the function that represents the histogram naming convention as a separate AliRoot class/function and use it from inside QA and AMORE.
 //...
  amore::core::String_t sourceName="TPCQA/", subscription; 
  subscription=sourceName+"TPCRAW";
  Subscribe(subscription.c_str());
  subscription=sourceName+"hist";
  Subscribe(subscription.c_str());
  
}

void UIQA::Update() { // This is executed after getting the updated contents of the subscribed MonitorObjects. Notice that the output of moInt[i] and moString[i] varies with time for a specific i because on the dqmAgent the "quality" check fails or succeeds. This is the essence of automatic data quality checks in AMORE. Try to use the moString[i] on a text widget to alert the shifter, or -depending of the value of moInt[i], 0 or 1- make part of the screen change color...
 std::ostringstream stringStream;
 // Example of accessing a normal TObject. The name is the name of the object in the QA framework


 //amore::core::MonitorObjectTObject* ptr=gSubscriber->At<amore::core::MOTObj>("TPCQA/TPCRAW");
 amore::core::MonitorObjectTObject* ptr=gSubscriber->At<amore::core::MOTObj>("TPCQA/hist");
 AliTPCdataQA *tpcqa=0;
 printf("Pointer - %p\n",ptr);
 if(ptr) {
   ptr->Object()->Print();
 }
 
 ptr=gSubscriber->At<amore::core::MOTObj>("TPCQA/TPCRAW");
 tpcqa=0;
 printf("Pointer - %p\n",ptr);
 if(ptr) {
   tpcqa=(AliTPCdataQA*)ptr->Object();
   printf("Pointertpcqa - %p\n",tpcqa);
   tpcqa->Print();
 }

 if (!tpcqa) return;
 //
 // Over threshold
 //

 TCanvas *canvas  = fEC[0]->GetCanvas();
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

 canvas  = fEC[1]->GetCanvas();
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
 if (tpcqa) MakeTree(tpcqa);

 // End of access example

 //amore::da::AmoreDA amoreDA;
 // hCalibCE->Dump();
 // TObject *temp=0;
 //amoreDA.Receive("CEDA/CE",temp);
 //temp->Dump();


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
    AliTPCCalPad * noise = GetNoise();
    if (noise)  preprocesor->AddComponent(new AliTPCCalPad(*noise));
    AliTPCCalPad * pedestal = GetPedestal();
    if (pedestal)  preprocesor->AddComponent(new AliTPCCalPad(*pedestal));
    char fname[10000];
    sprintf(fname,"QAtree%d.root",fCycle);
    preprocesor->DumpToFile(fname);
    fCycle++;
    delete noise;
    delete pedestal;
    //  /*CalibTree
    AliTPCCalibViewer *viewer = fViewerGUI->GetViewer();
    AliTPCCalibViewer *nviewer = new  AliTPCCalibViewer(fname, "calPads");
    fViewerGUI->Initialize(nviewer);
    //*/
    //
    //
    //
  delete preprocesor;
  }

  AliTPCCalPad *  UIQA::GetNoise(){
    //
    // GetNoise - if not in AmoreDB than from file
    //
    //amore::da::AmoreDA amoreDA;
    // TObject *temp=0;
    //amoreDA.Receive("PEDESTAL/NOISE",temp);
    //temp->Dump();
    // if (temp) return (AliTPCCalPad*) temp;
    TDirectory * dir = gDirectory;
    TFile *f = new TFile("$AMORE_SITE/PadNoise.root");
    AliCDBEntry * entry = (AliCDBEntry*)f->Get("AliCDBEntry");
    if (entry){
      AliTPCCalPad * pad = ((AliTPCCalPad *)entry->GetObject()->Clone());
      if (pad) {
	pad->SetName("Noise");
	pad->Print();
      	return pad;
      }
    }
    f->Close();
    dir->cd();
    return 0;
  }

  AliTPCCalPad *  UIQA::GetPedestal(){
    //
    //
    TFile f("$AMORE_SITE/Pedestal.root");
    AliCDBEntry * entry = (AliCDBEntry*)f.Get("AliCDBEntry");
    if (entry){
      AliTPCCalPad * pad = ((AliTPCCalPad *)entry->GetObject()->Clone());
      if (pad) {
	pad->SetName("Pedestals");
	pad->Print();
	return pad;
      }
    }
    return 0;
  }
  AliTPCCalPad *  UIQA::GetTime0(){
    return 0;
  }



};

};

};
