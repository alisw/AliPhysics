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
 
}


UIQA::~UIQA()
{
}

void UIQA::Construct() { // The custom GUI is constructed here. gRootFrame is the container of the custom widgets.

 fTab=new TGTab(amore::ui::gRootFrame);
 amore::ui::gRootFrame->AddFrame(fTab);
 TGCompositeFrame* tempFrame=fTab->AddTab("My Tab 1");
 fEC[0]=new TRootEmbeddedCanvas("fEC0", tempFrame, 800, 450);
 tempFrame->AddFrame(fEC[0]);
 fEC[0]->GetCanvas()->Divide(3, 3);
 tempFrame=fTab->AddTab("My Tab 2");
 amore::ui::gRootFrame->MapSubwindows();
 amore::ui::gRootFrame->Resize();
 amore::ui::gRootFrame->MapWindow();

}

void UIQA::SubscribeMonitorObjects() { // Before using any MonitorObject, a subscription should be made.

 std::ostringstream stringStream;
 amore::core::String_t sourceName="TPCQA", subscription; // The agent name acting as a source could be concatenated with all the objects it contains
 subscription=sourceName+"/moInt1";
 Subscribe(moInt1, subscription.c_str());
 subscription=sourceName+"/hHighPhosModules";
 Subscribe(subscription.c_str()); // Here you put a series of subscriptions where the string corresponds to the object name as published in the Publisher Module. As these names are internal to the QA framework, the recommended way of having consistency between AMORE and QA is to factor-out of QA the function that represents the histogram naming convention as a separate AliRoot class/function and use it from inside QA and AMORE.
 //...
 
}

void UIQA::Update() { // This is executed after getting the updated contents of the subscribed MonitorObjects. Notice that the output of moInt[i] and moString[i] varies with time for a specific i because on the dqmAgent the "quality" check fails or succeeds. This is the essence of automatic data quality checks in AMORE. Try to use the moString[i] on a text widget to alert the shifter, or -depending of the value of moInt[i], 0 or 1- make part of the screen change color...
 std::ostringstream stringStream;
 
 // Example of accessing a normal TObject. The name is the name of the object in the QA framework
 amore::core::MonitorObjectTObject* ptr=gSubscriber->At<amore::core::MOTObj>("TPCQA/hHighPhosModules");
 TH1F* hHighPhosModules=0;
 if(ptr) {
  hHighPhosModules=(TH1F*)ptr->Object();
 }
 // End of access example
 
 fEC[0]->GetCanvas()->cd(1);
 if(hHighPhosModules) hHighPhosModules->Draw();
 for(size_t i=1; i<=4; ++i) {
  //Fill the other pads similarly...
 }
 fEC[0]->GetCanvas()->Update();

}

void UIQA::Process() {

}

void UIQA::StartOfCycle() {

}

void UIQA::EndOfCycle() {

}

};

};

};
