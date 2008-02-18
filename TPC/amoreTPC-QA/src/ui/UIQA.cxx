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
#include "AliTPCCalibCE.h"
#include <AmoreDA.h>

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

  //std::ostringstream stringStream;
 //amore::core::String_t sourceName="CEDA/", subscription; // The agent name acting as a source could be concatenated with all the objects it contains
 //subscription=sourceName+"CE";
 //Subscribe(subscription.c_str()); // Here you put a series of subscriptions where the string corresponds to the object name as published in the Publisher Module. As these names are internal to the QA framework, the recommended way of having consistency between AMORE and QA is to factor-out of QA the function that represents the histogram naming convention as a separate AliRoot class/function and use it from inside QA and AMORE.
 //...
 
}

void UIQA::Update() { // This is executed after getting the updated contents of the subscribed MonitorObjects. Notice that the output of moInt[i] and moString[i] varies with time for a specific i because on the dqmAgent the "quality" check fails or succeeds. This is the essence of automatic data quality checks in AMORE. Try to use the moString[i] on a text widget to alert the shifter, or -depending of the value of moInt[i], 0 or 1- make part of the screen change color...
 std::ostringstream stringStream;
 amore::da::AmoreDA amoreDA;
 // Example of accessing a normal TObject. The name is the name of the object in the QA framework
 //amore::core::MonitorObjectTObject* ptr=gSubscriber->At<amore::core::MOTObj>("CEDA/CE");
 //AliTPCCalibCE* hCalibCE=0;
 //if(ptr) {
 //hCalibCE=(AliTPCCalibCE*)ptr->Object();
 //}
 // End of access example

 // hCalibCE->Dump();
 TObject *temp=0;
 amoreDA.Receive("CEDA/CE",temp);
 temp->Dump();


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
