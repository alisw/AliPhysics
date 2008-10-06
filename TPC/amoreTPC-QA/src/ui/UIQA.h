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
#ifndef AMORE_TPC_UIUIQA_H
#define AMORE_TPC_UIUIQA_H

#include <VisualModule.h>
#include <TRootEmbeddedCanvas.h>
#include "../common/Common.h" // Absolute path necessary to allow usage of ACLiC
#include <TGTab.h>
#include <TGNumberEntry.h>
#include <TGTextView.h>
#include <TString.h>
class TList;
class AliTPCCalibViewerGUI;
class AliTPCdataQA;
class AliTPCCalPad;
class AliTPCCalibPedestal;
class AliTPCCalibPulser;
class AliTPCCalibCE;
namespace amore {

namespace da {
  class AmoreDA;
};
  

namespace TPC {

namespace ui {

/**
@author Filimon Roukoutakis
*/

class UIQA : public amore::ui::VisualModule, public amore::TPC::common::Common { // VisualModule inheritance mandatory, Common inheritance optional for sharing the MonitorObject declarations with other modules

 public:
 
 UIQA();
 ~UIQA();
 
 virtual void Construct();
 virtual void Update();
 virtual void SubscribeMonitorObjects();
 virtual void Process();
 virtual void StartOfCycle();
 virtual void EndOfCycle();
 virtual void StartOfRun() {};
 virtual void EndOfRun() {};
 virtual void StartOfSession() {};
 virtual void EndOfSession() {};
 
 protected:
 void MakeTree(AliTPCdataQA * qa);
 void MergeCalPadPedestal();
 void RetrieveFromAmoreDB();
 void CollectFromLDCs(AliTPCCalPad *calPad, const TString &calName, const TString &daType);
 void CollectFromMon(AliTPCCalPad *calPad, const TString &calName, const TString &daType, const TString &mon);
 TList *fMapCalibObjects;
 TList *fListCalibObjInfo;
 amore::da::AmoreDA      *fAmoreDA;
 // gui
 TGTab* fTab;
 TGCompositeFrame* fExpert;
 TRootEmbeddedCanvas* fEC[10];
 TGNumberEntryField* fNEF[10];
 TGTextView* fTextView[10];
 AliTPCCalibViewerGUI *fViewerGUI;
 TList *fListGuiObjects;
 void SetupTabDACalib(TGCompositeFrame *frame);
 void UpdateAmoreDBValues();
 Int_t fCycle;
 ClassDef(UIQA, 1);

};

};

};

};

#endif
