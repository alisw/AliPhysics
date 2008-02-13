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
#include "PublisherQA.h"

ClassImp(amore::TPC::publisher::PublisherQA)

#include <AliPHOSQADataMakerRec.h>
#include <TObjArray.h>
#include <AliRawReaderDate.h>

namespace amore {

namespace TPC {

namespace publisher {

using amore::publisher::Publish;

PublisherQA::PublisherQA() : fqadm(new AliPHOSQADataMakerRec) {

}


PublisherQA::~PublisherQA() {

}

void PublisherQA::BookMonitorObjects() {
    
 Publish(moInt1, "moInt1", "My Integer MonitorObject 1");
 int cycleLength=0;
 fqadmList=fqadm->Init(AliQA::kRAWS, 0, cycleLength);
 TObjArrayIter* lIt=(TObjArrayIter*)fqadmList->MakeIterator();
 TNamed* obj;
 while((obj=(TNamed*)lIt->Next())) Publish(obj, obj->GetName());

}

void PublisherQA::StartOfCycle() {
 
 *moInt1=0;
 fqadm->StartOfCycle(AliQA::kRAWS);

}

void PublisherQA::EndOfCycle() {
 
 fqadm->EndOfCycle(AliQA::kRAWS);
 std::cerr << *moInt1 << std::endl;

}

void PublisherQA::MonitorEvent(amore::core::Event& event) {
 
 AliRawReaderDate* arr=new AliRawReaderDate(event.DATEEvent());
 fqadm->Exec(AliQA::kRAWS, arr);
 ++*moInt1;
 delete arr;

}

void PublisherQA::StartOfRun() {

}

void PublisherQA::EndOfRun() {

 fqadm->Finish();
 
}

};

};

};
