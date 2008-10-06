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
#include <AmoreDA.h>
#include <AliLog.h>

ClassImp(amore::TPC::publisher::PublisherQA)

#include <AliTPCQADataMakerRec.h>
#include <TObjArray.h>
#include <AliRawReaderDate.h>

namespace amore {

namespace TPC {

namespace publisher {

using amore::publisher::Publish;

PublisherQA::PublisherQA() : fqadm(new AliTPCQADataMakerRec) {
  // Constructor
  // make instance of the TPCQADataMakerRec
  //

  AliLog::SetClassDebugLevel("AliTPCRawStream",-5);
  AliLog::SetClassDebugLevel("AliRawReaderDate",-5);
  AliLog::SetClassDebugLevel("AliTPCAltroMapping",-5);
  AliLog::SetModuleDebugLevel("RAW",-5);
  printf("PublisherQA::Constructor\n");

}


PublisherQA::~PublisherQA() {

}

void PublisherQA::BookMonitorObjects() {
  //
  //  Called once at the beginning  after invocation of amoreAgent
  //  all Data which should be accessible to clents has to be published here
  //

 printf("PublisherQA::BookMonitorObject\n");  
 Publish(moInt1, "moInt1", "My Integer MonitorObject 1");
 int cycleLength=0;
 fqadmList=fqadm->Init(AliQA::kRAWS, cycleLength);
 TObjArrayIter* lIt=(TObjArrayIter*)fqadmList->MakeIterator();
 TNamed* obj;
 while((obj=(TNamed*)lIt->Next())) {
   obj->Dump();
   printf("%s\n",obj->GetName());
   Publish(obj, obj->GetName());
 }
}

void PublisherQA::StartOfCycle() {
 
 *moInt1=0;
 fqadm->StartOfCycle(AliQA::kRAWS);
 
}

void PublisherQA::EndOfCycle() {
  printf("PublisherQA::EndOfCycle\n");
 
 fqadm->EndOfCycle(AliQA::kRAWS);
 std::cerr << *moInt1 << std::endl;

}

void PublisherQA::MonitorEvent(amore::core::Event& event) {
   printf("MonitorEvent\n");

 AliRawReaderDate* arr=new AliRawReaderDate(event.DATEEvent());
 fqadm->Exec(AliQA::kRAWS, arr);
 ++*moInt1;
 delete arr;

}

void PublisherQA::StartOfRun(const amore::core::Run& /*run*/) {
  printf("PublisherQA::StartOfRun\n");

}

 void PublisherQA::EndOfRun(const amore::core::Run& /*run*/) {
  printf("PublisherQA::EndOfRun\n");

 fqadm->Finish();
 
}

};

};

};
