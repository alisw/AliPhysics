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
#ifndef AMORE_TPC_PUBLISHERMYMODULE2_H
#define AMORE_TPC_PUBLISHERMYMODULE2_H

#include <PublisherModule.h>
#include "../common/Common.h"

class AliTPCQADataMakerRec;
class TObjArray;

namespace amore {

namespace TPC {

namespace publisher {

/**
@author Filimon Roukoutakis
*/

class PublisherQA : public amore::publisher::PublisherModule, public amore::TPC::common::Common {

 public:
 
 PublisherQA();
 ~PublisherQA();
 virtual void BookMonitorObjects();
 virtual void MonitorEvent(amore::core::Event&);
 virtual void StartOfCycle();
 virtual void EndOfCycle();
 virtual void StartOfRun(const amore::core::Run& run);
 virtual void EndOfRun(const amore::core::Run& run);
 virtual void StartOfSession(const amore::core::Session& session) {};
 virtual void EndOfSession(const amore::core::Session& session) {};
 
 protected:
 
 AliTPCQADataMakerRec* fqadm;
 TObjArray* fqadmList;
 
 ClassDef(PublisherQA, 1);

};

};

};

};

#endif
