// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   testAliHLTDataBuffer.C
    @author Matthias Richter
    @date   
    @brief  Test program for the AliHLTDataBuffer class
 */

#ifndef __CINT__
#include "TDatime.h"
#include "TRandom.h"
#include "AliHLTDataTypes.h"
#include "algorithm"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "AliHLTDAQ.h"
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cerrno>
#include "AliHLTDataBuffer.h"
#endif

using namespace std;

int gVerbosity=0;
const int gPagesize=1024*1024*10;

int testAliHLTRawPage();
int GetRandom(int min, int max);

struct testProcessDescription {
  int fLevel;
  AliHLTUInt32_t fBufferSize;
  AliHLTDataBuffer::AliHLTRawBuffer* fRawBuffer;
};

int testAliHLTDataBuffer()
{
  int iResult=0;
  if ((iResult=testAliHLTRawPage())<0) return iResult;

  return 0;
}

int fillTestSample(int levels, int processes, vector<testProcessDescription>& descriptions)
{
  int availableProcesses=processes;
  for (int level=0; level<levels; level++) {
    int levelprocesses=GetRandom(availableProcesses/2, availableProcesses-(levels-level));
    for (int process=0; process<levelprocesses; process++) {
      testProcessDescription desc;
      desc.fLevel=level;
      desc.fBufferSize=GetRandom(1024, 3*gPagesize/processes);
      desc.fRawBuffer=NULL;
      descriptions.push_back(desc);
    }
    availableProcesses-=levelprocesses;
  }
  return 0;
}

int allocateBuffers(vector<testProcessDescription>& descriptions,
		    vector<AliHLTDataBuffer::AliHLTRawPage*>& pages,
		    int level=-1)
{
  for (unsigned i=0; i<descriptions.size(); i++) {
    if (level>=0 && descriptions[i].fLevel<level) continue;
    if (level>=0 && descriptions[i].fLevel>level) break;
    if (descriptions[i].fRawBuffer) {
      cerr << "warning: buffer already allocated" << endl;
      continue;
    }

    vector<AliHLTDataBuffer::AliHLTRawPage*>::iterator page=pages.begin();
    for (page=pages.begin();page!=pages.end(); page++) {
      if ((descriptions[i].fRawBuffer=(*page)->Alloc(GetRandom(descriptions[i].fBufferSize/2, descriptions[i].fBufferSize)))!=NULL) {
	if (gVerbosity>1) {
	  printf("allocated raw buffer %p from page %p\n", descriptions[i].fRawBuffer, *page);
	  descriptions[i].fRawBuffer->Print("min");
	}
	break;
      }
    }
    if (!descriptions[i].fRawBuffer) {
      AliHLTDataBuffer::AliHLTRawPage* rawpage=new AliHLTDataBuffer::AliHLTRawPage(gPagesize);
      if (!rawpage) {
	cerr << "can not create raw page" << endl;
	return -ENOMEM;
      }
      pages.push_back(rawpage);
      if ((descriptions[i].fRawBuffer=rawpage->Alloc(GetRandom(descriptions[i].fBufferSize/2, descriptions[i].fBufferSize)))!=NULL) {
	if (gVerbosity>1) {
	  printf("allocated raw buffer %p from page %p\n", descriptions[i].fRawBuffer, rawpage);
	  descriptions[i].fRawBuffer->Print("min");
	}
      }
    }
    if (!descriptions[i].fRawBuffer) {
      cerr << "failed to allocate buffer for process " << i << endl;
      return -EFAULT;
    }
  }
  return 0;
}
		   
int setBufferSizes(vector<testProcessDescription>& descriptions,
		   vector<AliHLTDataBuffer::AliHLTRawPage*>& pages,
		   int level=-1)
{
  // set buffer size for all processes of the specified level
  // buffer size is chosen randomly between 0 and the allocated size
  vector<unsigned> positions;
  for (unsigned i=0; i<descriptions.size(); i++) {
    if (!descriptions[i].fRawBuffer) continue;
    if (level>=0 && descriptions[i].fLevel<level) continue;
    if (level>=0 && descriptions[i].fLevel>level) break;
    positions.push_back(i);
  }

  random_shuffle(positions.begin(), positions.end());
  for (vector<unsigned>::iterator position=positions.begin();
       position!=positions.end(); position++) {
    vector<AliHLTDataBuffer::AliHLTRawPage*>::iterator page=pages.begin();
    for (; page!=pages.end(); page++) {
      if ((*page)->SetSize(descriptions[*position].fRawBuffer, GetRandom(0, descriptions[*position].fBufferSize))==0) {
	if (gVerbosity>1) {
	  cout << "setting size for raw buffer " << descriptions[*position].fRawBuffer << " of page " << *page << ", process level " << descriptions[*position].fLevel << endl;
	  descriptions[*position].fRawBuffer->Print("min");
	}
	break;
      }
    }
    if (page==pages.end()) {
      cerr << "failed to set size for raw buffer " << descriptions[*position].fRawBuffer << "of process " << *position << endl;
      return -EFAULT;
    }
  }
  return 0;
}
		   
int releaseBuffers(vector<testProcessDescription>& descriptions,
		   vector<AliHLTDataBuffer::AliHLTRawPage*>& pages,
		   int level=-1)
{
  // find the processes to be releases according to the specified level
  // shuffle the processes and then release in this random order
  vector<unsigned> positions;
  for (unsigned i=0; i<descriptions.size(); i++) {
    if (!descriptions[i].fRawBuffer) continue;
    if (level>=0 && descriptions[i].fLevel<level) continue;
    if (level>=0 && descriptions[i].fLevel>level) break;
    positions.push_back(i);
  }

  random_shuffle(positions.begin(), positions.end());
  for (vector<unsigned>::iterator position=positions.begin();
       position!=positions.end(); position++) {
    vector<AliHLTDataBuffer::AliHLTRawPage*>::iterator page=pages.begin();
    for (; page!=pages.end(); page++) {
      if ((*page)->Free(descriptions[*position].fRawBuffer)==0) {
	if (gVerbosity>1) cout << "released raw buffer " << descriptions[*position].fRawBuffer << " from page " << *page << ", process level " << descriptions[*position].fLevel << endl;
	descriptions[*position].fRawBuffer=NULL;
	break;
      }
    }
    if (page==pages.end()) {
      cerr << "failed to release raw buffer " << descriptions[*position].fRawBuffer << "of process " << *position << endl;
      return -EFAULT;
    }
  }
  return 0;
}
		   
int testAliHLTRawPage()
{
  int iResult=0;
  int nofLevels=GetRandom(3,5);
  int nofProcesses=GetRandom(nofLevels,10);
  int nofCycles=GetRandom(5,10);

  vector<testProcessDescription> descriptions;
  if ((iResult=fillTestSample(nofLevels, nofProcesses, descriptions))<0) {
    cerr << "failed to fill test sample" << endl;
    return iResult;
  }

  if (gVerbosity>1) {
    for (unsigned i=0; i<descriptions.size(); i++) {
      cout << "process " << i << ": level " << descriptions[i].fLevel << "   max size " << descriptions[i].fBufferSize << endl;
    }
  }

  vector<AliHLTDataBuffer::AliHLTRawPage*> pages;
  vector<AliHLTDataBuffer::AliHLTRawPage*>::iterator page=pages.end();
  for (int cycle=0; cycle<nofCycles; cycle++) {
    unsigned process=0;
    for (int level=0; level<=nofLevels; level++) {
      // allocate buffers
      if (level<nofLevels) {
	if ((iResult=allocateBuffers(descriptions, pages, level))<0) {
	  cerr << "failed to allocate buffers" << endl;
	  return iResult;
	}
      }

      if (gVerbosity>0) cout << "level " << level << " - status of pages: " << pages.size() << endl;
      for (page=pages.begin();page!=pages.end(); page++) {
	if (gVerbosity>0) (*page)->Print("");
      }

      if (level<nofLevels) {
	if ((iResult=setBufferSizes(descriptions, pages, level))<0) {
	  cerr << "failed setting  buffers" << endl;
	  return iResult;
	}
      }

      if (gVerbosity>0) cout << "level " << level << " - status of pages: " << pages.size() << endl;
      for (page=pages.begin();page!=pages.end(); page++) {
	if (gVerbosity>0) (*page)->Print("");
      }

      if (level>0) {
	if ((iResult=releaseBuffers(descriptions, pages, level-1))<0) {
	  cerr << "failed to release buffers" << endl;
	  return iResult;
	}
      }
    }

    if (gVerbosity>0) 
      cout << "status of released pages:" << endl;
    for (page=pages.begin();page!=pages.end(); page++) {
      if (gVerbosity>0) (*page)->Print("");
      if ((*page)->IsUsed()) {
	cerr << "page " << *page << " has used buffers" << endl;
	return -EFAULT;
      }
      if ((*page)->Size()!=(*page)->Capacity()) {
	cerr << "page " << *page << " not completely released" << endl;
	return -EFAULT;
      }
      if ((*page)->IsFragmented()) {
	cerr << "page " << *page << " is still fragmented" << endl;
	return -EFAULT;
      }
    }
  }

  return 0;
}

/**
 * Get a random number in the given range.
 */
int GetRandom(int min, int max)
{
  if (max-min<2) return min;
  static TRandom rand;
  static bool seedSet=false;
  if (!seedSet) {
    TDatime dt;
    rand.SetSeed(dt.Get());
    seedSet=true;
  }
  return min+rand.Integer(max-min);
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// main functions

int main(int /*argc*/, const char** /*argv*/)
{
  int iResult=0;
  iResult=testAliHLTDataBuffer();
  return iResult;
}
