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
const int gPagesize=1024*1024;

template<class T>
int testLoop(vector<T>& descriptions, int nofLevels, int nofCycles);
int GetRandom(int min, int max);

class AliHLTRandomBuffer {
public:
  AliHLTRandomBuffer(): fpBuffer(NULL), fBufferSize(0), fInstances() {}
  AliHLTRandomBuffer(AliHLTUInt32_t size): fpBuffer(NULL), fBufferSize(0), fInstances() {Init(size);}
  AliHLTRandomBuffer(const AliHLTRandomBuffer&): fpBuffer(NULL), fBufferSize(0), fInstances() {}
  AliHLTRandomBuffer& operator=(const AliHLTRandomBuffer& src) {if (&src!=this) {fpBuffer=NULL; fBufferSize=0;} return *this;}
  ~AliHLTRandomBuffer() {}

  int Init(AliHLTUInt32_t size) {
    if (fpBuffer) return -EINPROGRESS;
    
    fpBuffer=new AliHLTUInt8_t[size];
    if (!fpBuffer) return -ENOMEM;
    fBufferSize=size;

    int words=size/sizeof(AliHLTUInt32_t);
    AliHLTUInt32_t* pTgt=reinterpret_cast<AliHLTUInt32_t*>(fpBuffer);
    for (int i=0; i<words; i++, pTgt++) *pTgt=GetRandom(0, 0x7fffffff);
    return 0;
  }

  struct AliHLTRandomBufferInstance {
    int fId;
    AliHLTUInt32_t fOffset;
    AliHLTUInt32_t fSize;
  };

  int CreateInstance(AliHLTUInt32_t size) {
    if (size>fBufferSize) return -ENOSPC;

    AliHLTRandomBufferInstance instance;
    instance.fId=fInstances.size()>0?fInstances.back().fId+1:0;
    instance.fOffset=GetRandom(0, fBufferSize-size);
    instance.fSize=size;
    fInstances.push_back(instance);
    return instance.fId;
  }

  int ReleaseInstance(int id) {
    for (vector<AliHLTRandomBufferInstance>::iterator instance=fInstances.begin();
	 instance!=fInstances.end(); instance++) {
      if (instance->fId==id) {
	fInstances.erase(instance);
	return 0;
      }
    }
    return -ENOENT;
  }

  const AliHLTUInt8_t* GetBuffer(int id) {
    for (vector<AliHLTRandomBufferInstance>::iterator instance=fInstances.begin();
	 instance!=fInstances.end(); instance++) {
      if (instance->fId==id) {

	return fpBuffer+instance->fOffset;
      }
    }
    return NULL;
  }

  bool CompareBuffer(const AliHLTUInt8_t* buffer, AliHLTUInt32_t size, int id) {
    if (!buffer) return false;
    for (vector<AliHLTRandomBufferInstance>::iterator instance=fInstances.begin();
	 instance!=fInstances.end(); instance++) {
      if (instance->fId==id) {
	if (instance->fSize!=size) {
	  cout << "CompareBuffer ("<< instance->fId << "): size missmatch: " << size << ", expected " << instance->fSize << endl; 
	  return false;
	}
	return memcmp(buffer, fpBuffer+instance->fOffset, size)==0;
      }
    }
    return false;
  }

  void Print(const char* /*option*/) {
    cout << "AliHLTRandomBufferInstance " << this << ": buffer " << (void*)fpBuffer << " size " << fBufferSize << endl;
    if (!fpBuffer) return;

    int words=fBufferSize/sizeof(AliHLTUInt32_t);
    int lines=words/4;
    AliHLTUInt32_t* pTgt=reinterpret_cast<AliHLTUInt32_t*>(fpBuffer);
    for (int line=0; line<lines; line++) {
      for (int i=0; i<4; i++, pTgt++) cout << hex << " 0x" << (*pTgt);
      cout << endl;
    }
    if (words*sizeof(AliHLTUInt32_t)!=fBufferSize) cout << "some remaining bytes left due to alignment" << endl;
  }

private:
  AliHLTUInt8_t* fpBuffer;
  AliHLTUInt32_t fBufferSize;

  vector<AliHLTRandomBufferInstance> fInstances;
};

AliHLTRandomBuffer gData(2*gPagesize);

/**
 * @class AliHLTTestRawPage
 * Helper class for testing cycles of memory allocation, buffer modification, and
 * cleanup with the AliHLTRawPage within the function testLoop.
 */
class AliHLTTestRawPage {
public:
  AliHLTTestRawPage() : fBuffer(NULL) {}
  AliHLTTestRawPage(const AliHLTTestRawPage& src) : fBuffer(src.fBuffer) {}
  AliHLTTestRawPage& operator=(const AliHLTTestRawPage& src) {
    fBuffer=src.fBuffer;
    return *this;
  }

  ~AliHLTTestRawPage() {}

  bool IsAllocated() {return fBuffer!=NULL;}
  bool Alloc(AliHLTUInt32_t size, int verbosity) {
    fBuffer=AliHLTDataBuffer::AliHLTRawPage::GlobalAlloc(size);
    if (!fBuffer) return false;
    if (verbosity>1) {
      printf("allocated raw buffer %p from page %p\n", fBuffer, AliHLTDataBuffer::AliHLTRawPage::FindPage(fBuffer));
      fBuffer->Print("min");
    }
    return true;
  }

  AliHLTUInt8_t* GetPointer() {return fBuffer?fBuffer->GetPointer():NULL;}

  bool SetSize(AliHLTUInt32_t size, int verbosity) {
    AliHLTDataBuffer::AliHLTRawPage* rawpage=AliHLTDataBuffer::AliHLTRawPage::FindPage(fBuffer);
    if (rawpage) {
      if (rawpage->SetSize(fBuffer, size)==0) {
	if (verbosity>1) {
	  cout << "setting size for raw buffer " << fBuffer << " of page " << rawpage << endl;
	  fBuffer->Print("min");
	}
	return true;
      } else {
	cerr << "failed to set size for raw buffer " << fBuffer << endl; 
      }
    } else {
      cerr << "can not find raw page for buffer " << fBuffer << endl;
    }
    
    return false;
  }

  bool Free(int verbosity) {
    AliHLTDataBuffer::AliHLTRawPage* rawpage=AliHLTDataBuffer::AliHLTRawPage::FindPage(fBuffer);
    if (rawpage) {
      if (rawpage->Free(fBuffer)==0) {
	if (verbosity>1) cout << "released raw buffer " << fBuffer << " from page " << rawpage << endl;
	fBuffer=NULL;
	return true;
      } else {
	cerr << "failed to release raw buffer " << fBuffer  << endl;
      }
    } else {
      cerr << "can not find raw page for buffer " << fBuffer  << endl;
    }
    return false;
  }

  AliHLTUInt32_t GetTotalSize() {return fBuffer?fBuffer->GetTotalSize():0;}

  bool GlobalClean(int verbosity) {
    int nofPages=0;
    for (AliHLTDataBuffer::AliHLTRawPage* rawpage=AliHLTDataBuffer::AliHLTRawPage::NextPage(NULL);
	 rawpage!=NULL; 
	 rawpage=AliHLTDataBuffer::AliHLTRawPage::NextPage(rawpage)) {
      nofPages++;
    }
    if (verbosity>=1) {
      cout << "total number of pages: " << nofPages << endl;
    }

    AliHLTDataBuffer::AliHLTRawPage::GlobalClean();
    return true;
  }

  void Print(const char* options) {
    if (strcmp(options, "pages")==0 ||
	strcmp(options, "global")==0) {
      AliHLTDataBuffer::AliHLTRawPage::NextPage(NULL)->Print("global");
      return;
    }
  }
private:
  AliHLTDataBuffer::AliHLTRawBuffer* fBuffer;
};

struct testRawPageDescription {
  testRawPageDescription() : fLevel(0), fBufferSize(0), fDataId(-1), fTest() {};
  int fLevel;
  AliHLTUInt32_t fBufferSize;
  int fDataId;
  AliHLTTestRawPage fTest;
};

class AliHLTDataBufferWrapper : public AliHLTDataBuffer {
public:
  int ResetDataBuffer() {
    return AliHLTDataBuffer::ResetDataBuffer();
  }
};

/**
 * @class AliHLTTestDataBuffer
 * Helper class for testing cycles of memory allocation, buffer modification, and
 * cleanup with the AliHLTDataBuffer within the function testLoop.
 */
class AliHLTTestDataBuffer {
public:
  AliHLTTestDataBuffer() : fInstance(NULL), fBuffer(NULL), fSize(0) {}
  AliHLTTestDataBuffer(const AliHLTTestDataBuffer&)  : fInstance(NULL), fBuffer(NULL), fSize(0) {} ;
  AliHLTTestDataBuffer& operator=(const AliHLTTestDataBuffer&) {
    fInstance=NULL;
    fBuffer=NULL;
    fSize=0;
    return *this;
  }

  bool IsAllocated() {return fBuffer!=NULL;}
  bool Alloc(AliHLTUInt32_t size, int verbosity) {
    if (!fInstance) {
      fInstance=new AliHLTDataBufferWrapper;
      if (!fInstance) return false;
      fInstance->SetGlobalLoggingLevel(verbosity>1?kHLTLogAll:(AliHLTComponentLogSeverity)0x7c);
    }
    
    fBuffer=fInstance->GetTargetBuffer(size);
    if (!fBuffer) return false;
    fSize=size;
    if (verbosity>1) {
      cout << "allocated data buffer" << this << endl;
    }
    return true;
  }

  AliHLTUInt8_t* GetPointer() {return fBuffer;}

  bool SetSize(AliHLTUInt32_t size, int verbosity) {
    if (!fInstance) return false;
    AliHLTComponentBlockData bd;
    AliHLTComponent::FillBlockData(bd);
    bd.fSize=size;
    bd.fOffset=0;
    fInstance->SetSegments(fBuffer, &bd, 1);

    if (verbosity>1) {
      cout << "setting size for data buffer " << this << ": allocated size " << fSize << " -> " << size << endl;
    }
    fSize=size;
    return true;
  }

  bool Free(int verbosity) {
    if (!fInstance) return false;
    fSize=0;
    if (fInstance->ResetDataBuffer()==0) {
      if (verbosity>1) cout << "released data buffer " << this << endl;
      fBuffer=NULL;
      return true;
    } else {
      cerr << "failed to release raw buffer " << fBuffer  << endl;
    }

    return false;
  }

  AliHLTUInt32_t GetTotalSize() {return fSize;}

  bool GlobalClean(int verbosity) {
    if (!fInstance) return false;
    if (verbosity>=1) {
      fInstance->PrintStatistics();
    }
    delete fInstance;

    AliHLTDataBuffer::AliHLTRawPage::GlobalClean();
    return true;
  }

  void Print(const char* options) {
    if (strcmp(options, "pages")==0 ||
	strcmp(options, "global")==0) {
      fInstance->PrintStatistics();
      AliHLTDataBuffer::AliHLTRawPage::NextPage(NULL)->Print("global");
      return;
    }
  }
private:
  AliHLTDataBufferWrapper* fInstance;
  AliHLTUInt8_t* fBuffer;
  AliHLTUInt32_t fSize;
};

struct testDataBufferDescription {
  testDataBufferDescription() : fLevel(0), fBufferSize(0), fDataId(-1), fTest() {};
  int fLevel;
  AliHLTUInt32_t fBufferSize;
  int fDataId;
  AliHLTTestDataBuffer fTest;
};

template<class T,class S>
int CopyDescription(vector<T>& target, vector<S>& source)
{
  target.clear();
  for (unsigned i=0; i<source.size(); i++) {
    T clone;
    clone.fLevel=source[i].fLevel;
    clone.fBufferSize=source[i].fBufferSize;
    clone.fDataId=source[i].fDataId;
    target.push_back(clone);
  }
  return 0;
}

template<class T>
int fillTestSample(int levels, int processes, vector<T>& descriptions)
{
  int availableProcesses=processes;
  for (int level=0; level<levels; level++) {
    int levelprocesses=GetRandom(availableProcesses/2, availableProcesses-(levels-level));
    for (int process=0; process<levelprocesses; process++) {
      T desc;
      desc.fLevel=level;
      desc.fBufferSize=GetRandom(1024, 2*gPagesize);
      //desc.fTest.Clear();
      descriptions.push_back(desc);
    }
    availableProcesses-=levelprocesses;
  }
  return 0;
}

template<class T>
int allocateBuffers(vector<T>& descriptions,
		    int level=-1)
{
  for (unsigned i=0; i<descriptions.size(); i++) {
    if (level>=0 && descriptions[i].fLevel<level) continue;
    if (level>=0 && descriptions[i].fLevel>level) break;
    if (descriptions[i].fTest.IsAllocated()) {
      cerr << "warning: buffer already allocated" << endl;
      continue;
    }

    if (descriptions[i].fTest.Alloc(GetRandom(descriptions[i].fBufferSize/2, descriptions[i].fBufferSize), gVerbosity)) {
    }

    if (!descriptions[i].fTest.IsAllocated()) {
      cerr << "failed to allocate buffer for process " << i << endl;
      return -EFAULT;
    }
  }
  return 0;
}
		   
template<class T>
int setBufferSizes(vector<T>& descriptions,
		   int level=-1)
{
  // set buffer size for all processes of the specified level
  // buffer size is chosen randomly between 0 and the allocated size
  vector<unsigned> positions;
  for (unsigned i=0; i<descriptions.size(); i++) {
    if (!descriptions[i].fTest.IsAllocated()) continue;
    if (level>=0 && descriptions[i].fLevel<level) continue;
    if (level>=0 && descriptions[i].fLevel>level) break;
    positions.push_back(i);
  }

  random_shuffle(positions.begin(), positions.end());
  for (vector<unsigned>::iterator position=positions.begin();
       position!=positions.end(); position++) {
    AliHLTUInt32_t datasize=GetRandom(0, descriptions[*position].fTest.GetTotalSize());
    descriptions[*position].fDataId=gData.CreateInstance(datasize);
    const AliHLTUInt8_t* pSrc=gData.GetBuffer(descriptions[*position].fDataId);
    if (!pSrc) {
      cout << "error creating data instance" << endl;
      return -1;
    }
    AliHLTUInt8_t* pTgt=descriptions[*position].fTest.GetPointer();
    if (!pTgt) {
      cout << "error target buffer" << endl;
      return -1;
    }
    memcpy(pTgt, pSrc, datasize);
    if (descriptions[*position].fTest.SetSize(datasize, gVerbosity)) {
      if (gVerbosity>1) cout << "  process " << *position << " level " << descriptions[*position].fLevel << endl;
    } else {
      cout << "failed allocation for process " << *position << " level " << descriptions[*position].fLevel << endl;
      return -1;
    }
  }
  return 0;
}

template<class T>		   
int releaseBuffers(vector<T>& descriptions,
		   int level=-1)
{
  // find the processes to be releases according to the specified level
  // shuffle the processes and then release in this random order
  vector<unsigned> positions;
  for (unsigned i=0; i<descriptions.size(); i++) {
    if (!descriptions[i].fTest.IsAllocated()) continue;
    if (level>=0 && descriptions[i].fLevel<level) continue;
    if (level>=0 && descriptions[i].fLevel>level) break;
    positions.push_back(i);
  }

  random_shuffle(positions.begin(), positions.end());
  for (vector<unsigned>::iterator position=positions.begin();
       position!=positions.end(); position++) {
    if (!gData.CompareBuffer(descriptions[*position].fTest.GetPointer(), descriptions[*position].fTest.GetTotalSize(), descriptions[*position].fDataId)) {
      cout << "data comparison failed for process " << *position << " level " << descriptions[*position].fLevel << endl;
      return -1;
    }
    gData.ReleaseInstance(descriptions[*position].fDataId);
    if (descriptions[*position].fTest.Free(gVerbosity)) {
      if (gVerbosity>1) cout << "  process " << *position << " level " << descriptions[*position].fLevel << endl;
    } else {
      cout << "failed deallocation for process " << *position << " level " << descriptions[*position].fLevel << endl;
      return -1;
    }
  }
  return 0;
}

template<class T>
int testLoop(vector<T>& descriptions, int nofLevels, int nofCycles)
{
  // test the buffer allocation and deallocation from the AliHLTRawPage
  //
  // build a hierarchy of processes like in a real chain, in the first
  // implementation it is only a linear tree
  // the number of levels and processes are randomly chosen from a certain range
  // processes are distributed among the levels such that there is a decreasing
  // number of processes but at least one process
  //
  // The test cycle consists of the following steps
  // - loop over all levels
  //   - allocation of a maximum buffer for every process of the current level
  //   - set the buffer size to be used for every process of the current level
  //   - release buffers of the previous level
  // - check all created raw pages and throw error, if
  //   - not all buffers released
  //   - not all memory released
  //   - still fragmented
  // 
  int iResult=0;

  AliHLTDataBuffer::AliHLTRawPage* rawpage=NULL;
  for (int cycle=0; cycle<nofCycles && iResult>=0; cycle++) {
    for (int level=0; level<=nofLevels; level++) {
      // allocate buffers
      if (level<nofLevels) {
	if ((iResult=allocateBuffers(descriptions, level))<0) {
	  cerr << "failed to allocate buffers" << endl;
	  return iResult;
	}
      }

      if (gVerbosity>1) {
	cout << "finished allocation - level " << level << " ";
	descriptions[0].fTest.Print("pages");
      }

      if (level<nofLevels) {
	if ((iResult=setBufferSizes(descriptions, level))<0) {
	  cerr << "failed setting  buffers" << endl;
	  return iResult;
	}
      }

      if (gVerbosity>1) {
	cout << "finished size modification - level " << level << " ";
	descriptions[0].fTest.Print("pages");
      }

      if (level>0) {
	if ((iResult=releaseBuffers(descriptions, level-1))<0) {
	  cerr << "failed to release buffers" << endl;
	  return iResult;
	}
      }
    }

    if (gVerbosity>1) {
      cout << "status of released pages: " << endl << "  ";
      descriptions[0].fTest.Print("pages");
    }
    for (rawpage=AliHLTDataBuffer::AliHLTRawPage::NextPage(NULL);
	 rawpage!=NULL; 
	 rawpage=AliHLTDataBuffer::AliHLTRawPage::NextPage(rawpage)) {
      if (rawpage->IsUsed()) {
	cerr << "page " << rawpage << " has used buffers" << endl;
	iResult=-EFAULT;
      } else if (rawpage->Size()!=rawpage->Capacity()) {
	cerr << "page " << rawpage << " not completely released" << endl;
	iResult=-EFAULT;
      } else if (rawpage->IsFragmented()) {
	cerr << "page " << rawpage << " is still fragmented" << endl;
	iResult=-EFAULT;
      }
    }
    AliHLTDataBuffer::SetGlobalEventCount(cycle);
  }

  if (iResult<0) {
    descriptions[0].fTest.Print("pages");
  }

  descriptions[0].fTest.GlobalClean(gVerbosity);

  return iResult;
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

int testAliHLTDataBuffer()
{
  int iResult=0;
  const int iMaxLevels=10;
  int nofLevels=GetRandom(3,iMaxLevels);
  int nofProcesses=GetRandom(nofLevels,2*iMaxLevels);
  int nofCycles=GetRandom(100,1000);

  vector<testDataBufferDescription> descriptions;

  if ((iResult=fillTestSample(nofLevels, nofProcesses, descriptions))<0) {
    cerr << "failed to fill test sample" << endl;
    return iResult;
  }

  // adjust to the level of the last entry
  nofLevels=descriptions[descriptions.size()-1].fLevel;

  if (gVerbosity>0) {
    cout << "    Test setup: " << nofLevels << " levels/ " << nofProcesses << " processes/ " << nofCycles << " cycle(s)" << endl;
    cout << "    Pagesize: " << gPagesize << endl;
    for (unsigned i=0; i<descriptions.size(); i++) {
      cout << "      process " << i << ": level " << descriptions[i].fLevel << "   max size " << descriptions[i].fBufferSize << endl;
    }
  }

  if ((iResult=testLoop(descriptions, nofLevels, nofCycles))<0) return iResult;

  if (gVerbosity>0) {
    cout << "checking memory allocation using AliHLTRawPage directly ..." << endl;
  }
  vector<testRawPageDescription> rawpagedesc;
  CopyDescription(rawpagedesc, descriptions);
  if ((iResult=testLoop(rawpagedesc, nofLevels, nofCycles))<0) return iResult;

  return 0;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// main functions

int main(int /*argc*/, const char** /*argv*/)
{
  int iResult=0;
  const int nofCycles=10;
  for (int cycle=0; cycle<nofCycles && iResult>=0; cycle++) {
    iResult=testAliHLTDataBuffer();
    if (((cycle+1)*100/nofCycles)%10==0) {
      cout << (cycle+1)*100/nofCycles << "% complete" << endl;
    }
  }

  return iResult;
}
