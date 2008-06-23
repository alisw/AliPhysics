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

/** @file   testAliHLTEsdManager.C
    @author Matthias Richter
    @date   
    @brief  Test macro/program for the AliHLTEsdManager class
 */

#ifndef __CINT__
#include "TFile.h"
#include "TDatime.h"
#include "TRandom.h"
#include "TArrayI.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "AliESDEvent.h"
#include "AliHLTDataTypes.h"
#include "AliHLTComponent.h"
#include "AliHLTEsdManager.h"
#include "AliHLTMessage.h"
#include "AliHLTSystem.h"
#include <ostream>
#endif //__CINT__

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// configuration of the test program
//

// printouts or not
const bool bVerbose=false;


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// forward declarations
//
class AliHLTEsdManager;
int GetRandom(int min, int max);
int CreateAndWriteESD(AliHLTEsdManager* manager, int eventno, AliHLTComponentDataType dt, AliESDEvent* pTgt);
int CheckFields(const char* file, TArrayI* fields);
int CheckFields(TTree* pTree, TArrayI* fields, const char* file);

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// compiled version and root macro are not equivalent for this unit test
// The macro just tests with one data type
// This is mainly due to some restrictions in CINT which can not handle
// the array of data types correctly
int testAliHLTEsdManager()
{
  cout << "macro still not working with CINT, sorry!" << endl;
  return 0;

  int iResult=0;
#ifdef __CINT__
  if (gSystem->Load("libHLTrec.so")<0) {
    cerr << "error: error loading libHLTrec.so library" << endl;
    return -1;
  }
#endif

  AliHLTEsdManager*  pManager=AliHLTEsdManager::New();
  if (!pManager) {
    cerr << "error: can not create manager instance" << endl;
    return -1;
  }
  pManager->SetDirectory(gSystem->TempDirectory());

  int nofEvents=10;
  AliHLTComponentDataType tpcesd;
  AliHLTComponent::SetDataType(tpcesd, "ESD_TREE", "TPC ");
  cout << AliHLTComponent::DataType2Text(tpcesd).c_str() << endl;
  for (int event=0; event<nofEvents && iResult>=0; event++) {
    cout << AliHLTComponent::DataType2Text(tpcesd).c_str() << endl;
    CreateAndWriteESD(pManager, event, tpcesd, NULL);
  }

  AliHLTEsdManager::Delete(pManager);
  return iResult;
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
int main(int /*argc*/, const char** /*argv*/)
{
  int iResult=0;
  int nofEvents=10;
  AliHLTSystem gHLT;
  gHLT.SetGlobalLoggingLevel(kHLTLogDefault);

  AliHLTEsdManager*  pManager=AliHLTEsdManager::New();
  if (!pManager) {
    cerr << "error: can not create manager instance" << endl;
    return -1;
  }
  pManager->SetDirectory(gSystem->TempDirectory());

  AliHLTComponentDataType types[] = {
    // first entry is special, ESD is written to the global target ESD
    //kAliHLTDataTypeESDTree|kAliHLTDataOriginTPC,
    kAliHLTDataTypeESDTree|kAliHLTDataOriginTPC,
    kAliHLTDataTypeESDTree|kAliHLTDataOriginPHOS,
    kAliHLTDataTypeESDTree|kAliHLTDataOriginTRD,
    kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTRD,
    kAliHLTVoidDataType
  };

  TTree* pMasterTree=new TTree("esdTree", "Tree with HLT ESD objects");
  pMasterTree->SetDirectory(0);
  AliESDEvent* pMasterESD=new AliESDEvent;
  pMasterESD->CreateStdContent();
  pMasterESD->WriteToTree(pMasterTree);

  vector<TArrayI*> randomFields;
  for (int event=0; event<nofEvents && iResult>=0; event++) {
    pMasterESD->ResetStdContent();
    for (unsigned int type=0; types[type]!=kAliHLTVoidDataType && iResult>=0; type++) {
      if (randomFields.size()<=type) {
	randomFields.push_back(new TArrayI(nofEvents));
      }
      AliESDEvent* pTgt=NULL;
      //if (type==0) pTgt=pMasterESD;
      int field=CreateAndWriteESD(pManager, event, types[type], pTgt);
      if (field>=0) {
	(*randomFields[type])[event]=field;
      } else {
	iResult=-1;
	break;
      }
    }
    pMasterTree->Fill();
  }

  if (iResult>=0) {
    pManager->PadESDs(nofEvents);
  }

  vector<TString> filenames;
  for (int type=0; types[type]!=kAliHLTVoidDataType; type++) {
    TString filename=pManager->GetFileNames(types[type]);
    filenames.push_back(filename);
  }
  // delete the manager instance to make sure the files are closed
  AliHLTEsdManager::Delete(pManager);
  pManager=NULL;
  for (int type=0; types[type]!=kAliHLTVoidDataType; type++) {
    if (iResult>=0) {
      iResult=CheckFields(filenames[type], dynamic_cast<TArrayI*>(randomFields[type]));
    }
    TString shellcmd="rm -f ";
    shellcmd+=filenames[type];
    gSystem->Exec(shellcmd);    
  }

  vector<TArrayI*>::iterator element;
  while ((element=randomFields.end())!=randomFields.begin()) {
    element--;
    if (*element) delete *element;
    randomFields.pop_back();
  }

  delete pMasterESD;

  return iResult;
}

Bool_t seedSet=kFALSE;

/**
 * Get a random number in the given range.
 */
int GetRandom(int min, int max)
{
  if (max-min<2) return min;
  static TRandom rand;
  if (!seedSet) {
    TDatime dt;
    rand.SetSeed(dt.Get());
    seedSet=kTRUE;
  }
  return rand.Integer(max-min);
}

/**
 * Creates a dummy ESD object and sets the magnetic field to a random number.
 * The ESD is streamed via AliHLTMessage and processed by the AliHLTEsdmanager.
 * @return 0 no ESD created, >0 random number of the magnetic field, 
 *         neg error if failed
 */
int CreateAndWriteESD(AliHLTEsdManager* pManager, int eventno, AliHLTComponentDataType dt, AliESDEvent* pTgt)
{
  int iResult=0;
  int magField=0;
  if (!pManager) {
    cerr << "error: missing manager instance" << endl;
    return -1;
  }
  const char* message="";
  if ((GetRandom(0,10)%3)==0) {
    message=": adding ESD for block ";
    TTree* pTree=new TTree;
    AliESDEvent* pESD=new AliESDEvent;
    pESD->CreateStdContent();
    magField=GetRandom(1, 1000);
    pESD->SetMagneticField(magField);
    pESD->WriteToTree(pTree);
    pTree->Fill();
    pTree->GetUserInfo()->Add(pESD);
    AliHLTMessage msg(kMESS_OBJECT);
    msg.WriteObject(pTree);
    Int_t iMsgLength=msg.Length();
    if (iMsgLength>0) {
      msg.SetLength(); // sets the length to the first (reserved) word
      iResult=pManager->WriteESD((AliHLTUInt8_t*)msg.Buffer(), iMsgLength, dt, pTgt, eventno);
    }
    pTree->GetUserInfo()->Clear();
    delete pTree;
    delete pESD;
  } else {
    message=": omitting block       ";
  }
  if (iResult>=0) iResult=magField;
  if (bVerbose) cout << "event " << eventno << message << AliHLTComponent::DataType2Text(dt).c_str() << ": " << iResult << endl;
  return iResult;
}

/**
 * Read the ESD from the file and compare with the
 * random field values previously set to the ESDs
 */
int CheckFields(const char* file, TArrayI* fields)
{
  if (!file || !fields) {
    cerr << "error: invalid parameters" << endl;
    return 0;
  }
  TFile esdfile(file);
  if (!esdfile.IsZombie()) {
    TTree* pTree=NULL;
    esdfile.GetObject("esdTree", pTree);
    if (pTree) {
      int res=CheckFields(pTree, fields, file);
      if (res<0) return res;
    } else {
      cerr << "error: can not find esdTree in file " << file << endl;
    }
  } else {
    cerr << "error: can not open file " << file << endl;
    return -1;
  }
  cout << "checking: " << file << " ok" << endl;
  return 0;
}

/**
 * Compare ESD from tree with the
 * random field values previously set to the ESDs
 */
int CheckFields(TTree* pTree, TArrayI* fields, const char* file)
{
  if (fields->GetSize()!=pTree->GetEntries()) {
    cerr << "error: event number mismatch in file " << file << " : expected " << fields->GetSize() << "  found " << pTree->GetEntries() << endl;
    return -1;
  }
  AliESDEvent* pESD=new AliESDEvent;
  pESD->ReadFromTree(pTree);
  for (int event=0; event<pTree->GetEntries(); event++) {
    pTree->GetEvent(event);
    if (fields->At(event)!=pESD->GetMagneticField()) {
      cerr << "error: magnetic field mismatch in file " << file << " event " << event << ": expected " << fields->At(event) << "  found " << pESD->GetMagneticField() << endl;
      return -1;
    }
  }

  delete pESD;
  return 0;
}
