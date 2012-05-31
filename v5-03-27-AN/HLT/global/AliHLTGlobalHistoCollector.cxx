
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>          *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTGlobalHistoCollector.cxx
    @author Kalliopi Kanaki
    @date   
    @brief  The Histogram Handler component
*/

#if __GNUC__>= 3
using namespace std;
#endif
#include "AliHLTGlobalHistoCollector.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TH1.h"
#include "TTimeStamp.h"
#include "TSystem.h"

ClassImp(AliHLTGlobalHistoCollector) //ROOT macro for the implementation of ROOT specific class methods

  AliHLTGlobalHistoCollector::AliHLTGlobalHistoCollector()
    :    
    fUID(0),
    fStore(),
    fBenchmark("GlobalHistoCollector")
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTGlobalHistoCollector::~AliHLTGlobalHistoCollector() { 
  // see header file for class documentation
  Clear();
}



// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTGlobalHistoCollector::GetComponentID() 
{ 
  // see header file for class documentation
  return "GlobalHistoCollector";
}

void AliHLTGlobalHistoCollector::GetInputDataTypes( vector<AliHLTComponentDataType>& list) 
{ 
  // see header file for class documentation

  list.clear(); 
  list.push_back( kAliHLTAllDataTypes );
}

AliHLTComponentDataType AliHLTGlobalHistoCollector::GetOutputDataType() 
{ 
  // see header file for class documentation
  return kAliHLTAllDataTypes;
}


void AliHLTGlobalHistoCollector::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) 
{ 
  // see header file for class documentation
  constBase=0;
  inputMultiplier=1.0;
}

AliHLTComponent* AliHLTGlobalHistoCollector::Spawn() 
{ 
  // see header file for class documentation
  return new AliHLTGlobalHistoCollector();
}
	
int AliHLTGlobalHistoCollector::DoInit( int argc, const char** argv ) 
{ 
  // see header file for class documentation

  Clear(); 
  fBenchmark.Reset();
  fBenchmark.SetTimer(0,"total");
  fBenchmark.SetTimer(1,"merging");

  int iResult=0;
  
  TString configuration="";
  TString argument="";
  for (int j=0; j<argc && iResult>=0; j++) {
    
    argument=argv[j];
    if (!configuration.IsNull()) configuration+=" ";
    configuration+=argument;    
  }
   
  if (!configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  } else {
    iResult=Reconfigure(NULL, NULL);
  }
  fUID = 0;
  return iResult;
}


int AliHLTGlobalHistoCollector::DoDeinit() 
{ 
  // see header file for class documentation 
    
  Clear();
  fUID = 0;
  return 0;
}

int AliHLTGlobalHistoCollector::Configure(const char* arguments) 
{ 
  // see header file for class documentation
  
  int iResult=0;
  if (!arguments) return iResult;
  HLTInfo("parsing configuration string \'%s\'", arguments);

  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;
     
      //if (argument.CompareTo("-sum-noise-histograms")==0) {
      //fNoiseHistograms = kTRUE;
      //HLTInfo("got \'-sum-noise-histograms\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
      //else 
      {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    } // end for
    
    delete pTokens;
  
  } // end if pTokens
  
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  return iResult;
}


int AliHLTGlobalHistoCollector::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/) { 
  // see header file for class documentation
  
  return 0; // no CDB entry exist
  /*
  int iResult=0;  
  const char* path="HLT/ConfigGlobal/GlobalHistoCollector";
  const char* defaultNotify="";
  if (cdbEntry) {
    path=cdbEntry;
    defaultNotify=" (default)";
  }
  
  if (path) {
    HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path);
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
	HLTInfo("received configuration object string: \'%s\'", pString->GetString().Data());
	iResult=Configure(pString->GetString().Data());
      } else {
	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("cannot fetch object \"%s\" from CDB", path);
    }
  }
  return iResult;
*/
}


void AliHLTGlobalHistoCollector::Clear() 
{ 
  // reset the store

  for( unsigned int i=0; i<fStore.size(); i++ ){
    for( unsigned int j=0; j<fStore[i].fInstances.size(); j++ ){
      delete fStore[i].fInstances[j].fObject;
    }
    delete fStore[i].fMergedObject;
  }
  fStore.clear();
}




int AliHLTGlobalHistoCollector::DoEvent(const AliHLTComponentEventData & evtData, AliHLTComponentTriggerData& /*trigData*/)
{
  // see header file for class documentation
  //cout<<"\n\nDoEvent called"<<endl;

  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )) return 0;

  fBenchmark.StartNewEvent();
  fBenchmark.Start(0);

  if( fUID == 0 ){
    TTimeStamp t;
    fUID = ( gSystem->GetPid() + t.GetNanoSec())*10 + evtData.fEventID;
  }

  for( const AliHLTComponentBlockData *i= GetFirstInputBlock(); i!=NULL; i=GetNextInputBlock() ){
    fBenchmark.AddInput(i->fSize);
  }

  const TObject *iter = NULL;
  for(iter = GetFirstInputObject(); iter != NULL; iter = GetNextInputObject()){
            
    if( !iter->InheritsFrom(TH1::Class()) 
	&& !iter->InheritsFrom(TSeqCollection::Class()) ) continue;

    //cout<<"received object "<<iter->GetName()<<" with id="<<GetSpecification(iter)<<endl;

    //search for the base entry, if not exist then create a new entry   

    int iColl = -1;
    for( unsigned int i=0; i<fStore.size(); i++ ){
      if( fStore[i].fHLTDataType != GetDataType(iter) ) continue;
      if( fStore[i].fInstances.size()<1 ) continue; 
      TObject * obj = fStore[i].fInstances[0].fObject;
      if( !obj ) continue;
      if( TString(obj->GetName()).CompareTo(iter->GetName())==0){
	iColl = i;
	break;
      }
    }
    //cout<<"Collection found: "<<iColl<<endl;
    if( iColl<0 ){
      AliHLTGlobalHCCollection c;
      c.fHLTDataType = GetDataType(iter);
      c.fMergedObject = 0;
      c.fNeedToMerge = 1;
      fStore.push_back(c);
      iColl = fStore.size()-1;
    }else{
      fStore[iColl].fNeedToMerge = 1;
    }

    // search for the specific entry, if not exist then create a new one
    
    AliHLTGlobalHCCollection &c = fStore[iColl];
   
    int iSpec=-1;
    for( unsigned int i=0; i<c.fInstances.size(); i++ ){
      AliHLTGlobalHCInstance &inst = c.fInstances[i];
      if( inst.fHLTSpecification == GetSpecification(iter) ){
	iSpec = i;
	break;
      }
    }
    //cout<<"Instance found:"<<iSpec<<endl;
    if( iSpec<0 ){
      AliHLTGlobalHCInstance inst;
      inst.fHLTSpecification = GetSpecification(iter);
      inst.fObject = 0;
      c.fInstances.push_back(inst);
      iSpec = c.fInstances.size()-1;      
    }else{
      delete c.fInstances[iSpec].fObject;
    }

    c.fInstances[iSpec].fObject = iter->Clone();
    
    //cout<<"index = "<<iColl<<","<<iSpec<<endl;

  } // end for loop over input blocks

  fBenchmark.Start(1);
 
  
  // merge histos 

  for( unsigned int iColl = 0; iColl<fStore.size(); iColl++){
    AliHLTGlobalHCCollection &c = fStore[iColl];
    if( !c.fNeedToMerge && c.fMergedObject ) continue;
    if( c.fInstances.size() <1 ) continue;
    delete c.fMergedObject;
    c.fMergedObject = c.fInstances[0].fObject->Clone();
    TList l;
    for( unsigned int i=1; i<c.fInstances.size(); i++ ){
      l.Add(c.fInstances[i].fObject);
    }

    if( c.fMergedObject->InheritsFrom(TH1::Class()) ){
      TH1 *histo = dynamic_cast<TH1*>(c.fMergedObject);
      if( histo ) histo->Merge(&l);
    }
    else if( c.fMergedObject->InheritsFrom(TSeqCollection::Class()) ){
      TSeqCollection *list = dynamic_cast<TSeqCollection*>(c.fMergedObject);
      if( list ) list->Merge(&l);
    }    
    c.fNeedToMerge = 0;
  }
  fBenchmark.Stop(1);
 
  // send output 

  for( unsigned int i=0; i<fStore.size(); i++ ){
    if( fStore[i].fMergedObject ){
      PushBack((TObject*) fStore[i].fMergedObject, fStore[i].fHLTDataType, fUID );
      fBenchmark.AddOutput(GetLastObjectSize());
    }
  }

  fBenchmark.Stop(0);
  HLTInfo(fBenchmark.GetStatistics());

  return 0;
}
  
