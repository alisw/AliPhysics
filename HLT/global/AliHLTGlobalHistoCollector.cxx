
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
    fStore()
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
  list.push_back( kAliHLTDataTypeHistogram );
}

AliHLTComponentDataType AliHLTGlobalHistoCollector::GetOutputDataType() 
{ 
  // see header file for class documentation
  return kAliHLTDataTypeHistogram;
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


int AliHLTGlobalHistoCollector::Reconfigure(const char* cdbEntry, const char* chainId) { 
  // see header file for class documentation
  
  return 0; // no CDB entry exist

  int iResult=0;  
  const char* path="HLT/ConfigGlobal/GlobalHistoCollector";
  const char* defaultNotify="";
  if (cdbEntry) {
    path=cdbEntry;
    defaultNotify=" (default)";
  }
  
  if (path) {
    HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
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
}


void AliHLTGlobalHistoCollector::Clear() 
{ 
  // reset the store

  for( unsigned int i=0; i<fStore.size(); i++ ){
    for( unsigned int j=0; j<fStore[i].fHistos.size(); j++ ){
      delete fStore[i].fHistos[j].fHisto;
    }
    delete fStore[i].fMergedHisto;
  }
  fStore.clear();
}




int AliHLTGlobalHistoCollector::DoEvent(const AliHLTComponentEventData & evtData, AliHLTComponentTriggerData& /*trigData*/)
{
  // see header file for class documentation

  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )) return 0;

  if( fUID == 0 ){
    TTimeStamp t;
    fUID = ( gSystem->GetPid() + t.GetNanoSec())*10 + evtData.fEventID;
  }

 
  const TObject *iter = NULL;        
  for(iter = GetFirstInputObject(kAliHLTDataTypeHistogram); iter != NULL; iter = GetNextInputObject()){
  
    if( !( GetDataType(iter) == kAliHLTDataTypeHistogram )) continue;    
    
    const TH1 *h = dynamic_cast<TH1*>(const_cast<TObject*>( iter ) );
    if( !h ) continue;
    //cout<<"received histo "<<h->GetName()<<" with id="<<GetSpecification(iter)<<endl;
    // search for the base entry, if not exist then create a new entry   

    int iBase = -1;
    for( unsigned int i=0; i<fStore.size(); i++ ){
      if( fStore[i].fDataType == GetDataType(iter) &&
	  TString(fStore[i].fMergedHisto->GetName()).CompareTo(h->GetName())==0){
	iBase = i;
	break;
      }
    }

   if( iBase<0 ){
      AliHLTHistoBaseData b;
      b.fDataType = GetDataType(iter);
      b.fMergedHisto = (TH1*) iter->Clone();
      fStore.push_back(b);
      iBase = fStore.size()-1;
    }

    // search for the specific entry, if not exist then create a new one

    AliHLTHistoBaseData &b = fStore[iBase];
    
    int iSpec=-1;
    for( unsigned int i=0; i<fStore[iBase].fHistos.size(); i++ ){
      AliHLTHistoData &d = b.fHistos[i];
      if( d.fSpecification == GetSpecification(iter) ){
	iSpec = i;
	break;
      }
    }

    if( iSpec<0 ){
      AliHLTHistoData d;
      d.fSpecification = GetSpecification(iter);
      d.fHisto = (TH1*) iter->Clone();
      b.fHistos.push_back(d);
      iSpec = b.fHistos.size()-1;      
    }
    b.fHistos[iSpec].fHisto->Reset();
    b.fHistos[iSpec].fHisto->Add( (TH1*)iter, 1);
    //cout<<"index = "<<iBase<<","<<iSpec<<", nentr="<<b.fHistos[iSpec].fHisto->GetEntries()<<endl;
    
    // merge histos 
    b.fMergedHisto->Reset();

    for( unsigned int i=0; i<b.fHistos.size(); i++ ){
      b.fMergedHisto->Add(b.fHistos[i].fHisto,1);
    }
    //cout<<" merged="<<b.fMergedHisto->GetEntries()<<endl;

  } // end for loop over histogram blocks

  for( unsigned int i=0; i<fStore.size(); i++ ){
    PushBack((TObject*) fStore[i].fMergedHisto, fStore[i].fDataType, fUID );
  }
  return 0;
}
  
