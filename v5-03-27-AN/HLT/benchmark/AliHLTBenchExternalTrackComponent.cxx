// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/** @file   AliHLTBenchExternalTrackComponent.cxx
    @author Matthias Richter
    @date   2008-10-30
    @brief  Benchmark component for AliExternalTrackParam transportation.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cassert>
#include "AliHLTBenchExternalTrackComponent.h"
#include "AliExternalTrackParam.h"
#include "AliHLTExternalTrackParam.h"
#include "TString.h"
#include "TDatime.h"
#include "TRandom.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TClass.h"
#include "TList.h"

/** global object for component registration */
AliHLTBenchExternalTrackComponent gAliHLTBenchExternalTrackComponent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTBenchExternalTrackComponent)

AliHLTBenchExternalTrackComponent::AliHLTBenchExternalTrackComponent()
  : AliHLTProcessor()
  , fVerbosity(0)
  , fDisableRegistry(false)
  , fMode(kPublishingOff)
  , fMaxSize(10000)
  , fMinSize(100)
  , fEventModulo(-1)
  , fRangeOffset(0)
  , fRangeMultiplicator(1.0)
  , fpTcArray(NULL)
  , fpTObjArray(NULL)
  , fpDice(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTBenchExternalTrackComponent::~AliHLTBenchExternalTrackComponent()
{
  // see header file for class documentation
}

const char* AliHLTBenchExternalTrackComponent::GetComponentID()
{
  // see header file for class documentation
  return "BenchmarkAliExternalTrackParam";
}

void AliHLTBenchExternalTrackComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponentDataType AliHLTBenchExternalTrackComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTMultipleDataType;
}

int AliHLTBenchExternalTrackComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeTObject);
  return tgtList.size();
}

void AliHLTBenchExternalTrackComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // see header file for class documentation
  constBase=sizeof(AliHLTExternalTrackParam)*fMaxSize;
  inputMultiplier=0.0;
}

AliHLTComponent* AliHLTBenchExternalTrackComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTBenchExternalTrackComponent;
}

int AliHLTBenchExternalTrackComponent::DoInit(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=0;
  TString argument="";
  bool bMissingParam=0;
  char* cpErr=NULL;
  int i=0;
  for (; i<argc && iResult>=0; i++) {
    cpErr=NULL;
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -tobjarray
    if (argument.CompareTo("-tobjarray")==0) {
      fMode=ktobjarray;

    // -tclonesarray
    } else if (argument.CompareTo("-tclonesarray")==0) {
      fMode=ktclonesarray;

    // -carray
    } else if (argument.CompareTo("-carray")==0) {
      fMode=kcarray;

    // -maxsize
    } else if (argument.CompareTo("-maxsize")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fMaxSize=strtoul( argv[i], &cpErr ,0);
      if ( *cpErr ) break;

    // -minsize
    } else if (argument.CompareTo("-minsize")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fMinSize=strtoul( argv[i], &cpErr ,0);
      if ( *cpErr ) break;

    // -rangemodulo
    } else if (argument.CompareTo("-rangemodulo")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fEventModulo=strtoul( argv[i], &cpErr ,0);
      if ( *cpErr ) break;

    // -rangeoffset
    } else if (argument.CompareTo("-rangeoffset")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fRangeOffset=strtoul( argv[i], &cpErr ,0);
      if ( *cpErr ) break;

    // -rangefactor
    } else if (argument.CompareTo("-rangefactor")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fRangeMultiplicator=strtof( argv[i], &cpErr);
      if ( *cpErr ) break;

    // -verbosity
    } else if (argument.CompareTo("-verbosity")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fVerbosity=strtoul( argv[i], &cpErr ,0);
      if ( *cpErr ) break;

    // -nocheck
    } else if (argument.CompareTo("-nocheck")==0) {
      fDisableRegistry=true;

    } else {
      HLTError("unknown argument %s", argument.Data());
      iResult=-EINVAL;
    }
  }

  if (cpErr && *cpErr) {
    HLTError("Cannot convert specifier '%s' for argument '%s'", argv[i], argument.Data());
    iResult=-EINVAL;
  } else if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  // use a TClonesArray internally and for publishing in case of
  // ktclonesarray. The maximum size and all members are allocated
  fpTcArray=new TClonesArray("AliExternalTrackParam", fMaxSize);
  if (fpTcArray) {
    fpTcArray->ExpandCreate(fMaxSize);
    switch (fMode) {
    case kPublishingOff:
      // nothing to do
      break;
    case ktclonesarray:
    case kcarray:
      // nothing to do
      break;
    case ktobjarray:
      // use a TObjArray which owns the objects
      fpTObjArray=new TObjArray;
      if (!fpTObjArray) {
	iResult=-ENOMEM;
      }
      break;
    default:
      HLTError("unknown publishing mode %d", fMode);
      iResult=-EINVAL;
    }
  } else {
    iResult=-ENOMEM;
  }

  if (iResult>=0) {
    fpDice=new TRandom;
    if (fpDice) {
      TDatime dt;
      fpDice->SetSeed(dt.Get());
    } else {
       iResult=-ENOMEM;
    }
  }

  return iResult;
}

int AliHLTBenchExternalTrackComponent::DoDeinit()
{
  // see header file for class documentation
  if (fpTObjArray) delete fpTObjArray;
  fpTObjArray=NULL;

  if (fpTcArray) {
    fpTcArray->Delete();
    delete fpTcArray;
    fpTcArray=NULL;
  }

  if (fpDice) delete fpDice;
  fpDice=NULL;

  return 0;
}

int AliHLTBenchExternalTrackComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
							const AliHLTComponentBlockData* /*blocks*/,
							AliHLTComponentTriggerData& /*trigData*/,
							AliHLTUInt8_t* outputPtr, 
							AliHLTUInt32_t& size,
							AliHLTComponentBlockDataList& outputBlocks)
{
  // see header file for class documentation

  // scan the input data blocks for types kAliHLTDataTypeTrack and
  // kAliHLTDataTypeTObjArray with AliExternalTrackParam objects
  // 
  // If data is available, the array is extracted into the internal
  // TClonesArray fpTcArray. If no input data is available, the array
  // is filled randomly, both in size and content.
  //
  // The array is published according to the fMode member either as
  // TClonesArray, TObjArray or C-structure (AliHLTExternalTrackParam)

  int iResult=0;
  AliHLTUInt32_t capacity=size;
  size=0;
  if (!IsDataEvent()) return 0;
  if (!fpTcArray || !fpDice) return -ENODEV;

  unsigned int arraySize=0;
  const AliHLTComponentBlockData* pBlock=NULL;
  const TObject* pObject=NULL;
  const TObjArray* pObjArray=NULL;
  AliHLTUInt32_t spec=0;
  if ((pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack))!=NULL) {
    if (pBlock->fSize%sizeof(AliHLTExternalTrackParam)==0) {
      arraySize=pBlock->fSize/sizeof(AliHLTExternalTrackParam);
      if (fpTcArray) {
	spec=GetSpecification(pBlock);
	fpTcArray->ExpandCreate(arraySize);
	if ((iResult=ReadFromStruct(fpTcArray, reinterpret_cast<AliHLTExternalTrackParam*>(pBlock->fPtr), arraySize))<0) {
	  HLTError("can not convert struct to TObjArray, error %d", iResult);
	} else {
	  pObjArray=fpTcArray;
	  HLTDebug("created TObjArray from C struct, specification 0x%08x", spec);
	}
      } else {
	iResult=-ENOMEM;
      }
    } else {
      HLTWarning("data block does not match the AliHLTExternalTrackParam struct size");
    }
  } else if ((pObject=GetFirstInputObject(kAliHLTDataTypeTObjArray))!=NULL) {
    pObjArray=dynamic_cast<const TObjArray*>(pObject);
    if (pObjArray) {
      spec=GetSpecification(pObject);
      HLTDebug("extracted TObjArray, specification 0x%08x", spec);
    } else {
      HLTWarning("object of data type %s is not of type TObjArray, ignoring it",
		 DataType2Text(GetDataType()).c_str());
    }
  }
  
  if (pObjArray) {
    // we already have the object array, if publishing mode is >0 the
    // block is forwarded according to the mode
    // 
    // if specification is >0 the source object is searched in the list
    // and compared with the extracted one.
    if (spec>0) {
      TObject* pSrcObject=FindObject(spec);
      TObjArray* pSrcArray=NULL;
      if (!pSrcObject) {
	HLTError("can not find source object of id 0x%08x", spec);
      } else if ((pSrcArray=dynamic_cast<TObjArray*>(pSrcObject))==NULL) {
	HLTError("type cast failed for object");
      } else if (!Compare(pObjArray, pSrcArray)) {
	HLTError("extracted array differs from source array");
	if (fVerbosity>1) {
	  for (int i=0; i<pObjArray->GetEntries(); i++) {
	    pObjArray->At(i)->Print();
	  }
	  HLTError("dump of original array");
	  for (int i=0; i<pSrcArray->GetEntries(); i++) {
	    pSrcArray->At(i)->Print();
	  }
	}
      } else {
	if (fVerbosity>0) {
	  HLTInfo("extracted array of size %d", pObjArray->GetEntries());
	}
      }
    }
  } else {
    // no array to forward, create a new one
    arraySize=fpDice->Integer(fMaxSize-fMinSize)+fMinSize;
    fpTcArray->ExpandCreate(arraySize);
    for (unsigned int track=0; track<arraySize; track++) {
      FillRandom(reinterpret_cast<AliExternalTrackParam*>(fpTcArray->At(track)));
    }
    pObjArray=fpTcArray;
  }

  switch (fMode) {
  case kPublishingOff:
    break;
  case ktclonesarray:
      // register the object if it is not a forwarded array
    if (spec==0 && pObjArray==fpTcArray) spec=Register(fpTcArray);
    if ((iResult=PushBack(const_cast<TObjArray*>(pObjArray), kAliHLTDataTypeTObjArray, spec))>=0) {
      if (fVerbosity>0) {
	HLTInfo("publishing TClonesarray, %d elements, specification 0x%08x", pObjArray->GetEntries(), spec);
      }
    }
    break;
  case kcarray:
    iResult=SerializeToStruct(pObjArray, outputPtr, capacity);
    if (iResult>=0) {
      size=iResult;
      AliHLTComponentBlockData bd;
      FillBlockData(bd);
      bd.fPtr=NULL;
      bd.fOffset=0;
      bd.fSize=size;
      bd.fDataType=kAliHLTDataTypeTrack;
      // register the object if it is not a forwarded array
      if (spec==0 && pObjArray==fpTcArray) spec=Register(fpTcArray);
      bd.fSpecification=spec;
      outputBlocks.push_back(bd);
      if (fVerbosity>0) {
	HLTInfo("publishing C array, %d elements, specification 0x%08x", bd.fSize/sizeof(AliHLTExternalTrackParam), bd.fSpecification);
      }
    }
    break;
  case ktobjarray:
    // use a TObjArray which owns the objects
    if (fpTObjArray) {
      fpTObjArray->Clear();
      int entries=pObjArray->GetEntries();
      for (int i=0; i<entries; i++) {
	fpTObjArray->Add(pObjArray->At(i));
      }

      // register the object if it is not a forwarded array
      if (spec==0 && pObjArray==fpTcArray) spec=Register(fpTObjArray);
      if ((iResult=PushBack(fpTObjArray, kAliHLTDataTypeTObjArray, spec))>=0) {
	if (fVerbosity>0) {
	  HLTInfo("publishing TObjArray, %d elements, specification 0x%08x", pObjArray->GetEntries(), spec);
	}
      }
    } else {
      HLTError("object array not initialized");
      iResult=-EFAULT;
    }
    break;
  default:
    HLTError("unknown publishing mode %d", fMode);
    iResult=-EINVAL;
  }

  if (fEventModulo>0 &&
      ((GetEventCount()+1)%fEventModulo)==0 &&
      (fRangeMultiplicator!=1.0 || fRangeOffset!=0)) {
    fMaxSize+=fRangeOffset;
    fMaxSize=(int)(fMaxSize*fRangeMultiplicator);
    fMinSize+=fRangeOffset;
    fMinSize=(int)(fMinSize*fRangeMultiplicator);
    if (fMaxSize<1) fMaxSize=1;
    if (fMinSize<1) fMinSize=1;
  }
  return iResult;
}

int AliHLTBenchExternalTrackComponent::FillRandom(AliExternalTrackParam* track, int fillCov)
{
  // see header file for class documentation
  int iResult=0;
  if (!track) return -EINVAL;
  Double_t param[5];
  Double_t covar[15];

  static TRandom* rand=NULL;
  if (!rand) {
    rand=new TRandom;
    if (!rand) return -ENOMEM;
    TDatime dt;
    rand->SetSeed(dt.Get());
  }

  param[0]=(Double_t)rand->Integer(100);
  param[1]=(Double_t)rand->Integer(100);
  param[2]=(Double_t)rand->Integer(100)/100;
  param[3]=(Double_t)rand->Integer(100)/100;
  param[4]=(Double_t)rand->Integer(100);
  for (int i=0; i<15 && i<fillCov; i++) {
    covar[i]=(Double_t)rand->Integer(1000)/1000;
  }
  track->Set((Double_t)rand->Integer(1000),(Double_t)rand->Integer(100)/100, param, covar );
  return iResult;
}

int AliHLTBenchExternalTrackComponent::SerializeToStruct(const TObjArray* pArray, AliHLTUInt8_t* buffer, unsigned int size)
{
  // see header file for class documentation
  int iResult=0;
  if (!pArray || !buffer) return -EINVAL;
  int entries=pArray->GetEntries();
  if (entries*sizeof(AliHLTExternalTrackParam)>size) return -ENOSPC;

  AliHLTLogging log;
  AliHLTExternalTrackParam* tgtParam=reinterpret_cast<AliHLTExternalTrackParam*>(buffer);
  for (int track=0; track<entries; track++, tgtParam++) {
    TObject* pObj=pArray->At(track);
    if (!pObj) {
      log.LoggingVarargs(kHLTLogError, "AliHLTBenchExternalTrackComponent", "SerializeToStruct" , __FILE__ , __LINE__ ,
			  "internal mismatch: can not get object %d from array. aborting", track);
      iResult=-ENOENT;
      break;
    }
    AliExternalTrackParam* srcParam=dynamic_cast<AliExternalTrackParam*>(pObj);
    if (!srcParam) {
      log.LoggingVarargs(kHLTLogError, "AliHLTBenchExternalTrackComponent", "SerializeToStruct" , __FILE__ , __LINE__ ,
			  "object %d has wrong type %s, expecting %s", track, pObj->Class()->GetName(), "AliExternalTrackParam");
      iResult=-ENOENT;
      break;
    }
    
    //= srcParam->GetAlpha();
    tgtParam->fX = srcParam->GetX();
    tgtParam->fY = srcParam->GetY();
    tgtParam->fZ = srcParam->GetZ();
    tgtParam->fLastX=0;
    tgtParam->fLastY=0;
    tgtParam->fLastZ=0;
    tgtParam->fSinPsi = srcParam->GetSnp();
    tgtParam->fTgl = srcParam->GetTgl();
    tgtParam->fq1Pt = srcParam->GetSigned1Pt();
    const Double_t *cov=srcParam->GetCovariance();

    for (int i=0; i<15; i++) {
      tgtParam->fC[i]=cov[i];
    }
    tgtParam->fNPoints=0;
  }

  if (iResult>=0) {
    iResult=entries*sizeof(AliHLTExternalTrackParam);
  }

  return iResult;
}

int AliHLTBenchExternalTrackComponent::ReadFromStruct(TObjArray* pTgtArray, AliHLTExternalTrackParam* pArray, unsigned int arraySize)
{
  // see header file for class documentation
  AliHLTLogging log;
  if ((unsigned)pTgtArray->GetEntries()<arraySize) {
    log.LoggingVarargs(kHLTLogError, "AliHLTBenchExternalTrackComponent", "ReadFromStruct" , __FILE__ , __LINE__ ,
		       "not enough space in target array: %d, required %d", pTgtArray->GetEntries(), arraySize);
    return -ENOSPC;
  }

  Double_t param[5];
  Double_t covar[15];

  for (unsigned int track=0; track<arraySize; track++) {
    AliExternalTrackParam* tp=dynamic_cast<AliExternalTrackParam*>(pTgtArray->At(track));
    if (tp) {
      // enable if AliExternalTrackParam allows templates
      //tp->Set(pArray[track].fX, (Float_t)0.0, &pArray[track].fY, pArray[track].fC);
      param[0]=pArray[track].fY;
      param[1]=pArray[track].fZ;
      param[2]=pArray[track].fSinPsi;
      param[3]=pArray[track].fTgl;
      param[4]=pArray[track].fq1Pt;
      for (int i=0; i<15; i++) {
	covar[i]=pArray[track].fC[i];
      }
      tp->Set((Double_t)pArray[track].fX, (Double_t)0.0, param, covar);
    } else {
      log.LoggingVarargs(kHLTLogError, "AliHLTBenchExternalTrackComponent", "ReadFromStruct" , __FILE__ , __LINE__ ,
			 "invalid object type %s", pTgtArray->At(track)->Class()->GetName());
      return -EFAULT;
    }
  }
  return arraySize;
}

AliHLTUInt32_t AliHLTBenchExternalTrackComponent::CalcChecksum(const TObjArray* pArray)
{
  // see header file for class documentation
  AliHLTUInt32_t crc=0;
  int entries=0;
  if (pArray && (entries=pArray->GetEntries())>0) {
    AliHLTLogging log;
    unsigned int bufferSize=entries*sizeof(AliHLTExternalTrackParam);
    AliHLTUInt8_t* buffer=new AliHLTUInt8_t[bufferSize];
    if (buffer && SerializeToStruct(pArray, buffer, bufferSize)>0) {
      crc=CalculateChecksum(buffer, bufferSize);
    } else if (buffer) {
      log.LoggingVarargs(kHLTLogError, "AliHLTBenchExternalTrackComponent", "SerializeToStruct" , __FILE__ , __LINE__ ,
			 "failed to serialize TObjArray");
    }
  }
  return crc;
}

bool AliHLTBenchExternalTrackComponent::Compare(const TObjArray* array1, const TObjArray* array2) 
{
  // see header file for class documentation
  if (!array1 || !array2) return false;
  int entries=array1->GetEntries();
  if (entries!=array2->GetEntries()) {
    return false;
  }

  for (int i=0; i<entries; i++) {
    TObject* object1=array1->At(i);
    TObject* object2=array2->At(i);
    if (!object1 || !object2) return false;

    AliExternalTrackParam* param1=dynamic_cast<AliExternalTrackParam*>(object1);
    AliExternalTrackParam* param2=dynamic_cast<AliExternalTrackParam*>(object2);
    if (!param1 || !param2) return false;

    if (TMath::Abs(param1->GetX()-param2->GetX())>0.0001) return false;
    if (TMath::Abs(param1->GetY()-param2->GetY())>0.0001) return false;
    if (TMath::Abs(param1->GetZ()-param2->GetZ())>0.0001) return false;
    if (TMath::Abs(param1->GetSnp()-param2->GetSnp())>0.0001) return false;
    if (TMath::Abs(param1->GetTgl()-param2->GetTgl())>0.0001) return false;
    if (TMath::Abs(param1->GetSigned1Pt()-param2->GetSigned1Pt())>0.0001) return false;

    if (TMath::Abs(param1->GetSigmaY2()-param2->GetSigmaY2())>0.0001) return false;
    if (TMath::Abs(param1->GetSigmaZY()-param2->GetSigmaZY())>0.0001) return false;
    if (TMath::Abs(param1->GetSigmaZ2()-param2->GetSigmaZ2())>0.0001) return false;
    if (TMath::Abs(param1->GetSigmaSnpY()-param2->GetSigmaSnpY())>0.0001) return false;
    if (TMath::Abs(param1->GetSigmaSnpZ()-param2->GetSigmaSnpZ())>0.0001) return false;
    if (TMath::Abs(param1->GetSigmaSnp2()-param2->GetSigmaSnp2())>0.0001) return false;
    if (TMath::Abs(param1->GetSigmaTglY()-param2->GetSigmaTglY())>0.0001) return false;
    if (TMath::Abs(param1->GetSigmaTglZ()-param2->GetSigmaTglZ())>0.0001) return false;
    if (TMath::Abs(param1->GetSigmaTglSnp()-param2->GetSigmaTglSnp())>0.0001) return false;
    if (TMath::Abs(param1->GetSigmaTgl2()-param2->GetSigmaTgl2())>0.0001) return false;
    if (TMath::Abs(param1->GetSigma1PtY()-param2->GetSigma1PtY())>0.0001) return false;
    if (TMath::Abs(param1->GetSigma1PtZ()-param2->GetSigma1PtZ())>0.0001) return false;
    if (TMath::Abs(param1->GetSigma1PtSnp()-param2->GetSigma1PtSnp())>0.0001) return false;
    if (TMath::Abs(param1->GetSigma1PtTgl()-param2->GetSigma1PtTgl())>0.0001) return false;
    if (TMath::Abs(param1->GetSigma1Pt2()-param2->GetSigma1Pt2())>0.0001) return false;
  }

  return true;
}

TList* AliHLTBenchExternalTrackComponent::fgpRegistry=NULL;

TObject* AliHLTBenchExternalTrackComponent::FindObject(AliHLTUInt32_t id)
{
  // see header file for class documentation
  if (!fgpRegistry) return NULL;

  TIter iter(fgpRegistry);
  TObject* pObj=NULL;
  while ((pObj=iter.Next())) {
    if (pObj->GetUniqueID()==id) break;
  }
  return pObj;
}

AliHLTUInt32_t AliHLTBenchExternalTrackComponent::Register(TObject* pObject)
{
  // see header file for class documentation
  if (fDisableRegistry) return 0;

  if (!fgpRegistry) {
    fgpRegistry=new TList;
  }
  if (!fgpRegistry) return 0;

  TObject* pExist=fgpRegistry->FindObject(pObject);
  if (pExist) return pExist->GetUniqueID();

  AliHLTUInt32_t id=AliHLTComponent::CalculateChecksum(reinterpret_cast<AliHLTUInt8_t*>(&pObject), sizeof(TObject*));
  pExist=FindObject(id);
  assert(pExist==NULL);
  pObject->SetUniqueID(id);
  fgpRegistry->Add(pObject);
  HLTInfo("adding object %p with id 0x%08x to registry", pObject, id);
  return id;
}

int AliHLTBenchExternalTrackComponent::Unregister(TObject* pObject)
{
  // see header file for class documentation
  if (!fgpRegistry) return 0;

  if (fgpRegistry->Remove(pObject)==pObject) return 0;

  return -ENOENT;
}
