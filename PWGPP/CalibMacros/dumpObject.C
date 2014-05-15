/*
  //
  // Function to dump recursivally in human readable format, content of the objects + example use case
  // 
  //
  .L  $mcProd/dumpObject.C+
  ExampleUse(); > AliTPCClusterParam.dump
  AliTPCRecoParam param; 
  DumpObjectRecursive(&param); >> AliTPCRecoParam.dump

  Use case examples 
  1.) compare oontent of alignent OCDB files for differnt yers
  2.) compare ClusterParam for different periods
   
=================================================================================================================
  // 1.)
  // Compare alignment example:
  // Compare TPC alignemnt 2013 and 2010
  //
  DumpOCDBFile("/cvmfs/alice.gsi.de/alice/data/2013/OCDB/TPC/Align/Data/Run0_999999999_v1_s0.root","TPCalign2013.dump",1);
  DumpOCDBFile("/cvmfs/alice.gsi.de/alice/data/2010/OCDB/TPC/Align/Data/Run0_999999999_v1_s0.root","TPCalign2010.dump",1);
  diff  TPCalign2013.dump TPCalign2010.dump > TPCalign2013_TPCalign2010.diff
  //
  //    
=================================================================================================================
//  2.) 
  // Compare CluterParam OCDB etry
  //
  DumpOCDBFile("/cvmfs/alice.gsi.de/alice/data/2010/OCDB/TPC/Calib/ClusterParam/Run131541_999999999_v2_s0.root","2010_TPC_Calib_ClusterParam_Run131541_999999999_v2_s0.dump",1);
  DumpOCDBFile("/cvmfs/alice.gsi.de/alice/data/2010/OCDB/TPC/Calib/ClusterParam/Run0_999999999_v1_s0.root","2010_TPC_Calib_ClusterParam_Run0_999999999_v1_s0.dump",1);
  DumpOCDBFile("/cvmfs/alice.gsi.de/alice/data/2013/OCDB/TPC/Calib/ClusterParam/Run0_999999999_v1_s0.root","2013_TPC_Calib_ClusterParam_Run0_999999999_v1_s0.dump",1);
  
  diff 2010_TPC_Calib_ClusterParam_Run131541_999999999_v2_s0.dump 2010_TPC_Calib_ClusterParam_Run0_999999999_v1_s0.dump


*/
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include <TVectorF.h>
#include <TLinearFitter.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TProfile2D.h>
#include <TVectorD.h>
#include <TObjArray.h>
#include "THnBase.h"
#include "TRealData.h"
#include "TDataMember.h"
#include "TClass.h"
#include "AliCDBEntry.h"
#include "TObjArray.h"
#include "TNamed.h"
#include "TROOT.h"
#endif

void DumpObjectRecursive(TObject *obj);
void DumpObjectRecursive(TObject *obj, TString prefix, Int_t &counterRec);

void ExampleUse(){
  //
  //
  //
  TFile *f = TFile::Open("/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual/TPC/Calib/ClusterParam/Run127712_130850_v3_s0.root");
  //TFile *f = TFile::Open("./OCDB_NoClustersBelowThreshold/TPC/Calib/RecoParam/Run127712_130850_v0_s0.root");
  AliCDBEntry * entry = (AliCDBEntry*)f->Get("AliCDBEntry");
  TObject  *obj = ( TObject*)entry->GetObject();
  DumpObjectRecursive(obj);
  //
  //
}

void DumpOCDBFile(const char *finput , const char *foutput, Bool_t dumpMetaData, Bool_t xml){
  //
  //  
  //  DumpOCDBFile("$ALICE_ROOT/OCDB/ITS/Align/Data/Run0_999999999_v0_s0.root", "ITS_Align_Data_Run0_999999999_v0_s0.dump")
  //
  if (finput==0) return ;
  TFile *falignITS  = TFile::Open(finput);
  AliCDBEntry *entry  = (AliCDBEntry*)falignITS->Get("AliCDBEntry");
  if (!entry) return; 
  TObject *obj = ((AliCDBEntry*)falignITS->Get("AliCDBEntry"))->GetObject();  

  //
  if (!xml){
    if (dumpMetaData) gROOT->ProcessLine(TString::Format("((TObject*)%p)->Dump(); >%s",entry, foutput).Data());
    if (!obj) return;
    gROOT->ProcessLine(TString::Format("DumpObjectRecursive((TObject*)%p); >>%s",obj, foutput).Data());
  }
  if (xml){
    TFile * f = TFile::Open(TString::Format("%s.xml",foutput).Data(),"recreate");
    if (dumpMetaData) entry->Write("AliCDBEntry");
    else obj->Write("AliCDBEntry");
    f->Close();
  }
}



void DumpObjectRecursive(TObject *obj){
  //
  //
  //
  Int_t counterRec=0;
  printf("==> Dumping object at: %p, name=%s, class=%s)\n", obj, obj->GetName(), (obj->IsA()->GetName()));
  DumpObjectRecursive(obj, TString(obj->IsA()->GetName())+".",counterRec);
}
 
//
//
//
void DumpObjectRecursive(TObject *obj, TString prefix, Int_t &counterRec){
  //
  // Recursive dump of the TObject
  // Dump all basic types and follow pointers to the objects
  // current limitation:
  //    a.) clases and structures not derived from TObject not followed (to fix)
  //    b.) dynamic arrays not followed
  //    c.) std maps,array ....  not followed
  //    
  //
  if (!obj) return;
  //
  // Special case of Collection classes
  //
  if (obj->IsA()->InheritsFrom(TCollection::Class())) {
    TIter myiter((TCollection*)obj);
    TObject  *arObject=0;
    Int_t counter=0;
    while ((arObject = (TObject*)myiter.Next())) {
      TString prefixArr = TString::Format("%s[%d]",prefix.Data(),counter);
      DumpObjectRecursive(arObject,prefixArr,counterRec);
      counter++;
    } 
    counterRec++;
    return;
  }

  TClass * cl = obj->IsA();
  if (!(cl->GetListOfRealData())) cl->BuildRealData();
  TRealData* rd = 0;
  TIter next(cl->GetListOfRealData());  
  while ((rd = (TRealData*) next())) {
    counterRec++;
    TDataMember* dm = rd->GetDataMember();
    TDataType* dtype = dm->GetDataType();
    Int_t offset = rd->GetThisOffset();
    char* pointer = ((char*) obj) + offset;
    
    if (dm->IsaPointer()) {
      // We have a pointer to an object or a pointer to an array of basic types.
      TClass* clobj = 0;
      if (!dm->IsBasic()) {
	clobj = TClass::GetClass(dm->GetTypeName());
      }
      if (clobj) {
	// We have a pointer to an object.
	//
	if (!clobj->InheritsFrom(TObject::Class())) {
	  // It must be a TObject object.
	  continue; 
	}
	char** apointer = (char**) pointer;
	TObject* robj = (TObject*) *apointer;
	//	
	if(!robj)
	  printf("M:%s%s\n",prefix.Data(),dm->GetName()); // Missing - 0 pointer
	else{
	  printf("T:%s\t%s%s\n", clobj->GetName(),prefix.Data(), dm->GetName());
	  TString prefixNew=prefix;
	  prefixNew+=dm->GetName();
	  prefixNew+=".";
	  if (robj!=obj) DumpObjectRecursive(robj,prefixNew,counterRec);  // trivial check 
	  if (robj==obj){
	    printf("R:%s\t%s%s\n",clobj->GetName(),prefix.Data(), dm->GetName());
	  }
	}
      }
    } else if (dm->IsBasic()) {
      //
      // Basic data type
      //
      const char* index = dm->GetArrayIndex();
      if (dm->GetArrayDim()==0){
	printf("B:\t%s%s\t%s\n", prefix.Data(),rd->GetName(), dtype->AsString(pointer));
      }
      //
      // Basic array - fixed length
      //
      //      if (dm->GetArrayDim()>0 && strlen(index) != 0){
      if (dm->GetArrayDim()>0 ){
	printf("A:\t%s%s\t",prefix.Data(),rd->GetName());
	Int_t counter=0;
	for  (Int_t idim=0; idim<dm->GetArrayDim(); idim++){
	  //printf("A:%d\t%d\n", dm->GetArrayDim(),dm->GetMaxIndex(idim));
	  for (Int_t j=0; j<dm->GetMaxIndex(idim); j++){
	    printf("%s\t",dtype->AsString(pointer+dm->GetUnitSize()*counter));
	    counter++;
	    if (counter%5==0) printf("\nA:\t%s%s\t",prefix.Data(),rd->GetName());
	  }
	}
	printf("\n");
      }
      //
      // Basic array - dynamic length
      //
      if (dm->GetArrayDim()>0 && strlen(index) != 0){
	//
	// Dump first only for the moment
	//  
	printf("B:\t%s%s\t%s\n",prefix.Data(),rd->GetName(), dtype->AsString(pointer));
      }
    } else {
    }
  }
}  

//
// Small checks to test the TRealData and TDataType
//



void DumpDataSimple(){
  //
  // Dump example for elenatr data types 
  //
  TObject *obj = new TVectorD(20);
  TClass * cl = obj->IsA();
  if (!cl->GetListOfRealData()) cl->BuildRealData();
  //
  TRealData* rd = 0;
  rd = (TRealData*)(cl->GetListOfRealData()->FindObject("fNrows"));
  TDataMember* dm = rd->GetDataMember();
  TDataType* dtype = dm->GetDataType();
  //
  Int_t offset = rd->GetThisOffset();
  char* pointer = ((char*) obj) + offset;
  printf("%s\n",dtype->AsString(pointer));
  
}

void DumpDataArray(){
  //
  // print array example
  // 
  TObject *obj = new TVectorD(20);
  TClass * cl = obj->IsA();
  if (!cl->GetListOfRealData()) cl->BuildRealData();
  TRealData* rd = 0;
  rd = (TRealData*)(cl->GetListOfRealData()->FindObject("*fElements"));
  TDataMember* dm = rd->GetDataMember();
  TDataType* dtype = dm->GetDataType();
  dtype->Print();
  //
  Int_t offset = rd->GetThisOffset();
  char* pointer = ((char*) obj) + offset; 
  printf("%s\n",dtype->AsString(pointer));
}

void DumpTObjectArray(){
  //
  //
  //
  TObjArray *array = new TObjArray(10);
  for (Int_t i=0; i<10; i++) array->AddLast(new TNamed(Form("n%d",i), Form("n%d",i)));  
  DumpObjectRecursive(array);
  //
  //
  TObject *obj = array;
  TClass * cl = obj->IsA();
  if (!cl->GetListOfRealData()) cl->BuildRealData();
  TRealData* rd = 0;
  rd = (TRealData*)(cl->GetListOfRealData()->FindObject("*fCont"));
  TDataMember* dm = rd->GetDataMember();
  TDataType* dtype = dm->GetDataType();
  //
  Int_t offset = rd->GetThisOffset();
  char* pointer = ((char*) obj) + offset;
  char** apointer = (char**) pointer;
  //we have pointer to pointer here
  TObject** ppobj = (TObject**) *apointer;
  (*ppobj)->Print();
  //
  TIter myiter(array);
  TObject  *arObject; 
  dtype->Print();
  while ((arObject = (TObject*)myiter.Next())) {
    DumpObjectRecursive(arObject);
  } 


}
