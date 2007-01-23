/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
//  marian.ivanov@cern.ch
//
//  ------------------------------------------------------------------------------------------------
//  TTreeStream
//  Standard stream (cout) like input for the tree
//  Run and see TTreeStreamer::Test() - to see TTreeStreamer functionality
//  ------------------------------------------------------------------------------------------------  
//
//  -------------------------------------------------------------------------------------------------
//  TTreeSRedirector
//  Redirect file to  different TTreeStreams  
//  Run and see   TTreeSRedirector::Test() as an example of TTreeSRedirectorer functionality 
// 

#include <TClass.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TTreeStream.h>

ClassImp(TTreeDataElement)
ClassImp(TTreeStream)
ClassImp(TTreeSRedirector)



void TTreeStream::Test()
{
  //
  // 
  TFile *ftest = new TFile("teststreamer.root","recreate");
  if (!ftest) ftest = new TFile("teststreamer.root","new");
  //
  //create to streems Tree1 and Tree2
  TTreeStream stream1("Tree1");
  TTreeStream stream2("Tree2");
  //
  Char_t ch='s';
  Float_t f=3.;
  Float_t f2=1;
  TObject *po  = new TObject;
  TObject *po2 = new TObject;
  for (Int_t i=0;i<100000;i++) {
    f=i*100;
    po->SetUniqueID(i);
    po2->SetUniqueID(i*100);
    ch=i%120;
    //
    //    Stream the data
    //    The data layout of stream is defined during first invocation of streamer.
    //    Endl is the trigger which define the end of structure.
    // 
    //    The name of branch can be specified using strings with = at the the end
    //    if string is not specified automatic convention is u (sed B0, B1, ...Bn)
    stream1<<"i="<<i<<"ch="<<ch<<"f="<<f<<"po="<<po<<"\n";
    f  = 1./(100.1+i);
    f2 = -f;     
    //3.) just another example - we can fill the same tree with different objects
    //
    stream2<<f<<po<<"\n";
    stream2<<f2<<po2<<"\n";
  }
  //
  //4.) Close the streeamers (Write the streamed tree's to the file) and close the corresponding file.
  //
  stream1.Close();
  stream2.Close();
  ftest->Close();
  delete ftest;
  //
  //5.) and now see results  in file tteststreamer.root
}



void TTreeSRedirector::Test()
{
  //
  //Example test function to show functionality of TTreeSRedirector
  //
  //
  //1.)create the  redirector associated with file (testredirector.root)
  //
  //
  TTreeSRedirector *pmistream= new TTreeSRedirector("testredirector.root");
  TTreeSRedirector &mistream = *pmistream;
  Char_t ch='s';
  Float_t f=3.;
  Float_t f2=1;
  TObject *po  = new TObject;
  TObject *po2 = new TObject;
  for (Int_t i=0;i<100000;i++) {
    f=i*100;
    po->SetUniqueID(i);
    po2->SetUniqueID(i*100);
    ch=i%120;
    //
    //2.) create the tree with identifier specified by first argument
    //                                layout specified by sequence of arguments
    //                                Tree identifier has to be specified as first argument !!! 
    //    if the tree and layout was already defined the consistency if layout is checked
    //                                if the data are consisten fill given tree 
    //    the name of branch can be specified using strings with = at the the end
    //    if string is not specified use automatic convention  B0, B1, ...Bn
    mistream<<"TreeIdentifier"<<"i="<<i<<"ch="<<ch<<"f="<<f<<"po="<<po<<"\n";
    f  = 1./(100.1+i);
    f2 = -f; 
    
    //3.) just another example - we can fill the same tree with different objects
    //
    mistream<<"TreeK"<<f<<po<<"\n";
    mistream<<"TreeK"<<f2<<po2<<"\n";
  }
  //
  //4.) write the streamed tree's to the file and close the corresponding file in destructor
  //
  delete pmistream;
  //
  //5.) and now see results in file testredirector.root 
}


TTreeSRedirector::TTreeSRedirector(const char *fname) :
  fFile(new TFile(fname,"recreate")),
  fDataLayouts(0)
{
  //
  // Constructor
  //
  if (!fFile){
    fFile = new TFile(fname,"new");
  }
}

TTreeSRedirector::~TTreeSRedirector()
{
  //
  // Destructor
  //
  Close();       //write the tree to the selected file
  fFile->Close();
  delete fFile;
}
void TTreeSRedirector::StoreObject(TObject* object){
  //
  //
  //
  TFile * backup = gFile;
  fFile->cd();
  object->Write();
  if (backup) backup->cd();
}



TTreeStream  & TTreeSRedirector::operator<<(Int_t id)
{
  //
  // return reference to the data layout with given identifier
  // if not existing - creates new
  if (!fDataLayouts) fDataLayouts = new TObjArray(10000);
  TTreeStream *clayout=0;
  Int_t entries = fDataLayouts->GetEntriesFast();
  for (Int_t i=0;i<entries;i++){
    TTreeStream * layout = (TTreeStream*)fDataLayouts->At(i);
    if (!layout) continue;
    if (layout->fId==id) {
      clayout = layout;
      break;
    }
  }
  if (!clayout){
    fFile->cd();
    char chname[100];
    sprintf(chname,"Tree%d",id);
    clayout = new TTreeStream(chname);
    clayout->fId=id;
    fDataLayouts->AddAt(clayout,entries);
  }
  return *clayout;
}


TTreeStream  & TTreeSRedirector::operator<<(const char* name)
{
  //
  // return reference to the data layout with given identifier
  // if not existing - creates new
  if (!fDataLayouts) fDataLayouts = new TObjArray(10000);
  TTreeStream *clayout=(TTreeStream*)fDataLayouts->FindObject(name);
  Int_t entries = fDataLayouts->GetEntriesFast();

  if (!clayout){
    fFile->cd();
    clayout = new TTreeStream(name);
    clayout->fId=-1;
    clayout->SetName(name);
    fDataLayouts->AddAt(clayout,entries);    
  }
  return *clayout;
}




void TTreeSRedirector::Close(){
  //
  //
  TFile * backup = gFile;
  fFile->cd();
  if (fDataLayouts){
    Int_t entries = fDataLayouts->GetEntriesFast();
    for (Int_t i=0;i<entries;i++){
      TTreeStream * layout = (TTreeStream*)fDataLayouts->At(i);
      if (layout){
	if (layout->fTree) layout->fTree->Write(layout->GetName());
      }
    }
    delete fDataLayouts;
    fDataLayouts=0;
  }
  if (backup) backup->cd();
}



//-------------------------------------------------------------
TTreeDataElement:: TTreeDataElement(Char_t type) :
  TNamed(),
  fType(type),
  fDType(0),
  fClass(0),
  fPointer(0)
{
  //
  //
  //
}

TTreeDataElement:: TTreeDataElement(TDataType* type) :
  TNamed(),
  fType(0),
  fDType(type),
  fClass(0),
  fPointer(0)
{
  //
  //
  //
}

TTreeDataElement:: TTreeDataElement(TClass* cl) :
  TNamed(),
  fType(0),
  fDType(0),
  fClass(cl),
  fPointer(0)
{
  //
  //
  //
}

//-------------------------------------------------------------------
TTreeStream::TTreeStream(const char *treename):
  TNamed(treename,treename),
  fElements(0),
  fBranches(0),
  fTree(new TTree(treename, treename)),
  fCurrentIndex(0),
  fId(0),
  fNextName(),
  fNextNameCounter(),
  fStatus(0)
{
  //
  // Standard ctor
  //
}

TTreeStream::~TTreeStream()
{
  //
  // Class dtor
  //
  fElements->Delete();
  fBranches->Clear();
  delete fElements;
  delete fBranches;
}

void TTreeStream::Close()
{
  //
  // Flush data to disk and close
  //
  fTree->Write();
}

Int_t TTreeStream::CheckIn(Char_t type, void *pointer)
{
  //
  // Insert object of given type
  //
  if (!fElements) fElements = new TObjArray(1000);
  TTreeDataElement* element = (TTreeDataElement*)fElements->At(fCurrentIndex);
  if (!element) {
    element = new TTreeDataElement(type);
    //
    char name[1000];
    if (fNextName.Length()>0){
      if (fNextNameCounter==0){
	sprintf(name,"%s",(const char*)fNextName);
      }
      if (fNextNameCounter>0){
	sprintf(name,"%s%d",(const char*)fNextName,fNextNameCounter);
      }      
    }
    else{
      sprintf(name,"B%d.",fCurrentIndex);
    }
    element->SetName(name);
    //
    element->SetPointer(pointer);
    fElements->AddAt(element,fCurrentIndex);
    fCurrentIndex++;
    return 0; //new element added
  }
  if (element->GetType()!=type){
    fStatus++;
    return 1; //mismatched data element
  }
  element->SetPointer(pointer);
  fCurrentIndex++;
  return 0;
}

Int_t TTreeStream::CheckIn(TObject *o){
  //
  // Insert TObject
  //
  if (!o) return 0;
  if (!fElements) fElements = new TObjArray(1000);
  TTreeDataElement* element = (TTreeDataElement*)fElements->At(fCurrentIndex);
  if (!element) {
    element = new TTreeDataElement(o->IsA());
    //
    char name[1000];
    if (fNextName.Length()>0){
      if (fNextNameCounter==0){
	sprintf(name,"%s",(const char*)fNextName);
      }
      if (fNextNameCounter>0){
	sprintf(name,"%s%d",(const char*)fNextName,fNextNameCounter);
      }      
    }
    else{
      sprintf(name,"B%d",fCurrentIndex);
    }
    element->SetName(name);

    element->SetPointer(o);
    fElements->AddAt(element,fCurrentIndex);
    fCurrentIndex++;
    return 0; //new element added
  }
  if (element->fClass!=o->IsA()){
    fStatus++;
    return 1; //mismatched data element
  }
  element->SetPointer(o);
  fCurrentIndex++;
  return 0;  
}

void TTreeStream::BuildTree(){
  //
  // Build the Tree
  //
  if (fTree->GetEntries()>0) return;
  fTree = new TTree(GetName(),GetName());
  Int_t entries = fElements->GetEntriesFast();  
  fBranches = new TObjArray(entries);
  
  for (Int_t i=0;i<entries;i++){
    //
    TTreeDataElement* element = (TTreeDataElement*)fElements->At(i);
    char bname1[1000];
    if (element->GetName()[0]==0){
      sprintf(bname1,"B%d",i);
    }
    else{
      sprintf(bname1,element->GetName());
    }
    if (element->fClass){
      if (element->fClass->GetBaseClass("TClonesArray")){
	TBranch * br = fTree->Branch(bname1,element->fClass->GetName(),&(element->fPointer));
	fBranches->AddAt(br,i);
      }else
	{
	  TBranch * br = fTree->Branch(bname1,element->fClass->GetName(),&(element->fPointer));
	  fBranches->AddAt(br,i);
	}
    }
    if (element->GetType()>0){
      char bname2[1000];
      sprintf(bname2,"B%d/%c",i,element->GetType());
      TBranch * br = fTree->Branch(bname1,element->fPointer,bname2);
      fBranches->AddAt(br,i);
    }
  }
}

void TTreeStream::Fill(){
  //
  // Fill the tree
  //
  if (fTree) { 
    Int_t entries=fElements->GetEntriesFast();
    if (entries>fTree->GetNbranches()) BuildTree();
    for (Int_t i=0;i<entries;i++){    
      TTreeDataElement* el  = (TTreeDataElement*)fElements->At(i);
      if (!el) continue;
      if (!el->GetType()) continue;
      TBranch      * br  = (TBranch*)fBranches->At(i);
      if (br &&el){
	if (el->GetType())  br->SetAddress(el->fPointer);
      }
    }
    if (fStatus==0) fTree->Fill(); //fill only in case of non conflicts
    fStatus=0;
  }
}

TTreeStream & TTreeStream::Endl()
{
  //
  // Perform pseudo endl operation
  //
  if (fTree->GetNbranches()==0) BuildTree();
  Fill();
  fStatus =0;
  fCurrentIndex=0;
  return *this;
}


TTreeStream  &TTreeStream::operator<<(Char_t *name)
{
  //
  // Endl 
  //
  if (name[0]=='\n'){
    return Endl();
  }
  //
  //if tree was already defined ignore
  if (fTree->GetEntries()>0) return *this;
  //check branch name if tree was not 
  //
  Int_t last=0;
  for (last=0;;last++){
    if (name[last]==0) break;    
  }
  
  if (last>0&&name[last-1]=='='){
    fNextName = name;
    fNextName[last-1]=0;
    fNextNameCounter=0;
  }
  return *this;
}

