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

/*
$Log$
*/
#include "TROOT.h"
#include "AliArrayBranch.h"
#include "TFile.h"
#include "TTree.h" 
#include "TBasket.h"
#include "TClass.h"
#include "TRealData.h"
#include "TDataType.h"
#include "TDataMember.h"

#include "TBranch.h"
#include "TBranchClones.h"
#include "TLeaf.h"
#include "TLeafB.h"
#include "TLeafC.h"
#include "TLeafF.h"
#include "TLeafD.h"
#include "TLeafI.h"
#include "TLeafS.h"
#include "TLeafObject.h"

#include "AliObjectArray.h"
#include "AliDataType.h"

//-----------------------------------------------------
// A Branch for the case of an array of clone objects. 
//-----------------------------------------------------

//*KEND.

R__EXTERN TTree *gTree;

ClassImp(AliArraySubBranch)
ClassImp(AliArrayBranch)
ClassImp(AliObjectBranch)
ClassImp(AliTree)



Int_t AliArraySubBranch::GetEntryExport(Int_t entry, Int_t getall, AliObjectArray *list, Int_t nentries)
{
//*-*-*-*-*-*Read all leaves of entry and return total number of bytes*-*-*
//*-* export buffers to real objects in the AliObjectArray list.
//*-*

   if (TestBit(kDoNotProcess)) return 0;
   if (fReadEntry == entry) return 1;
   if (entry < 0 || entry >= fEntryNumber) return 0;
   Int_t nbytes;
   Int_t first  = fBasketEntry[fReadBasket];
   Int_t last;
   if (fReadBasket == fWriteBasket) last = fEntryNumber - 1;
   else                             last = fBasketEntry[fReadBasket+1] - 1;
//
//      Are we still in the same ReadBasket?
   if (entry < first || entry > last) {
      fReadBasket = TMath::BinarySearch(fWriteBasket+1, fBasketEntry, entry);
      first       = fBasketEntry[fReadBasket];
   }

//     We have found the basket containing this entry.
//     make sure basket buffers are in memory.
   TBasket *basket = GetBasket(fReadBasket);
   if (!basket) return 0;
   TBuffer *buf    = basket->GetBufferRef();
//     Set entry offset in buffer and read data from all leaves
   if (!buf->IsReading()) {
      basket->SetReadMode();
   }
//   Int_t bufbegin = basket->GetEntryPointer(entry-first);
   Int_t bufbegin;
   Int_t *entryOffset = basket->GetEntryOffset();
   if (entryOffset) bufbegin = entryOffset[entry-first];
   else             bufbegin = basket->GetKeylen() + (entry-first)*basket->GetNevBufSize();
   buf->SetBufferOffset(bufbegin);

   TLeaf *leaf = (TLeaf*)fLeaves.UncheckedAt(0);
   //   leaf->ReadBasketExport(*buf,list,nentries);  //!!! MI
   ReadBasketExport(*buf,leaf, list,nentries);
   nbytes = buf->Length() - bufbegin;
   fReadEntry = entry;

   return nbytes;
}

void AliArraySubBranch::ReadBasketExport(TBuffer &b, TLeaf *leaf, AliObjectArray *list, Int_t n)
{
  //
  // 
  Int_t len    = leaf->GetLenStatic();
  Int_t offset = leaf->GetOffset();
  void *value  = leaf->GetValuePointer();

  //8bit integer
  if (leaf->IsA()==TLeafB::Class()){   
    Int_t j = 0;
    for (Int_t i=0;i<n;i++) {
      memcpy((char*)list->UncheckedAt(i) + offset,&((Char_t*)value)[j], len);
      j += len;
    } 
  } 
  //variable length string.
  if (leaf->IsA()==TLeafC::Class()){  
    UChar_t len;
    b >> len;
    if (len) {
      if (len >= len) len = len-1;
      b.ReadFastArray((Char_t*)value,len);
      ((Char_t*)value)[len] = 0;
    } else {
      value = 0;
    }    
    Int_t j = 0;
    for (Int_t i=0;i<n;i++) {
      memcpy((char*)list->UncheckedAt(i) + offset,&((Char_t*)value)[j], 1);
      j += len;
    }
  }
  //double
  if (leaf->IsA()==TLeafD::Class()){   
    b.ReadFastArray(((Double_t*)value),n*len);    
    Int_t j = 0;
    for (Int_t i=0;i<n;i++) {
      memcpy((char*)list->UncheckedAt(i) + offset,&((Double_t*)value)[j], 8*len);
      j += len;
   }
  }
  //float
  if (leaf->IsA()==TLeafF::Class()){   
    if (n*len == 1) {
      b >> ((Float_t*)value)[0];
    } else {
      b.ReadFastArray(((Float_t*)value),n*len);
    }
    
    Float_t *val = (Float_t*)value;
    for (Int_t i=0;i<n;i++) {
      char *first = (char*)list->UncheckedAt(i);
      Float_t *ff = (Float_t*)&first[offset];
      for (Int_t j=0;j<len;j++) {
	ff[j] = val[j];
      }
      val += len;
    }
    return;
  }
  //int2
  if (leaf->IsA()==TLeafS::Class()){       
    if (n*len == 1) {
      b >> ((Short_t*)value)[0];
    } else {
      b.ReadFastArray(((Short_t*)value),n*len);
    }
    Short_t *val = (Short_t*)value;
    for (Int_t i=0;i<n;i++) {
      char *first = (char*)list->UncheckedAt(i);
      Short_t *ii = (Short_t*)&first[offset];
      for (Int_t j=0;j<len;j++) {
	ii[j] = val[j];
      }
      val += len;
    }    
    return;
  }     
  //int4
  if (leaf->IsA()==TLeafI::Class()){       
    if (n*len == 1) {
      b >> ((Int_t*)value)[0];
    } else {
      b.ReadFastArray(((Int_t*)value),n*len);
    }
    Int_t *val = (Int_t*)value;
    for (Int_t i=0;i<n;i++) {
      char *first = (char*)list->UncheckedAt(i);
      Int_t *ii = (Int_t*)&first[offset];
      for (Int_t j=0;j<len;j++) {
	ii[j] = val[j];
      }
      val += len;
    }    
    return;
  }      
}


//______________________________________________________________________________
AliArrayBranch::AliArrayBranch(): TBranch()
{
//*-*-*-*-*-*Default constructor for BranchClones*-*-*-*-*-*-*-*-*-*
//*-*        ====================================

   fList        = 0;
   fRead        = 0;
   fN           = 0;
   fNdataMax    = 0;
   fBranchCount = 0;
}


//______________________________________________________________________________
AliArrayBranch::AliArrayBranch(const Text_t *name, void *pointer, TTree * tree,  Int_t basketsize, Int_t compress)
    :TBranch()
{
//*-*-*-*-*-*-*-*-*-*-*-*-*Create a BranchClones*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                      =====================
//
   char leaflist[80];
   char branchname[80];
   char branchcount[64];
   fTree       = tree;
   gTree       = tree; // MI because some bug in ROOT I didn't obtain proper gTree
   // it is necesary to set gTree because oder subranchces  defined below need it
   SetName(name);
   if (compress == -1) {
      TFile *bfile = fTree->GetDirectory()->GetFile();
      if (bfile) compress = bfile->GetCompressionLevel();
   }
   char *cpointer  = (char*)pointer;
   char **ppointer = (char**)(cpointer);
   fList     = (AliObjectArray*)(*ppointer);
   fAddress  = cpointer;
   fRead     = 0;
   fN        = 0;
   fNdataMax = 0;

   AliClassInfo *clinfo = fList->GetClassInfo();
   if (!clinfo) return;
   fClassName = clinfo->GetName();   
   
//*-*- Create a branch to store the array count
   if (basketsize < 100) basketsize = 100;
   sprintf(leaflist,"%s_/I",name);
   sprintf(branchcount,"%s_",name);
   fBranchCount = new TBranch(branchcount,&fN,leaflist,basketsize);
   fBranchCount->SetBit(kIsClone);
   TLeaf *leafcount = (TLeaf*)fBranchCount->GetListOfLeaves()->UncheckedAt(0);

//*-*-  Create the first basket

   fDirectory  = fTree->GetDirectory();
   fFileName   = "";

   TBasket *basket = new TBasket(branchcount,fTree->GetName(),this);
   fBaskets.Add(basket);

//*-*- Loop on all public data members of the class and its base classes

   TClass *cl = fList->GetClass();
   if (cl){
     if (!cl->GetListOfRealData())  cl->BuildRealData();
     
     const char *itype = 0;
     TRealData *rd;
     TIter      next(cl->GetListOfRealData());
     while ((rd = (TRealData *) next())) {
       TDataMember *member = rd->GetDataMember();
       if (!member->IsBasic()) {
         Warning("BranchClones","Cannot process member:%s",member->GetName());
         continue;
       }
       if (!member->IsPersistent()) continue; //do not process members with a ! as the first
       // character in the comment field
       TDataType *membertype = member->GetDataType();
       Int_t type = membertype->GetType();
       if (type == 0) {
         Warning("BranchClones","Cannot process member:%s",member->GetName());
         continue;
       }
       if (type == 1)  itype = "B";
       if (type == 11) itype = "b";
       if (type == 3)  itype = "I";
       if (type == 5)  itype = "F";
       if (type == 8)  itype = "D";
       if (type == 13) itype = "i";
       if (type == 2)  itype = "S";
       if (type == 12) itype = "s";
       
       
       Int_t arraydim = member->GetArrayDim();
       if (arraydim!=1){
	 //   OLD Version 
	 sprintf(leaflist,"%s[%s]/%s",member->GetName(),branchcount,itype);
	 Int_t comp = compress;
	 if (type == 5) comp--;
	 sprintf(branchname,"%s.%s",name,rd->GetName());
	 TBranch *branch  = new AliArraySubBranch(branchname,this,leaflist,basketsize,comp);
	 branch->SetBit(kIsClone);
	 TObjArray *leaves = branch->GetListOfLeaves();
	 TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(0);
	 leaf->SetOffset(rd->GetThisOffset());
	 leaf->SetLeafCount(leafcount);
	 Int_t arraydim = member->GetArrayDim();
	 if (arraydim) {
	   Int_t maxindex = member->GetMaxIndex(arraydim-1);
	   leaf->SetLen(maxindex);
	 }
	 fBranches.Add(branch);
       }
       else
	 for (Int_t i=0;i< member->GetMaxIndex(0);i++){
	   const char * dmname = member->GetName() ;
	   char  bname[200];
	   Int_t j=0;
	   while ( (dmname[j]!='[') && (dmname[j]!=0) ){
	     bname[j]=dmname[j];
	     j++;
	   }
	   bname[j]=0;
	   sprintf(leaflist,"%s(%d)[%s]/%s",bname,i,branchcount,itype);
	   Int_t comp = compress;
	   if (type == 5) comp--;
	   sprintf(branchname,"%s.%s(%d)",name,bname,i);
	   TBranch *branch  = new AliArraySubBranch(branchname,this,leaflist,basketsize,comp);
	   branch->SetBit(kIsClone);
	   TObjArray *leaves = branch->GetListOfLeaves();
	   TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(0);
	   leaf->SetOffset(rd->GetThisOffset()+membertype->Size()*i);
	   leaf->SetLeafCount(leafcount);
	   fBranches.Add(branch);
	 }                
       
     }     
   }
   else if (clinfo->IsA()->InheritsFrom("AliDataType")){ //branch for basic type
     Int_t type = (((AliDataType*)clinfo)->GetDataType())->GetType();
     char *itype = 0;
     if (type <=0) 
       Warning("BranchClones","Cannot process member:%s",clinfo->GetName());       
     else{
     
       if (type == 1)  itype = "B";
       if (type == 11) itype = "b";
       if (type == 3)  itype = "I";
       if (type == 5)  itype = "F";
       if (type == 8)  itype = "D";
       if (type == 13) itype = "i";
       if (type == 2)  itype = "S";
       if (type == 12) itype = "s";
       sprintf(leaflist,"%s[%s]/%s",name,branchcount,itype);
       Int_t comp = compress;
       if (type == 5) comp--;
       sprintf(branchname,"%s",clinfo->GetName());       
       TBranch *branch  = new AliArraySubBranch(branchname,this,leaflist,basketsize,comp);
       branch->SetBit(kIsClone);
       TObjArray *leaves = branch->GetListOfLeaves();
       TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(0);
       leaf->SetOffset(0);
       leaf->SetLeafCount(leafcount);
       fBranches.Add(branch);  
     }
   }

}


//______________________________________________________________________________
AliArrayBranch::~AliArrayBranch()
{
//*-*-*-*-*-*Default destructor for a BranchClones*-*-*-*-*-*-*-*-*-*-*-*
//*-*        =====================================

   delete fBranchCount;
   fBranchCount = 0;
   fBranches.Delete();
   fList = 0;
}


//______________________________________________________________________________
void AliArrayBranch::Browse(TBrowser *b)
{
   fBranches.Browse( b );
}

//______________________________________________________________________________
Int_t AliArrayBranch::Fill()
{
//*-*-*-*-*Loop on all Branches of this BranchClones to fill Basket buffer*-*
//*-*      ===============================================================

   Int_t i;
   Int_t nbytes = 0;
   Int_t nbranches = fBranches.GetEntriesFast();
   char **ppointer = (char**)(fAddress);
   if (ppointer == 0) return 0;
   fList = (AliObjectArray*)(*ppointer);
   //   fN    = fList->GetEntriesFast();
   fN    = fList->GetSize();
   fEntries++;

   if (fN > fNdataMax) {
      fNdataMax = fList->GetSize();
      char branchcount[64];
      sprintf(branchcount,"%s_",GetName());
      TLeafI *leafi = (TLeafI*)fBranchCount->GetLeaf(branchcount);
      leafi->SetMaximum(fNdataMax);
      for (i=0;i<nbranches;i++)  {
         TBranch *branch = (TBranch*)fBranches.UncheckedAt(i);
         TObjArray *leaves = branch->GetListOfLeaves();
         TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(0);
         leaf->SetAddress();
      }
   }
   nbytes += fBranchCount->Fill();
   for (i=0;i<nbranches;i++)  {
      TBranch *branch = (TBranch*)fBranches.UncheckedAt(i);
      TObjArray *leaves = branch->GetListOfLeaves();
      TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(0);
      // leaf->Import(fList, fN);   // MI
      Import(leaf,fN);              // MI change
      nbytes += branch->Fill();
   }
   return nbytes;
}

//______________________________________________________________________________
Int_t AliArrayBranch::GetEntry(Int_t entry, Int_t getall)
{
//*-*-*-*-*Read all branches of a BranchClones and return total number of bytes
//*-*      ====================================================================

   if (TestBit(kDoNotProcess) && !getall) return 0;
   Int_t nbytes = fBranchCount->GetEntry(entry);
   TLeaf *leafcount = (TLeaf*)fBranchCount->GetListOfLeaves()->UncheckedAt(0);
   fN = Int_t(leafcount->GetValue());
   if (fN <= 0) return 0;

   TBranch *branch;
   Int_t nbranches = fBranches.GetEntriesFast();

     // if fList exists, create clonesarray objects
   if (fList) {
     //fList->ExpandCreateFast(fN);   //MI  
     fList->Resize(fN);    //MI change 
     for (Int_t i=0;i<nbranches;i++)  {
       branch = (TBranch*)fBranches.UncheckedAt(i);  
       nbytes += ((AliArraySubBranch*)branch)->GetEntryExport(entry, getall, fList, fN);  // !!!MI
      }
   } else {
      for (Int_t i=0;i<nbranches;i++)  {
         branch = (TBranch*)fBranches.UncheckedAt(i);
         nbytes += branch->GetEntry(entry, getall);
      }
   }
  return nbytes;
}

//______________________________________________________________________________
void AliArrayBranch::Print(Option_t *option)
{
//*-*-*-*-*-*-*-*-*-*-*-*Print TBranch parameters*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                    ========================

   fBranchCount->Print(option);
   Int_t i;
   Int_t nbranches = fBranches.GetEntriesFast();
   for (i=0;i<nbranches;i++)  {
      TBranch *branch = (TBranch*)fBranches[i];
      branch->Print(option);
   }
}

//______________________________________________________________________________
void AliArrayBranch::Reset(Option_t *option)
{
//*-*-*-*-*-*-*-*Reset a Branch*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*            ====================
//
//    Existing buffers are deleted
//    Entries, max and min are reset
//

   fEntries        = 0;
   fTotBytes       = 0;
   fZipBytes       = 0;
   Int_t i;
   Int_t nbranches = fBranches.GetEntriesFast();
   for (i=0;i<nbranches;i++)  {
      TBranch *branch = (TBranch*)fBranches[i];
      branch->Reset(option);
   }
   fBranchCount->Reset();
} 

//______________________________________________________________________________
void  AliArrayBranch::SetAddress(void *add)  
{
//*-*-*-*-*-*-*-*Set address of this branch*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*            ====================
//*-*

   fReadEntry = -1;
   fAddress = (char*)add;
   char **ppointer = (char**)(fAddress);
   if ( (*ppointer)==0 ) {  //MI change 
     *ppointer = (char*) new AliObjectArray(fClassName);
     fAddress = (char*)ppointer;
   }
   fList = (AliObjectArray*)(*ppointer);
   fBranchCount->SetAddress(&fN);

}

//______________________________________________________________________________
void AliArrayBranch::SetBasketSize(Int_t buffsize)
{
//*-*-*-*-*-*-*-*Reset basket size for all subbranches of this branchclones
//*-*            ==========================================================
//

   fBasketSize = buffsize;
   Int_t i;
   Int_t nbranches = fBranches.GetEntriesFast();
   for (i=0;i<nbranches;i++)  {
      TBranch *branch = (TBranch*)fBranches[i];
      branch->SetBasketSize(buffsize);
   }
}

//_______________________________________________________________________
void AliArrayBranch::Streamer(TBuffer &b)
{
//*-*-*-*-*-*-*-*-*Stream a class object*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*              =========================================
   if (b.IsReading()) {
      b.ReadVersion();  //Version_t v = b.ReadVersion();
      TNamed::Streamer(b);
      b >> fCompress;
      b >> fBasketSize;
      b >> fEntryOffsetLen;
      b >> fMaxBaskets;
      b >> fWriteBasket;
      b >> fEntryNumber;
      b >> fEntries;
      b >> fTotBytes;
      b >> fZipBytes;
      b >> fOffset;
      b >> fBranchCount;
      fClassName.Streamer(b);
      fBranches.Streamer(b);
      fTree = gTree;
      TBranch *branch;
      TLeaf *leaf;
      Int_t nbranches = fBranches.GetEntriesFast();
      for (Int_t i=0;i<nbranches;i++)  {
         branch = (TBranch*)fBranches[i];
         branch->SetBit(kIsClone);
         leaf = (TLeaf*)branch->GetListOfLeaves()->UncheckedAt(0);
         leaf->SetOffset(0);
      }
      fRead = 1;

      AliClassInfo *clinfo = AliClassInfo::FindClassInfo(fClassName);
      if (!clinfo) {
	AliObjectArray tmp(fClassName); 
	//MI change - object array to construct class description
	clinfo = AliClassInfo::FindClassInfo(fClassName);
      }
      if (!clinfo) return;
      
      TClass *cl = clinfo->GetClass();
      //      if (!cl) {
      //   Warning("Streamer","Unknow class: %s. Cannot read BranchClones: %s",
      //      fClassName.Data(),GetName());
      //   return;
      //}
      if (cl){

	if (!cl->GetListOfRealData())  cl->BuildRealData();
	char branchname[80];
	TRealData *rd;
	TIter      next(cl->GetListOfRealData());
	while ((rd = (TRealData *) next())) {
	  TDataMember *member = rd->GetDataMember();
	  if (!member->IsBasic())      continue;
	  if (!member->IsPersistent()) continue;
	  TDataType *membertype = member->GetDataType();
	  if (membertype->GetType() == 0) continue; 
	 //MI change - for array spliting
	  Int_t arraydim = member->GetArrayDim();
	  if (arraydim==1){
	    for (Int_t i=0;i< member->GetMaxIndex(0);i++){
	      const char * dmname = member->GetName() ;
	      char  bname[200];
	      Int_t j=0;
	      while ( (dmname[j]!='[') && (dmname[j]!=0) ){
		bname[j]=dmname[j];
		j++;
	      }
	      bname[j]=0;
	      sprintf(branchname,"%s.%s(%d)",GetName(),bname,i);
	      branch  = (TBranch*)fBranches.FindObject(branchname);
	      if (!branch) continue;
	      TObjArray *leaves = branch->GetListOfLeaves();
	      leaf = (TLeaf*)leaves->UncheckedAt(0);
	      leaf->SetOffset(rd->GetThisOffset()+membertype->Size()*i);
	    }
	  }
	  sprintf(branchname,"%s.%s",GetName(),rd->GetName());
	  branch  = (TBranch*)fBranches.FindObject(branchname);
	  if (!branch) continue;
	  TObjArray *leaves = branch->GetListOfLeaves();
	  leaf = (TLeaf*)leaves->UncheckedAt(0);
	  leaf->SetOffset(rd->GetThisOffset());     	
	}

      }
      else if (clinfo->IsA()->InheritsFrom("AliDataType")){ //branch for basic type
	char branchname[100];
	sprintf(branchname,"%s",clinfo->GetName());  
	branch  = (TBranch*)fBranches.FindObject(branchname);
	if (branch){
	  TObjArray *leaves = branch->GetListOfLeaves();
	  leaf = (TLeaf*)leaves->UncheckedAt(0);
	  leaf->SetOffset(0);
	}   
      }
   }
   else{
	
      b.WriteVersion(AliArrayBranch::IsA());
      TNamed::Streamer(b);
      b << fCompress;
      b << fBasketSize;
      b << fEntryOffsetLen;
      b << fMaxBaskets;
      b << fWriteBasket;
      b << fEntryNumber;
      b << fEntries;
      b << fTotBytes;
      b << fZipBytes;
      b << fOffset;
      b << fBranchCount;
      fClassName.Streamer(b);
      fBranches.Streamer(b);
   }
}



void AliArrayBranch::Import(TLeaf * leaf, Int_t n)
{

  const Int_t kIntUndefined = -9999;
  Int_t j = 0;
  char *clone;
  Int_t len    = leaf->GetLenStatic();
  Int_t fOffset = leaf->GetOffset();
  void *value  = leaf->GetValuePointer();
  //
  for (Int_t i=0;i<n;i++) {
    clone = (char*)fList->UncheckedAt(i);
    //8bit int
    if (leaf->IsA()==TLeafB::Class()){  
      memcpy(&((Char_t*)value)[j],clone + fOffset, len);
    }
    //var size
    if (leaf->IsA()==TLeafC::Class()){  
      memcpy(&((Char_t*)value)[j],clone + fOffset, 1);
    }
    //double
    if (leaf->IsA()==TLeafD::Class()){  
      if (clone) memcpy(&((Double_t*)value)[j],clone + fOffset, 8*len);
      else       memcpy(&((Double_t*)value)[j],&kIntUndefined,  8*len);      
    }
    //float
    if (leaf->IsA()==TLeafF::Class()){   
      if (clone) memcpy(&((Float_t*)value)[j],clone + fOffset, 4*len);
      else       memcpy(&((Float_t*)value)[j],&kIntUndefined,  4*len);     
    }    
    //int
    if (leaf->IsA()==TLeafI::Class()){  
      if (clone) memcpy(&((Int_t*)value)[j],clone + fOffset, 4*len);
      else       memcpy(&((Int_t*)value)[j],&kIntUndefined,  4*len);
    }
   //short
    if (leaf->IsA()==TLeafS::Class()){  
      if (clone) memcpy(&((Short_t*)value)[j],clone + fOffset, 2*len);
      else       memcpy(&((Short_t*)value)[j],&kIntUndefined,  2*len);
    }     
    j += len;  
  }
  //
}

AliObjectBranch::AliObjectBranch(const Text_t *name, const Text_t *classname, void *addobj,
				 TTree * tree, 
				 Int_t basketsize, Int_t splitlevel, Int_t compress): TBranchObject()
{
//*-*-*-*-*-*-*-*-*-*-*-*-*Create a BranchObject*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                      =====================
//
   TClass *cl      = gROOT->GetClass(classname);
   fTree =tree;  //MI change
   if (!cl) {
      Error("TBranchObject","Cannot find class:%s",classname);
      return;
   }
   if (!cl->GetListOfRealData())  cl->BuildRealData();
   Int_t bufsize=basketsize; //MI ?
   SetName(name);
   SetTitle(name);
   fCompress = compress;
   if (compress == -1) {
     TFile *bfile = fTree->GetDirectory()->GetFile(); //MI chnge fTrre - gTree
      if (bfile) fCompress = bfile->GetCompressionLevel();
   }
   if (basketsize < 100) basketsize = 100;
   fBasketSize     = basketsize;
   fAddress        = (char*)addobj;
   fClassName      = classname;
   fBasketEntry    = new Int_t[fMaxBaskets];
   fBasketBytes    = new Int_t[fMaxBaskets];
   fBasketSeek     = new Seek_t[fMaxBaskets];
   fOldObject      = 0;

   fBasketEntry[0] = fEntryNumber;
   fBasketBytes[0] = 0;

   TLeaf *leaf     = new TLeafObject(name,classname);
   leaf->SetBranch(this);
   leaf->SetAddress(addobj);
   fNleaves = 1;
   fLeaves.Add(leaf);
   fTree->GetListOfLeaves()->Add(leaf);  //MI change fTree-gTree

// Set the bit kAutoDelete to specify that when reading
// in TLeafObject::ReadBasket, the object should be deleted
// before calling Streamer.
// It is foreseen to not set this bit in a future version.
   SetAutoDelete(kTRUE);

//*-*-  Create the first basket
   //   fTree       = gTree;  //MI change - no need anymore 
   fDirectory  = fTree->GetDirectory();
   fFileName   = "";

   if (!splitlevel){
     TBasket *basket = new TBasket(name,fTree->GetName(),this);
     fBaskets.Add(basket);
     return;
   }

   //
   // 
   TBranch * branch =this;
   TObjArray *blist = branch->GetListOfBranches();
   const char *rdname;
   const char *dname;
   char branchname[64];
   if (!cl->GetListOfRealData()) cl->BuildRealData();
   char **apointer = (char**)(addobj);
   TObject *obj = (TObject*)(*apointer);
   Bool_t delobj = kFALSE;
   if (!obj) {
      obj = (TObject*)cl->New();
      delobj = kTRUE;
   }
//*-*- Loop on all public data members of the class and its base classes
   Int_t lenName = strlen(name);
   Int_t isDot = 0;
   if (name[lenName-1] == '.') isDot = 1;
   TBranch *branch1 = 0;
   TRealData *rd;
   TIter      next(cl->GetListOfRealData());
   while ((rd = (TRealData *) next())) {
      TDataMember *dm = rd->GetDataMember();
      if (!dm->IsPersistent()) continue; //do not process members with a ! as the first
                                         // character in the comment field
      rdname = rd->GetName();
      dname  = dm->GetName();

  //  Next line now commented, functionality to process arrays is now implemented
  //  the statement is left to show how to use Property() and kIsArray
  //     if (dm->Property() & kIsArray) continue;

      TDataType *dtype = dm->GetDataType();
      Int_t code = 0;
      if (dtype) code = dm->GetDataType()->GetType();

//*-*- Encode branch name. Use real data member name
      sprintf(branchname,"%s",rdname);
      if (isDot) {
         if (dm->IsaPointer()) sprintf(branchname,"%s%s",name,&rdname[1]);
         else                  sprintf(branchname,"%s%s",name,&rdname[0]);
      }
      char leaflist[64];
      Int_t offset    = rd->GetThisOffset();
      char *pointer   = (char*)obj + offset;
      if (dm->IsaPointer()) {
         TClass *clobj = 0;
         if (!dm->IsBasic()) clobj = gROOT->GetClass(dm->GetTypeName());
         if (clobj && !strcmp("TClonesArray",clobj->GetName())) {
            char *cpointer  =(char*)pointer;
            char **ppointer =(char**)cpointer;
            TClonesArray *list = (TClonesArray*)(*ppointer);
            if (splitlevel != 2) {
               if (isDot) branch1 = new TBranchClones(&branchname[0],pointer,bufsize);
               else       branch1 = new TBranchClones(&branchname[1],pointer,bufsize);
               blist->Add(branch1);
            } else {
               if (isDot) branch1 = new TBranchObject(&branchname[0],list->ClassName(),pointer,bufsize);
               else       branch1 = new TBranchObject(&branchname[1],list->ClassName(),pointer,bufsize);
               blist->Add(branch1);
            }
         }
	 else
	   if (clobj && !strcmp("AliObjectArray",clobj->GetName())) {
	     char *cpointer  =(char*)pointer;
	     char **ppointer =(char**)cpointer;
	     TClonesArray *list = (TClonesArray*)(*ppointer);
	     if (splitlevel != 2) {
               if (isDot) branch1 = new AliArrayBranch(&branchname[0],pointer,fTree,bufsize,compress);
               else       branch1 = new AliArrayBranch(&branchname[1],pointer,fTree,bufsize,compress);
               blist->Add(branch1);
	     } else {
               if (isDot) branch1 = new AliObjectBranch(&branchname[0],list->ClassName(),pointer,fTree,bufsize);
               else       branch1 = new AliObjectBranch(&branchname[1],list->ClassName(),pointer,fTree,bufsize);
               blist->Add(branch1);
	     }
	   }
	   else {
	     if (!clobj) {
               if (code != 1) continue;
               sprintf(leaflist,"%s/%s",dname,"C");
               branch1 = new TBranch(branchname,pointer,leaflist,bufsize);
               branch1->SetTitle(dname);
               blist->Add(branch1);
	     } else {
               if (!clobj->InheritsFrom(TObject::Class())) continue;
               //branch1 = new TBranchObject(dname,clobj->GetName(),pointer,bufsize,0); //MI change
	       branch1 = new AliObjectBranch(dname,clobj->GetName(),pointer,fTree,bufsize,splitlevel);
               if (isDot) branch1->SetName(&branchname[0]);
               else       branch1->SetName(&branchname[1]);  //do not use the first character (*)
               blist->Add(branch1);
	     }
	   }
      }else {
//*-*-------------Data Member is a basic data type----------
	if (dm->IsBasic()) {
	  if      (code ==  1) sprintf(leaflist,"%s/%s",rdname,"B");
	  else if (code == 11) sprintf(leaflist,"%s/%s",rdname,"b");
	  else if (code ==  2) sprintf(leaflist,"%s/%s",rdname,"S");
	  else if (code == 12) sprintf(leaflist,"%s/%s",rdname,"s");
	  else if (code ==  3) sprintf(leaflist,"%s/%s",rdname,"I");
	  else if (code == 13) sprintf(leaflist,"%s/%s",rdname,"i");
	  else if (code ==  5) sprintf(leaflist,"%s/%s",rdname,"F");
	  else if (code ==  8) sprintf(leaflist,"%s/%s",rdname,"D");
	  else {
	    printf("Cannot create branch for rdname=%s, code=%d\n",branchname, code);
	    leaflist[0] = 0;
	  }
	  branch1 = new TBranch(branchname,pointer,leaflist,bufsize);
	  branch1->SetTitle(rdname);
	  blist->Add(branch1);
	}
      }
      if (branch1) branch1->SetOffset(offset);
      else Warning("Branch","Cannot process member:%s",rdname);      
   }
   if (delobj) delete obj;
}



//______________________________________________________________________________
void AliObjectBranch::SetAddress(void *add)
{
//*-*-*-*-*-*-*-*Set address of this branch*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*            ====================
//

   //special case when called from code generated by TTree::MakeClass
   if (Long_t(add) == -1) {
      SetBit(kWarn);
      return;
   }
   fReadEntry = -1;
   Int_t nbranches = fBranches.GetEntriesFast();
   TLeaf *leaf = (TLeaf*)fLeaves.UncheckedAt(0);
   if (leaf) leaf->SetAddress(add);
   TBranch *branch;
   fAddress = (char*)add;
   char *pointer   = fAddress;
   void **ppointer = (void**)add;
   TObject *obj = (TObject*)(*ppointer);
   TClass *cl = gROOT->GetClass(fClassName.Data());
   if (!obj && cl) {
      obj = (TObject*)cl->New();
      *ppointer = (void*)obj;
   }
   Int_t i, offset;
   if (!cl) {
      for (i=0;i<nbranches;i++)  {
         branch  = (TBranch*)fBranches[i];
         pointer = (char*)obj;
         branch->SetAddress(pointer);
      }
      return;
   }
   if (!cl->GetListOfRealData())  cl->BuildRealData();
   char *fullname = new char[200];
   const char *bname = GetName();
   Int_t lenName = strlen(bname);
   Int_t isDot = 0;
   if (bname[lenName-1] == '.') isDot = 1;
   const char *rdname;
   TRealData *rd;
   TIter      next(cl->GetListOfRealData());
   while ((rd = (TRealData *) next())) {
      TDataMember *dm = rd->GetDataMember();
      if (!dm->IsPersistent()) continue;
      rdname = rd->GetName();
      TDataType *dtype = dm->GetDataType();
      Int_t code = 0;
      if (dtype) code = dm->GetDataType()->GetType();
      offset  = rd->GetThisOffset();
      pointer = (char*)obj + offset;
      branch  = 0;
      if (dm->IsaPointer()) {
         TClass *clobj = 0;
         if (!dm->IsBasic()) clobj = gROOT->GetClass(dm->GetTypeName());
         if (clobj && !strcmp("TClonesArray",clobj->GetName())) {
            if (isDot) sprintf(fullname,"%s%s",bname,&rdname[1]);
            else       sprintf(fullname,"%s",&rdname[1]);
            branch = (TBranch*)fBranches.FindObject(fullname);
         }
	 else
	   if (clobj && !strcmp("AliObjectArray",clobj->GetName())) {
	     if (isDot) sprintf(fullname,"%s%s",bname,&rdname[1]);
	     else       sprintf(fullname,"%s",&rdname[1]);
	     branch = (TBranch*)fBranches.FindObject(fullname);
	   }
	 else {
            if (!clobj) {
               if (code != 1) continue;
               if (isDot) sprintf(fullname,"%s%s",bname,&rdname[0]);
               else       sprintf(fullname,"%s",&rdname[0]);
               branch = (TBranch*)fBranches.FindObject(fullname);
            } else {
               if (!clobj->InheritsFrom(TObject::Class())) continue;
               if (isDot) sprintf(fullname,"%s%s",bname,&rdname[1]);
               else       sprintf(fullname,"%s",&rdname[1]);
               branch = (TBranch*)fBranches.FindObject(fullname);
            }
         }
      } else {
         if (dm->IsBasic()) {
            if (isDot) sprintf(fullname,"%s%s",bname,&rdname[0]);
            else       sprintf(fullname,"%s",&rdname[0]);
            branch = (TBranch*)fBranches.FindObject(fullname);
         }
      }
      if(branch) branch->SetAddress(pointer);
   }
   delete [] fullname;
}





AliTree::AliTree(const char *name,const char *title, Int_t maxvirtualsize):
  TTree(name,title,maxvirtualsize)
{
  //
  //default constructor for AliTree
  gTree =this;
}

TBranch * AliTree::AliBranch(const char *name, void *clonesaddress, Int_t bufsize, Int_t splitlevel,
			     Int_t compres)
{
  if (clonesaddress == 0) return 0;
  char *cpointer =(char*)clonesaddress;
  char **ppointer =(char**)cpointer;
  AliObjectArray *list = (AliObjectArray*)(*ppointer);
  if (list == 0) return 0;
  gTree = this;
  if (splitlevel) {
    TBranch *branch = new AliArrayBranch(name,clonesaddress,this,bufsize, compres);
    fBranches.Add(branch);
    return branch;
  } else {
    TBranchObject *branch = new TBranchObject(name,list->ClassName(),clonesaddress,bufsize,0);
    fBranches.Add(branch);
    return branch;
  }
}

TBranch* AliTree::AliBranch(const char *name, const char *classname, void *addobj, 
		     Int_t bufsize, Int_t splitlevel)
{
  gTree = this;
  TClass *cl = gROOT->GetClass(classname);
  if (!cl) {
    Error("BranchObject","Cannot find class:%s",classname);
    return 0;
  }
  TBranch * branch = new AliObjectBranch(name,classname,addobj, this, bufsize,splitlevel);
  fBranches.Add(branch);
  return branch;
}


