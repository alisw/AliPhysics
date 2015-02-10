/// \file MakeTreeStat.C
///
/// Macro to get the size of the Tree
/// As a improvment to  tree->Print() function, this algorithm
/// gives the size of all of the branchces and in addition print them 
/// sorted according total tree size (MEMORY USAGE if one event in tree)
/// or zip size (THE storage size on disk)
///
/// Printed statistic:
/// 1. Order
/// 2. TotSize (in memory) + fraction of total size
/// 3. ZipSize (in memory) + fraction of zip size
/// 4. Compression ratio 
/// 
/// Usage:
///
/// 1. Enable macro
///
/// ~~~ 
/// .L $ALICE_ROOT/TPC/macros/MakeTreeStat.C+
/// ~~~
///
/// 2. Open the tree (eg.)
///
/// ~~~{.cpp}
/// TFile f("AliESDs.root");
/// TTree * tree = (TTree*)f.Get("esdTree");
/// ~~~
///
/// 3. Print statistic (sorting according secon argument - either zip Bytes (kTRUE or TotSize (kFALSE)
///
/// ~~~{.cpp}
/// MakeStat(tree, kTRUE);
/// ~~~
/// 
/// \author M.Ivanov, GSI, m.ivanov@gsi.de







#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TArrayF.h"
#include "TObjArray.h"
#include "TObjString.h"




TTree     *fTree;   // tree of interest
TObjArray aReport;  // array with branch statistic
TArrayF   totSize;  // total size for branch
TArrayF   zipSize;  // zip Size for branch
TArrayF   zipRatio; // zip Ratio for branch

void MakeStat(TTree *tree, Bool_t zipSort);



void PrintSorted(Bool_t zipSort){
  /// print statistic

  Int_t entries   = aReport.GetEntries();
  Int_t* indexes  = new Int_t[entries];
  if (zipSort) TMath::Sort(entries,zipSize.GetArray(),indexes,kTRUE);
  else{
    TMath::Sort(entries,totSize.GetArray(),indexes,kTRUE);
  }
  Float_t zipBytes = zipSize[indexes[0]];
  Float_t totBytes =  totSize[indexes[0]];
  Float_t ratioT = 100.*zipBytes/totBytes;
  for (Int_t i=entries-1; i>=0; i--){
    Int_t ib = indexes[i];
    Float_t ratio0= 100.*totSize[ib]/totBytes;
    Float_t ratio1= 100.*zipSize[ib]/zipBytes;
    if (i==0) {
      printf("\n------------------------------------------------------------\n");
      printf("%d  \t\t%5.f(%.2f\%)  \t\t%5.f(%.2f\%)  \t%.2f  \t%s\n",i, 
	     totSize[ib],100., zipSize[ib],100., 100.*zipRatio[ib], aReport.At(ib)->GetName());
    }else{
    printf("%d  \t\t%5.f(%.2f\%)  \t\t%5.f(%.2f\%)  \t%.2f  \t%s\n",i, 
	   totSize[ib],ratio0, zipSize[ib],ratio1, 100.*zipRatio[ib], aReport.At(ib)->GetName());
    }
  }
    

}


void AddToReport(const char *prefix,const char * name, Float_t size[2], Float_t ratio){
  /// add branch info to array

  char fullname[10000];
  sprintf(fullname,"%s.%s",prefix,name);
  aReport.AddLast(new TObjString(fullname));
  Int_t entry = aReport.GetEntries();
  if (totSize.GetSize()<entry){
    totSize.Set(entry*2);
    zipSize.Set(entry*2);
    zipRatio.Set(entry*2);
  }
  totSize[entry-1]=Float_t(size[0]);
  zipSize[entry-1]=Float_t(size[1]);
  zipRatio[entry-1]=Float_t(ratio);
}





void MakeStat(const char *prefix, TBranch * branch, Float_t* size, Float_t mratio);



void MakeStat(TTree *tree, Bool_t zipSort){
  /// make recursve loop over tree branches

  fTree= tree;
  aReport.Clear();
  TObjArray * array = tree->GetListOfBranches(); 
  Float_t size[2]={0,0};
  char * prefix ="";
  Float_t mratio=tree->GetZipBytes()/float(tree->GetTotBytes());
  for (Int_t i=0; i<array->GetEntries(); i++){
    MakeStat(prefix,(TBranch*)array->At(i),size, mratio);
  }  
  Float_t ratio= (size[0]>0) ? size[1]/Float_t(size[0]): 0;
  //  printf("Sum :\t%f\t%f\t%f\t%s.%s\t\n", float(size[0]), float(size[1]),ratio, prefix,tree->GetName());

  AddToReport(prefix, tree->GetName(),size,ratio);
  PrintSorted(zipSort);
}


void MakeStat(const char *prefix, TBranch * branch, Float_t *size, Float_t mratio){
  /// Recursive function to get size of the branches
  /// and ratios

  TObjArray * array = branch->GetListOfBranches();
  Float_t bsizeSum[2]={0,0};

  if (!array || array->GetEntries()==0){
    Float_t bsize[2] = {0,0};
    bsize[0]=branch->GetTotalSize();
    bsize[1]=branch->GetZipBytes();
    if (bsize[1]>0){
      Float_t ratio= (bsize[0]>0) ? bsize[1]/Float_t(bsize[0]): 0;
      //      printf("Term:\t%f\t%f\t%f\t%s.%s\t\n",float(bsize[0]), float(bsize[1]),ratio, prefix,branch->GetName());
      AddToReport(prefix, branch->GetName(),bsize,ratio);	
      //branch->Print();
      size[0]+=bsize[0];
      size[1]+=bsize[1];
    }else{
      Float_t ratio= mratio;
      //printf("Ter?:\t%f\t%f\t%f\t%s.%s\t\n",float(bsize[0]), float(-1),ratio, prefix,branch->GetName());
      AddToReport(prefix, branch->GetName(),bsize,ratio);
      //branch->Print();
      size[0]+=bsize[0];
      size[1]+=TMath::Nint(bsize[0]*mratio);
    }

    return;
  }
  for (Int_t i=0; i<array->GetEntries(); i++){
    Float_t bsize[2] = {0,0};
    TString str=prefix;
    str+= branch->GetName();
    MakeStat(str.Data(),(TBranch*)array->At(i), bsize, mratio);
    bsizeSum[0]+=bsize[0];
    bsizeSum[1]+=bsize[1];
    size[0]+=bsize[0];
    size[1]+=bsize[1];
  } 
  Float_t ratio= (size[0]>0) ? size[1]/Float_t(size[0]): 0;
  //printf("Sum :\t%f\t%f\t%f\t%s.%s\t\n", float(bsizeSum[0]), float(bsizeSum[1]),ratio, prefix,branch->GetName()); 
  AddToReport(prefix,branch->GetName(),bsizeSum,ratio);
}
