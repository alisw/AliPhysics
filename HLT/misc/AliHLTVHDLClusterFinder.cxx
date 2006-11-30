// @(#) $Id$

// Author: Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright & copy ALICE HLT Group
/** \class AliHLTVHDLClusterFinder
<pre>
//____________________________________________________
// AliHLTVHDLClusterFinder
//
// The current VHDL cluster finder for HLT
// Based on STAR L3
//
// Most important parameters:
// fThreshold - threshold for noise clusters
// fMatch - length in time for overlapping sequences
</pre> 
*/

#include "AliHLTStandardIncludes.h"
#include "AliHLTRootTypes.h"
#include "AliHLTLogging.h"
#include "AliHLTAltroMemHandler.h"

//#define VHDLDEBUG
#include "AliHLTVHDLClusterFinder.h"

#if __GNUC__ >= 3
using namespace std;
#endif


ClassImp(AliHLTVHDLClusterFinder)

AliHLTVHDLClusterFinder::AliHLTVHDLClusterFinder()
{
  // default constructor
  fMatch = 4;
  fThreshold = 10;
  fMinMerge = 1;
  fNClusters=0;

  fXYErr = 0.2;
  fZErr = 0.3;
  fDeconvPad = kTRUE;
  fDeconvTime = kTRUE;
  fstdout = kFALSE;
  fcalcerr = kTRUE;

  Clear();
#ifdef VHDLDEBUG
  fdeb=fopen("vhdlclusterfinder.debug","w");
  //fdeb=stderr;
#endif
}

AliHLTVHDLClusterFinder::~AliHLTVHDLClusterFinder()
{
  // destructor
#ifdef VHDLDEBUG
  fclose(fdeb);
#endif
}

void AliHLTVHDLClusterFinder::ProcessDigits()
{
  //Loop over data like the analyzer of the VHDL code
  const UChar_t kn=255;
  UShort_t rrow=0,rtime=0;
  UChar_t rpad=0,i=kn;
  UShort_t *charges=new UShort_t[kn];

  fNClusters=0;
  fRow=0;
  fPad=0;
  Clear();

  //loop over input data
  while(fAltromem.ReadSequence(rrow,rpad,rtime,i,&charges)){
#if 0
    cout << "Padrow " << (int)rrow << " pad " << (int)rpad << " time " <<(int)rtime << " charges ";
    for(UChar_t ii=0;ii<i;ii++) cout << (int)charges[ii] << " ";
    cout << endl;
#endif
#ifdef VHDLDEBUG
    fprintf(fdeb,"ProcessDigits: Input Data: %d %d %d charges:",(int)rrow,(int)rpad,(int)rtime);
    for(UChar_t ii=0;ii<i;ii++) fprintf(fdeb," %d",(int)charges[ii]);
    fprintf(fdeb,"\n");
#endif

    fNRow=rrow;
    fNPad=rpad;
    fTC=0,fMT=0,fST=0;

    //calculate sequence values 
    if(fDeconvTime){
      UChar_t ii=0;
      Int_t sl=0;
      Int_t charge=0,lcharge=0;
      Bool_t falling=kFALSE;
      while(ii<i){
	charge=charges[ii];
	if((falling)&&(charge>lcharge)){
	  fSM=rtime-sl/2;
	  MakeSequence();
	  ProcessSequence();
	  OutputMemory();

	  rtime=rtime-sl;      
	  falling=kFALSE;
	  sl=0;
	  fTC=0,fMT=0,fST=0;
	} else if(charge<lcharge) falling=kTRUE;

	fTC+=charge;
	fMT+=(rtime-sl)*charge;
	fST+=(rtime-sl)*(rtime-sl)*charge;
	sl++;
	ii++;
	lcharge=charge;
      } //end loop over sequence
      fSM=rtime-sl/2;
    } else { /* no deconvolution */
      for(UChar_t ii=0;ii<i;ii++){
	fTC+=charges[ii];
	fMT+=(rtime-ii)*charges[ii];
	fST+=(rtime-ii)*(rtime-ii)*charges[ii];
      }
      fSM=rtime-i/2;
    }

    MakeSequence();
    ProcessSequence();
    OutputMemory();      

#ifdef VHDLDEBUG
    fflush(fdeb);
#endif
    i=kn; //store size of charges array
  } //loop over data

  //flush everything left
  FlushMemory();
  while(fOP!=fFP) OutputMemory();
}

void AliHLTVHDLClusterFinder::MakeSequence(){
  // makes the sequence
  if(!fTC) return;

  Int_t mp=fNPad*fTC;
  Int_t sp=fNPad*fNPad*fTC;

  fSeq.fTotalCharge=fTC;
  fSeq.fPad=mp;
  fSeq.fTime=fMT;
  fSeq.fPad2=sp;
  fSeq.fTime2=fST; 
  fSeq.fMean=0;
  fSeq.fMean=fSM;
  fSeq.fMerge=0;
  fSeq.fRow=fNRow;
  fSeq.fLastPad=fNPad;
  fSeq.fChargeFalling=0;
  if(fDeconvPad) fSeq.fLastCharge=fTC;
  else fSeq.fLastCharge=0;
}

void AliHLTVHDLClusterFinder::ProcessSequence()
{
  // processes the sequence
  if(fNRow!=fRow) FlushMemory();
  else if(fNPad==fPad+1) PrepareMemory();
  else if(fNPad!=fPad) FlushMemory();

  //store new row and pad values
  fRow=fNRow;
  fPad=fNPad;

#ifdef VHDLDEBUG
  fprintf(fdeb,"ProcessSequence: Mean=%d TC=%d ",fSeq.fMean,fSeq.fTotalCharge);
#endif

  CompareSeq(); //merge or insert
}

void AliHLTVHDLClusterFinder::PrepareMemory()
{
  // prepares the memory
#ifdef VHDLDEBUG
  fprintf(fdeb,"PrepareMemory %d %d %d\n",fRP,fEP,fWP);
#endif
  fRP=fEP;
  fEP=fWP;
}

void AliHLTVHDLClusterFinder::FlushMemory()
{
  // flushes the memory
#ifdef VHDLDEBUG
  fprintf(fdeb,"FlushMemory %d %d %d %d\n",fFP,fRP,fEP,fWP);
#endif
  fFP=fWP;
  fRP=fWP;
  fEP=fWP;
}

void AliHLTVHDLClusterFinder::CompareSeq()
{
  // compares sequences
  while(fRP!=fEP){
    Int_t diff=fSeqs[fPList[fRP]].fMean-fSeq.fMean;

    if(diff>fMatch){ //no match
      IncRPointer(); //cluster finished
      continue;
    } else if(diff<-fMatch){ //no match
      break;                 //insert new cluster       
    } else { //match found, merge it
      MergeSeq(); 
      return;
    }
  }

  InsertSeq(); //start new cluster
}

void AliHLTVHDLClusterFinder::MergeSeq()
{
  // merges sequences
#ifdef VHDLDEBUG
  fprintf(fdeb,"merged with Mean=%d TC=%d (new Merge=%d) %d %d\n",fSeqs[fPList[fRP]].fMean,fSeqs[fPList[fRP]].fTotalCharge,fSeqs[fPList[fRP]].fMerge+1,fRow,fPad);
#endif
  if(fSeqs[fPList[fRP]].fRow!=fSeq.fRow){
    LOG(AliHLTLog::kWarning,"AliHLTVHDLClusterFinder::","Memory Check")
      <<"Sequences can be merged on the same rows only."<<ENDLOG;
  }
  if(fSeqs[fPList[fRP]].fLastPad+1!=fSeq.fLastPad){
    LOG(AliHLTLog::kWarning,"AliHLTVHDLClusterFinder::","Memory Check")
      <<"Sequences can be merged on consecutive pads only."<<ENDLOG;
  }
  
  if(fDeconvPad){
    if(fSeq.fTotalCharge > fSeqs[fPList[fRP]].fLastCharge){
      if(fSeqs[fPList[fRP]].fChargeFalling){ //The previous pad was falling
	IncRPointer();
	InsertSeq(); //start a new cluster
	return;
      }		    
    }
    else fSeqs[fPList[fRP]].fChargeFalling = 1;

    fSeqs[fPList[fRP]].fLastCharge = fSeq.fTotalCharge;
  }
  
  fSeqs[fPList[fRP]].fMean=fSeq.fMean; //take the new mean
  fSeqs[fPList[fRP]].fLastPad=fSeq.fLastPad;
  fSeqs[fPList[fRP]].fTotalCharge+=fSeq.fTotalCharge;
  fSeqs[fPList[fRP]].fPad+=fSeq.fPad;
  fSeqs[fPList[fRP]].fTime+=fSeq.fTime;
  fSeqs[fPList[fRP]].fPad2+=fSeq.fPad2;
  fSeqs[fPList[fRP]].fTime2+=fSeq.fTime2;
  fSeqs[fPList[fRP]].fMerge+=1;

  fPList[fWP]=fPList[fRP];
  fPList[fRP]=N_clmem; //mark it to be free

  IncRPointer();
  IncWPointer();
}

void AliHLTVHDLClusterFinder::InsertSeq()
{
  // inserts sequence
#ifdef VHDLDEBUG
  fprintf(fdeb,"inserted %d %d\n",fRow,fPad);
#endif
  NextFreeIndex();    //get next index
  fSeqs[fFirst]=fSeq; //store data
  fPList[fWP]=fFirst; //store index
  IncWPointer();
}

void AliHLTVHDLClusterFinder::OutputMemory()
{
  // output memory?
  Float_t mtime=0,mpad=0;
  Float_t mtime2=0,mpad2=0;
  UInt_t tc,row,mno;

  if(fOP!=fFP){
  //while(fOP!=fFP){

    UInt_t index=fPList[fOP];
    fPList[fOP]=N_clmem;
    //cout << fOP << " - " << fFP << " - " << index << endl;
    IncPointer(fOP);
    if(index>=N_clmem) return; //nothing to do
    //if(index>=N_clmem) continue; //nothing to do
    
    tc=fSeqs[index].fTotalCharge;
    mno=fSeqs[index].fMerge;
    row=fSeqs[index].fRow;
    if(tc!=0){
      mtime=(Float_t)(fSeqs[index].fTime)/tc;
      mtime2=sqrt((Float_t)(fSeqs[index].fTime2)/tc-mtime*mtime);
    }
    if(tc!=0){
      mpad=(Float_t)(fSeqs[index].fPad)/tc;
      mpad2=sqrt((Float_t)(fSeqs[index].fPad2)/tc-mpad*mpad);
    }

    if(mno<fMinMerge){ //noise cluster
#ifdef VHDLDEBUG
    fprintf(fdeb,"OutputMemory %d %d: noise cluster (merge): Row %d Pad %.3f Time %.3f TC %d Merge %d\n",fOP,fFP,row,mpad,mtime,tc,mno);
#endif          
      FreeSeq(index);  
      return;
      //continue;
    }

    if(tc<fThreshold){ //total charge below threshold
#ifdef VHDLDEBUG
    fprintf(fdeb,"OutputMemory %d %d: noise cluster (threshold): Row %d Pad %.3f Time %.3f TC %d Merge %d\n",fOP,fFP,row,mpad,mtime,tc,mno);
#endif          
      FreeSeq(index);
      return;
      //continue;
    }

#ifdef VHDLDEBUG
    fprintf(fdeb,"OutputMemory %d: Row %d Pad %.3f Time %.3f TC %d Merge %d\n",fNClusters,row,mpad,mtime,tc,mno);
#endif    

    //    if(stdout)
    //      cout<<"WriteCluster: padrow "<<row<<" pad "<<mpad<< " +- "<<mpad2<<" time "<<mtime<<" +- "<<mtime2<<" charge "<<tc<<endl;
    
    fNClusters++;
    FreeSeq(index);
  }
}

void AliHLTVHDLClusterFinder::FreeSeq(UShort_t i)
{
  // frees the sequence
  ClearSeq(i);
  IncPointer(fLast,1,N_clmem);
}

void AliHLTVHDLClusterFinder::Clear()
{
  // clears everything
  fFirst=0; //first list pointer
  fLast=0;  //last  list pointer

  fWP=0; //write pointer
  fRP=0; //read pointer
  fEP=0; //end read pointer
  fFP=0; //end flush pointer
  fOP=0; //output pointer

  fSeq.fRow=0;
  fSeq.fLastPad=0;
  fSeq.fTotalCharge=0;
  fSeq.fPad=0;
  fSeq.fTime=0;
  fSeq.fPad2=0;
  fSeq.fTime2=0; 
  fSeq.fMean=0;
  fSeq.fMerge=0;

  for(UInt_t i=0;i<N_mem;i++) fPList[i]=N_clmem;
  for(UInt_t i=0;i<N_clmem;i++) ClearSeq(i);
}

void AliHLTVHDLClusterFinder::ClearSeq(UShort_t i){
  // clears a sequence
  fSeqs[i].fRow=0;
  fSeqs[i].fLastPad=0;
  fSeqs[i].fTotalCharge=0;
  fSeqs[i].fPad=0;
  fSeqs[i].fTime=0;
  fSeqs[i].fPad2=0;
  fSeqs[i].fTime2=0; 
  fSeqs[i].fMean=0;
  fSeqs[i].fMerge=0;
}

void AliHLTVHDLClusterFinder::IncPointer(UShort_t &p, Short_t add, UShort_t N){
  // increments pointer by 'add' in a circular buffer 
  Short_t pp=p;
  pp+=add;  
  if(pp>=N) p=UShort_t(pp-N);
  else if(pp<0) p=UShort_t(pp+N);
  else p=UShort_t(pp);
}

void AliHLTVHDLClusterFinder::IncRPointer(){
  // increments pointer fRP 
  IncPointer(fRP);
}

void AliHLTVHDLClusterFinder::IncWPointer(){
  // increments pointer fWP
  IncPointer(fWP);

  if(fWP==fOP){
    LOG(AliHLTLog::kWarning,"AliHLTVHDLClusterFinder::IncWPointer","Memory Check")
      <<"Write pointer overwrites output pointer."<<ENDLOG;
  }
}

void AliHLTVHDLClusterFinder::NextFreeIndex(){
  // finds next free index
  IncPointer(fFirst,1,N_clmem);
  if(fFirst==fLast) {
    LOG(AliHLTLog::kFatal,"AliHLTVHDLClusterFinder::GetFreeIndex","Memory Check")
      <<"No space left in sequence list: "<<fFirst<<"=="<<fLast<<ENDLOG;
  }
}

#if 0
void AliHLTClustFinderNew::WriteClusters(Int_t n_clusters,ClusterData *list)
{
  // writes clusters
  Int_t thisrow,thissector;
  UInt_t counter = fNClusters;
  
  for(int j=0; j<n_clusters; j++)
    {
      if(!list[j].fFlags) continue; //discard 1 pad clusters
      if(list[j].fTotalCharge < fThreshold) continue; //noise cluster

      Float_t xyz[3];      
      Float_t fpad =(Float_t)list[j].fPad /(Float_t)list[j].fTotalCharge;
      Float_t fpad2=fXYErr;
      Float_t ftime =(Float_t)list[j].fTime /(Float_t)list[j].fTotalCharge;
      Float_t ftime2=fZErr;

      if(fcalcerr) {
	fpad2=(Float_t)list[j].fPad2/(Float_t)list[j].fTotalCharge - fpad*fpad;
	fpad2 = sqrt(fpad2);
	ftime2=(Float_t)list[j].fTime2/(Float_t)list[j].fTotalCharge - ftime*ftime;
	ftime2 = sqrt(ftime2); 
      }
       
      if(fstdout==kTRUE)
	cout<<"WriteCluster: padrow "<<fCurrentRow<<" pad "<<fpad << "+-"<<fpad2<<" time "<<ftime<<"+-"<<ftime2<<" charge "<<list[j].fTotalCharge<<endl;

      AliHLTTransform::Slice2Sector(fCurrentSlice,fCurrentRow,thissector,thisrow);
      AliHLTTransform::Raw2Local(xyz,thissector,thisrow,fpad,ftime);
      if(xyz[0]==0) LOG(AliHLTLog::kError,"AliHLTClustFinder","Cluster Finder")
	<<AliHLTLog::kDec<<"Zero cluster"<<ENDLOG;
      if(fNClusters >= fMaxNClusters)
	{
	  LOG(AliHLTLog::kError,"AliHLTClustFinder::WriteClusters","Cluster Finder")
	    <<AliHLTLog::kDec<<"Too many clusters"<<ENDLOG;
	  return;
	}  
      fSpacePointData[counter].fCharge = list[j].fTotalCharge;
      fSpacePointData[counter].fX = xyz[0];
      fSpacePointData[counter].fY = xyz[1];
      fSpacePointData[counter].fZ = xyz[2];
      fSpacePointData[counter].fPadRow = fCurrentRow;
      fSpacePointData[counter].fXYErr = fpad2;
      fSpacePointData[counter].fZErr  = ftime2;
      fSpacePointData[counter].fID = counter
	+((fCurrentSlice&0x7f)<<25)+((fCurrentPatch&0x7)<<22);//Uli

      fNClusters++;
      counter++;
    }
}
#endif

