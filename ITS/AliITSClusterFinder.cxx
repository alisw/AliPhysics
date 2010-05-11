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
////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Base Class used to find                                                //
// the reconstructed points for ITS                                       //
// See also AliITSClusterFinderSPD, AliITSClusterFinderSDD,               //
// AliITSClusterFinderSDD  AliITSClusterFinderV2                          //
////////////////////////////////////////////////////////////////////////////

#include "AliRun.h"
#include "AliITSClusterFinder.h"
#include "AliITSRecPoint.h"
#include "AliITSdigit.h"
#include "AliITSDetTypeRec.h"
#include "AliITSMap.h"
#include "AliITSgeomTGeo.h"
#include <TParticle.h>
#include "AliMC.h"

ClassImp(AliITSClusterFinder)

extern AliRun *gAlice;

//----------------------------------------------------------------------
AliITSClusterFinder::AliITSClusterFinder():
TObject(),
fModule(0),
fDigits(0),
fNdigits(0),
fDetTypeRec(0),
fClusters(0),
fMap(0),
fNPeaks(-1),
fNModules(AliITSgeomTGeo::GetNModules()),
fEvent(0),
fZmin(0),
fZmax(0),
fXmin(0),
fXmax(0){
    // default cluster finder
    // Input:
    //   none.
    // Output:
    //   none.
    // Return:
    //   A default constructed AliITSCulsterFinder
}
//----------------------------------------------------------------------
AliITSClusterFinder::AliITSClusterFinder(AliITSDetTypeRec* dettyp):
TObject(),
fModule(0),
fDigits(0),
fNdigits(0),
fDetTypeRec(dettyp),
fClusters(0),
fMap(0),
fNPeaks(-1),
fNModules(AliITSgeomTGeo::GetNModules()),
fEvent(0),
fZmin(0),
fZmax(0),
fXmin(0),
fXmax(0){
    // default cluster finder
    // Standard constructor for cluster finder
    // Input:
    //   AliITSsegmentation *seg  The segmentation class to be used
    //   AliITSresponse     *res  The response class to be used
    // Output:
    //   none.
    // Return:
    //   A Standard constructed AliITSCulsterFinder

}
//----------------------------------------------------------------------
AliITSClusterFinder::AliITSClusterFinder(AliITSDetTypeRec* dettyp,
					 TClonesArray *digits):
TObject(),
fModule(0),
fDigits(digits),
fNdigits(0),
fDetTypeRec(dettyp),
fClusters(0),
fMap(0),
fNPeaks(-1),
fNModules(AliITSgeomTGeo::GetNModules()),
fEvent(0),
fZmin(0),
fZmax(0),
fXmin(0),
fXmax(0){
    // default cluster finder
    // Standard + cluster finder constructor
    // Input:
    //   AliITSsegmentation *seg  The segmentation class to be used
    //   AliITSresponse     *res  The response class to be used
    //   TClonesArray    *digits  Array of digits to be used
    // Output:
    //   none.
    // Return:
    //   A Standard constructed AliITSCulsterFinder

    fNdigits = fDigits->GetEntriesFast();
}

//______________________________________________________________________
AliITSClusterFinder::AliITSClusterFinder(const AliITSClusterFinder &source) : 
  TObject(source),
  fModule(source.fModule),
  fDigits(),
  fNdigits(source.fNdigits),
  fDetTypeRec(),
  fClusters(),
  fMap(),
  fNPeaks(source.fNPeaks),
  fNModules(source.fNModules),
  fEvent(source.fEvent),
  fZmin(source.fZmin),
  fZmax(source.fZmax),
  fXmin(source.fXmin),
  fXmax(source.fXmax) 
{
  // Copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  AliError("Copy constructor not allowed\n");
}


//______________________________________________________________________
//AliITSClusterFinder& AliITSClusterFinder::operator=(const AliITSClusterFinder& /* source */){
  // Assignment operator
  // Assignment is not allowed. The method is protected to avoid misuse.
//  Fatal("= operator","Assignment operator not allowed\n");
//  return *this;
//}

//----------------------------------------------------------------------
AliITSClusterFinder::~AliITSClusterFinder(){
    // destructor cluster finder
    // Input:
    //   none.
    // Output:
    //   none.
    // Return:
    //   none.

    if(fMap) {delete fMap;}
    // Zero local pointers. Other classes own these pointers.
    fMap          = 0;
    fDigits       = 0;
    fNdigits      = 0;
    fNPeaks       = 0;
    fDetTypeRec   = 0;

}
//__________________________________________________________________________
void AliITSClusterFinder::InitGeometry(){
 //
 // Initialisation of ITS geometry
 //
  Int_t mmax=AliITSgeomTGeo::GetNModules();
  for (Int_t m=0; m<mmax; m++) {
    Int_t lay,lad,det; AliITSgeomTGeo::GetModuleId(m,lay,lad,det);
    fNdet[m] = (lad-1)*AliITSgeomTGeo::GetNDetectors(lay) + (det-1);
    fNlayer[m] = lay-1;
  }
}




//______________________________________________________________________
Bool_t AliITSClusterFinder::IsNeighbor(TObjArray *digs,Int_t i,Int_t n[])const{
    // Locagical function which checks to see if digit i has a neighbor.
    // If so, then it returns kTRUE and its neighbor index j.
    // This routine checks if the digits are side by side or one before the
    // other. Requires that the array of digits be in proper order.
    // Returns kTRUE in the following cases.
    //                 ji   0j   if kdiagonal  j0    0i
    //                 00   0i   if kdiagonal  0i    j0
    // Inputs:
    //    TObjArray *digs   Array to search for neighbors in
    //    Int_t      i      Index of digit for which we are searching for
    //                      a neighbor of.
    // Output:
    //    Int_t      j[4]   Index of one or more of the digits which is a
    //                      neighbor of digit a index i.
    // Return:
    //    Bool_t            kTRUE if a neighbor was found kFALSE otherwise.
    Int_t ix,jx,iz,jz,j;
    const Bool_t kdiagonal=kFALSE;
    Bool_t nei[4];

    // No neighbors found if array empty.
    if(digs->GetEntriesFast()<=0) return kFALSE;
    // can not be a digit with first element or elements out or range
    if(i<=0 || i>=digs->GetEntriesFast()) return kFALSE;

    for(j=0;j<4;j++){n[j] = -1;nei[j]=kFALSE;}
    ix = ((AliITSdigit*)(digs->At(i)))->GetCoord1();
    iz = ((AliITSdigit*)(digs->At(i)))->GetCoord2();
    for(j=0;j<i;j++){
        jx = ((AliITSdigit*)(digs->At(j)))->GetCoord1();
        jz = ((AliITSdigit*)(digs->At(j)))->GetCoord2();
        if(jx+1==ix && jz  ==iz){n[0] = j;nei[0] = kTRUE;}
        if(jx  ==ix && jz+1==iz){n[1] = j;nei[1] = kTRUE;}
        if(jx+1==ix && jz+1==iz){n[2] = j;nei[2] = kTRUE;}
        if(jx+1==ix && jz-1==iz){n[3] = j;nei[3] = kTRUE;}
    } // end for k
    if(nei[0]||nei[1]) return kTRUE;
    if(kdiagonal&&(nei[2]||nei[3])) return kTRUE;
    // no Neighbors found.
    return kFALSE;
}

//______________________________________________________________________
void AliITSClusterFinder::Print(ostream *os) const{
    //Standard output format for this class
    // Inputs:
    //    ostream *os   Output stream
    // Output:
    //    ostream *os   Output stream
    // Return:
    //    none.

    *os << fModule<<",";
    *os << fNdigits<<",";
    *os << fNPeaks<<endl;
}
//______________________________________________________________________
void AliITSClusterFinder::Read(istream *is)  {
    //Standard input for this class
    // Inputs:
    //    istream *is   Input stream
    // Output:
    //    istream *is   Input stream
    // Return:
    //    none.

    *is >> fModule;
    *is >> fNdigits;
    *is >> fNPeaks;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSClusterFinder &source){
    // Standard output streaming function.
    // Inputs:
    //    ostream             *os     Output stream
    //    AliITSClusterFinder &source Class to be printed
    // Output:
    //    ostream             *os     Output stream
    // Return:
    //    none.

    source.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSClusterFinder &source){
    // Standard output streaming function.
    // Inputs:
    //    istream              *is      Input stream
    //     AliITSClusterFinder &source  Class to be read in.
    // Output:
    //    istream              *is      Input stream
    // Return:
    //    none.

    source.Read(&is);
    return is;
}
//______________________________________________________________________
void AliITSClusterFinder::CheckLabels2(Int_t lab[10]) {
  //------------------------------------------------------------
  // Tries to find mother's labels
  //------------------------------------------------------------
  AliRunLoader *rl = AliRunLoader::Instance();
  if(!rl) return;
  TTree *trK=(TTree*)rl->TreeK();

  if(trK){
    Int_t nlabels =0; 
    for (Int_t i=0;i<10;i++) if (lab[i]>=0) nlabels++;
    if(nlabels == 0) return; // In case of no labels just exit


    Int_t ntracks = gAlice->GetMCApp()->GetNtrack();

    for (Int_t i=0;i<10;i++){
      Int_t label = lab[i];
      if (label>=0 && label<ntracks) {
	TParticle *part=(TParticle*)gAlice->GetMCApp()->Particle(label);

	if (part->P() < 0.02) {
	  Int_t m=part->GetFirstMother();
	  if (m<0) {	
	    continue;
	  }
	  if (part->GetStatusCode()>0) {
	    continue;
	  }
	  lab[i]=m;       
	}
	else
	  if (part->P() < 0.12 && nlabels>3) {
	    lab[i]=-2;
	    nlabels--;
	  } 
      }
      else{
	if ( (label>ntracks||label <0) && nlabels>3) {
	  lab[i]=-2;
	  nlabels--;
	} 
      }
    }  
    if (nlabels>3){
      for (Int_t i=0;i<10;i++){
	if (nlabels>3){
	  Int_t label = lab[i];
	  if (label>=0 && label<ntracks) {
	    TParticle *part=(TParticle*)gAlice->GetMCApp()->Particle(label);
	    if (part->P() < 0.1) {
	      lab[i]=-2;
	      nlabels--;
	    }
	  }
	}
      }
    }

    //compress labels -- if multi-times the same
    Int_t lab2[10];
    for (Int_t i=0;i<10;i++) lab2[i]=-2;
    for (Int_t i=0;i<10  ;i++){
      if (lab[i]<0) continue;
      for (Int_t j=0;j<10 &&lab2[j]!=lab[i];j++){
	if (lab2[j]<0) {
	  lab2[j]= lab[i];
	  break;
	}
      }
    }
    for (Int_t j=0;j<10;j++) lab[j]=lab2[j];
  
  }
}

//______________________________________________________________________
void AliITSClusterFinder::AddLabel(Int_t lab[10], Int_t label) {
  //add label to the cluster
  AliRunLoader *rl = AliRunLoader::Instance();
  TTree *trK=(TTree*)rl->TreeK();
  if(trK){
    if(label<0) return; // In case of no label just exit

    Int_t ntracks = gAlice->GetMCApp()->GetNtrack();
    if (label>ntracks) return;
    for (Int_t i=0;i<10;i++){
      //    if (label<0) break;
      if (lab[i]==label) break;
      if (lab[i]<0) {
	lab[i]= label;
	break;
      }
    }
  }
}


//______________________________________________________________________
void AliITSClusterFinder:: 
FindCluster(Int_t k,Int_t maxz,AliBin *bins,Int_t &n,Int_t *idx) {
  //------------------------------------------------------------
  // returns an array of indices of digits belonging to the cluster
  // (needed when the segmentation is not regular) 
  //------------------------------------------------------------
  if (n<200) idx[n++]=bins[k].GetIndex();
  bins[k].Use();

  if (bins[k-maxz].IsNotUsed()) FindCluster(k-maxz,maxz,bins,n,idx);
  if (bins[k-1   ].IsNotUsed()) FindCluster(k-1   ,maxz,bins,n,idx);
  if (bins[k+maxz].IsNotUsed()) FindCluster(k+maxz,maxz,bins,n,idx);
  if (bins[k+1   ].IsNotUsed()) FindCluster(k+1   ,maxz,bins,n,idx);
  /*
  if (bins[k-maxz-1].IsNotUsed()) FindCluster(k-maxz-1,maxz,bins,n,idx);
  if (bins[k-maxz+1].IsNotUsed()) FindCluster(k-maxz+1,maxz,bins,n,idx);
  if (bins[k+maxz-1].IsNotUsed()) FindCluster(k+maxz-1,maxz,bins,n,idx);
  if (bins[k+maxz+1].IsNotUsed()) FindCluster(k+maxz+1,maxz,bins,n,idx);
  */
}

//______________________________________________________________________
Bool_t AliITSClusterFinder::IsMaximum(Int_t k,Int_t max,const AliBin *bins) {
  //------------------------------------------------------------
  //is this a local maximum ?
  //------------------------------------------------------------
  UShort_t q=bins[k].GetQ();
  if (q==1023) return kFALSE;
  if (bins[k-max].GetQ() > q) return kFALSE;
  if (bins[k-1  ].GetQ() > q) return kFALSE; 
  if (bins[k+max].GetQ() > q) return kFALSE; 
  if (bins[k+1  ].GetQ() > q) return kFALSE; 
  if (bins[k-max-1].GetQ() > q) return kFALSE;
  if (bins[k+max-1].GetQ() > q) return kFALSE; 
  if (bins[k+max+1].GetQ() > q) return kFALSE; 
  if (bins[k-max+1].GetQ() > q) return kFALSE;
  return kTRUE; 
}

//______________________________________________________________________
void AliITSClusterFinder::
FindPeaks(Int_t k,Int_t max,AliBin *b,Int_t *idx,UInt_t *msk,Int_t& n) {
  //------------------------------------------------------------
  //find local maxima
  //------------------------------------------------------------
  if (n<31)
  if (IsMaximum(k,max,b)) {
    idx[n]=k; msk[n]=(2<<n);
    n++;
  }
  b[k].SetMask(0);
  if (b[k-max].GetMask()&1) FindPeaks(k-max,max,b,idx,msk,n);
  if (b[k-1  ].GetMask()&1) FindPeaks(k-1  ,max,b,idx,msk,n);
  if (b[k+max].GetMask()&1) FindPeaks(k+max,max,b,idx,msk,n);
  if (b[k+1  ].GetMask()&1) FindPeaks(k+1  ,max,b,idx,msk,n);
}

//______________________________________________________________________
void AliITSClusterFinder::
MarkPeak(Int_t k, Int_t max, AliBin *bins, UInt_t m) {
  //------------------------------------------------------------
  //mark this peak
  //------------------------------------------------------------
  UShort_t q=bins[k].GetQ();

  bins[k].SetMask(bins[k].GetMask()|m); 

  if (bins[k-max].GetQ() <= q)
     if ((bins[k-max].GetMask()&m) == 0) MarkPeak(k-max,max,bins,m);
  if (bins[k-1  ].GetQ() <= q)
     if ((bins[k-1  ].GetMask()&m) == 0) MarkPeak(k-1  ,max,bins,m);
  if (bins[k+max].GetQ() <= q)
     if ((bins[k+max].GetMask()&m) == 0) MarkPeak(k+max,max,bins,m);
  if (bins[k+1  ].GetQ() <= q)
     if ((bins[k+1  ].GetMask()&m) == 0) MarkPeak(k+1  ,max,bins,m);
}

//______________________________________________________________________
void AliITSClusterFinder::
MakeCluster(Int_t k,Int_t max,AliBin *bins,UInt_t m,AliITSRecPoint &c) {
  //------------------------------------------------------------
  //make cluster using digits of this peak
  //------------------------------------------------------------
  Float_t q=(Float_t)bins[k].GetQ();
  Int_t i=k/max, j=k-i*max;
  if(c.GetQ()<0.01){ // first entry in cluster
    fXmin=i;
    fXmax=i;
    fZmin=j;
    fZmax=j;
  }else{  // check cluster extension
    if(i<fXmin) fXmin=i;
    if(i>fXmax) fXmax=i;
    if(j<fZmin) fZmin=j;
    if(j>fZmax) fZmax=j;
  }
  c.SetQ(c.GetQ()+q);
  c.SetY(c.GetY()+i*q);
  c.SetZ(c.GetZ()+j*q);
  c.SetSigmaY2(c.GetSigmaY2()+i*i*q);
  c.SetSigmaZ2(c.GetSigmaZ2()+j*j*q);

  bins[k].SetMask(0xFFFFFFFE);
  
  if (bins[k-max].GetMask() == m) MakeCluster(k-max,max,bins,m,c);
  if (bins[k-1  ].GetMask() == m) MakeCluster(k-1  ,max,bins,m,c);
  if (bins[k+max].GetMask() == m) MakeCluster(k+max,max,bins,m,c);
  if (bins[k+1  ].GetMask() == m) MakeCluster(k+1  ,max,bins,m,c);
}
