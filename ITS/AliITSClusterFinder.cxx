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
//
// Base Class used to find
// the reconstructed points for ITS
// See also AliITSClusterFinderSPD, AliITSClusterFinderSDD, 
// AliITSClusterFinderSDD
//

#include "AliITSClusterFinder.h"
#include "AliITSdigitSPD.h"
#include "AliITSdigitSDD.h"
#include "AliITSdigitSSD.h"
#include "AliITSMap.h"
#include "AliRun.h"
#include "AliITS.h"

ClassImp(AliITSClusterFinder)

//----------------------------------------------------------------------
AliITSClusterFinder::AliITSClusterFinder():
TObject(),
fDebug(0),
fModule(0),
fDigits(0),
fNdigits(0),
fResponse(0),
fSegmentation(0),
fClusters(0),
fNRawClusters(0),
fMap(0),
fNperMax(0),
fDeclusterFlag(0),
fClusterSize(0),
fNPeaks(-1){
    // default cluster finder
    // Input:
    //   none.
    // Output:
    //   none.
    // Return:
    //   A default constructed AliITSCulsterFinder
}
//----------------------------------------------------------------------
AliITSClusterFinder::AliITSClusterFinder(AliITSsegmentation *seg, 
                                         AliITSresponse *res):
TObject(),
fDebug(0),
fModule(0),
fDigits(0),
fNdigits(0),
fResponse(res),
fSegmentation(seg),
fClusters(0),
fNRawClusters(0),
fMap(0),
fNperMax(0),
fDeclusterFlag(0),
fClusterSize(0),
fNPeaks(-1){
    // Standard constructor for cluster finder
    // Input:
    //   AliITSsegmentation *seg  The segmentation class to be used
    //   AliITSresponse     *res  The response class to be used
    // Output:
    //   none.
    // Return:
    //   A Standard constructed AliITSCulsterFinder

    SetNperMax();
    SetClusterSize();
    SetDeclusterFlag();
}
//----------------------------------------------------------------------
AliITSClusterFinder::AliITSClusterFinder(AliITSsegmentation *seg, 
                                         AliITSresponse *response, 
                                         TClonesArray *digits):
TObject(),
fDebug(0),
fModule(0),
fDigits(digits),
fNdigits(0),
fResponse(response),
fSegmentation(seg),
fClusters(0),
fNRawClusters(0),
fMap(0),
fNperMax(0),
fDeclusterFlag(0),
fClusterSize(0),
fNPeaks(-1){
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
    SetNperMax();
    SetClusterSize();
    SetDeclusterFlag();
}
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
    fSegmentation = 0;
    fResponse     = 0;
    fMap          = 0;
    fDigits       = 0;
    fNdigits      = 0;
    fNRawClusters = 0;
    fNperMax      = 0;
    fDeclusterFlag= 0;
    fClusterSize  = 0;
    fNPeaks       = 0;
}
//__________________________________________________________________________
AliITSClusterFinder::AliITSClusterFinder(const AliITSClusterFinder &source) :
    TObject(source){
    //     Copy Constructor
    // Input:
    //   AliITSClusterFinder &source  The class which will become a copy of 
    //                                this class
    // Output:
    //   none.
    // Return:
    //   A copy of this class

    if(&source == this) return;
    *this = source;
    return;
}
//______________________________________________________________________
AliITSClusterFinder& AliITSClusterFinder::operator=(const AliITSClusterFinder 
                                                    &source) {
    //    Assignment operator
    // Input:
    //   AliITSClusterFinder &source  The class which will become a copy of 
    //                                this class
    // Output:
    //   none.
    // Return:
    //   A copy of this class

    if(&source == this) return *this;
    this->fDigits = source.fDigits;
    this->fNdigits = source.fNdigits;
    this->fResponse = source.fResponse;
    this->fSegmentation = source.fSegmentation;
    this->fNRawClusters = source.fNRawClusters;
    this->fMap = source.fMap;
    this->fNperMax = source.fNperMax;
    this->fDeclusterFlag = source.fDeclusterFlag;
    this->fClusterSize = source.fClusterSize;
    this->fNPeaks = source.fNPeaks;
    return *this;
}
//----------------------------------------------------------------------
void AliITSClusterFinder::AddCluster(Int_t branch, AliITSRawCluster *c){
    // Add a raw cluster copy to the list
    // Input:
    //   Int_t       branch  The branch to which the cluster is to be added to
    //   AliITSRawCluster *c The cluster to be added to the array of clusters
    // Output:
    //   none.
    // Return:
    //   none.

    AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");
    iTS->AddCluster(branch,c); 
    fNRawClusters++;
}
//----------------------------------------------------------------------
void AliITSClusterFinder::AddCluster(Int_t branch, AliITSRawCluster *c, 
				     AliITSRecPoint &rp){
    // Add a raw cluster copy to the list and the RecPoint
    // Input:
    //   Int_t       branch  The branch to which the cluster is to be added to
    //   AliITSRawCluster *c The cluster to be added to the array of clusters
    //   AliITSRecPoint  &rp The RecPoint to be added to the array of RecPoints
    // Output:
    //   none.
    // Return:
    //   none.

    AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");
    iTS->AddCluster(branch,c); 
    fNRawClusters++;
    iTS->AddRecPoint(rp); 
}
//______________________________________________________________________
void AliITSClusterFinder::FindRawClusters(Int_t module){
    // Default Cluster finder.
    // Input:
    //   Int_t module   Module number for which culster are to be found.
    // Output:
    //   none.
    // Return:
    //   none.
    const Int_t kelms = 10;
    Int_t ndigits = fDigits->GetEntriesFast();
    TObjArray *digs = new TObjArray(ndigits);
    TObjArray *clusts = new TObjArray(ndigits); // max # cluster
    TObjArray *clust0=0; // A spacific cluster of digits
    TObjArray *clust1=0; // A spacific cluster of digits
    AliITSdigit *dig=0; // locat pointer to a digit
    Int_t i=0,nc=0,j[4],k,k2=0;

    // Copy all digits for this module into a local TObjArray.
    for(i=0;i<ndigits;i++) digs->AddAt(new AliITSdigit(*(GetDigit(i))),i);
    digs->Sort();
    // First digit is a cluster.
    i  = 0;
    nc = 0;
    clusts->AddAt(new TObjArray(kelms),nc);
    clust0 = (TObjArray*)(clusts->At(nc));
    clust0->AddAtFree(digs->At(i)); // move owner ship from digs to clusts
    nc++;
    for(i=1;i<ndigits;i++){
        if(IsNeighbor(digs,i,j)){
            dig = (AliITSdigit*)(digs->At(j[0]));
            // Add to existing cluster. Find which cluster this digis 
            for(k=0;k<nc;k++){
                clust0 = ((TObjArray*)(clusts->At(k)));
                if(clust0->IndexOf(dig)>=0) break;
            } // end for k
            if(k>=nc){
                Fatal("FindRawClusters","Digit not found as expected");
            } // end if
            if(j[1]>=0){
                dig = (AliITSdigit*)(digs->At(j[1]));
                // Add to existing cluster. Find which cluster this digis 
                for(k2=0;k2<nc;k2++){
                    clust1 = ((TObjArray*)(clusts->At(k2)));
                    if(clust1->IndexOf(dig)>=0) break;
                } // end for k2
                if(k2>=nc){
                    Fatal("FindRawClusters","Digit not found as expected");
                } // end if
            } // end if j[1]>=0
            // Found cluster with neighboring digits add this one to it.
            if(clust0==clust1){ // same cluster
                clust0->AddAtFree(digs->At(i));
                clust0 = 0; // finished with cluster. zero for safty
                clust1 = 0; // finished wit hcluster. zero for safty
            }else{ // two different clusters which need to be merged.
                clust0->AddAtFree(digs->At(i)); // Add digit to this cluster.
                for(k=0;k<clust1->GetEntriesFast();k++){
                    // move clust1 into clust0
                    //move digit to this cluster
                    clust0->AddAtFree(clust1->At(k));
                    clust1->AddAt(0,k); // zero this one
                } // end for k
                delete clust1;
                clusts->AddAt(0,k2); // zero array of clusters element clust1
                clust0 = 0; // finished with cluster. zero for safty
                clust1 = 0; // finished wit hcluster. zero for safty
            } // end if clust0==clust1
        }else{// New cluster
            clusts->AddAt(new TObjArray(kelms),nc);
            clust0 = ((TObjArray*)(clusts->At(nc)));
            // move owner ship from digs to clusts
            clust0->AddAtFree(digs->At(i));
            clust0 = 0; // finished with cluster. zero for safty
            nc++;
        } // End if IsNeighbor
    } // end for i
    // There are now nc clusters in clusts. Each element of clust is an
    // array of digits which are clustered together.

    // For each cluster call detector specific CreateRecPoints
    for(i=0;i<nc;i++) CreateRecPoints((TObjArray*)(clusts->At(i)),module);

    // clean up at the end.
    for(i=0;i<nc;i++){ 
        clust0 =(TObjArray*)(clusts->At(i));
        // Digits deleted below, so zero this TObjArray
        for(k=0;k<clust0->GetEntriesFast();k++) clust0->AddAt(0,k);
        delete clust0; // Delete this TObjArray
        clusts->AddAt(0,i); // Contents deleted above, so zero it.
    } // end for i
    delete clusts; // Delete this TObjArray/
    // Delete the digits then the TObjArray which containted them.
    for(i=0;i<ndigits;i++) delete ((AliITSdigit*)(digs->At(i)));
    delete digs;
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
void AliITSClusterFinder::Print(ostream *os){
    //Standard output format for this class
    // Inputs:
    //    ostream *os   Output stream
    // Output:
    //    ostream *os   Output stream
    // Return:
    //    none.

    *os << fDebug<<",";
    *os << fModule<<",";
    *os << fNdigits<<",";
    *os << fNRawClusters<<",";
    *os << fNperMax<<",";
    *os << fDeclusterFlag<<",";
    *os << fClusterSize<<",";
    *os << fNPeaks<<endl;
}
//______________________________________________________________________
void AliITSClusterFinder::Read(istream *is){
    //Standard input for this class
    // Inputs:
    //    istream *is   Input stream
    // Output:
    //    istream *is   Input stream
    // Return:
    //    none.

    *is >> fDebug;
    *is >> fModule;
    *is >> fNdigits;
    *is >> fNRawClusters;
    *is >> fNperMax;
    *is >> fDeclusterFlag;
    *is >> fClusterSize;
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
