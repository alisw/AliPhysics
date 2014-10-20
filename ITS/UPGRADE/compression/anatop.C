#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TObjArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TBits.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "../ITS/UPGRADE/AliITSUClusterPix.h"
#include "../ITS/UPGRADE/AliITSURecoLayer.h"
#include "../ITS/UPGRADE/AliITSURecoDet.h"
#include "../ITS/UPGRADE/AliITSUHit.h"
#include "../ITS/UPGRADE/AliITSUGeomTGeo.h"
#include "AliITSsegmentation.h"
#include "AliGeomManager.h"

#endif

TFile* fl=0;
TTree* tr=0;
TObjArray patterns;
TArrayF pattCount;
TArrayI sortIndex,sortColRow;
TClonesArray *modClStore;
Int_t nPatterns;

TH1F* hPattFreq = 0;
TH1F* hVolDiffFreq = 0;
TH2F* hColDiffFreq = 0;
TH2F* hColFreq = 0;
TH2F* hRowDiffFreq = 0;
TH1F* hDataVol = 0;
int lastVol=9999999,volDif=0,volID=0;
TBits* bitTop=0;


const int    kNRow   = 650*2; 
const int    kNCol   = 750*2;

const int    kNBinClVol = 6; // number of module multiplicity bins
const double kBinsClVol[kNBinClVol+1] = {0.5, 5.5, 15.5, 31.5, 65.5, 125.5, 500.5};
//
void ProcessModule(int nc);


void anatop(const char* topTreeFile="clusterTopology.root")
{
  // analyze tree with cluster topologies
  AliGeomManager::LoadGeometry("geometry.root");
  AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE);
  Int_t nLayers = gm->GetNLayers();
  //
  fl = TFile::Open(topTreeFile);
  tr = (TTree*)fl->Get("clTop");
  UShort_t mnCol,mnRow;
  tr->SetBranchAddress("tp",&bitTop);
  tr->SetBranchAddress("volID",&volID);
  tr->SetBranchAddress("mnCol",&mnCol);
  tr->SetBranchAddress("mnRow",&mnRow);
  //
  int nCl = tr->GetEntries();
  //
  hDataVol     = new TH1F("hDataVol","Data Volume", kNBinClVol,kBinsClVol);
  hVolDiffFreq = new TH1F("hVolDiffFreq","Delta Vol ID",gm->GetNChips(),0,gm->GetNChips());
  hColDiffFreq = new TH2F("hColDiffFreq","Delta Col",kNBinClVol,kBinsClVol,2*kNCol,-kNCol,kNCol);
  hColFreq     = new TH2F("hColFreq","Col"          ,kNBinClVol,kBinsClVol,kNCol,0,kNCol);
  hRowDiffFreq = new TH2F("hRowDiffFreq","Delta Row",kNBinClVol,kBinsClVol,kNRow,0,kNRow);
  //
  printf("%d clusters in the tree\n",nCl);
  modClStore = new TClonesArray("TBits");
  modClStore->ExpandCreate(5000);
  sortIndex.Set(5000);
  sortColRow.Set(5000);
  //
  int lastEv=-1,nClVol=0;
  Bool_t newEv;
  nPatterns  = 0;
  //
  if (nCl<1) {
    printf("No clusters\n");
    return;
  }
  for (int icl=0;icl<nCl;icl++) {
    if (icl%(nCl/10)==0) {
      printf("Done %d clusters of %d, seen %d patterns\n",icl,nCl,nPatterns);
    }
    newEv = kFALSE;
    tr->GetEntry(icl);
    //
    if (volID != lastVol) {
      if (volID<lastVol) {
	newEv = kTRUE;
	lastVol = 0;
	lastEv++;
      }
      // store previous module info in the frequency histos
      if (nClVol>0) ProcessModule(nClVol);
      nClVol = 0;
      //
      //
    }
    //-------------------------
    //
    *((TBits*)modClStore->At(nClVol)) = *bitTop;
    sortColRow[nClVol] = (mnRow<<16) + mnCol;
    nClVol++;

    //-------------------------
  }
  //
  printf("Sorting %d patterns\n",nPatterns);
  // sort patterns
  sortIndex.Set(nPatterns);
  TMath::Sort(nPatterns,pattCount.GetArray(),sortIndex.GetArray());
  //
  hPattFreq = new TH1F("hPattFreq","Patterns Frequency",nPatterns,0,nPatterns);
  for (int i=0;i<nPatterns;i++) {
    int ind = sortIndex[i];
    TBits* patt = (TBits*)patterns[ind];
    pattCount[ind] /= nCl;
    Double_t freq = pattCount[ind];
    hPattFreq->SetBinContent(i+1,freq);
    //
    int nrow = (patt->GetUniqueID()>>8) & 0xff;
    int ncol = patt->GetUniqueID() & 0xff;
    printf("\nPattern %5d, Freq: %e Ncol=%d Nrow=%d\n",i,freq,ncol,nrow);
    for (int ir=0;ir<nrow;ir++) {
      for (int ic=0;ic<ncol;ic++) {
	int bt = ir*ncol + ic;
	printf("%c",patt->TestBitNumber(bt) ? '*':'.');
      }
      printf("\n");
    }
  }

}


//___________________________________________
void ProcessModule(int nc)
{
  // sort clusters
  int lastCol = 0;
  int lastRow = 0;
  TMath::Sort(nc,sortColRow.GetArray(),sortIndex.GetArray(),kFALSE);
  hVolDiffFreq->Fill(volID - lastVol);
  for (int iclv=0;iclv<nc;iclv++) {
    hDataVol->Fill(nc);
    int ind = sortIndex[iclv];
    int col = sortColRow[ind]&0xffff;
    int row = (sortColRow[ind]>>16)&0xffff;
    TBits *refPat=(TBits*)modClStore->At(ind);
    //    printf("%5d %5d | %4d %4d  | %4d %4d\n",volID,lastVol, col,row, lastCol,lastRow);
    hColDiffFreq->Fill(nc, col - lastCol);
    hColFreq->Fill(nc, col);
    hRowDiffFreq->Fill(nc, row - lastRow);
    lastCol = col;
    lastRow = row;
    //
    // check if the pattern is already encountered
    Bool_t newPatt = kTRUE;
    for (int ip=0;ip<nPatterns;ip++) {
      TBits* patt = (TBits*)patterns.At(ip);
      if ( *patt == *refPat && patt->GetUniqueID()==refPat->GetUniqueID()) {
	// require that the matrix size and the bit patterns are equal
	newPatt = kFALSE;
	pattCount[ip]++;
	break;
      }
    }
    //
    if (newPatt) {
      TBits* pt = new TBits(*refPat);
      patterns.AddLast(pt);
      if (pattCount.GetSize()<nPatterns+100) {
	pattCount.Set(nPatterns+100);
      }
      pattCount[nPatterns++] = 1;
    }
    //-------------------------

  }  
  lastVol = volID;
  //
}
