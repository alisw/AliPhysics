#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TMath.h"
#include "TRandom.h"
#include "TH2.h"
#include "TAxis.h"
#include "TBits.h"
#endif

const int    kNRow   = 650; 
const int    kNCol   = 750*2;
TBits bmap(kNRow*kNCol);

TH2 *hCl,*hRw,*hInd;

void cmbHufN(int nhit);
Int_t GetIndex(int col,int row) {return kNCol*row+col;}
void  GetColRow(int index, int &col, int &row) {
  col = index%kNCol;
  row = index/kNCol;
}


void cmbHuf(int ntst=100,int nhitMn=1,int nhitMx=0.001*kNRow*kNCol, int nbn=20)
{
  hCl  = new TH2F("hCl" ,"hCl" ,nbn,nhitMn-0.5,nhitMx+0.5,kNCol,0,kNCol);
  hRw  = new TH2F("hRw" ,"hRw" ,nbn,nhitMn-0.5,nhitMx+0.5,kNRow,0,kNRow);
  hInd = new TH2F("hInd","hInd",nbn,nhitMn-0.5,nhitMx+0.5,kNCol*kNRow,0,kNCol*kNRow);
  //
  TAxis* xax = hCl->GetXaxis();
  for (int ib=1;ib<=nbn;ib++) {    
    int nh = xax->GetBinCenter(ib);
    printf("Bin%d: nh=%d\n",ib,nh);
    for (int it=0;it<ntst;it++) cmbHufN(nh);
  }
  //
}

void cmbHufN(int nhit)
{
  bmap.ResetAllBits();
  for (int i=nhit;i--;) {
    int rw = gRandom->Integer(kNRow);
    int cl = gRandom->Integer(kNCol);
    bmap.SetBitNumber(GetIndex(cl,rw));
  }
  //
  int clp=0,rwp=0,indp=0;
  for (int ind=0;ind<kNRow*kNCol;ind++) {
    if (!bmap.TestBitNumber(ind)) continue;
    int cl,rw;
    GetColRow(ind,cl,rw);
    int dcl = cl-clp;
    int drw = rw-rwp;
    int dind= ind-indp;
    //
    hCl->Fill(nhit,dcl);
    hRw->Fill(nhit,drw);
    hInd->Fill(nhit,dind);
    //
    clp = cl;
    rwp = rw;
    indp = ind;
  }

}


