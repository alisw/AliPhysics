// $Header$

//__________________________________________________________________________
// TPCDigitsInfo
//
//


#include <TFile.h>
#include <TStopwatch.h>
#include <Reve/TTreeTools.h>

#include "TPCDigitsInfo.h"


using namespace Reve;
using namespace Alieve;
using namespace std;

void TPCSeg::Dump() const
{
  printf("TPCSeg: pw %f, pl %f, R %f, nRows %d, nMaxPads %d \n",
	 fPadWidth, fPadLength,fRlow,fNRows,fNMaxPads);
}


ClassImp(TPCDigitsInfo)

/**************************************************************************/

void TPCDigitsInfo::Init()
{
  fTree = 0;
  fParameter= 0;
}


TPCDigitsInfo::~TPCDigitsInfo()
{
  delete fParameter;
  delete fTree;
}

/**************************************************************************/

void TPCDigitsInfo::SetData(AliTPCParam* par, TTree* digits)
{ 
 
  static const Exc_t eH("TPCDigitsInfo::SetData");

  fParameter = par;
  fTree = digits;

  TStopwatch* sw = new  TStopwatch();
  sw->Start();
  fTree->LoadBaskets();
  sw->Stop();
  // printf("TPCDigitsInfo::SetData timer %f\n", sw->RealTime());
  // sw->Dump();

  fSegEnt.assign(72,-1);
  AliSimDigits *digit = &fSimDigits;
  fTree->GetBranch("Segment")->SetAddress(&digit);
  
  Int_t sbr=(Int_t)fTree->GetEntries();
  for (Int_t ent=0; ent<sbr; ent++) {
    fTree->GetEntry(ent);
    Int_t s, row;
    par->AdjustSectorRow(digit->GetID(),s,row);
    // printf("found entry %d in sec %d row %d \n",ent, s, row);

    if(row == 0) fSegEnt[s] = ent;
  }


  // read TPC Seg data
  fInnSeg.fPadWidth   = par->GetInnerPadPitchWidth();
  fInnSeg.fPadLength  = par->GetInnerPadPitchLength();
  fInnSeg.fRlow       = par->GetPadRowRadiiLow(0);
  fInnSeg.fNRows      = par->GetNRowLow();
  fInnSeg.fNMaxPads   = par->GetNPadsLow(fInnSeg.fNRows - 1);


  fOut1Seg.fPadWidth   = par->GetOuterPadPitchWidth();
  fOut1Seg.fPadLength  = par->GetOuter1PadPitchLength();
  fOut1Seg.fRlow       = par->GetPadRowRadiiUp(0);
  fOut1Seg.fNRows      = par->GetNRowUp1();
  fOut1Seg.fNMaxPads   = par->GetNPadsUp(fOut1Seg.fNRows-1);
 

  fOut2Seg.fPadWidth   = par->GetOuterPadPitchWidth();
  fOut2Seg.fPadLength  = par->GetOuter2PadPitchLength();
  fOut2Seg.fRlow       = par->GetPadRowRadiiUp(fOut1Seg.fNRows);
  fOut2Seg.fNRows      = par->GetNRowUp() - fOut1Seg.fNRows;
  fOut2Seg.fNMaxPads   = par->GetNPadsUp(par->GetNRowUp()-1);


  // set stepsize array
  Int_t k, npads;
  // Inn
  k=0, npads = par->GetNPadsLow(0);
  for (int row = 0; row < fInnSeg.fNRows ;row++ ){
    if (par->GetNPadsLow(row) > npads){
      npads = par->GetNPadsLow(row);
      fInnSeg.fStepY[k]=row* fInnSeg.fPadLength +fInnSeg.fRlow;
      k++;
    }
  }
  fInnSeg.fNsteps=k;
  // Out1 seg 
  k=0; npads = par->GetNPadsUp(0);
  for (int row = 0; row < fOut1Seg.fNRows ;row++ ){
    if (par->GetNPadsUp(row) > npads){
      npads = par->GetNPadsUp(row);
      fOut1Seg.fStepY[k]=row* fOut1Seg.fPadLength + fOut1Seg.fRlow ;
      k++;
    }
  }
  fOut1Seg.fNsteps=k;
  // Out2 seg
  k=0; npads = par->GetNPadsUp(fOut1Seg.fNRows);
  for (int row = fOut1Seg.fNRows; row < par->GetNRowUp() ;row++ ){
    if (par->GetNPadsUp(row) > npads){
      npads = par->GetNPadsUp(row);
      fOut2Seg.fStepY[k]=(row - fOut1Seg.fNRows)* fOut2Seg.fPadLength + fOut2Seg.fRlow ;
      k++;
    }
  }
  fOut2Seg.fNsteps=k;
}

/**************************************************************************/

void TPCDigitsInfo::Print(Option_t* ) const
{
  fInnSeg.Dump();
  fOut1Seg.Dump();
  fOut2Seg.Dump();
}
