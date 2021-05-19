/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// class for inherited from AliRDHFBDTD0toKpi with additional selection based on AliRDHFBDT class objects 
//
// Author: M.Cai, cai.mengke@cern.ch
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>

#include "AliRDHFCutsD0toKpiBDT.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAODVertex.h"
#include "TObjArray.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliRDHFCutsD0toKpiBDT);
/// \endcond

//--------------------------------------------------------------------------
AliRDHFCutsD0toKpiBDT::AliRDHFCutsD0toKpiBDT(const char* name) :
  AliRDHFCutsD0toKpi(name),
  fNBDTOpt(0),
  fnPtBinsBDT(1),
  fnPtBinLimitsBDT(1),
  fPtBinLimitsBDT(0),
  fBDTNames(0),
  fBDTCutGlobalIndex(1),
  fBDTCuts(0),
  fIsUpperCutBDT(0),
  fRejFraction(0)
{
  //
  // Default Constructor
  //
  fListOfBDT = new TList();
}
//--------------------------------------------------------------------------
AliRDHFCutsD0toKpiBDT::AliRDHFCutsD0toKpiBDT(const AliRDHFCutsD0toKpi &source) :
  AliRDHFCutsD0toKpi(source),
  fNBDTOpt(0),
  fnPtBinsBDT(source.GetNPtBins()),
  fnPtBinLimitsBDT(source.GetNPtBins()+1),
  fPtBinLimitsBDT(0),
  fBDTNames(0),
  fBDTCutGlobalIndex(1),
  fBDTCuts(0),
  fIsUpperCutBDT(0),
  fRejFraction(0)
{
  //
  // Copy constructor
  //
  fListOfBDT = new TList();
  if(source.GetPtBinLimits()) SetPtBinsBDT(source.GetNPtBins()+1,source.GetPtBinLimits());
  printf("Copying from a AliRDHFCutsD0toKpi, please remember to add the list of input BDT\n");
}
//--------------------------------------------------------------------------
AliRDHFCutsD0toKpiBDT::AliRDHFCutsD0toKpiBDT(const AliRDHFCutsD0toKpiBDT &source) :
  AliRDHFCutsD0toKpi(source),
  fNBDTOpt(source.fNBDTOpt),
  fnPtBinsBDT(source.fnPtBinsBDT),
  fnPtBinLimitsBDT(source.fnPtBinLimitsBDT),
  fPtBinLimitsBDT(0),
  fBDTNames(0),
  fBDTCutGlobalIndex(source.fBDTCutGlobalIndex),
  fBDTCuts(0),
  fIsUpperCutBDT(0),
  fRejFraction(0)
{
  //
  // Copy constructor
  //
  fListOfBDT = new TList();
  if(source.fPtBinLimitsBDT) SetPtBinsBDT(source.fnPtBinLimitsBDT,source.fPtBinLimitsBDT);
  if(source.fBDTCuts) SetBDTCuts(source.fBDTCutGlobalIndex,source.fBDTCuts);
  if(source.fBDTNames) SetBDTNames(source.fNBDTOpt,source.fBDTNames,source.fIsUpperCutBDT);
  if(source.fListOfBDT) SetListOfBDT(source.fListOfBDT);
  if(source.fRejFraction) SetRejFraction(source.fRejFraction);
}
//--------------------------------------------------------------------------
AliRDHFCutsD0toKpiBDT &AliRDHFCutsD0toKpiBDT::operator=(const AliRDHFCutsD0toKpiBDT &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliRDHFCutsD0toKpi::operator=(source); 
  fnPtBinsBDT=source.fnPtBinsBDT;
  fnPtBinLimitsBDT=source.fnPtBinLimitsBDT;
  if(source.fPtBinLimitsBDT) SetPtBinsBDT(source.fnPtBinLimitsBDT,source.fPtBinLimitsBDT);
  fBDTCutGlobalIndex=source.fBDTCutGlobalIndex;
  if(source.fBDTCuts) SetBDTCuts(source.fBDTCutGlobalIndex,source.fBDTCuts);
  if(source.fBDTNames) SetBDTNames(source.fNBDTOpt,source.fBDTNames,source.fIsUpperCutBDT);
  fNBDTOpt=source.fNBDTOpt;
  fListOfBDT = new TList(); if(source.fListOfBDT) SetListOfBDT(source.fListOfBDT);
  fRejFraction = new Float_t[source.fnPtBinsBDT]; SetRejFraction(source.fRejFraction);
  return *this;
}

//---------------------------------------------------------------------------
AliRDHFCutsD0toKpiBDT::~AliRDHFCutsD0toKpiBDT()
{
  //
  // Destructor
  //
  if (fPtBinLimitsBDT) {
    delete[] fPtBinLimitsBDT;
    fPtBinLimitsBDT = 0;
  }
  if (fBDTCuts) {
    delete[] fBDTCuts;
    fBDTCuts = 0;
  }
  if (fBDTNames){
    delete [] fBDTNames;
    fBDTNames =0;
  }
  if (fIsUpperCutBDT){
    delete [] fIsUpperCutBDT;
    fIsUpperCutBDT =0;
  }
  if (fListOfBDT) fListOfBDT->Clear();
}


void AliRDHFCutsD0toKpiBDT::SetBDTNames(Int_t nOpt,TString *BDTNames,Bool_t *isUpperCut)
{
  // Set the BDT names
  if(fBDTNames) { delete [] fBDTNames; fBDTNames = NULL; }
  if(nOpt!=fNBDTOpt){
    printf("WARNING: Wrong number of NBDTOpt: it has to be %d\n",fNBDTOpt);
    return;
  }
  //fnVars=nVars;
  fBDTNames = new TString[nOpt];
  fIsUpperCutBDT = new Bool_t[nOpt];
  for(Int_t i=0; i<fNBDTOpt; i++) {
    fBDTNames[i] = BDTNames[i];
    fIsUpperCutBDT[i] = isUpperCut[i];
  }
  return;
}

//---------------------------------------------------------------------------
void AliRDHFCutsD0toKpiBDT::SetPtBinsBDT(Int_t nPtBinLimits,Float_t *ptBinLimits)
{
  // Set the pt bins
  if(fPtBinLimitsBDT) { delete [] fPtBinLimitsBDT; fPtBinLimitsBDT = NULL; printf("Changing the pt bins\n"); }
  if(fnPtBinLimitsBDT != fnPtBinsBDT+1){
    cout<<"Warning: nPtBinLimitsBDT dimention "<<nPtBinLimits<<" != nPtBinsBDT+1 ("<<fnPtBinsBDT+1<<")\nSetting nPtBinsBDT to "<<nPtBinLimits-1<<endl;
    SetNPtBinsBDT(nPtBinLimits-1);
  }
  fnPtBinLimitsBDT = nPtBinLimits;
  SetBDTCutGlobalIndex();
  //cout<<"Changing also Global Index -> "<<fGlobalIndex<<endl;
  fPtBinLimitsBDT = new Float_t[fnPtBinLimitsBDT];
  for(Int_t ib=0; ib<nPtBinLimits; ib++) fPtBinLimitsBDT[ib]=ptBinLimits[ib];
  //~ printf("INFO: Reseting random-rejected factor to 0...");
  //~ for(Int_t ib=0; ib<nPtBinLimits-1; ib++) fRejFraction[ib]=1.;
  return;
}

//---------------------------------------------------------------------------
void AliRDHFCutsD0toKpiBDT::SetBDTCuts(Int_t nOpt,Int_t nPtBins,Float_t** cuts)
{
  // the cuts on BDT outputs
  if(nOpt!=fNBDTOpt) { printf("FATAL: Wrong number of NBDT: it has to be %d\n",fNBDTOpt); AliFatal("exiting"); }
  if(nPtBins!=fnPtBinsBDT) { printf("FATAL: Wrong number of pt bins: it has to be %d\n",fnPtBinsBDT); AliFatal("exiting"); }

  if(!fBDTCuts)  fBDTCuts = new Float_t[fBDTCutGlobalIndex];
  for(Int_t iv=0; iv<fNBDTOpt; iv++){
    for(Int_t ib=0; ib<fnPtBinsBDT; ib++){
      //check
      if(GetBDTCutGlobalIndex(iv,ib)>=fBDTCutGlobalIndex) { cout<<"ERROR: Overflow, exit..."<<endl; return; }
      fBDTCuts[GetBDTCutGlobalIndex(iv,ib)] = cuts[iv][ib];
    }
  }
  return;
}

//---------------------------------------------------------------------------
void AliRDHFCutsD0toKpiBDT::SetBDTCuts(Int_t glIndex,Float_t* cutsRDGlob)
{
  // the cuts on BDT outputs
  if(glIndex != fBDTCutGlobalIndex){ cout<<"FATAL: Wrong array size: it has to be "<<fBDTCutGlobalIndex<<endl; AliFatal("exiting"); }
  if(!fBDTCuts)  fBDTCuts = new Float_t[fBDTCutGlobalIndex];

  for(Int_t iGl=0;iGl<fBDTCutGlobalIndex;iGl++) fBDTCuts[iGl] = cutsRDGlob[iGl];
  return;
}

//---------------------------------------------------------------------------
void AliRDHFCutsD0toKpiBDT::SetRejFraction(Float_t *rej)
{
  fRejFraction = new Float_t[fnPtBinsBDT];
  for(Int_t i=0;i<fnPtBinsBDT;i++)
	  fRejFraction[i] = rej[i];
  return;
}

//---------------------------------------------------------------------------
void AliRDHFCutsD0toKpiBDT::SetListOfBDT(TList *l) const
{
  // the list of BDT
  if(fListOfBDT->GetEntries()>0) fListOfBDT->Clear();
  for(Int_t i=0;i<l->GetEntries();i++)
	  fListOfBDT->Add(l->At(i));
  fListOfBDT->SetOwner(1);
  
  if(fListOfBDT->GetEntries()!=fBDTCutGlobalIndex) { cout<<"FATAL: Wrong N_BDT: it has to be "<<fBDTCutGlobalIndex<<endl; AliFatal("exiting"); }
  return;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpiBDT::GetBDTCutGlobalIndex(Int_t iOpt,Int_t iPtBin) const
{
  // give the global index from variable and pt bin
  return iPtBin*fNBDTOpt+iOpt;
}

//---------------------------------------------------------------------------
Double_t AliRDHFCutsD0toKpiBDT::GetRDHFVarsForSel(AliAODRecoDecayHF2Prong *d, AliAODEvent *aod, TString VarName, Int_t isD0bar) const
{
  if(VarName=="topo1"){
    Double_t diffIP, errdiffIP;
	d->Getd0MeasMinusExpProng(0,aod->GetMagneticField(),diffIP,errdiffIP);
	return diffIP/errdiffIP;
  }
  if(VarName=="topo2"){
    Double_t diffIP, errdiffIP;
	d->Getd0MeasMinusExpProng(1,aod->GetMagneticField(),diffIP,errdiffIP);
	return diffIP/errdiffIP;
  }
  if(VarName=="lxy")	return d->AliAODRecoDecay::DecayLengthXY(d->GetPrimaryVtx());
  if(VarName=="nlxy")	return d->AliAODRecoDecay::DecayLengthXY(d->GetPrimaryVtx())/TMath::Sqrt(d->GetSecondaryVtx()->Error2DistanceXYToVertex(d->GetPrimaryVtx()));
  if(VarName=="d0d0")	return d->Getd0Prong(0)*d->Getd0Prong(1);
  if(VarName=="cosp")	return d->CosPointingAngle();
  if(VarName=="dca")	return d->GetDCA();
  if(VarName=="cospxy")	return d->GetDCA();
  if(VarName=="d0k")	return d->Getd0Prong(0);
  if(VarName=="d0pi")	return d->Getd0Prong(1);
  if(VarName=="cosstar"){
	  if(!isD0bar)	return d->CosThetaStarD0();
	  else			return d->CosThetaStarD0bar();
  }
  cout<<"WARNING: Unknown variable name "<< VarName <<" used, plz check your input.."<<endl;
  return -99999999;
  
}

//---------------------------------------------------------------------------
Double_t AliRDHFCutsD0toKpiBDT::GetBDTResponse(AliAODRecoDecayHF2Prong *d, AliAODEvent *aod, Int_t iOpt, Int_t isD0bar) const
{
  // Read output of the BDT
  Int_t ptbin = PtBinBDT(d->Pt());
  
  AliRDHFBDT *thisBDT = (AliRDHFBDT*)fListOfBDT->At(GetBDTCutGlobalIndex(iOpt,ptbin));
  if(!thisBDT) return -1; // rejected failed reading BDT
  
  Int_t Nvar = thisBDT->GetNVars();
  std::vector<Double_t> BDTClsVar; BDTClsVar.resize(Nvar);
  for(Int_t ivar=0;ivar<Nvar;ivar++) BDTClsVar[ivar]=GetRDHFVarsForSel(d,aod,thisBDT->GetVarName(ivar),isD0bar);
  return thisBDT->GetResponse(BDTClsVar);
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpiBDT::IsSelectedBDT(AliAODRecoDecayHF2Prong *d, AliAODEvent *aod) const
{
  //
  // Apply selection
  //

  //~ fIsSelectedCuts=0;
  //~ fIsSelectedPID=0;

  if(!fBDTCuts){
    cout<<"BDT Cut matrix not inizialized. Exit..."<<endl;
    return 0;
  }
  //PrintAll();
  if(!d) { cout<<"AliAODRecoDecayHF2Prong null"<<endl; return 0; }

  Int_t ptbin = PtBinBDT(d->Pt());
  if(ptbin==-1) return 0; // rejected pT out of range
  
  Int_t okD0=1,okD0bar=1;
  
  AliAODVertex *origownvtx=0x0;
  // data recalculate vertex w/o daughters (for ImpPar)
  if(fRemoveDaughtersFromPrimary && !fUseMCVertex) {
	if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());
	if(!RecalcOwnPrimaryVtx(d,aod)) {
	  CleanOwnPrimaryVtx(d,aod,origownvtx);
	  return 0;
    }
  }
  // MC
  if(fUseMCVertex) {
	if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());
	if(!SetMCPrimaryVtx(d,aod)) {
	  CleanOwnPrimaryVtx(d,aod,origownvtx);
	  return 0;
	}
  }
  for(Int_t i=0;i<fNBDTOpt;i++){
	if(fIsUpperCutBDT[i]){
	  if(GetBDTResponse(d,aod,i,0)>fBDTCuts[GetBDTCutGlobalIndex(i,ptbin)]) okD0=0;
	  if(GetBDTResponse(d,aod,i,1)>fBDTCuts[GetBDTCutGlobalIndex(i,ptbin)]) okD0bar=0;
	}
	else{
	  if(GetBDTResponse(d,aod,i,0)<=fBDTCuts[GetBDTCutGlobalIndex(i,ptbin)]) okD0=0;
	  if(GetBDTResponse(d,aod,i,1)<=fBDTCuts[GetBDTCutGlobalIndex(i,ptbin)]) okD0bar=0;
	}
  }
  // unset recalculated primary vertex when not needed any more
  CleanOwnPrimaryVtx(d,aod,origownvtx);
  
  
  if(okD0&&!okD0bar) return 1;
  if(!okD0&&okD0bar) return 2;
  if(okD0&&okD0bar)	return 3;
  return 0;
  
}

//___________________________________________________________
void AliRDHFCutsD0toKpiBDT::PrintAll() const {
  //
  // print all cuts values
  // 
  AliRDHFCutsD0toKpi::PrintAll();
  for(Int_t i=0;i<fnPtBinsBDT;i++)
	printf("Sampling fraction: %f	",fRejFraction[i]);
  if(fListOfBDT->GetEntries()>0){
	printf("\nList of BDT:\n");
	for(Int_t i=0;i<fnPtBinsBDT;i++){
		for(Int_t j=0;j<fNBDTOpt;j++){
			printf("%s		",fListOfBDT->At(GetBDTCutGlobalIndex(j,i))->GetName());
			printf("cut: Resp %c %.4f\n",!fIsUpperCutBDT[j]?'>':'<',fBDTCuts[GetBDTCutGlobalIndex(j,i)]);
		}
	}
  }
  else printf("No BDT added...\n");
  return;
}
