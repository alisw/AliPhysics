/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//                Dielectron HF                                  //
//                                                                       //
//                                                                       //
/*
Detailed description


*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TVectorD.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <THnSparse.h>
#include <TAxis.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>

#include <AliVParticle.h>
#include <AliLog.h>

#include "AliDielectron.h"
#include "AliDielectronHelper.h"
#include "AliDielectronMC.h"
#include "AliDielectronPair.h"
#include "AliDielectronSignalMC.h"

#include "AliDielectronHistos.h"
#include "AliDielectronHF.h"

ClassImp(AliDielectronHF)

AliDielectronHF::AliDielectronHF() :
  TNamed(),
  fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
  fArrPairType(),
  fPairType(kSeOnlyOS),
  fSignalsMC(0x0),
  fVarCutType(new TBits(kMaxCuts)),
  fAxes(kMaxCuts),
  fHasMC(kFALSE),
  fStepGenerated(kFALSE),
  fEventArray(kFALSE),
  fRefObj(1)
{
  //
  // Default Constructor
  //
  for (Int_t i=0; i<kMaxCuts; ++i){
    fVarCuts[i]=0;
    //    fVarCutType[i]=0;
    fBinType[i]=kStdBin;
  }
  fAxes.SetOwner(kTRUE);
  fRefObj.SetOwner(kTRUE);
  fArrPairType.SetOwner(kTRUE);
}

//______________________________________________
AliDielectronHF::AliDielectronHF(const char* name, const char* title) :
  TNamed(name, title),
  fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
  fArrPairType(),
  fPairType(kSeOnlyOS),
  fSignalsMC(0x0),
  fVarCutType(new TBits(kMaxCuts)),
  fAxes(kMaxCuts),
  fHasMC(kFALSE),
  fStepGenerated(kFALSE),
  fEventArray(kFALSE),
  fRefObj(1)
{
  //
  // Named Constructor
  //
  for (Int_t i=0; i<kMaxCuts; ++i){
    fVarCuts[i]=0;
    //    fVarCutType[i]=0;
    fBinType[i]=kStdBin;
  }
  fAxes.SetOwner(kTRUE);
  fRefObj.SetOwner(kTRUE);
  fArrPairType.SetOwner(kTRUE);
}

//______________________________________________
AliDielectronHF::~AliDielectronHF()
{
  //
  // Default Destructor
  //
  if(fUsedVars)   delete fUsedVars;
  if(fVarCutType) delete fVarCutType;
  fAxes.Delete();
  fRefObj.Delete();
  fArrPairType.Delete();
}

//_____________________________________________________________________________
void AliDielectronHF::UserProfile(const char* histClass, UInt_t valTypeP,
				      const TVectorD * const binsX,
				      UInt_t valTypeX, TString option, UInt_t valTypeW)
{
  //
  // Histogram creation 1D case with arbitraty binning X
  // the TVectorD is assumed to be surplus after the creation and will be deleted!!!
  //

  TH1 *hist=0x0;
  if(valTypeP==AliDielectronHistos::kNoProfile)
    hist=new TH1F("","",binsX->GetNrows()-1,binsX->GetMatrixArray());
  else {
    TString opt=""; Double_t pmin=0., pmax=0.;
    if(!option.IsNull()) {
      TObjArray *arr=option.Tokenize(";");
      arr->SetOwner();
      opt=((TObjString*)arr->At(0))->GetString();
      if(arr->GetEntriesFast()>1) pmin=(((TObjString*)arr->At(1))->GetString()).Atof();
      if(arr->GetEntriesFast()>2) pmax=(((TObjString*)arr->At(2))->GetString()).Atof();
      delete arr;
    }
    hist=new TProfile("","",binsX->GetNrows()-1,binsX->GetMatrixArray(),pmin,pmax,opt.Data());
    //      printf(" name %s PROFILE options: pmin %.1f pmax %.1f err %s \n",name,((TProfile*)hist)->GetYmin(),((TProfile*)hist)->GetYmax(),((TProfile*)hist)->GetErrorOption() );
  }

  // store variales in axes
  UInt_t valType[4] = {0};
  valType[0]=valTypeX;     valType[1]=valTypeP;
  AliDielectronHistos::StoreVariables(hist, valType);
  hist->SetUniqueID(valTypeW); // store weighting variable

  for(Int_t i=0; i<4; i++)   fUsedVars->SetBitNumber(valType[i],kTRUE);
  fUsedVars->SetBitNumber(valTypeW,kTRUE);

  // adapt the name and title of the histogram in case they are empty
  AliDielectronHistos::AdaptNameTitle(hist, histClass);
  hist->SetName(Form("HF_%s",hist->GetName()));

  fRefObj.AddLast(hist);
  delete binsX;
}

//_____________________________________________________________________________
void AliDielectronHF::UserProfile(const char* histClass, UInt_t valTypeP,
				      const TVectorD * const binsX, const TVectorD * const binsY,
				      UInt_t valTypeX, UInt_t valTypeY, TString option, UInt_t valTypeW)
{
  //
  // Histogram creation 2D case with arbitraty binning X and Y
  // the TVectorD is assumed to be surplus after the creation and will be deleted!!!
  //

  TH1 *hist=0x0;
  if(valTypeP==AliDielectronHistos::kNoProfile) {
    hist=new TH2F("","",
		  binsX->GetNrows()-1,binsX->GetMatrixArray(),
		  binsY->GetNrows()-1,binsY->GetMatrixArray());
  }
  else  {
    TString opt=""; Double_t pmin=0., pmax=0.;
    if(!option.IsNull()) {
      TObjArray *arr=option.Tokenize(";");
      arr->SetOwner();
      opt=((TObjString*)arr->At(0))->GetString();
      if(arr->GetEntriesFast()>1) pmin=(((TObjString*)arr->At(1))->GetString()).Atof();
      if(arr->GetEntriesFast()>2) pmax=(((TObjString*)arr->At(2))->GetString()).Atof();
      delete arr;
    }
    hist=new TProfile2D("","",
			binsX->GetNrows()-1,binsX->GetMatrixArray(),
			binsY->GetNrows()-1,binsY->GetMatrixArray());
    ((TProfile2D*)hist)->BuildOptions(pmin,pmax,opt.Data());
  }

  // store variales in axes
  UInt_t valType[4] = {0};
  valType[0]=valTypeX;     valType[1]=valTypeY; valType[2]=valTypeP;
  AliDielectronHistos::StoreVariables(hist, valType);
  hist->SetUniqueID(valTypeW); // store weighting variable

  for(Int_t i=0; i<4; i++)   fUsedVars->SetBitNumber(valType[i],kTRUE);
  fUsedVars->SetBitNumber(valTypeW,kTRUE);

  // adapt the name and title of the histogram in case they are empty
  AliDielectronHistos::AdaptNameTitle(hist, histClass);
  hist->SetName(Form("HF_%s",hist->GetName()));

  fRefObj.AddLast(hist);
  delete binsX;
  delete binsY;
}

//_____________________________________________________________________________
void AliDielectronHF::UserProfile(const char* histClass, UInt_t valTypeP,
				      const TVectorD * const binsX, const TVectorD * const binsY, const TVectorD * const binsZ,
				      UInt_t valTypeX, UInt_t valTypeY, UInt_t valTypeZ, TString option, UInt_t valTypeW)
{
  //
  // Histogram creation 3D case with arbitraty binning X, Y, Z
  // the TVectorD is assumed to be surplus after the creation and will be deleted!!!
  //
  TH1 *hist=0x0;
  if(valTypeP==AliDielectronHistos::kNoProfile) {
    hist=new TH3F("","",
		  binsX->GetNrows()-1,binsX->GetMatrixArray(),
		  binsY->GetNrows()-1,binsY->GetMatrixArray(),
		  binsZ->GetNrows()-1,binsZ->GetMatrixArray());
  }
  else {
    TString opt=""; Double_t pmin=0., pmax=0.;
    if(!option.IsNull()) {
      TObjArray *arr=option.Tokenize(";");
      arr->SetOwner();
      opt=((TObjString*)arr->At(0))->GetString();
      if(arr->GetEntriesFast()>1) pmin=(((TObjString*)arr->At(1))->GetString()).Atof();
      if(arr->GetEntriesFast()>2) pmax=(((TObjString*)arr->At(2))->GetString()).Atof();
      delete arr;
    }
    hist=new TProfile3D("","",
			binsX->GetNrows()-1,binsX->GetMatrixArray(),
			binsY->GetNrows()-1,binsY->GetMatrixArray(),
			binsZ->GetNrows()-1,binsZ->GetMatrixArray());
    ((TProfile3D*)hist)->BuildOptions(pmin,pmax,opt.Data());
  }

  // store variales in axes
  UInt_t valType[4] = {0};
  valType[0]=valTypeX;     valType[1]=valTypeY;     valType[2]=valTypeZ;     valType[3]=valTypeP;
  AliDielectronHistos::StoreVariables(hist, valType);
  hist->SetUniqueID(valTypeW); // store weighting variable

  for(Int_t i=0; i<4; i++)   fUsedVars->SetBitNumber(valType[i],kTRUE);
  fUsedVars->SetBitNumber(valTypeW,kTRUE);

  // adapt the name and title of the histogram in case they are empty
  AliDielectronHistos::AdaptNameTitle(hist, histClass);
  hist->SetName(Form("HF_%s",hist->GetName()));

  fRefObj.AddLast(hist);
  delete binsX;
  delete binsY;
  delete binsZ;
}

//_____________________________________________________________________________
void AliDielectronHF::UserSparse(const char* histClass, Int_t ndim, TObjArray *limits, UInt_t *vars, UInt_t valTypeW)
{
  //
  // THnSparse creation with non-linear binning
  //

  THnSparseF *hist=0;
  Int_t *bins=new Int_t[ndim];
  // get number of bins
  for(Int_t idim=0 ;idim<ndim; idim++) {
    TVectorD *vec = (TVectorD*) limits->At(idim);
    bins[idim]=vec->GetNrows()-1;
  }

  hist=new THnSparseF("",histClass, ndim, bins, 0x0, 0x0);
  delete [] bins;

  // set binning
  for(Int_t idim=0 ;idim<ndim; idim++) {
    TVectorD *vec = (TVectorD*) limits->At(idim);
    hist->SetBinEdges(idim,vec->GetMatrixArray());
  }

  // store variales in axes
  AliDielectronHistos::StoreVariables(hist, vars);
  hist->SetUniqueID(valTypeW); // store weighting variable

  // store which variables are used
  for(Int_t i=0; i<20; i++)   fUsedVars->SetBitNumber(vars[i],kTRUE);
  fUsedVars->SetBitNumber(valTypeW,kTRUE);

  // adapt the name and title of the histogram in case they are empty
  TString name;
  for(Int_t iv=0; iv < ndim; iv++) name+=Form("%s_",AliDielectronVarManager::GetValueName(vars[iv]));
  name.Resize(name.Length()-1);
  hist->SetName(Form("HF_%s",name.Data()));

  fRefObj.AddLast(hist);
  delete limits;

}

//________________________________________________________________
void AliDielectronHF::AddCutVariable(AliDielectronVarManager::ValueTypes type,
				     Int_t nbins, Double_t min, Double_t max, Bool_t log, Bool_t leg, EBinType btype)
{
  //
  // Add a variable to the mixing handler
  //

  // limit number of variables to kMaxCuts
  if (fAxes.GetEntriesFast()>=kMaxCuts) return;

  TVectorD *binLimits=0x0;
  if (!log) binLimits=AliDielectronHelper::MakeLinBinning(nbins,min,max);
  else binLimits=AliDielectronHelper::MakeLogBinning(nbins,min,max);
  if (!binLimits) return;

  Int_t size=fAxes.GetEntriesFast();
  fVarCuts[size]=(UShort_t)type;
  //  fVarCutType[size]=leg;
  fVarCutType->SetBitNumber(size,leg);
  fAxes.Add(binLimits->Clone());
  fBinType[size]=btype;
  fUsedVars->SetBitNumber(type,kTRUE);
}

//________________________________________________________________
void AliDielectronHF::AddCutVariable(AliDielectronVarManager::ValueTypes type,
                                             const char* binLimitStr, Bool_t leg, EBinType btype)
{
  //
  // Add a variable to the mixing handler with arbitrary binning
  //

  // limit number of variables to kMaxCuts
  if (fAxes.GetEntriesFast()>=kMaxCuts) return;

  TVectorD *binLimits=AliDielectronHelper::MakeArbitraryBinning(binLimitStr);
  if (!binLimits) return;

  Int_t size=fAxes.GetEntriesFast();
  fVarCuts[size]=(UShort_t)type;
  //  fVarCutType[size]=leg;
  fVarCutType->SetBitNumber(size,leg);
  fAxes.Add(binLimits);
  fBinType[size]=btype;
  fUsedVars->SetBitNumber(type,kTRUE);
}

//________________________________________________________________
void AliDielectronHF::AddCutVariable(AliDielectronVarManager::ValueTypes type,
                                             TVectorD * binLimits, Bool_t leg, EBinType btype)
{
  //
  // Add a variable to the mixing handler with a vector
  // the TVectorD is assumed to be surplus after the creation and will be deleted!!!
  //

  // limit number of variables to kMaxCuts
  if (fAxes.GetEntriesFast()>=kMaxCuts) return;

  if (!binLimits) return;

  Int_t size=fAxes.GetEntriesFast();
  fVarCuts[size]=(UShort_t)type;
  //  fVarCutType[size]=leg;
  fVarCutType->SetBitNumber(size,leg);
  fAxes.Add(binLimits);
  fBinType[size]=btype;
  fUsedVars->SetBitNumber(type,kTRUE);
}

//______________________________________________
void AliDielectronHF::Fill(Int_t label1, Int_t label2, Int_t nSignal)
{
  //
  // fill the pure MC part of the container starting from a pair of 2 particles (part1 and part2 are legs)
  //
  // fill only if we have asked for these steps
  if(!fStepGenerated || fEventArray) return;

  AliVParticle* part1 = AliDielectronMC::Instance()->GetMCTrackFromMCEvent(label1);
  AliVParticle* part2 = AliDielectronMC::Instance()->GetMCTrackFromMCEvent(label2);
  if(!part1 || !part2) return;

  AliDielectronMC* dieMC = AliDielectronMC::Instance();

  Int_t mLabel1 = dieMC->GetMothersLabel(label1);    // should work for both ESD and AOD
  Int_t mLabel2 = dieMC->GetMothersLabel(label2);

  // check the same mother option
  AliDielectronSignalMC* sigMC = (AliDielectronSignalMC*)fSignalsMC->At(nSignal);
  if(sigMC->GetMothersRelation()==AliDielectronSignalMC::kSame && mLabel1!=mLabel2) return;
  if(sigMC->GetMothersRelation()==AliDielectronSignalMC::kDifferent && mLabel1==mLabel2) return;

  AliDielectronVarManager::SetFillMap(fUsedVars);
  // fill the leg variables
  Double_t valuesLeg1[AliDielectronVarManager::kNMaxValues];
  Double_t valuesLeg2[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::Fill(part1,valuesLeg1);
  AliDielectronVarManager::Fill(part2,valuesLeg2);

  // fill the pair and event variables
  Double_t valuesPair[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::Fill(dieMC->GetMCEvent(), valuesPair);
  AliDielectronVarManager::FillVarMCParticle2(part1,part2,valuesPair);

  // if pair types are filled, fill mc sources at the end
  Int_t istep=0;
  if(fPairType!=kMConly) istep=AliDielectron::kEv1PMRot+1;

  // only OS at the moment
  if(part1->Charge()*part2->Charge()<0) {
    Fill(istep+nSignal+fSignalsMC->GetEntries(), valuesPair,  valuesLeg1, valuesLeg2);
  }

  return;
}
//______________________________________________
void AliDielectronHF::Fill(Int_t pairIndex, const AliDielectronPair *particle)
{
  //
  // fill histograms for event, pair and daughter cuts and pair types
  //

  // only OS pairs in case of MC
  //////////////////////////////  if(fHasMC && pairIndex!=AliDielectron::kEv1PM) return;

  // only selected pair types in case of data
  if(!IsPairTypeSelected(pairIndex) || fEventArray) return;

  // get event and pair variables
  Double_t valuesPair[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::SetFillMap(fUsedVars);
  AliDielectronVarManager::Fill(particle,valuesPair);

  // get leg variables (TODO: do not fill for the moment since leg cuts are not opened)
  Double_t valuesLeg1[AliDielectronVarManager::kNMaxValues]={0};
  if(fVarCutType->CountBits())  AliDielectronVarManager::Fill(particle->GetFirstDaughterP(),valuesLeg1);
  Double_t valuesLeg2[AliDielectronVarManager::kNMaxValues]={0};
  if(fVarCutType->CountBits())  AliDielectronVarManager::Fill(particle->GetSecondDaughterP(),valuesLeg2);

  // fill

  // if pair types are filled, fill mc sources at the end
  Int_t istep = 0;
  if(fPairType!=kMConly) istep=AliDielectron::kEv1PMRot+1;

  // mc source steps (only OS SE pairs)
  if(fHasMC && fSignalsMC && pairIndex==AliDielectron::kEv1PM) {
    for(Int_t i=0; i<fSignalsMC->GetEntries(); i++) {
      if(AliDielectronMC::Instance()->IsMCTruth(particle, (AliDielectronSignalMC*)fSignalsMC->At(i)))
	Fill(istep+i, valuesPair,  valuesLeg1, valuesLeg2);
    }
  }

  // all pair types w/o use of mc information
  if(fPairType==kMConly) return;

  // remove comments
  //// select correct step if we are looking at signals too
  ////  if(fHasMC && fSignalsMC) pairIndex += ( fSignalsMC->GetEntries() * (fStepGenerated ? 2 : 1) );
  Fill(pairIndex, valuesPair,  valuesLeg1, valuesLeg2);

  return;
}

//______________________________________________
void AliDielectronHF::Fill(Int_t index, Double_t * const valuesPair, Double_t * const valuesLeg1, Double_t * const valuesLeg2)
{
  //
  // main fill function using index and values as input
  //

  TObjArray *histArr = static_cast<TObjArray*>(fArrPairType.At(index));
  if(!histArr) return;

  Int_t size  = GetNumberOfBins();
  // loop over all histograms
  for(Int_t ihist=0; ihist<size; ihist++) {

    Int_t sizeAdd   = 1;
    Bool_t selected = kTRUE;

    // loop over all cut variables
    Int_t nvars = fAxes.GetEntriesFast();
    for(Int_t ivar=0; ivar<nvars; ivar++) {

      // get bin limits
      TVectorD *bins = static_cast<TVectorD*>(fAxes.At(ivar));
      Int_t nbins    = bins->GetNrows()-1;

      // bin limits for current ivar bin
      Int_t ibin   = (ihist/sizeAdd)%nbins;
      Double_t lowEdge = (*bins)[ibin];
      Double_t upEdge  = (*bins)[ibin+1];
      switch(fBinType[ivar]) {
      case kStdBin:     upEdge=(*bins)[ibin+1];     break;
      case kBinToMax:   upEdge=(*bins)[nbins];      break;
      case kBinFromMin: lowEdge=(*bins)[0];         break;
      case kSymBin:     upEdge=(*bins)[nbins-ibin];
	if(ibin>=((Double_t)(nbins+1))/2) upEdge=(*bins)[nbins]; // to avoid low>up
	break;
      }

      // leg variable
      if(fVarCutType->TestBitNumber(ivar)) {
	if( (valuesLeg1[fVarCuts[ivar]] < lowEdge || valuesLeg1[fVarCuts[ivar]] >= upEdge) ||
	    (valuesLeg2[fVarCuts[ivar]] < lowEdge || valuesLeg2[fVarCuts[ivar]] >= upEdge) ) {
	  selected=kFALSE;
	  break;
	}
      }
      else { // pair and event variables
	if( (valuesPair[fVarCuts[ivar]] < lowEdge || valuesPair[fVarCuts[ivar]] >= upEdge) ) {
	  selected=kFALSE;
	  break;
	}
      }

      sizeAdd*=nbins;
    } //end of var cut loop

    // do not fill the histogram
    if(!selected) continue;

    // fill the object with Pair and event values
    TObjArray *tmp = (TObjArray*) histArr->At(ihist);
    TString title = tmp->GetName();
    AliDebug(10,title.Data());
    for(Int_t i=0; i<tmp->GetEntriesFast(); i++) {
      AliDielectronHistos::FillValues(tmp->At(i), valuesPair);
    }
    //    AliDebug(10,Form("Fill var %d %s value %f in %s \n",fVar,AliDielectronVarManager::GetValueName(fVar),valuesPair[fVar],tmp->GetName()));
  } //end of hist loop

}

//______________________________________________
void AliDielectronHF::Init()
{
  //
  // initialise event buffers
  //

  // has MC signals
  fHasMC=AliDielectronMC::Instance()->HasMC();

  Int_t steps = 0;
  if(!fSignalsMC && fHasMC){
    AliFatal("Running on MC data but AliDielectronHF::fSignalsMC is a Nullptr - Exiting AliDielectronHF::Init() now!");
    return;
  }
  if(fHasMC) steps=fSignalsMC->GetEntries();
  if(fStepGenerated) steps*=2;
  if(fEventArray) steps=1;

  // init pair type array
  fArrPairType.SetName(Form("%s_HF",GetName()));
  if( (fHasMC && fPairType==kMConly) || fEventArray) fArrPairType.Expand(steps);
  else fArrPairType.Expand(AliDielectron::kEv1PMRot+1+steps);

  Int_t size  = GetNumberOfBins();
  AliDebug(10,Form("Creating a histo array with size %d \n",size));

  Int_t sizeAdd  = 1;

  // fill object array with the array of bin cells
  TObjArray *histArr = new TObjArray(0);
  if(!histArr) return;
  histArr->SetOwner(kTRUE);
  histArr->Expand(size);

  //  printf("fRefObj %p \n",fRefObj);
  // array of histograms to each bin cell
  for(Int_t ihist=0; ihist<size; ihist++) {
    histArr->AddAt(fRefObj.Clone(""), ihist);
    //histArr->AddAt(fRefObj.Clone(Form("h%04d",ihist)), ihist);
  }

  // loop over all cut variables and do the naming according to its bin cell
  Int_t nvars = fAxes.GetEntriesFast();
  for(Int_t ivar=0; ivar<nvars; ivar++) {

    // get bin limits
    TVectorD *bins = static_cast<TVectorD*>(fAxes.At(ivar));
    Int_t nbins    = bins->GetNrows()-1;


    // loop over all bin cells an set unique titles
    for(Int_t ihist=0; ihist<size; ihist++) {

      // get the lower limit for current ivar bin
      Int_t ibin   = (ihist/sizeAdd)%nbins;
      Double_t lowEdge = (*bins)[ibin];
      Double_t upEdge  = (*bins)[ibin+1];
      switch(fBinType[ivar]) {
      case kStdBin:     upEdge=(*bins)[ibin+1];     break;
      case kBinToMax:   upEdge=(*bins)[nbins];      break;
      case kBinFromMin: lowEdge=(*bins)[0];         break;
      case kSymBin:     upEdge=(*bins)[nbins-ibin];
	if(ibin>=((Double_t)(nbins+1))/2) upEdge=(*bins)[nbins]; // to avoid low>up
	break;
      }

      TObjArray *tmp= (TObjArray*) histArr->At(ihist);
      TString title = tmp->GetName();
      if(!ivar)             title ="";
      if( ivar)             title+=":";
      if(fVarCutType->TestBitNumber(ivar)) title+="Leg";
      title+=AliDielectronVarManager::GetValueName(fVarCuts[ivar]);
      title+=Form("#%.2f#%.2f",lowEdge,upEdge);
      tmp->SetName(title.Data());
      AliDebug(10,title.Data());
    } // end: array of bin cell
    sizeAdd*=nbins;
  } //end: cut loop

  // copy array to the selected event,  pair types/ MC sources
  Int_t istep=0;

  ////////////////// only event array
  if(fEventArray) {
    // add a deep copy of the array
    fArrPairType[istep]=(TObjArray*)histArr->Clone("Event");
    ((TObjArray*)fArrPairType[istep])->SetOwner();
  }
  else {
    /////////////// pair types
    if(fPairType != kMConly) {
      for(istep=0; istep<AliDielectron::kEv1PMRot+1; istep++) {

	// pair type should be filled
	if(IsPairTypeSelected(istep)) {
	  // add a deep copy of the array
	  fArrPairType[istep]=(TObjArray*)histArr->Clone(AliDielectron::PairClassName(istep));
	  ((TObjArray*)fArrPairType[istep])->SetOwner();
	}
	else { //empty array
	  fArrPairType[istep]=new TObjArray(0);
	  ((TObjArray*)fArrPairType[istep])->SetOwner();
	  ((TObjArray*)fArrPairType[istep])->SetName(AliDielectron::PairClassName(istep));
	}
      } //end: loop over pair types
    }

    // mc sources
    if(fHasMC) {
      for(Int_t i=0; i<fSignalsMC->GetEntries(); i++) {
	TString title = Form("(Signal: %s)",fSignalsMC->At(i)->GetTitle());
	fArrPairType[istep+i]=(TObjArray*)histArr->Clone(title.Data());
	if(fStepGenerated)  {
	  title+=" MC truth";
	  fArrPairType[istep+i+fSignalsMC->GetEntries()]=(TObjArray*)histArr->Clone(title.Data());
	}
      } // end: loop over sources
    } //end: hasMC
  } //end: pair type array

  // clean up
  if(histArr) {
    delete histArr;
    histArr=0;
  }
}

//______________________________________________
Int_t AliDielectronHF::GetNumberOfBins() const
{
  //
  // return the number of bins this histogram grid has
  //
  Int_t size=1;
  for (Int_t i=0; i<fAxes.GetEntriesFast(); ++i)
    size*=((static_cast<TVectorD*>(fAxes.At(i)))->GetNrows()-1);
  return size;
}

//______________________________________________
Bool_t AliDielectronHF::IsPairTypeSelected(Int_t itype)
{
  //
  // check whether a pair type was selected
  // TODO: cross check or replace by mixinghandlers processsing

  Bool_t selected = kFALSE;

  // fill all
  if(fPairType==kAll) return kTRUE;

  switch(itype) {
  case AliDielectron::kEv1PP:
  case AliDielectron::kEv1MM:
    if(fPairType==kSeAll || fPairType==kSeMeAll || fPairType==kSeReAll )   selected = kTRUE;
    break;
  case AliDielectron::kEv1PM:
    if(fPairType!=kMeOnlyOS)  selected = kTRUE;
    break;
  case AliDielectron::kEv1PEv2P:
  case AliDielectron::kEv1MEv2M:
    if(fPairType==kMeAll || fPairType==kSeMeAll)   selected = kTRUE;
    break;
  case AliDielectron::kEv1PEv2M:
  case AliDielectron::kEv1MEv2P:
    if(fPairType==kMeAll || fPairType==kSeMeAll || fPairType==kMeOnlyOS || fPairType==kSeMeOnlyOS)   selected = kTRUE;
    break;
  case AliDielectron::kEv2PP:
  case AliDielectron::kEv2PM:
  case AliDielectron::kEv2MM:
    selected = kFALSE;
    break;
  case AliDielectron::kEv1PMRot:
    if(fPairType==kSeReAll || fPairType==kSeReOnlyOS)   selected = kTRUE;
    break;
  }

  return selected;

}
