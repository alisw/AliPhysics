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
  fArrPairType(0x0),
  fPairType(kOSonly),
  fSignalsMC(0x0),
  fAxes(kMaxCuts),
  fHasMC(kFALSE),
  fStepGenerated(kFALSE),
  fRefObj(1)
{
  //
  // Default Constructor
  //
  for (Int_t i=0; i<kMaxCuts; ++i){
    fVarCuts[i]=0;
    fVarCutType[i]=0;
    fBinType[i]=kStdBin;
  }
  fAxes.SetOwner(kTRUE);
  fRefObj.SetOwner(kTRUE);
}

//______________________________________________
AliDielectronHF::AliDielectronHF(const char* name, const char* title) :
  TNamed(name, title),
  fArrPairType(0x0),
  fPairType(kOSonly),
  fSignalsMC(0x0),
  fAxes(kMaxCuts),
  fHasMC(kFALSE),
  fStepGenerated(kFALSE),
  fRefObj(1)
{
  //
  // Named Constructor
  //
  for (Int_t i=0; i<kMaxCuts; ++i){
    fVarCuts[i]=0;
    fVarCutType[i]=0;
    fBinType[i]=kStdBin;
  }
  fAxes.SetOwner(kTRUE);
  fRefObj.SetOwner(kTRUE);
}

//______________________________________________
AliDielectronHF::~AliDielectronHF()
{
  //
  // Default Destructor
  //
  fAxes.Delete();
}

//_____________________________________________________________________________
void AliDielectronHF::UserProfile(const char* histClass, UInt_t valTypeP,
				      const TVectorD * const binsX,
				      UInt_t valTypeX, TString option)
{
  //
  // Histogram creation 1D case with arbitraty binning X
  // the TVectorD is assumed to be surplus after the creation and will be deleted!!!
  //

  TH1 *hist=0x0;
  if(valTypeP==999)
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
    hist=new TProfile("","",binsX->GetNrows()-1,binsX->GetMatrixArray());
    ((TProfile*)hist)->BuildOptions(pmin,pmax,opt.Data());
    //      printf(" name %s PROFILE options: pmin %.1f pmax %.1f err %s \n",name,((TProfile*)hist)->GetYmin(),((TProfile*)hist)->GetYmax(),((TProfile*)hist)->GetErrorOption() );
  }

  // store variales in axes
  UInt_t valType[4] = {0};
  valType[0]=valTypeX;     valType[1]=valTypeP;
  AliDielectronHistos::StoreVariables(hist, valType);

  // adapt the name and title of the histogram in case they are empty
  AliDielectronHistos::AdaptNameTitle(hist, histClass);
  hist->SetName(Form("HF_%s",hist->GetName()));

  fRefObj.AddLast(hist);
  delete binsX;
}

//_____________________________________________________________________________
void AliDielectronHF::UserProfile(const char* histClass, UInt_t valTypeP,
				      const TVectorD * const binsX, const TVectorD * const binsY,
				      UInt_t valTypeX, UInt_t valTypeY, TString option)
{
  //
  // Histogram creation 2D case with arbitraty binning X and Y
  // the TVectorD is assumed to be surplus after the creation and will be deleted!!!
  //

  TH1 *hist=0x0;
  if(valTypeP==999) {
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
  valType[0]=valTypeX;     valType[1]=valTypeY; valType[3]=valTypeP;
  AliDielectronHistos::StoreVariables(hist, valType);

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
				      UInt_t valTypeX, UInt_t valTypeY, UInt_t valTypeZ, TString option)
{
  //
  // Histogram creation 3D case with arbitraty binning X, Y, Z
  // the TVectorD is assumed to be surplus after the creation and will be deleted!!!
  //
  TH1 *hist=0x0;
  if(valTypeP==999) {
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

  // adapt the name and title of the histogram in case they are empty
  AliDielectronHistos::AdaptNameTitle(hist, histClass);
  hist->SetName(Form("HF_%s",hist->GetName()));

  fRefObj.AddLast(hist);
  delete binsX;
  delete binsY;
  delete binsZ;
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
  fVarCutType[size]=leg;
  fAxes.Add(binLimits->Clone());
  fBinType[size]=btype;
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
  fVarCutType[size]=leg;
  fAxes.Add(binLimits);
  fBinType[size]=btype;
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
  fVarCutType[size]=leg;
  fAxes.Add(binLimits);
  fBinType[size]=btype;
}

//______________________________________________
void AliDielectronHF::Fill(Int_t label1, Int_t label2, Int_t nSignal) 
{
  //
  // fill the pure MC part of the container starting from a pair of 2 particles (part1 and part2 are legs)
  //
  // fill only if we have asked for these steps
  if(!fStepGenerated) return;

  AliVParticle* part1 = AliDielectronMC::Instance()->GetMCTrackFromMCEvent(label1);
  AliVParticle* part2 = AliDielectronMC::Instance()->GetMCTrackFromMCEvent(label2);

  AliDielectronMC* dieMC = AliDielectronMC::Instance();
  
  Int_t mLabel1 = dieMC->GetMothersLabel(label1);    // should work for both ESD and AOD
  Int_t mLabel2 = dieMC->GetMothersLabel(label2);

  // check the same mother option
  AliDielectronSignalMC* sigMC = (AliDielectronSignalMC*)fSignalsMC->At(nSignal);
  if(sigMC->GetMothersRelation()==AliDielectronSignalMC::kSame && mLabel1!=mLabel2) return;
  if(sigMC->GetMothersRelation()==AliDielectronSignalMC::kDifferent && mLabel1==mLabel2) return;
    
  // fill the leg variables
  Double_t valuesLeg1[AliDielectronVarManager::kNMaxValues];
  Double_t valuesLeg2[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::Fill(part1,valuesLeg1);
  AliDielectronVarManager::Fill(part2,valuesLeg2);
    
  // fill the pair and event variables
  Double_t valuesPair[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::Fill(dieMC->GetMCEvent(), valuesPair);
  AliDielectronVarManager::FillVarMCParticle2(part1,part2,valuesPair);

  if(part1->Charge()*part2->Charge()<0) {
    //valuesPair[AliDielectronVarManager::kPairType]=1;
//     for(Int_t i=0; i<fAxes.GetEntriesFast(); i++) {
//       printf("[D]: %s %f %f %f \n",
// 	     AliDielectronVarManager::GetValueName(fVarCuts[i]),
// 	     valuesLeg1[fVarCuts[i]], valuesLeg2[fVarCuts[i]], valuesPair[fVarCuts[i]]); 
//     }
    Fill(nSignal+fSignalsMC->GetEntries(), valuesPair,  valuesLeg1, valuesLeg2);
  }
  // only OS at the moment
  // else if(part1->Charge()>0)
  //   valuesPair[AliDielectronVarManager::kPairType]=0;
  // else
  //   valuesPair[AliDielectronVarManager::kPairType]=2; // if one of the two particles is neutral, the pair will go here

  return;
}
//______________________________________________
void AliDielectronHF::Fill(Int_t pairIndex, const AliDielectronPair *particle)
{
  //
  // fill histograms for event, pair and daughter cuts and pair types
  //
  
  // only OS pairs in case of MC
  if(fHasMC && pairIndex!=AliDielectron::kEv1PM) return;

  // only selected pair types in case of data
  if(!IsPairTypeSelected(pairIndex)) return;

  // get event and pair variables
  Double_t valuesPair[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::Fill(particle,valuesPair);

  // get leg variables
  Double_t valuesLeg1[AliDielectronVarManager::kNMaxValues]={0};
  AliDielectronVarManager::Fill(particle->GetFirstDaughter(),valuesLeg1);
  Double_t valuesLeg2[AliDielectronVarManager::kNMaxValues]={0};
  AliDielectronVarManager::Fill(particle->GetSecondDaughter(),valuesLeg2);

  // fill
  if(!fHasMC) { Fill(pairIndex, valuesPair,  valuesLeg1, valuesLeg2); }
  if(fHasMC && fSignalsMC) {
    for(Int_t i=0; i<fSignalsMC->GetEntries(); i++) {
      if(AliDielectronMC::Instance()->IsMCTruth(particle, (AliDielectronSignalMC*)fSignalsMC->At(i))) 
	Fill(i, valuesPair,  valuesLeg1, valuesLeg2);
    }
  }
  
  return;
  
}

//______________________________________________
void AliDielectronHF::Fill(Int_t Index, Double_t * const valuesPair, Double_t * const valuesLeg1, Double_t * const valuesLeg2)
{
  //
  // main fill function using index and values as input
  //

  TObjArray *histArr = static_cast<TObjArray*>(fArrPairType.At(Index));
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
      if(fVarCutType[ivar]) {
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

    // fill the object with Pair and event values (TODO: this needs to be changed)
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
  if(fHasMC) steps=fSignalsMC->GetEntries();
  if(fStepGenerated) steps*=2; 

  // init pair type array
  fArrPairType.SetName(Form("%s_HF",GetName()));
  if(fHasMC) fArrPairType.Expand(steps);
  else fArrPairType.Expand(AliDielectron::kEv1PMRot+1);

  Int_t size  = GetNumberOfBins();
  AliDebug(10,Form("Creating a histo array with size %d \n",size));

  Int_t sizeAdd  = 1; 

  // fill object array with the histograms
  TObjArray *histArr = new TObjArray();
  histArr->Expand(size);

  //  printf("fRefObj %p \n",fRefObj);
  for(Int_t ihist=0; ihist<size; ihist++) {
    histArr->AddAt(fRefObj.Clone(""), ihist);
    //histArr->AddAt(fRefObj.Clone(Form("h%04d",ihist)), ihist);
  }

  // loop over all cut variables
  Int_t nvars = fAxes.GetEntriesFast();
  for(Int_t ivar=0; ivar<nvars; ivar++) {
    
    // get bin limits
    TVectorD *bins = static_cast<TVectorD*>(fAxes.At(ivar));
    Int_t nbins    = bins->GetNrows()-1;


    // loop over all histograms an set unique titles
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
      if(fVarCutType[ivar]) title+="Leg";
      title+=AliDielectronVarManager::GetValueName(fVarCuts[ivar]);
      title+=Form("#%.2f#%.2f",lowEdge,upEdge);
      tmp->SetName(title.Data());
      AliDebug(10,title.Data());
    }
    sizeAdd*=nbins;
  }

  // copy array to the selected pair types/ MC sources
  if(fHasMC) {
    for(Int_t i=0; i<fSignalsMC->GetEntries(); i++) {
	TString title = Form("(Signal: %s)",fSignalsMC->At(i)->GetTitle());
	fArrPairType[i]=(TObjArray*)histArr->Clone(title.Data());
	if(fStepGenerated)  {
	  title+=" MC truth";
	  fArrPairType[i+fSignalsMC->GetEntries()]=(TObjArray*)histArr->Clone(title.Data());
	}
      }
   }
   else {
    for(Int_t i=0; i<AliDielectron::kEv1PMRot+1; i++) {
      if(IsPairTypeSelected(i)) fArrPairType[i]=(TObjArray*)histArr->Clone(Form("%s",AliDielectron::PairClassName(i)));
      else fArrPairType[i]=0x0;
    }
  }
  
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
  // return the number of bins this mixing handler has
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
  //
  
  Bool_t selected = kFALSE;

  // fill all
  if(fPairType==kAll) selected = kTRUE;

  switch(itype) {
  case AliDielectron::kEv1PP: 
  case AliDielectron::kEv1MM:
    if(fPairType==kOSandLS)   selected = kTRUE;
    break;
  case AliDielectron::kEv1PM: selected = kTRUE;
    break;
  case AliDielectron::kEv1PEv2P:
  case AliDielectron::kEv1MEv2P:
  case AliDielectron::kEv1PEv2M:
  case AliDielectron::kEv1MEv2M:
    if(fPairType==kOSandMIX)  selected = kTRUE;
    break;
  case AliDielectron::kEv2PP:
  case AliDielectron::kEv2PM: 
  case AliDielectron::kEv2MM:
    selected = kFALSE;
    break;
  case AliDielectron::kEv1PMRot:
    if(fPairType==kOSandROT)  selected = kTRUE;
    break;
  }
  
  return selected;

}


