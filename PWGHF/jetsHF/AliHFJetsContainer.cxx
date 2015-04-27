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

/* mailto: svallero@to.infn.it */

/* Class defining containers for the HF b-jets analysis */

#include "AliHFJetsContainer.h"
#include "AliCFContainer.h"
#include "AliCFGridSparse.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TList.h"
#include "TObjString.h"
#include "TCollection.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "AliLog.h"
#include "TCanvas.h"
#include "TF1.h"
#include "AliTHn.h"
#include "THn.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

ClassImp(AliHFJetsContainer)

	//const Int_t AliHFJetsContainer::fgkCFSteps = 6;
	//const Int_t AliHFJetsContainer::fgkCFVars = 4;

	//AliHFJetsContainer::AliHFJetsContainer(): 
	//  TNamed("",""),
	//  fCustomVarNames(0x0),
	//  fContainer(0)
	//{
	//  // Dummy constructor for daughter classes
	//  AliInfo("Creating dummy container."); 
	//}

AliHFJetsContainer::AliHFJetsContainer(const char* name, Bool_t dummy): 
	TNamed(name,name),
	fCustomVarNames(0x0),
	fContainer(0){
  
      fBinning=new Double_t*[fgkCFVars];
      fAxisTitle=new const char*[fgkCFVars];
	for (Int_t k=0; k<fgkCFVars; k++){
		fBinning[k] = new Double_t[1000];
	}

	CreateDefaultBinning();

	if (!dummy){

		// Constructor
		AliInfo(MAG"Creating default container." Bee); 


		CreateContainer("fContainerStandard", "Standard container for corrections", fgkCFVars, fNbins, fBinning, fAxisTitle);  
	}

}

AliHFJetsContainer::AliHFJetsContainer(const char* name, const Int_t nvars, const char* varnames[], Int_t *nbins, Double_t **binning, const char*  axistitle[]): 
	TNamed(name,name),
	fCustomVarNames(0x0),
	fContainer(0)
{
      fBinning=new Double_t*[fgkCFVars];
      fAxisTitle=new const char*[fgkCFVars];
	for (Int_t k=0; k<fgkCFVars; k++){
		fBinning[k] = new Double_t[1000];
	}
	CreateCustomContainer(nvars, varnames, nbins, binning, axistitle);

}

//----------------------------------------------------------------
AliHFJetsContainer::AliHFJetsContainer(const AliHFJetsContainer &c) :
	TNamed(),
	fContainer(0)
{

	// AliHFJetsContainer copy constructor

	((AliHFJetsContainer &) c).Copy(*this);
}

//----------------------------------------------------------------
AliHFJetsContainer::~AliHFJetsContainer()
{
	// Destructor

	if (fContainer)
	{
		delete fContainer;
		fContainer = 0;
	}

	// Delete arrays 
	for (Int_t k=0; k<fgkCFVars; k++){
		delete fBinning[k];
                fBinning[k]=0;
	}
	delete [] fBinning;
}


//----------------------------------------------------------------
AliHFJetsContainer &AliHFJetsContainer::operator=(const AliHFJetsContainer &c)
{
	// assigment operator

	if (this != &c)
		((AliHFJetsContainer &) c).Copy(*this);

	return *this;
}

//----------------------------------------------------------------
void AliHFJetsContainer::Copy(TObject& c) const
{
	// copy function

	AliHFJetsContainer& target = (AliHFJetsContainer &) c;

	if (fContainer)
		target.fContainer = dynamic_cast<AliCFContainer*> (fContainer->Clone());
	if (fCustomVarNames)
		target.fCustomVarNames = dynamic_cast<TList*> (fCustomVarNames->Clone());
	if (fCustomVarNames)
		target.fCustomVarNames = dynamic_cast<TList*> (fCustomVarNames->Clone());

}

//----------------------------------------------------------------
void AliHFJetsContainer::Add(const AliHFJetsContainer* aContainerToAdd, Double_t c)
{

	//add the content of container aContainerToAdd to the current one
	// uses the Merge method of AliCFContainer
	// there the consistency between number of steps, variables and bins 
	// is checked

        fContainer->Add(aContainerToAdd->fContainer,c);
}


//----------------------------------------------------------------
Long64_t AliHFJetsContainer::Merge(TCollection* list)
{
	// Merge a list of AliCorrection objects with this (needed for
	// PROOF). 
	// Returns the number of merged objects (including this).

	if (!list)
		return 0;

	if (list->IsEmpty())
		return 1;

	TIter iter(list);
	TObject* obj;

	Int_t count = 0;
	while ((obj = iter())) {
		AliHFJetsContainer* entry = dynamic_cast<AliHFJetsContainer*> (obj);
		if (entry == 0)
			continue;
		this->Add(entry);
		count++;
	}

	return count+1;
}

//----------------------------------------------------------------
void AliHFJetsContainer::CreateContainer(TString name, TString title, Int_t nvars, Int_t *nbins, Double_t **binning, const char *axistitle[])
{

	fContainer = new AliCFContainer(name, title, fgkCFSteps,nvars,nbins);
	for (Int_t j=0; j<nvars; j++){
		fContainer->SetBinLimits(j, binning[j]);
		fContainer->SetVarTitle(j, axistitle[j]);
	}

	SetStepNames(fContainer);

}


//----------------------------------------------------------------
void AliHFJetsContainer::CreateCustomContainer(const Int_t nvars, const char* varnames[], Int_t *nbins, Double_t **binning, const char*  axistitle[]){

	AliInfo(MAG"Creating custom container: standard variables will be added at positions 0,1,2,3!" Bee);

	CreateDefaultBinning();

	const Int_t totnvars = nvars + fgkCFVars; 
	Int_t totnbins[nvars + fgkCFVars];
	Double_t *totbinning[totnvars+1];  
	const char *totaxistitle[totnvars]; 
	for (Int_t i=0; i<fgkCFVars; i++){
		totnbins[i]=fNbins[i];
		totbinning[i]=fBinning[i];
		totaxistitle[i]=fAxisTitle[i];
		AliDebug(AliLog::kDebug,Form(mage"Standard vars: ID \"%d\" NAME \"%s\" TITLE \"%s\"" Bee, i, GetVarName((CFVars)i), totaxistitle[i])) ;
	}

	// Write custom variables names in a global list,
	// to ease SetAxis and Project functions
	fCustomVarNames = new TList();
	fCustomVarNames->SetOwner(kTRUE);
	for (Int_t j=fgkCFVars; j<totnvars; j++){
		totnbins[j]=nbins[j-fgkCFVars];
		totbinning[j]=binning[j-fgkCFVars];
		totaxistitle[j]=axistitle[j-fgkCFVars];
		TObjString sobj(varnames[j-fgkCFVars]);
		fCustomVarNames->Add((TObjString*)sobj.Clone());
		AliDebug(AliLog::kDebug,Form(cy"Custom vars: ID \"%d\" NAME \"%s\" TITLE \"%s\"" Bee, j, GetVarName((CFVars) j), totaxistitle[j]));
	}
	CreateContainer("fContainerCustom", "Custom container for corrections", totnvars, totnbins, totbinning, totaxistitle);  
	//delete fCustomVarNames;
}
//----------------------------------------------------------------
//Double_t* AliHFJetsContainer::GetBinning(TString var, Int_t& nBins, const char*& axistitle)
void AliHFJetsContainer::GetBinning(TString var, Int_t& nBins,Double_t* bins, const char*& axistitle)
{
	// Assigns variable-specific bin edges 
	// (you can define array of bins "by hand" to allow for non uniform bin width) 
	//if (var.Contains("DP") || var.Contains("BP") || var.Contains("BH") ){
	//  if (var.Contains("part")) var.Form("idPart");
	//  else if (var.Contains("pt")) var.Form("ptPart");
	//}

	Float_t binmin=0., binmax=0.;    
	if (var.EqualTo("mult")){
		axistitle="Multiplicity";
		nBins = 1000; binmin= 0.5; binmax= 1000.5;
	} else if (var.EqualTo("jetPt")){
		axistitle="p_{T,jet} (GeV/c)";
		nBins = 100; binmin= 0.; binmax= 100.;
		// Double_t *bins = {5.,10.,15., ...};
		// return bins;
	} else if (var.EqualTo("jetEta")){
		axistitle="#eta_{jet}";
		nBins = 20; binmin= -1.; binmax= 1.;
	} else if (var.EqualTo("jetPhi")){
		axistitle="#phi_{jet} (rad)";
		nBins = 20; binmin= -TMath::Pi(); binmax= TMath::Pi();
	} else {
		AliError(Form(RED"Variable %s not defined!" Bee, var.Data()));
	}


	// Define regular binning
	Double_t binwidth = (binmax-binmin)/(1.*nBins);
	//Double_t* bins = new Double_t[nBins+1];
	//bins = new Double_t[nBins+1];
	//Double_t bins[nBins+1];
	for (Int_t j=0; j<nBins+1; j++){
		if (j==0) bins[j]= binmin;
		else bins[j] = bins[j-1]+binwidth;
		//Printf(RED"*** Bin %d value %f" Bee,j, bins[j]);
	}
	//return bins;
}

void AliHFJetsContainer::SetAxisRangeStep(const char* axisname, Double_t min, Double_t max, CFSteps step, Bool_t overflow)
{

	AliInfo(Form(MAG"Setting range for axis: \"%s\" step: \"%s\"" Bee, axisname, GetStepName(step)));
	Int_t axis = GetVarAxis(axisname); 
	//AliInfo(Form("Resetting axis %d", axis));
        Int_t startbin= fContainer->GetAxis(axis, step)->FindBin(min);
        Int_t stopbin= fContainer->GetAxis(axis, step)->FindBin(max);
        Int_t lastbin= fContainer->GetAxis(axis, step)->GetLast();
	//AliInfo(Form(RED"startbin: %d stopbin: %d lastbin: %d" Bee, startbin, stopbin, lastbin));
        //if (stopbin < startbin || !stopbin){
        if (max < min || !stopbin){
		AliInfo("Invalid axis range! Setting maximum to last bin!");
 		stopbin= lastbin;
        } 
        if (!overflow){ 
	   //fContainer->GetAxis(axis, step)->SetRangeUser(min,max);
	   fContainer->GetAxis(axis, step)->SetRange(startbin,stopbin);
        } else {
           if (stopbin != lastbin) {
             AliError(RED"You requested overflow, but your max bin is not the last!!! I'm setting it to last bin!" Bee);
             stopbin = lastbin;
             //return;
             }
           Double_t ofw=fContainer->GetOverFlows(axis,step);
           Printf("OVERFLOW: %f", ofw);
           if (ofw > 0.0)stopbin=stopbin+1;
           else AliInfo("No overflow, projecting up to last bin only!");
	   fContainer->GetAxis(axis, step)->SetRange(startbin,stopbin);
        }
}

void AliHFJetsContainer::SetAxisRangeAllSteps(const char* axisname, Double_t min, Double_t max, Bool_t overflow)
{

	for (Int_t i=0; i<fgkCFSteps; i++){
		SetAxisRangeStep(axisname, min, max, (CFSteps) i, overflow);
	}  

}  

void AliHFJetsContainer::PrintVars()
{
	Int_t nvars=fContainer->GetNVar();
	for (Int_t i=0; i<nvars; i++){
		Printf(cy"Var %d: %s -> %s" Bee,i, GetVarName((CFVars) i), fContainer->GetVarTitle(i));
	}
}  

void AliHFJetsContainer::PrintSteps()
{
	Int_t nstep=fContainer->GetNStep();
	for (Int_t i=0; i<nstep; i++){
		Printf(mage"Step %d: %s -> %s" Bee,i, GetStepName((CFSteps) i), GetStepTitle((CFSteps) i));
	}
}

void AliHFJetsContainer::ResetAxisStep(const char* axisname, CFSteps step)
{

	AliInfo(Form(MAG"Resetting range for axis: \"%s\" step: \"%s\"" Bee, axisname, GetStepName(step)));
	Int_t axis = GetVarAxis(axisname); 
	fContainer->GetAxis(axis, step)->SetRange(0, -1);

}

void AliHFJetsContainer::ResetAxisAllSteps(const char* axisname)
{

	for (Int_t i=0; i<fgkCFSteps; i++){
		ResetAxisStep(axisname,(CFSteps) i);
	}  

}  

TH1D *AliHFJetsContainer::Project1D(CFSteps step, const char* varname)
{

	AliInfo(Form(MAG"Projecting axis: \"%s\" step: \"%s\"" Bee, varname, GetStepName(step)));
	Int_t var = GetVarAxis(varname);  
	TH1D* h1=(TH1D*)fContainer->Project(step, var);

	return h1;
}

TH2D *AliHFJetsContainer::Project2D(CFSteps step, const char* varname1, const char* varname2)
{

	AliInfo(Form(MAG"Projecting axis: \"%s\" and \"%s\"  step: \"%s\"" Bee, varname1, varname2, GetStepName(step)));

	Int_t var1 = GetVarAxis(varname1);  
	Int_t var2 = GetVarAxis(varname2);

	TH2D* h2=(TH2D*)fContainer->Project(step, var1, var2);

	return h2;
}

void AliHFJetsContainer::ScaleStep(Double_t factor, CFSteps step)
{
	Double_t fact[2] = {factor,0};

	AliInfo(Form(MAG"Scaling container at step %d..." Bee, step));
	fContainer->GetGrid(step)->Scale(fact);

}


void AliHFJetsContainer::ScaleAllSteps(Double_t factor)
{
	for (Int_t i=0; i<fgkCFSteps; i++){
		ScaleStep(factor,(CFSteps) i);
	}
}  

void AliHFJetsContainer::SetStepNames(AliCFContainer* container)
{
	// sets the names of the correction steps

	for (Int_t i=0; i<fgkCFSteps; i++)
		container->SetStepTitle(i, GetStepTitle((CFSteps) i));
}

const char* AliHFJetsContainer::GetStepTitle(CFSteps step)
{
	// returns the name of the given step
	switch (step){
		case kCFStepAll:
			return "All events";
		case kCFStepTriggered:
			return "Triggered";
		case kCFStepVertex:
			return "Primary vertex";
		case kCFStepMatchedB:
			return "True B-jets (reco B-jet matches true B-jet)";
		case kCFStepMatchedAny:
			return "True B-jets (reco B-jet matches any true jet)";
		case kCFStepRecoB:
			return "Reconstructed B-jets";
		case kCFStepMatchedC:
			return "True jets (reco jet matches true C-jet)";
		case kCFStepMatchedLight:
			return "True jets (reco jet matches true light-jet)";
		case kCFStepMatchedGluon:
			return "True jets (reco jet matches true gluon-jet)";
	}

	return 0;
}

const char* AliHFJetsContainer::GetStepName(CFSteps step)
{
	switch (step){
		case kCFStepAll:
			return "kCFStepAll";
		case kCFStepTriggered:
			return "kCFStepTriggered";
		case kCFStepVertex:
			return "kCFStepVertex";
		case kCFStepMatchedB:
			return "kCFStepMatchedB";
		case kCFStepMatchedAny:
			return "kCFStepMatchedAny";
		case kCFStepRecoB:
			return "kCFStepRecoB";
		case kCFStepMatchedC:
			return "kCFStepMatchedC";
		case kCFStepMatchedLight:
			return "kCFStepMatchedLight";
		case kCFStepMatchedGluon:
			return "kCFStepMatchedGluon";
	}
	return 0;
}

const char* AliHFJetsContainer::GetVarName(CFVars var)
{ 

	switch (var){
		case kCFMult:
			return "kCFMult";
		case kCFJetPt:
			return "kCFJetPt";
		case kCFJetEta:
			return "kCFJetEta";
		case kCFJetPhi:
			return "kCFJetPhi";
		default:
			//const char* value=Form("var%d",var);
			TObjString *sname = (TObjString*)fCustomVarNames->At(var-fgkCFVars);
			const char* value = sname->GetName();
			return value;
	}
	return 0;
}  

Int_t AliHFJetsContainer::GetVarAxis(const char* varname){
	if (!strncmp(varname,"kCFMult",50))
		return kCFMult;
	else if (!strncmp(varname,"kCFJetPt",50))
		return kCFJetPt;
	else if (!strncmp(varname,"kCFJetEta",50))
		return kCFJetEta;
	else if (!strncmp(varname,"kCFJetPhi",50))
		return kCFJetPhi;
	else {
		//const char* value=Form("var%d",var);
		TObjString *obj = (TObjString*)fCustomVarNames->FindObject(varname);
		if (!obj) {
			AliError(Form(RED"Variable \"%s\" does not exist!" Bee,varname));
			exit(0);
		}
		Int_t value = fCustomVarNames->IndexOf((TObjString*)obj);
		//Printf("%d", value);
		value+=fgkCFVars;
		return value;
		//return 1;
	}

}

void AliHFJetsContainer::CreateDefaultBinning()
{
	TString vars("mult;jetPt;jetEta;jetPhi");
	// Get binning for each variable
	TObjArray *arr;
	TObjString *objstr;
	arr = vars.Tokenize(";");
	//static const Int_t nvars = arr->GetEntriesFast();
	Int_t nvars = arr->GetEntriesFast();
	if (nvars != fgkCFVars)AliError(RED"Number of initialized variables is not fgkCFVars!" Bee);
	//Int_t nbins[4];       // number of bins for each variable
	//const char* axistitle[4]; // axis title for each variable
	//Double_t *binning[4]; // array of bins for each variable
	TIter next(arr);
	Int_t i = 0;
	while ((objstr=(TObjString*)next())){
		//fBinning[i] = new Double_t[1000];
		//fBinning[i]=GetBinning(objstr->GetString(), fNbins[i], fAxisTitle[i]);
		GetBinning(objstr->GetString(), fNbins[i], fBinning[i], fAxisTitle[i]);
		i++;
	}
	delete arr;
}


TH1* AliHFJetsContainer::StepsRatio(CFSteps num, CFSteps denom, Int_t var1, Int_t var2)
{
	TH1 *hnum=0;
	TH1 *hdenom=0;
	if (var2 >= 0){
		hnum=fContainer->Project(num, var1, var2);
		hdenom=fContainer->Project(denom, var1, var2);
	} else {
		hnum=fContainer->Project(num, var1);
		hdenom=fContainer->Project(denom, var1);
	}
	hnum->Divide(hnum, hdenom, 1, 1, "B"); // "B" means binomial error
	//hnum->Divide(hnum, hdenom);
	delete hdenom;
	return hnum;
}

TH1D* AliHFJetsContainer::GetBEfficiencyPt(const char* method)
{
   TH1D *h = (TH1D*)GetEfficiencyPt(method, 4);
   return h;
}

TH1D* AliHFJetsContainer::GetCEfficiencyPt(const char* method)
{
   TH1D *h = (TH1D*)GetEfficiencyPt(method, 3);
   return h;
}

TH1D* AliHFJetsContainer::GetLightEfficiencyPt(const char* method)
{
   TH1D *h = (TH1D*)GetEfficiencyPt(method, 2);
   return h;
}

TH1D* AliHFJetsContainer::GetGluonEfficiencyPt(const char* method)
{
   TH1D *h = (TH1D*)GetEfficiencyPt(method, 1);
   return h;
}

TH1D* AliHFJetsContainer::GetMatchingEfficiencyPt(const char* method)
{
   TH1D *h = (TH1D*)GetEfficiencyPt(method, 0);
   return h;
}

TH1D* AliHFJetsContainer::GetEfficiencyPt(const char* method, Int_t flavour)
{
        // Determine step "matched" according to chosen flavour
        CFSteps step_matched;
        Float_t fmin= flavour*1.;
        Float_t fmax = flavour*1.;
        switch (flavour){
          case 0:
	    step_matched = kCFStepMatchedAny;
            fmin=0; fmax=4;
	    break;
	  case 1:
	    step_matched = kCFStepMatchedGluon;
            break;
	  case 2:
            step_matched = kCFStepMatchedLight;		
            break;
	  case 3:
            step_matched = kCFStepMatchedC;		
            break;
	  case 4:
            step_matched = kCFStepMatchedB;		
            break;
	default:
            AliError(RED"Wrong flavour!" Bee);
            break;
        } 
         
        // restric flavour of MC jet (denominator)
        ResetAxisStep(method,kCFStepVertex);
        SetAxisRangeStep(method,fmin,fmax, kCFStepVertex);
        // restric flavour of RECO jet (numerator)
        ResetAxisStep(method,step_matched);
        SetAxisRangeStep(method,fmin,fmax,step_matched);
	TH1 *h = StepsRatio(step_matched, kCFStepVertex, kCFJetPt);
	h->SetTitle(Form("%s/kCFStepVertex",GetStepName(step_matched)));
	return dynamic_cast<TH1D*>(h);
}

TH1D* AliHFJetsContainer::GetEfficiencyEta()
{
	TH1 *h = StepsRatio(kCFStepMatchedB, kCFStepVertex, kCFJetEta);
	h->SetTitle("kCFStepMatchedB/kCFStepVertex");
	return dynamic_cast<TH1D*>(h);
}

TH2D* AliHFJetsContainer::GetEfficiency2D()
{
	TH1 *h = StepsRatio(kCFStepMatchedB, kCFStepVertex, kCFJetPt, kCFJetEta);
	h->SetTitle("kCFStepMatchedB/kCFStepVertex");
	return dynamic_cast<TH2D*>(h);
}

TH1D* AliHFJetsContainer::GetPurityPt()
{
	TH1 *h = StepsRatio(kCFStepMatchedAny, kCFStepMatchedB, kCFJetPt);
	h->SetTitle("kCFStepMatchedAny/kCFStepMatchedB");
	return dynamic_cast<TH1D*>(h);
}

TH1D* AliHFJetsContainer::GetPurityEta()
{
	TH1 *h = StepsRatio(kCFStepMatchedAny, kCFStepMatchedB, kCFJetEta);
	h->SetTitle("kCFStepMatchedAny/kCFStepMatchedB");
	return dynamic_cast<TH1D*>(h);
}

TH2D* AliHFJetsContainer::GetPurity2D()
{
	TH1 *h = StepsRatio(kCFStepMatchedAny, kCFStepMatchedB, kCFJetPt, kCFJetEta);
	h->SetTitle("kCFStepMatchedAny/kCFStepMatchedB");
	return dynamic_cast<TH2D*>(h);
}

void AliHFJetsContainer::FillStep(CFSteps step,const TArrayD *point){

	Int_t expectedvars=fCustomVarNames->GetEntries()+4;
	Int_t givenvars=point->GetSize();

	if (expectedvars!=givenvars){
		AliError(Form(RED"Wrong number of values: expected %d, provided %d!" Bee, expectedvars, givenvars));
	}

	fContainer->Fill(point->GetArray(), step, 1.);
}

void AliHFJetsContainer::CloneStep(CFSteps step1, CFSteps step2){
  
  // This method copies the AliCFContainer grid of step1 to that of step2
  
  AliCFGridSparse *grid = (AliCFGridSparse*) fContainer->GetGrid(step1);
  fContainer->SetGrid(step2, (AliCFGridSparse*)grid->Clone());


}
