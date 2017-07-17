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
	fContainer(0)
{
  fBinning=new Double_t*[fgkCFVars];
  fAxisTitle=new const char*[fgkCFVars];
  for (Int_t k=0; k<fgkCFVars; k++) {
    fBinning[k] = new Double_t[1000];
  }

	CreateDefaultBinning();

	if (!dummy) {

		// Constructor
		AliInfo(MSGINFO("Creating default container."));
    
		CreateContainer("fContainerStandard", "Standard container for corrections", fgkCFVars, fNbins, fBinning, fAxisTitle);  
	}
}

AliHFJetsContainer::AliHFJetsContainer(const char *name, const Int_t nvars, const char *varnames[], Int_t *nbins, Double_t **binning, const char *axistitle[]):
	TNamed(name, name),
	fCustomVarNames(0x0),
	fContainer(0)
{
  fBinning=new Double_t*[fgkCFVars];
  fAxisTitle=new const char*[fgkCFVars];
  for (Int_t k = 0; k<fgkCFVars; k++) {
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
	for (Int_t k=0; k<fgkCFVars; k++) {
    delete fBinning[k];
    fBinning[k] = 0;
  }
  
  delete [] fBinning;
}

//----------------------------------------------------------------
AliHFJetsContainer &AliHFJetsContainer::operator=(const AliHFJetsContainer &c)
{
	// assigment operato
	if (this != &c)
		((AliHFJetsContainer &) c).Copy(*this);

	return *this;
}

//----------------------------------------------------------------
void AliHFJetsContainer::Copy(TObject& c) const
{
	// copy function

	AliHFJetsContainer &target = (AliHFJetsContainer &) c;

  if (fContainer)
    target.fContainer = dynamic_cast<AliCFContainer *> (fContainer->Clone());
	if (fCustomVarNames)
		target.fCustomVarNames = dynamic_cast<TList *> (fCustomVarNames->Clone());
}

//----------------------------------------------------------------
void AliHFJetsContainer::Add(const AliHFJetsContainer *aContainerToAdd, Double_t c)
{
	//add the content of container aContainerToAdd to the current one
	// uses the Merge method of AliCFContainer
	// there the consistency between number of steps, variables and bins 
	// is checked

        fContainer->Add(aContainerToAdd->fContainer,c);
}

//----------------------------------------------------------------
Long64_t AliHFJetsContainer::Merge(TCollection *list)
{
	// Merge a list of AliCorrection objects with this (needed for
	// PROOF). 
	// Returns the number of merged objects (including this).

	if (!list)
		return 0;

	if (list->IsEmpty())
		return 1;

	TIter iter(list);
	TObject *obj;

	Int_t count = 0;
	while ((obj = iter())) {
		AliHFJetsContainer *entry = dynamic_cast<AliHFJetsContainer *> (obj);
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
	for (Int_t j = 0; j<nvars; j++) {
		fContainer->SetBinLimits(j, binning[j]);
		fContainer->SetVarTitle(j, axistitle[j]);
	}

	SetStepNames(fContainer);
  
}

//----------------------------------------------------------------
void AliHFJetsContainer::CreateCustomContainer(const Int_t nvars, const char *varnames[], Int_t *nbins, Double_t **binning, const char *axistitle[]){

	AliInfo(MSGINFO("Creating custom container: standard variables will be added at positions 0,1,2,3!"));

	CreateDefaultBinning();

	const Int_t totnvars = nvars + fgkCFVars; 
	Int_t totnbins[nvars + fgkCFVars];
	Double_t *totbinning[totnvars+1];  
	const char *totaxistitle[totnvars]; 
	for (Int_t i = 0; i<fgkCFVars; i++) {
		totnbins[i] = fNbins[i];
		totbinning[i] = fBinning[i];
		totaxistitle[i] = fAxisTitle[i];
		AliDebug(AliLog::kDebug, Form(MSGDEBUG("Standard vars: ID \"%d\" NAME \"%s\" TITLE \"%s\""), i, GetVarName((CFVars)i), totaxistitle[i])) ;
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
		AliDebug(AliLog::kDebug,Form(MSGDEBUG("Custom vars: ID \"%d\" NAME \"%s\" TITLE \"%s\""), j, GetVarName((CFVars) j), totaxistitle[j]));
	}
	CreateContainer("fContainerCustom", "Custom container for corrections", totnvars, totnbins, totbinning, totaxistitle);  
	//delete fCustomVarNames;
  
}

//----------------------------------------------------------------
//Double_t* AliHFJetsContainer::GetBinning(TString var, Int_t& nBins, const char*& axistitle)
void AliHFJetsContainer::GetBinning(TString var, Int_t &nBins, Double_t *bins, const char *&axistitle)
{
	// Assigns variable-specific bin edges 
	// (you can define array of bins "by hand" to allow for non uniform bin width) 
	//if (var.Contains("DP") || var.Contains("BP") || var.Contains("BH") ){
	//  if (var.Contains("part")) var.Form("idPart");
	//  else if (var.Contains("pt")) var.Form("ptPart");
	//}

	Float_t binmin=0., binmax=0.;    
	if (var.EqualTo("cent")) {
		axistitle="Multiplicity percentile";
		nBins = 100; binmin= 0.5; binmax= 100.5;
	}
  else if (var.EqualTo("jetPt")) {
		axistitle="p_{T,jet} (GeV/c)";
		nBins = 100; binmin= 0.; binmax= 100.;
		// Double_t *bins = {5.,10.,15., ...};
		// return bins;
	}
  else if (var.EqualTo("jetEta")) {
		axistitle="#eta_{jet}";
		nBins = 20; binmin= -1.; binmax= 1.;
	}
  else if (var.EqualTo("jetPhi")) {
		axistitle="#phi_{jet} (rad)";
		nBins = 20; binmin= -TMath::Pi(); binmax= TMath::Pi();
	}
  else {
		AliError(Form(MSGERROR("Variable %s not defined!"), var.Data()));
	}

	// Define regular binning
	Double_t binwidth = (binmax-binmin)/(1.*nBins);
	//Double_t* bins = new Double_t[nBins+1];
	//bins = new Double_t[nBins+1];
	//Double_t bins[nBins+1];
	for (Int_t j = 0; j<nBins+1; j++) {
		if (j==0) bins[j] = binmin;
		else bins[j] = bins[j-1]+binwidth;
	}
	//return bins;
}

//----------------------------------------------------------------
void AliHFJetsContainer::SetAxisRangeStep(const char *axisname, Double_t min, Double_t max, CFSteps step, Bool_t overflow)
{

	AliInfo(Form(MSGINFO("Setting range for axis: \"%s\" step: \"%s\""), axisname, GetStepName(step)));
	Int_t axis = GetVarAxis(axisname); 
  Int_t startbin= fContainer->GetAxis(axis, step)->FindBin(min);
  Int_t stopbin= fContainer->GetAxis(axis, step)->FindBin(max);
  Int_t lastbin= fContainer->GetAxis(axis, step)->GetLast();
  if (max < min || !stopbin){
    AliInfo(MSGINFO("Invalid axis range! Setting maximum to last bin!"));
    stopbin = lastbin;
  }
  if (!overflow) {
	   fContainer->GetAxis(axis, step)->SetRange(startbin,stopbin);
  }
  else {
    if (stopbin != lastbin) {
      AliError(MSGERROR("You requested overflow, but your max bin is not the last!!! I'm setting it to last bin!"));
      stopbin = lastbin;
      //return;
    }
    Double_t ofw = fContainer->GetOverFlows(axis,step);
    Printf("OVERFLOW: %f", ofw);
    if (ofw > 0.0) stopbin=stopbin+1;
    else AliInfo(MSGINFO("No overflow, projecting up to last bin only!"));
	   
    fContainer->GetAxis(axis, step)->SetRange(startbin,stopbin);
  }
  
}

//----------------------------------------------------------------
void AliHFJetsContainer::SetAxisRangeAllSteps(const char *axisname, Double_t min, Double_t max, Bool_t overflow)
{

	for (Int_t i = 0; i<fgkCFSteps; i++) {
    SetAxisRangeStep(axisname, min, max, (CFSteps) i, overflow);
	}
  
}  

//----------------------------------------------------------------
void AliHFJetsContainer::PrintVars()
{
	Int_t nvars=fContainer->GetNVar();
	for (Int_t i=0; i<nvars; i++){
		Printf("Var %d: %s -> %s",i, GetVarName((CFVars) i), fContainer->GetVarTitle(i));
	}
  
}  

//----------------------------------------------------------------
void AliHFJetsContainer::PrintSteps()
{
	Int_t nstep=fContainer->GetNStep();
	for (Int_t i=0; i<nstep; i++){
		Printf("Step %d: %s -> %s", i, GetStepName((CFSteps) i), GetStepTitle((CFSteps) i));
	}
  
}

//----------------------------------------------------------------
void AliHFJetsContainer::ResetAxisStep(const char* axisname, CFSteps step)
{
	AliInfo(Form(MSGINFO("Resetting range for axis: \"%s\" step: \"%s\""), axisname, GetStepName(step)));
	Int_t axis = GetVarAxis(axisname); 
	fContainer->GetAxis(axis, step)->SetRange(0, -1);

}

//----------------------------------------------------------------
void AliHFJetsContainer::ResetAxisAllSteps(const char* axisname)
{

	for (Int_t i=0; i<fgkCFSteps; i++){
		ResetAxisStep(axisname, (CFSteps) i);
	}  

}

//----------------------------------------------------------------
TH1D *AliHFJetsContainer::Project1D(CFSteps step, const char* varname)
{

	AliInfo(Form(MSGINFO("Projecting axis: \"%s\" step: \"%s\""), varname, GetStepName(step)));
	Int_t var = GetVarAxis(varname);  
	TH1D *h1  = (TH1D *) fContainer->Project(step, var);

	return h1;
}

//----------------------------------------------------------------
TH2D *AliHFJetsContainer::Project2D(CFSteps step, const char* varname1, const char* varname2)
{

	AliInfo(Form(MSGINFO("Projecting axis: \"%s\" and \"%s\"  step: \"%s\""), varname1, varname2, GetStepName(step)));

	Int_t var1 = GetVarAxis(varname1);  
	Int_t var2 = GetVarAxis(varname2);

	TH2D *h2 = (TH2D *) fContainer->Project(step, var1, var2);

	return h2;
}

//----------------------------------------------------------------
void AliHFJetsContainer::ScaleStep(Double_t factor, CFSteps step)
{
	Double_t fact[2] = {factor,0};

	AliInfo(Form(MSGINFO("Scaling container at step %d..."), step));
	fContainer->GetGrid(step)->Scale(fact);

}

//----------------------------------------------------------------
void AliHFJetsContainer::ScaleAllSteps(Double_t factor)
{
	for (Int_t i=0; i<fgkCFSteps; i++){
		ScaleStep(factor, (CFSteps) i);
	}
}  

//----------------------------------------------------------------
void AliHFJetsContainer::SetStepNames(AliCFContainer *container)
{
	// sets the names of the correction steps
	for (Int_t i=0; i<fgkCFSteps; i++)
		container->SetStepTitle(i, GetStepTitle((CFSteps) i));
}

//----------------------------------------------------------------
const char *AliHFJetsContainer::GetStepTitle(CFSteps step)
{
	// returns the name of the given step
	switch (step) {
    case kCFStepEventSelected:
      return "All events";
    case kCFStepMatchedAny:
      return "True jets (reco jet matches any true jet)";
    case kCFStepReco:
      return "Reconstructed jets";
  }

  return 0;
}

//----------------------------------------------------------------
const char *AliHFJetsContainer::GetStepName(CFSteps step)
{
  switch (step){
    case kCFStepEventSelected:
      return "kCFStepEventSelected";
    case kCFStepMatchedAny:
      return "kCFStepMatchedAny";
    case kCFStepReco:
      return "kCFStepReco";
  }
  
  return 0;
}

//----------------------------------------------------------------
const char *AliHFJetsContainer::GetVarName(CFVars var)
{ 

	switch (var) {
    case kCFCent:
      return "kCFCent";
    case kCFJetPt:
      return "kCFJetPt";
    case kCFJetEta:
      return "kCFJetEta";
    case kCFJetPhi:
      return "kCFJetPhi";
    default:
      //const char* value=Form("var%d",var);
      TObjString *sname = (TObjString *) fCustomVarNames->At(var-fgkCFVars);
      const char *value = sname->GetName();
      return value;
  }
  
  return 0;
}

//----------------------------------------------------------------
Int_t AliHFJetsContainer::GetVarAxis(const char* varname){
	if (!strncmp(varname,"kCFCent",50))
    return kCFCent;
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
      AliError(Form(MSGERROR("Variable \"%s\" does not exist!"),varname));
      exit(0);
    }
    Int_t value = fCustomVarNames->IndexOf((TObjString*)obj);
    //Printf("%d", value);
    value+=fgkCFVars;
    return value;
    //return 1;
  }
  
}

//----------------------------------------------------------------
void AliHFJetsContainer::CreateDefaultBinning()
{
	TString vars("cent;jetPt;jetEta;jetPhi");
	// Get binning for each variable
	TObjArray *arr;
	TObjString *objstr;
	arr = vars.Tokenize(";");
	//static const Int_t nvars = arr->GetEntriesFast();
	Int_t nvars = arr->GetEntriesFast();
	if (nvars != fgkCFVars)AliError(MSGERROR("Number of initialized variables is not fgkCFVars!"));
	TIter next(arr);
	Int_t i = 0;
	while ((objstr=(TObjString*)next())){
		GetBinning(objstr->GetString(), fNbins[i], fBinning[i], fAxisTitle[i]);
		i++;
	}
	delete arr;
}

//----------------------------------------------------------------
TH1 *AliHFJetsContainer::StepsRatio(CFSteps num, CFSteps denom, Int_t var1, Int_t var2)
{
	TH1 *hnum = NULL;
	TH1 *hdenom = NULL;
	if (var2 >= 0) {
		hnum=fContainer->Project(num, var1, var2);
		hdenom=fContainer->Project(denom, var1, var2);
	}
  else {
		hnum   = fContainer->Project(num, var1);
		hdenom = fContainer->Project(denom, var1);
	}
	hnum->Divide(hnum, hdenom, 1, 1, "B"); // "B" means binomial error
	delete hdenom;
	return hnum;
}

//----------------------------------------------------------------
TH1D *AliHFJetsContainer::GetMatchingEfficiencyPt(const char* method, Int_t flavour1, Int_t flavour2)
{
   SetAxisRangeStep(method, flavour1, flavour2, kCFStepMatchedAny);
   SetAxisRangeStep(method, flavour1, flavour2, kCFStepEventSelected);
   
   TH1 *h = StepsRatio(kCFStepMatchedAny, kCFStepEventSelected, kCFJetPt);
   h->SetTitle("Matching");
   return dynamic_cast<TH1D *> (h);
}

//----------------------------------------------------------------
TH1D *AliHFJetsContainer::GetPurityVariable(const char* method, Int_t flavour, Int_t variable)
{
  SetAxisRangeStep(method, flavour, flavour, kCFStepMatchedAny);
  TH1 *hnum = (TH1 *) fContainer->Project(kCFStepMatchedAny, variable);
  ResetAxisStep(method,kCFStepMatchedAny);
  TH1 *hdenom = (TH1 *) fContainer->Project(kCFStepMatchedAny, variable);
  hnum->Divide(hnum,hdenom,1,1,"B");
  hnum->SetTitle("Purity");
 return dynamic_cast<TH1D *> (hnum);
}

//----------------------------------------------------------------
// TH2D* AliHFJetsContainer::GetPurity2D()
// {
// 	TH1 *h = StepsRatio(kCFStepMatchedAny, kCFStepMatchedB, kCFJetPt, kCFJetEta);
// 	h->SetTitle("kCFStepMatchedAny/kCFStepMatchedB");
// 	return dynamic_cast<TH2D*>(h);
// }

//----------------------------------------------------------------
void AliHFJetsContainer::FillStep(CFSteps step, const TArrayD *point, Double_t weight){

	Int_t expectedvars=fCustomVarNames->GetEntries()+4;
	Int_t givenvars=point->GetSize();

	if (expectedvars != givenvars){
		AliError(Form(MSGERROR("Wrong number of values: expected %d, provided %d!"), expectedvars, givenvars));
	}

	fContainer->Fill(point->GetArray(), step, weight);
}

//----------------------------------------------------------------
void AliHFJetsContainer::CloneStep(CFSteps step1, CFSteps step2){
  
  // This method copies the AliCFContainer grid of step1 to that of step2
  
  AliCFGridSparse *grid = (AliCFGridSparse*) fContainer->GetGrid(step1);
  fContainer->SetGrid(step2, (AliCFGridSparse*)grid->Clone());


}
