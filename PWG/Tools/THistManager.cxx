/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
#include <cfloat>
#include <cstring>
#include <iostream>   // for unit tests
#include <sstream>
#include <string>
#include <exception>
#include <vector>
#include <TArrayD.h>
#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <THashList.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TProfile.h>
#include <TString.h>

#include "TBinning.h"
#include "THistManager.h"

/// \cond CLASSIMP
ClassImp(THistManager)
/// \endcond

THistManager::THistManager():
		TNamed(),
		fHistos(NULL),
		fIsOwner(true)
{
}

THistManager::THistManager(const char *name):
		TNamed(name, Form("Histogram container %s", name)),
		fHistos(NULL),
		fIsOwner(true)
{
	fHistos = new THashList();
	fHistos->SetName(Form("histos%s", name));
	fHistos->SetOwner();
}

THistManager::~THistManager(){
	if(fHistos && fIsOwner) delete fHistos;
}

THashList* THistManager::CreateHistoGroup(const char *groupname) {
  // At first step check whether the group already exists.
  THashList *foundgroup = FindGroup(groupname);
  if(foundgroup){
    // already existing, don't create it again
    return foundgroup;
  }
  TString parentname = basename(groupname);
	THashList *parentgroup = FindGroup(parentname);
	if(!parentgroup) parentgroup = CreateHistoGroup(parentname);
	THashList *childgroup = new THashList();
	childgroup->SetName(histname(groupname));
	childgroup->SetOwner();
	parentgroup->Add(childgroup);
	return childgroup;
}

TH1* THistManager::CreateTH1(const char *name, const char *title, int nbins, double xmin, double xmax, Option_t *opt){
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent) parent = CreateHistoGroup(dirname);
	if(parent->FindObject(hname)){
		Fatal("THistManager::CreateTH1", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TH1* h = new TH1D(hname, title, nbins, xmin, xmax);
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    h->Sumw2();
	parent->Add(h);
	return h;
}

TH1* THistManager::CreateTH1(const char *name, const char *title, int nbins, const double *xbins, Option_t *opt){
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent) parent = CreateHistoGroup(dirname);
	if(parent->FindObject(hname)){
		Fatal("THistManager::CreateTH1", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TH1* h = new TH1D(hname, title, nbins, xbins);
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    h->Sumw2();
	parent->Add(h);
	return h;
}

TH1* THistManager::CreateTH1(const char *name, const char *title, const TArrayD &xbins, Option_t *opt){
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent) parent = CreateHistoGroup(dirname);
	if(parent->FindObject(hname)){
		Fatal("THistManager::CreateTH1", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TH1* h = new TH1D(hname, title, xbins.GetSize()-1, xbins.GetArray());
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    h->Sumw2();
	parent->Add(h);
	return h;
}

TH1* THistManager::CreateTH1(const char *name, const char *title, const TBinning &xbin, Option_t *opt){
  TArrayD myxbins;
  try{
    xbin.CreateBinEdges(myxbins);
  } catch(std::exception &e){
    Fatal("THistManager::CreateTH1", "Exception raised: %s", e.what());
  }
  return CreateTH1(name, title, myxbins, opt);
}

TH2* THistManager::CreateTH2(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, Option_t *opt){
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent) parent = CreateHistoGroup(dirname);
	if(parent->FindObject(hname)){
		Fatal("THistManager::CreateTH2", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TH2* h = new TH2D(hname, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    h->Sumw2();
	parent->Add(h);
	return h;
}

TH2* THistManager::CreateTH2(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy, const double *ybins, Option_t *opt){
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent) parent = CreateHistoGroup(dirname);
	if(parent->FindObject(hname)){
		Fatal("THistManager::CreateTH2", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TH2* h = new TH2D(hname, title, nbinsx, xbins, nbinsy, ybins);
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    h->Sumw2();
	parent->Add(h);
	return h;
}

TH2* THistManager::CreateTH2(const char *name, const char *title, const TArrayD &xbins, const TArrayD &ybins, Option_t *opt){
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent) parent = CreateHistoGroup(dirname);
	if(parent->FindObject(hname)){
		Fatal("THistManager::CreateTH2", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TH2* h = new TH2D(hname, title, xbins.GetSize() - 1, xbins.GetArray(), ybins.GetSize() - 1, ybins.GetArray());
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    h->Sumw2();
	parent->Add(h);
	return h;
}

TH2* THistManager::CreateTH2(const char *name, const char *title, const TBinning &xbins, const TBinning &ybins, Option_t *opt){
  TArrayD myxbins, myybins;
  try{
    xbins.CreateBinEdges(myxbins);
  } catch (std::exception &e) {
    Fatal("THistManager::CreateTH2 (x-dir)", "Exception raised: %s", e.what());
  }
  try{
    ybins.CreateBinEdges(myybins);
  } catch (std::exception &e) {
    Fatal("THistManager::CreateTH2 (y-dir)", "Exception raised: %s", e.what());
  }

  return CreateTH2(name, title, myxbins, myybins, opt);
}

TH3* THistManager::CreateTH3(const char* name, const char* title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, int nbinsz, double zmin, double zmax, Option_t *opt) {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent) parent = CreateHistoGroup(dirname);
	if(parent->FindObject(hname)){
		Fatal("THistManager::CreateTH3", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TH3* h = new TH3D(hname, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax);
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    h->Sumw2();
	parent->Add(h);
	return h;
}

TH3* THistManager::CreateTH3(const char* name, const char* title, int nbinsx, const double* xbins, int nbinsy, const double* ybins, int nbinsz, const double* zbins, Option_t *opt) {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent) parent = CreateHistoGroup(dirname);
	if(parent->FindObject(hname)){
		Fatal("THistManager::CreateTH3", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
  }
	TH3* h = new TH3D(hname, title, nbinsx, xbins, nbinsy, ybins, nbinsz, zbins);
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    h->Sumw2();
	parent->Add(h);
	return h;
}

TH3* THistManager::CreateTH3(const char* name, const char* title, const TArrayD& xbins, const TArrayD& ybins, const TArrayD& zbins, Option_t *opt) {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent) parent = CreateHistoGroup(dirname);
	if(parent->FindObject(hname)){
		Fatal("THistManager::CreateTH3", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TH3* h = new TH3D(hname, title, xbins.GetSize()-1, xbins.GetArray(), ybins.GetSize()-1, ybins.GetArray(), zbins.GetSize()-1, zbins.GetArray());
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    h->Sumw2();
	parent->Add(h);
	return h;
}

TH3* THistManager::CreateTH3(const char *name, const char *title, const TBinning &xbins, const TBinning &ybins, const TBinning &zbins, Option_t *opt){
  TArrayD myxbins, myybins, myzbins;
  try{
    xbins.CreateBinEdges(myxbins);
  } catch (std::exception &e) {
    Fatal("THistManager::CreateTH2 (x-dir)", "Exception raised: %s",  e.what());
  }
  try{
    ybins.CreateBinEdges(myybins);
  } catch (std::exception &e) {
    Fatal("THistManager::CreateTH2 (y-dir)", "Exception raised: %s", e.what());
  }
  try{
    zbins.CreateBinEdges(myzbins);
  } catch (std::exception &e) {
    Fatal("THistManager::CreateTH2 (z-dir)", "Exception raised: %s", e.what());
  }

  return CreateTH3(name, title, myxbins, myybins, myzbins);
}

THnSparse* THistManager::CreateTHnSparse(const char *name, const char *title, int ndim, const int *nbins, const double *min, const double *max, Option_t *opt) {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent) parent = CreateHistoGroup(dirname);
	if(parent->FindObject(hname)){
		Fatal("THistManager::CreateTHnSparse", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	THnSparse* h = new THnSparseD(hname, title, ndim, nbins, min, max);
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    h->Sumw2();
	parent->Add(h);
	return h;
}

THnSparse* THistManager::CreateTHnSparse(const char *name, const char *title, int ndim, const TAxis **axes, Option_t *opt) {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent) parent = CreateHistoGroup(dirname);
	if(parent->FindObject(hname)){
		Fatal("THistManager::CreateTHnSparse", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TArrayD xmin(ndim), xmax(ndim);
	TArrayI nbins(ndim);
	for(int idim = 0; idim < ndim; ++idim){
		const TAxis &myaxis = *(axes[idim]);
		nbins[idim] = myaxis.GetNbins();
		xmin[idim] = myaxis.GetXmin();
		xmax[idim] = myaxis.GetXmax();
	}
	THnSparseD *hsparse = new THnSparseD(hname, title, ndim, nbins.GetArray(), xmin.GetArray(), xmax.GetArray());
	for(int id = 0; id < ndim; ++id)
		*(hsparse->GetAxis(id)) = *(axes[id]);
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    hsparse->Sumw2();
	parent->Add(hsparse);
	return hsparse;
}

THnSparse* THistManager::CreateTHnSparse(const char *name, const char *title, int ndim, const TBinning **axes, Option_t *opt){
  TString dirname(basename(name)), hname(histname(name));
  THashList *parent(FindGroup(dirname));
  if(!parent) parent = CreateHistoGroup(dirname);
  if(parent->FindObject(hname)){
    Fatal("THistManager::CreateTHnSparse", "Object %s already exists in group %s", hname.Data(), dirname.Data());
    return 0;
  }
  TArrayD xmin(ndim), xmax(ndim);
  TArrayI nbins(ndim);
  std::vector<TArrayD> binnings;
  for(int idim = 0; idim < ndim; ++idim){
    const TBinning &myaxis = *(axes[idim]);
    TArrayD binEdges;
    myaxis.CreateBinEdges(binEdges);
    nbins[idim] = binEdges.GetSize()-1;
    xmin[idim] = binEdges[0];
    xmax[idim] = binEdges[binEdges.GetSize()-1];
    binnings.push_back(binEdges);
  }
  THnSparseD *hsparse = new THnSparseD(hname, title, ndim, nbins.GetArray(), xmin.GetArray(), xmax.GetArray());
  for(int id = 0; id < ndim; ++id){
    TArrayD &binEdges = binnings[id];
    hsparse->GetAxis(id)->Set(binEdges.GetSize()-1, binEdges.GetArray());
  }
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    hsparse->Sumw2();
  parent->Add(hsparse);
  return hsparse;
}

void THistManager::CreateTProfile(const char* name, const char* title, int nbinsX, double xmin, double xmax, Option_t *opt) {
  TString dirname(basename(name)), hname(histname(name));
  THashList *parent(FindGroup(dirname));
  if(!parent) parent = CreateHistoGroup(dirname);
  if(parent->FindObject(hname))
		Fatal("THistManager::CreateTProfile", "Object %s already exists in group %s", hname.Data(), dirname.Data());
  TProfile *hist = new TProfile(hname, title, nbinsX, xmin, xmax, opt);
  parent->Add(hist);
}

void THistManager::CreateTProfile(const char* name, const char* title, int nbinsX, const double* xbins, Option_t *opt) {
  TString dirname(basename(name)), hname(histname(name));
  THashList *parent(FindGroup(dirname));
  if(!parent) parent = CreateHistoGroup(dirname);
  if(parent->FindObject(hname))
		Fatal("THistManager::CreateTHnSparse", "Object %s already exists in group %s", hname.Data(), dirname.Data());
  TProfile *hist = new TProfile(hname, title, nbinsX, xbins, opt);
  parent->Add(hist);
}

void THistManager::CreateTProfile(const char* name, const char* title, const TArrayD& xbins, Option_t *opt){
  TString dirname(basename(name)), hname(histname(name));
  THashList *parent(FindGroup(dirname));
  if(!parent) parent = CreateHistoGroup(dirname);
  if(parent->FindObject(hname))
		Fatal("THistManager::CreateTHnSparse", "Object %s already exists in group %s", hname.Data(), dirname.Data());
  TProfile *hist = new TProfile(hname.Data(), title, xbins.GetSize()-1, xbins.GetArray(), opt);
  parent->Add(hist);
}

void THistManager::CreateTProfile(const char *name, const char *title, const TBinning &xbins, Option_t *opt){
  TArrayD myxbins;
  try{
    xbins.CreateBinEdges(myxbins);
  } catch (std::exception &e){
    Fatal("THistManager::CreateProfile", "Exception raised: %s", e.what());
  }
  CreateTProfile(name, title, myxbins, opt);
}

void THistManager::SetObject(TObject * const o, const char *group) {
	THashList *parent(FindGroup(group));
	if(!parent) CreateHistoGroup(group);
	if(parent->FindObject(o->GetName())){
		Fatal("THistManager::SetObject", "Parent %s does not exist", strcmp(group, "/") ? group : "");
		return;
	}
	if(!(dynamic_cast<THnBase *>(o) || dynamic_cast<TH1 *>(o))){
		Fatal("THistManager::SetObject",  "Object %s is not of a histogram type",o->GetName());
		return;
	}
	fHistos->Add(o);
}

void THistManager::FillTH1(const char *name, double x, double weight, Option_t *opt) {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent){
		Fatal("THistManager::FillTH1", "Parent group %s does not exist", dirname.Data());
		return;
	}
	TH1 *hist = dynamic_cast<TH1 *>(parent->FindObject(hname));
	if(!hist){
		Fatal("THistManager::FillTH1", "Histogram %s not found in parent group %s", hname.Data(), dirname.Data());
		return;
	}
	TString optionstring(opt);
	if(optionstring.Contains("w")){
	  // use bin width as weight
	  Int_t bin = hist->GetXaxis()->FindBin(x);
	  // check if not overflow or underflow bin
	  if(bin != 0 && bin != hist->GetXaxis()->GetNbins())
	    weight = 1./hist->GetXaxis()->GetBinWidth(bin);
	}
	hist->Fill(x, weight);
}

void THistManager::FillTH1(const char *name, const char *label, double weight, Option_t *opt) {
  TString dirname(basename(name)), hname(histname(name));
  THashList *parent(FindGroup(dirname));
  if(!parent){
    Fatal("THistManager::FillTH1", "Parent group %s does not exist", dirname.Data());
    return;
  }
  TH1 *hist = dynamic_cast<TH1 *>(parent->FindObject(hname));
  if(!hist){
    Fatal("THistManager::FillTH1", "Histogram %s not found in parent group %s", hname.Data(), dirname.Data());
    return;
  }
	TString optionstring(opt);
	if(optionstring.Contains("w")){
	  // use bin width as weight
	  // get bin for label
	  Int_t bin = hist->GetXaxis()->FindBin(label);
	  // check if not overflow or underflow bin
	  if(bin != 0 && bin != hist->GetXaxis()->GetNbins())
	    weight = 1./hist->GetXaxis()->GetBinWidth(bin);
	}
  hist->Fill(label, weight);
}

void THistManager::FillTH2(const char *name, double x, double y, double weight, Option_t *opt) {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent){
		Fatal("THistManager::FillTH2", "Parent group %s does not exist", dirname.Data());
		return;
	}
	TH2 *hist = dynamic_cast<TH2 *>(parent->FindObject(hname));
	if(!hist){
		Fatal("THistManager::FillTH2", "Histogram %s not found in parent group %s", hname.Data(), dirname.Data());
		return;
	}
	TString optstring(opt);
	Double_t myweight = optstring.Contains("w") ? 1. : weight;
	if(optstring.Contains("wx")){
	  Int_t binx = hist->GetXaxis()->FindBin(x);
	  if(binx != 0 && binx != hist->GetXaxis()->GetNbins()) myweight *= 1./hist->GetXaxis()->GetBinWidth(binx);
	}
	if(optstring.Contains("wy")){
	  Int_t biny = hist->GetYaxis()->FindBin(y);
	  if(biny != 0 && biny != hist->GetYaxis()->GetNbins()) myweight *= 1./hist->GetYaxis()->GetBinWidth(biny);
	}
	hist->Fill(x, y, myweight);
}

void THistManager::FillTH2(const char *name, double *point, double weight, Option_t *opt) {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent){
		Fatal("THistManager::FillTH2", "Parent group %s does not exist", dirname.Data());
		return;
	}
	TH2 *hist = dynamic_cast<TH2 *>(parent->FindObject(hname));
	if(!hist){
		Fatal("THistManager::FillTH2", "Histogram %s not found in parent group %s", hname.Data(), dirname.Data());
		return;
	}
	TString optstring(opt);
	Double_t myweight = optstring.Contains("w") ? 1. : weight;
	if(optstring.Contains("wx")){
	  Int_t binx = hist->GetXaxis()->FindBin(point[0]);
	  if(binx != 0 && binx != hist->GetXaxis()->GetNbins()) myweight *= 1./hist->GetXaxis()->GetBinWidth(binx);
	}
	if(optstring.Contains("wy")){
	  Int_t biny = hist->GetYaxis()->FindBin(point[1]);
	  if(biny != 0 && biny != hist->GetYaxis()->GetNbins()) myweight *= 1./hist->GetYaxis()->GetBinWidth(biny);
	}
	hist->Fill(point[0], point[1], weight);
}

void THistManager::FillTH2(const char *name, const char *labelX, const char *labelY, double weight, Option_t *opt) {
  TString dirname(basename(name)), hname(histname(name));
  THashList *parent(FindGroup(dirname));
  if(!parent){
    Fatal("THistManager::FillTH2", "Parent group %s does not exist", dirname.Data());
    return;
  }
  TH2 *hist = dynamic_cast<TH2 *>(parent->FindObject(hname));
  if(!hist){
    Fatal("THistManager::FillTH2", "Histogram %s not found in parent group %s", hname.Data(), dirname.Data());
    return;
  }
  TString optstring(opt);
  Double_t myweight = optstring.Contains("w") ? 1. : weight;
  if(optstring.Contains("wx")){
    Int_t binx = hist->GetXaxis()->FindBin(labelY);
    if(binx != 0 && binx != hist->GetXaxis()->GetNbins()) myweight *= 1./hist->GetXaxis()->GetBinWidth(binx);
  }
  if(optstring.Contains("wy")){
    Int_t biny = hist->GetYaxis()->FindBin(labelX);
    if(biny != 0 && biny != hist->GetYaxis()->GetNbins()) myweight *= 1./hist->GetYaxis()->GetBinWidth(biny);
  }
  hist->Fill(labelX, labelY, weight);
}

void THistManager::FillTH3(const char* name, double x, double y, double z, double weight, Option_t *opt) {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent){
		Fatal("THistManager::FillTH3", "Parent group %s does not exist", dirname.Data());
		return;
	}
	TH3 *hist = dynamic_cast<TH3 *>(parent->FindObject(hname));
	if(!hist){
		Fatal("THistManager::FillTH3", "Histogram %s not found in parent group %s", hname.Data(), dirname.Data());
		return;
	}
	TString optstring(opt);
	Double_t myweight = optstring.Contains("w") ? 1. : weight;
	if(optstring.Contains("wx")){
	  Int_t binx = hist->GetXaxis()->FindBin(x);
	  if(binx != 0 && binx != hist->GetXaxis()->GetNbins()) myweight *= 1./hist->GetXaxis()->GetBinWidth(binx);
	}
	if(optstring.Contains("wy")){
	  Int_t biny = hist->GetYaxis()->FindBin(y);
	  if(biny != 0 && biny != hist->GetYaxis()->GetNbins()) myweight *= 1./hist->GetYaxis()->GetBinWidth(biny);
	}
	if(optstring.Contains("wz")){
	  Int_t binz = hist->GetZaxis()->FindBin(z);
	  if(binz != 0 && binz != hist->GetZaxis()->GetNbins()) myweight *= 1./hist->GetZaxis()->GetBinWidth(binz);
	}
	hist->Fill(x, y, z, weight);
}

void THistManager::FillTH3(const char* name, const double* point, double weight, Option_t *opt) {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent){
		Fatal("THistManager::FillTH3", "Parent group %s does not exist", dirname.Data());
		return;
	}
	TH3 *hist = dynamic_cast<TH3 *>(parent->FindObject(hname));
	if(!hist){
		Fatal("THistManager::FillTH3", "Histogram %s not found in parent group %s", hname.Data(), dirname.Data());
		return;
	}
	TString optstring(opt);
	Double_t myweight = optstring.Contains("w") ? 1. : weight;
	if(optstring.Contains("wx")){
	  Int_t binx = hist->GetXaxis()->FindBin(point[0]);
	  if(binx != 0 && binx != hist->GetXaxis()->GetNbins()) myweight *= 1./hist->GetXaxis()->GetBinWidth(binx);
	}
	if(optstring.Contains("wy")){
	  Int_t biny = hist->GetYaxis()->FindBin(point[1]);
	  if(biny != 0 && biny != hist->GetYaxis()->GetNbins()) myweight *= 1./hist->GetYaxis()->GetBinWidth(biny);
	}
	if(optstring.Contains("wz")){
	  Int_t binz = hist->GetZaxis()->FindBin(point[2]);
	  if(binz != 0 && binz != hist->GetZaxis()->GetNbins()) myweight *= 1./hist->GetZaxis()->GetBinWidth(binz);
	}
	hist->Fill(point[0], point[1], point[2], weight);
}

void THistManager::FillTHnSparse(const char *name, const double *x, double weight, Option_t *opt) {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent){
		Fatal("THistManager::FillTHnSparse", "Parent group %s does not exist", dirname.Data());
		return;
	}
	THnSparseD *hist = dynamic_cast<THnSparseD *>(parent->FindObject(hname));
	if(!hist){
		Fatal("THistManager::FillTHnSparse", "Histogram %s not found in parent group %s", hname.Data(), dirname.Data());
		return;
	}
	TString optstring(opt);
	Double_t myweight = optstring.Contains("w") ? 1. : weight;
	for(Int_t iaxis = 0; iaxis < hist->GetNdimensions(); iaxis++){
	  std::stringstream weighthandler;
	  weighthandler << "w" << iaxis;
	  if(optstring.Contains(weighthandler.str().c_str())){
	    Int_t bin = hist->GetAxis(iaxis)->FindBin(x[iaxis]);
	    if(bin != 0 && bin != hist->GetAxis(iaxis)->GetNbins()) myweight *= hist->GetAxis(iaxis)->GetBinWidth(bin);
	  }
	}

	hist->Fill(x, weight);
}

void THistManager::FillProfile(const char* name, double x, double y, double weight){
  TString dirname(basename(name)), hname(histname(name));
  THashList *parent(FindGroup(dirname));
  if(!parent)
		Fatal("THistManager::FillTProfile", "Parent group %s does not exist", dirname.Data());
  TProfile *hist = dynamic_cast<TProfile *>(parent->FindObject(hname));
  if(!hist)
		Fatal("THistManager::FillTProfile", "Histogram %s not found in parent group %s", hname.Data(), dirname.Data());
  hist->Fill(x, y, weight);
}

TObject *THistManager::FindObject(const char *name) const {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname));
	if(!parent) return NULL;
	return parent->FindObject(hname);
}

TObject* THistManager::FindObject(const TObject* obj) const {
	TString dirname(basename(obj->GetName())), hname(histname(obj->GetName()));
	THashList *parent(FindGroup(dirname));
	if(!parent) return NULL;
	return parent->FindObject(hname);
}

THashList *THistManager::FindGroup(const char *dirname) const {
	if(!strlen(dirname) || !strcmp(dirname, "/")) return fHistos;
	// recursive find - avoids tokenizing filename
	THashList *parentlist = FindGroup(basename(dirname));
	if(parentlist) return static_cast<THashList *>(parentlist->FindObject(histname(dirname)));
	return nullptr;
}

TString THistManager::basename(const TString &path) const {
	int index = path.Last('/');
	if(index < 0) return "";  // no directory structure
	return TString(path(0, index));
}

TString THistManager::histname(const TString &path) const {
	int index = path.Last('/');
	if(index < 0) return path;    // no directory structure
	return TString(path(index+1, path.Length() - (index+1)));
}

//////////////////////////////////////////////////////////
///                                                    ///
/// Implementation of THistManager::iterator           ///
///                                                    ///
//////////////////////////////////////////////////////////

THistManager::iterator::iterator(const THistManager * hmgr, Int_t currentpos, THMIDirection_t dir):
    fkArray(hmgr),
    fCurrentPos(),
    fNext(),
    fDirection(dir)
{}

THistManager::iterator::iterator(const THistManager::iterator &ref):
    fkArray(ref.fkArray),
    fCurrentPos(ref.fCurrentPos),
    fNext(ref.fNext),
    fDirection(ref.fDirection)
{}

THistManager::iterator &THistManager::iterator::operator=(const THistManager::iterator &ref){
  if(this != &ref){
    fkArray = ref.fkArray;
    fCurrentPos = ref.fCurrentPos;
    fNext = ref.fNext;
    fDirection = ref.fDirection;
  }
  return *this;
}

bool THistManager::iterator::operator!=(const THistManager::iterator &other) const{
  return fCurrentPos == other.fCurrentPos;
}

THistManager::iterator &THistManager::iterator::operator++(){
  if(fDirection == kTHMIforward)
    fCurrentPos++;
  else
    fCurrentPos--;
  return *this;
}

THistManager::iterator THistManager::iterator::operator++(int){
  iterator tmp(*this);
  operator++();
  return tmp;
}

THistManager::iterator &THistManager::iterator::operator--(){
  if(fDirection == kTHMIforward)
    fCurrentPos--;
  else
    fCurrentPos++;
  return *this;
};

THistManager::iterator THistManager::iterator::operator--(int){
  iterator tmp(*this);
  operator--();
  return tmp;
};

TObject *THistManager::iterator::operator*() const{
  if(fCurrentPos >=0 && fCurrentPos < fkArray->GetListOfHistograms()->GetEntries())
    return fkArray->GetListOfHistograms()->At(fCurrentPos);
  return NULL;
}


//////////////////////////////////////////////////////////////////////////////////////////////
///
///  Unit tests
///
//////////////////////////////////////////////////////////////////////////////////////////////

namespace TestTHistManager {

  int THistManagerTestSuite::TestBuildSimpleHistograms(){
    THistManager testmgr("testmgr");

    // Create histogram of each type
    testmgr.CreateTH1("Test1D", "Test Histogram 1D", 1, 0., 1.);
    testmgr.CreateTH2("Test2D", "Test Histogram 2D", 2, 0., 2., 10., 0., 10);
    testmgr.CreateTH3("Test3D", "Test Histogram 3D", 3, 2, 6., 10., 0., 10., 50., 0., 50.);
    int nbins[3] = {3, 3, 3}; double min[3] = {0., 0., 0}, max[3] = {6, 9, 12};   // dimensions for THnSparseTest
    testmgr.CreateTHnSparse("TestNSparse", "Test Histogram NSparse", 3, nbins, min, max);
    testmgr.CreateTProfile("TestProfile", "Test TProfile", 1, 0., 1);

    // Find histograms in the list (evaluate test)
    // Tell user why test has failed
    Bool_t found(true);
    if(!dynamic_cast<TH1 *>(testmgr.GetListOfHistograms()->FindObject("Test1D"))){
      std::cout << "Not found: Test1D" << std::endl;
      found = false;
    }
    if(!dynamic_cast<TH2 *>(testmgr.GetListOfHistograms()->FindObject("Test2D"))){
      std::cout << "Not found: Test2D" << std::endl;
      found = false;
    }
    if(!dynamic_cast<TH3 *>(testmgr.GetListOfHistograms()->FindObject("Test3D"))){
      std::cout << "Not found: Test3D" << std::endl;
      found = false;
    }
    if(!dynamic_cast<THnSparse *>(testmgr.GetListOfHistograms()->FindObject("TestNSparse"))){
      std::cout << "Not found: TestNSparse" << std::endl;
      found = false;
    }
    if(!dynamic_cast<TProfile *>(testmgr.GetListOfHistograms()->FindObject("TestProfile"))){
      std::cout << "Not found: TestProfile" << std::endl;
      found = false;
    }

    return found ? 0 : 1;
  }

  int THistManagerTestSuite::TestBuildGroupedHistograms(){
    THistManager testmgr("testmgr");

    // Creating 3 groups
    testmgr.CreateTH1("Group1/Test1", "Test Histogram 1 in group 1", 1, 0., 1.);
    testmgr.CreateTH1("Group1/Test2", "Test Histogram 2 in group 1", 1, 0., 1.);
    testmgr.CreateTH2("Group2/Test1", "Test Histogram 1 in group 2", 1, 0., 1., 1, 0., 1.);
    testmgr.CreateTH2("Group2/Test2", "Test Histogram 2 in group 2", 1, 0., 1., 1, 0., 1.);
    testmgr.CreateTH3("Group3/Test1", "Test Histogram 1 in group 3", 1, 0., 1., 1, 0., 1., 1, 0., 1.);
    testmgr.CreateTH3("Group3/Test2", "Test Histogram 2 in group 3", 1, 0., 1., 1, 0., 1., 1, 0., 1.);
    // Creating histogram in a subgroup
    testmgr.CreateTProfile("Group4/Subgroup1/Test1", "Test histogram for subgroup handling", 1, 0., 1);

    // Evalutate test
    // Tell user why test has failed
    Bool_t found(true);
    THashList *currentdir(nullptr);
    if(!(currentdir = dynamic_cast<THashList *>(testmgr.GetListOfHistograms()->FindObject("Group1")))){
      std::cout << "Not found: Group1" << std::endl;
      found = false;
    } else {
      if(!dynamic_cast<TH1 *>(currentdir->FindObject("Test1"))){
        std::cout << "Not found in Group1: Test1" << std::endl;
        found = false;
      }
      if(!dynamic_cast<TH1 *>(currentdir->FindObject("Test2"))){
        std::cout << "Not found in Group1: Test2" << std::endl;
        found = false;
      }
    }
    if(!(currentdir = dynamic_cast<THashList *>(testmgr.GetListOfHistograms()->FindObject("Group2")))){
      std::cout << "Not found: Group2" << std::endl;
      found = false;
    } else {
      if(!dynamic_cast<TH2 *>(currentdir->FindObject("Test1"))){
        std::cout << "Not found in Group2: Test1" << std::endl;
        found = false;
      }
      if(!dynamic_cast<TH2 *>(currentdir->FindObject("Test2"))){
        std::cout << "Not found in Group2: Test2" << std::endl;
        found = false;
      }
    }
    if(!(currentdir = dynamic_cast<THashList *>(testmgr.GetListOfHistograms()->FindObject("Group3")))){
      std::cout << "Not found: Group3" << std::endl;
      found = false;
    } else {
      if(!static_cast<TH3 *>(currentdir->FindObject("Test1"))){
        std::cout << "Not found in Group3: Test1" << std::endl;
        found = false;
      }
      if(!static_cast<TH3 *>(currentdir->FindObject("Test2"))){
        std::cout << "Not found in Group3: Test2" << std::endl;
        found = false;
      }
    }
    if(!(currentdir = dynamic_cast<THashList *>(testmgr.GetListOfHistograms()->FindObject("Group4")))){
      std::cout << "Not found: Group4" << std::endl;
      found = false;
    } else {
      if(!(currentdir = dynamic_cast<THashList *>(currentdir->FindObject("Subgroup1")))){
        std::cout << "Not found in Group4: Subgroup1" << std::endl;
        found = false;
      } else {
        if(!dynamic_cast<TH1 *>(currentdir->FindObject("Test1"))){
          std::cout << "Not found in Subgroup1: Test1" << std::endl;
          found = false;
        }
      }
    }
    return found ? 0 : 1;
  }

  int THistManagerTestSuite::TestFillSimpleHistograms(){
    THistManager testmgr("testmgr");

    testmgr.CreateTH1("Test1", "Test fill 1D histogram", 1, 0., 1.);
    testmgr.CreateTH2("Test2", "Test fill 2D histogram", 1, 0., 1., 1, 0., 1.);
    testmgr.CreateTH3("Test3", "Test fill 3D histogram", 1, 0., 1., 1, 0., 1., 1, 0., 1.);
    int nbins[4] = {1,1,1,1}; double min[4] = {0.,0.,0.,0.}, max[4] = {1.,1.,1.,1.};
    testmgr.CreateTHnSparse("TestN", "Test Fill THnSparse", 4, nbins, min, max);
    testmgr.CreateTProfile("TestProfile", "Test fill Profile histogram", 1, 0., 1.);

    double point[4] = {0.5, 0.5, 0.5, 0.5};
    for(int i = 0; i < 100; i++){
      testmgr.FillTH1("Test1", 0.5);
      testmgr.FillTH2("Test2", 0.5, 0.5);
      testmgr.FillTH3("Test3", 0.5, 0.5, 0.5);
      testmgr.FillProfile("TestProfile", 0.5, 1.);
      testmgr.FillTHnSparse("TestN", point);
    }

    // Evalutate test
    // tell user why test has failed
    bool success(true);

    TH1 *test1 = dynamic_cast<TH1 *>(testmgr.GetListOfHistograms()->FindObject("Test1"));
    if(test1){
      if(TMath::Abs(test1->GetBinContent(1) - 100) > DBL_EPSILON){
        std::cout << "Test1: Mismatch in values, expected 100, found " <<  test1->GetBinContent(1) << std::endl;
        success = false;
      }
    } else {
      std::cout << "Not found: Test1" << std::endl;
      success = false;
    }

    TH2 *test2 = dynamic_cast<TH2 *>(testmgr.GetListOfHistograms()->FindObject("Test2"));
    if(test2){
      if(TMath::Abs(test2->GetBinContent(1, 1) - 100) > DBL_EPSILON){
        std::cout << "Test2: Mismatch in values, expected 100, found " <<  test2->GetBinContent(1,1) << std::endl;
        success = false;
      }
    } else {
      std::cout << "Not found: Test2" << std::endl;
      success = false;
    }

    TH3 *test3 = dynamic_cast<TH3 *>(testmgr.GetListOfHistograms()->FindObject("Test3"));
    if(test3){
      if(TMath::Abs(test3->GetBinContent(1, 1, 1) - 100) > DBL_EPSILON){
        std::cout << "Test3: Mismatch in values, expected 100, found " <<  test3->GetBinContent(1,1,1) << std::endl;
        success = false;
      }
    } else {
      std::cout << "Not found: Test3" << std::endl;
      success = false;
    }

    THnSparse *testN = dynamic_cast<THnSparse *>(testmgr.GetListOfHistograms()->FindObject("TestN"));
    if(testN){
      int index[4] = {1,1,1,1};
      if(TMath::Abs(testN->GetBinContent(index) - 100) > DBL_EPSILON){
        std::cout << "TestN: Mismatch in values, expected 100, found " <<  testN->GetBinContent(index) << std::endl;
        success = false;
      }
    } else {
      std::cout << "Not found: TestN" << std::endl;
      success = false;
    }

    TProfile *testProfile = dynamic_cast<TProfile *>(testmgr.GetListOfHistograms()->FindObject("TestProfile"));
    if(testProfile){
      if(TMath::Abs(testProfile->GetBinContent(1) - 1) > DBL_EPSILON){
        std::cout << "TestProfile: Mismatch in values, expected 1, found " <<  testProfile->GetBinContent(1) << std::endl;
        success = false;
      }
    } else {
      std::cout << "Not found: TestProfile" << std::endl;
      success = false;
    }

    return success ? 0 : 1;
  }


  int THistManagerTestSuite::TestFillGroupedHistograms(){
    THistManager testmgr("testmgr");

    // Creating 3 groups, 1 with 1D and 1 with 2D histograms, and a third with a subgroup and a TProfile
    testmgr.CreateTH1("Group1/Test1", "Test 1 Group 1D", 1, 0., 1.);
    testmgr.CreateTH1("Group1/Test2", "Test 2 Group 1D", 1, 0., 1.);
    testmgr.CreateTH2("Group2/Test1", "Test 1 Group 2D", 1, 0., 1., 1, 0., 1.);
    testmgr.CreateTH2("Group2/Test2", "Test 2 Group 2D", 1, 0., 1., 1, 0., 1.);
    testmgr.CreateTProfile("Group3/Subgroup1/Test1", "Test 1 with subgroup", 1, 0., 1.);

    for(int i = 0; i < 100; i++){
      testmgr.FillTH1("Group1/Test1", 0.5);
      testmgr.FillTH1("Group1/Test2", 0.5);
      testmgr.FillTH2("Group2/Test1", 0.5, 0.5);
      testmgr.FillTH2("Group2/Test2", 0.5, 0.5);
      testmgr.FillProfile("Group3/Subgroup1/Test1", 0.5, 1);
    }

    // Evaluate test
    bool success(true);

    THashList *currentdir = dynamic_cast<THashList *>(testmgr.GetListOfHistograms()->FindObject("Group1"));
    if(currentdir){
      TH1 *test1 = dynamic_cast<TH1 *>(currentdir->FindObject("Test1"));
      if(test1){
        if(TMath::Abs(test1->GetBinContent(1) - 100) > DBL_EPSILON){
          std::cout << "Group1/Test1: Value mismatch: expected 100, found " << test1->GetBinContent(1) << std::endl;
          success = false;
        }
      } else {
        std::cout << "Not found in Group1: Test1" << std::endl;
        success = false;
      }
      test1 = dynamic_cast<TH1 *>(currentdir->FindObject("Test2"));
      if(test1){
        if(TMath::Abs(test1->GetBinContent(1) - 100) > DBL_EPSILON){
          std::cout << "Group1/Test2: Value mismatch: expected 100, found " << test1->GetBinContent(1) << std::endl;
          success = false;
        }
      } else {
        std::cout << "Not found in Group1: Test2" << std::endl;
        success = false;
      }
    } else {
      std::cout << "Not found: Group1" << std::endl;
      success = false;
    }

    currentdir = dynamic_cast<THashList *>(testmgr.GetListOfHistograms()->FindObject("Group2"));
    if(currentdir){
      TH2 *test2 = dynamic_cast<TH2 *>(currentdir->FindObject("Test1"));
      if(test2){
        if(TMath::Abs(test2->GetBinContent(1,1) - 100) > DBL_EPSILON){
          std::cout << "Group2/Test1: Value mismatch: expected 100, found " << test2->GetBinContent(1,1) << std::endl;
          success = false;
        }
      } else {
        std::cout << "Not found in Group2: Test1" << std::endl;
        success = false;
      }
      test2 = dynamic_cast<TH2 *>(currentdir->FindObject("Test2"));
      if(test2){
        if(TMath::Abs(test2->GetBinContent(1,1) - 100) > DBL_EPSILON){
          std::cout << "Group2/Test2: Value mismatch: expected 100, found " << test2->GetBinContent(1,1) << std::endl;
          success = false;
        }
      } else {
        std::cout << "Not found in Group2: Test2" << std::endl;
        success = false;
      }
    } else {
      std::cout << "Not found: Group2" << std::endl;
      success = false;
    }

    currentdir = dynamic_cast<THashList *>(testmgr.GetListOfHistograms()->FindObject("Group3"));
    if(currentdir){
      currentdir = dynamic_cast<THashList *>(currentdir->FindObject("Subgroup1"));
      if(currentdir){
        TProfile *testprofile = dynamic_cast<TProfile *>(currentdir->FindObject("Test1"));
        if(testprofile){
          if(TMath::Abs(testprofile->GetBinContent(1) - 1) > DBL_EPSILON){
            std::cout << "Group3/Subgroup1/Test1: Value mismatch: expected 1, found " << testprofile->GetBinContent(1) << std::endl;
            success = false;
          }
        } else {
          std::cout << "Not found in Group3/Subgroup1: Test1" << std::endl;
          success = false;
        }
      } else {
        std::cout << "Not found in Group3: Subgroup1" << std::endl;
        success = false;
      }
    } else {
      std::cout << "Not found: Group3" << std::endl;
      success = false;
    }
    return success ? 0 : 1;
  }

  int TestRunAll(){
    int testresult(0);
    THistManagerTestSuite testsuite;

    std::cout << "Running test: Build simple" << std::endl;
    testresult += testsuite.TestBuildSimpleHistograms();
    std::cout << "Result after test: " << testresult << std::endl;

    std::cout << "Running test: Build grouped" << std::endl;
    testresult += testsuite.TestBuildGroupedHistograms();
    std::cout << "Result after test: " << testresult << std::endl;

    std::cout << "Running test: Fill Simple" << std::endl;
    testresult += testsuite.TestFillSimpleHistograms();
    std::cout << "Result after test: " << testresult << std::endl;

    std::cout << "Running test: Fill Grouped" << std::endl;
    testresult += testsuite.TestFillGroupedHistograms();
    std::cout << "Result after test: " << testresult << std::endl;

    return testresult;
  }

  int TestRunBuildSimple(){
    THistManagerTestSuite testsuite;
    return testsuite.TestBuildSimpleHistograms();
  }

  int TestRunBuildGrouped(){
    THistManagerTestSuite testsuite;
    return testsuite.TestBuildGroupedHistograms();
  }

  int TestRunFillSimple(){
    THistManagerTestSuite testsuite;
    return testsuite.TestFillSimpleHistograms();
  }

  int TestRunFillGrouped(){
    THistManagerTestSuite testsuite;
    return testsuite.TestFillGroupedHistograms();
  }
}
