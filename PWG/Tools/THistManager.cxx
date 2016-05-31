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
#include <cstring>
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

THashList* THistManager::CreateHistoGroup(const char *groupname, const char *parent) {
	THashList *parentgroup = FindGroup(parent);
	if(!parentgroup){
		Fatal("THistManager::CreateHistoGroup", "Parent group %s does not exist", parent);
		return 0;
	}
	THashList *childgroup = new THashList();
	childgroup->SetName(groupname);
	parentgroup->Add(childgroup);
	return childgroup;
}

TH1* THistManager::CreateTH1(const char *name, const char *title, int nbins, double xmin, double xmax, Option_t *opt){
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent){
		Fatal("THistManager::CreateTH1", "Parent %s does not exist", dirname.Data());
		return 0 ;
	}
	if(parent->FindObject(hname.Data())){
		Fatal("THistManager::CreateTH1", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TH1* h = new TH1D(hname.Data(), title, nbins, xmin, xmax);
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
	if(!parent){
		Fatal("THistManager::CreateTH1", "Parent %s does not exist", dirname.Data());
		return 0;
	}
	if(parent->FindObject(hname.Data())){
		Fatal("THistManager::CreateTH1", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TH1* h = new TH1D(hname.Data(), title, nbins, xbins);
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
	if(!parent){
		Fatal("THistManager::CreateTH1", "Parent %s does not exist", dirname.Data());
		return 0;
	}
	if(parent->FindObject(hname.Data())){
		Fatal("THistManager::CreateTH1", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TH1* h = new TH1D(hname.Data(), title, xbins.GetSize()-1, xbins.GetArray());
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
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent){
		Fatal("THistManager::CreateTH2", "Parent %s does not exist", dirname.Data());
		return 0;
	}
	if(parent->FindObject(hname.Data())){
		Fatal("THistManager::CreateTH2", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TH2* h = new TH2D(hname.Data(), title, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    h->Sumw2();
	parent->Add(h);
	return h;
}

TH2* THistManager::CreateTH2(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy, const double *ybins, Option_t *opt){
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent){
		Fatal("THistManager::CreateTH2", "Parent %s does not exist", dirname.Data());
		return 0;
	}
	if(parent->FindObject(hname.Data())){
		Fatal("THistManager::CreateTH2", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TH2* h = new TH2D(hname.Data(), title, nbinsx, xbins, nbinsy, ybins);
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    h->Sumw2();
	parent->Add(h);
	return h;
}

TH2* THistManager::CreateTH2(const char *name, const char *title, const TArrayD &xbins, const TArrayD &ybins, Option_t *opt){
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent){
		Fatal("THistManager::CreateTH2", "Parent %s does not exist", dirname.Data());
		return 0;
	}
	if(parent->FindObject(hname.Data())){
		Fatal("THistManager::CreateTH2", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TH2* h = new TH2D(hname.Data(), title, xbins.GetSize() - 1, xbins.GetArray(), ybins.GetSize() - 1, ybins.GetArray());
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
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent){
		Fatal("THistManager::CreateTH3", "Parent %s does not exist", dirname.Data());
		return 0;
	}
	if(parent->FindObject(hname.Data())){
		Fatal("THistManager::CreateTH3", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TH3* h = new TH3D(hname.Data(), title, nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax);
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    h->Sumw2();
	parent->Add(h);
	return h;
}

TH3* THistManager::CreateTH3(const char* name, const char* title, int nbinsx, const double* xbins, int nbinsy, const double* ybins, int nbinsz, const double* zbins, Option_t *opt) {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent){
		Fatal("THistManager::CreateTH3", "Parent %s does not exist", dirname.Data());
		return 0;
	}
	if(parent->FindObject(hname.Data())){
		Fatal("THistManager::CreateTH3", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
  }
	TH3* h = new TH3D(hname.Data(), title, nbinsx, xbins, nbinsy, ybins, nbinsz, zbins);
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    h->Sumw2();
	parent->Add(h);
	return h;
}

TH3* THistManager::CreateTH3(const char* name, const char* title, const TArrayD& xbins, const TArrayD& ybins, const TArrayD& zbins, Option_t *opt) {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent){
		Fatal("THistManager::CreateTH3", "Parent %s does not exist", dirname.Data());
		return 0;
	}
	if(parent->FindObject(hname.Data())){
		Fatal("THistManager::CreateTH3", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	TH3* h = new TH3D(hname.Data(), title, xbins.GetSize()-1, xbins.GetArray(), ybins.GetSize()-1, ybins.GetArray(), zbins.GetSize()-1, zbins.GetArray());
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
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent){
		Fatal("THistManager::CreateTHnSparse", "Parent %s does not exist", dirname.Data());
		return 0;
	}
	if(parent->FindObject(hname.Data())){
		Fatal("THistManager::CreateTHnSparse", "Object %s already exists in group %s", hname.Data(), dirname.Data());
		return 0;
	}
	THnSparse* h = new THnSparseD(hname.Data(), title, ndim, nbins, min, max);
  TString optionstring(opt);
  optionstring.ToLower();
  if(optionstring.Contains("s"))
    h->Sumw2();
	parent->Add(h);
	return h;
}

THnSparse* THistManager::CreateTHnSparse(const char *name, const char *title, int ndim, const TAxis **axes, Option_t *opt) {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent){
		Fatal("THistManager::CreateTHnSparse", "Parent %s does not exist", dirname.Data());
		return 0;
	}
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
	THnSparseD *hsparse = new THnSparseD(hname.Data(), title, ndim, nbins.GetArray(), xmin.GetArray(), xmax.GetArray());
	for(int id = 0; id < ndim; ++id)
		*(hsparse->GetAxis(id)) = *(axes[id]);
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
  if(!parent)
		Fatal("THistManager::CreateTProfile", "Parent %s does not exist", dirname.Data());
  if(parent->FindObject(hname.Data()))
		Fatal("THistManager::CreateTProfile", "Object %s already exists in group %s", hname.Data(), dirname.Data());
  TProfile *hist = new TProfile(hname.Data(), title, nbinsX, xmin, xmax, opt);
  parent->Add(hist);
}

void THistManager::CreateTProfile(const char* name, const char* title, int nbinsX, const double* xbins, Option_t *opt) {
  TString dirname(basename(name)), hname(histname(name));
  THashList *parent(FindGroup(dirname));
  if(!parent)
		Fatal("THistManager::CreateTHnSparse", "Parent %s does not exist", dirname.Data());
  if(parent->FindObject(hname.Data()))
		Fatal("THistManager::CreateTHnSparse", "Object %s already exists in group %s", hname.Data(), dirname.Data());
  TProfile *hist = new TProfile(hname.Data(), title, nbinsX, xbins, opt);
  parent->Add(hist);
}

void THistManager::CreateTProfile(const char* name, const char* title, const TArrayD& xbins, Option_t *opt){
  TString dirname(basename(name)), hname(histname(name));
  THashList *parent(FindGroup(dirname));
  if(!parent)
		Fatal("THistManager::CreateTHnSparse", "Parent %s does not exist", dirname.Data());
  if(parent->FindObject(hname.Data()))
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
	if(!parent){
		Fatal("THistManager::SetObject", "Parent %s does not exist", strcmp(group, "/") ? group : "");
		return;
	}
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
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent){
		Fatal("THistManager::FillTH1", "Parnt group %s does not exist", dirname.Data());
		return;
	}
	TH1 *hist = dynamic_cast<TH1 *>(parent->FindObject(hname.Data()));
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
  THashList *parent(FindGroup(dirname.Data()));
  if(!parent){
    Fatal("THistManager::FillTH1", "Parnt group %s does not exist", dirname.Data());
    return;
  }
  TH1 *hist = dynamic_cast<TH1 *>(parent->FindObject(hname.Data()));
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
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent){
		Fatal("THistManager::FillTH2", "Parent group %s does not exist", dirname.Data());
		return;
	}
	TH2 *hist = dynamic_cast<TH2 *>(parent->FindObject(hname.Data()));
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
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent){
		Fatal("THistManager::FillTH2", "Parent group %s does not exist", dirname.Data());
		return;
	}
	TH2 *hist = dynamic_cast<TH2 *>(parent->FindObject(hname.Data()));
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

void THistManager::FillTH3(const char* name, double x, double y, double z, double weight, Option_t *opt) {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent){
		Fatal("THistManager::FillTH3", "Parent group %s does not exist", dirname.Data());
		return;
	}
	TH3 *hist = dynamic_cast<TH3 *>(parent->FindObject(hname.Data()));
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
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent){
		Fatal("THistManager::FillTH3", "Parent group %s does not exist", dirname.Data());
		return;
	}
	TH3 *hist = dynamic_cast<TH3 *>(parent->FindObject(hname.Data()));
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
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent){
		Fatal("THistManager::FillTHnSparse", "Parent group %s does not exist", dirname.Data());
		return;
	}
	THnSparseD *hist = dynamic_cast<THnSparseD *>(parent->FindObject(hname.Data()));
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
  THashList *parent(FindGroup(dirname.Data()));
  if(!parent)
		Fatal("THistManager::FillTProfile", "Parent group %s does not exist", dirname.Data());
  TProfile *hist = dynamic_cast<TProfile *>(parent->FindObject(hname.Data()));
  if(!hist)
		Fatal("THistManager::FillTProfile", "Histogram %s not found in parent group %s", hname.Data(), dirname.Data());
  hist->Fill(x, y, weight);
}

TObject *THistManager::FindObject(const char *name) const {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent) return NULL;
	return parent->FindObject(hname);
}

TObject* THistManager::FindObject(const TObject* obj) const {
	TString dirname(basename(obj->GetName())), hname(histname(obj->GetName()));
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent) return NULL;
	return parent->FindObject(hname);
}

THashList *THistManager::FindGroup(const char *dirname) const {
	if(!strlen(dirname) || !strcmp(dirname, "/")) return fHistos;
	std::vector<std::string> tokens;
	TokenizeFilename(dirname, "/", tokens);
	THashList *currentdir(fHistos);
	for(std::vector<std::string>::iterator it = tokens.begin(); it != tokens.end(); ++it){
		currentdir = dynamic_cast<THashList *>(currentdir->FindObject(it->c_str()));
		if(!currentdir) break;
	}
	return currentdir;
}

void THistManager::TokenizeFilename(const char *name, const char *delim, std::vector<std::string> &listoftokens) const {
	TString s(name);
	TObjArray *arr = s.Tokenize(delim);
	TObjString *ostr(NULL);
	TIter toks(arr);
	while((ostr = dynamic_cast<TObjString *>(toks()))){
		listoftokens.push_back(std::string(ostr->String().Data()));
	}
	delete arr;
}

TString THistManager::basename(const char *path) const {
	TString s(path);
	int index = s.Last('/');
	if(index < 0) return "";  // no directory structure
	return TString(s(0, index)).Data();
}

TString THistManager::histname(const char *path) const {
	TString s(path);
	int index = s.Last('/');
	if(index < 0) return path;    // no directory structure
	return TString(s(index+1, s.Length() - (index+1))).Data();
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
