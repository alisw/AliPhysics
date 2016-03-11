/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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
#include <string>
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

#include "THistManager.h"

/// \cond CLASSIMP
ClassImp(THistManager)
/// \endcond

/**
 * Default constructor, only initialising pointers with 0
 */
THistManager::THistManager():
		TNamed(),
		fHistos(NULL),
		fIsOwner(true)
{
}

/**
 * Main constructor, creating also a list for the histograms
 *
 * @param name Name of the object (list named accordingly)
 */
THistManager::THistManager(const char *name):
		TNamed(name, Form("Histogram container %s", name)),
		fHistos(NULL),
		fIsOwner(true)
{
	fHistos = new THashList();
	fHistos->SetName(Form("histos%s", name));
	fHistos->SetOwner();
}

/**
 * Destructor, deletes the list of histograms if it is the owner
 */
THistManager::~THistManager(){
	if(fHistos && fIsOwner) delete fHistos;
}

/**
 * Create a new group of histograms within a parent group. Groups are represented as list. The default parent is
 * always the top list. List name structure accouding to unix paths (i.e. top list /, hirarchies separated by /).
 *
 * \param groupname Name of the new group
 * \param parent (default "/") Name of the parent group
 * \throw HistoContainerContentException
 */
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

/**
 * Create a new TH1 within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the histogram
 * \param title Title of the histogram
 * \param nbins number of bins
 * \param xmin min. value of the range
 * \param xmax max. value of the range
 * \param opt Additonal options (s for sumw2)
 */
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

/**
 * Create a new TH1 within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the histogram
 * \param title Title of the histogram
 * \param nbins number of bins
 * \param xbins array of bin limits
 * \param opt Additonal options (s for sumw2)
 */
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

/**
 * Create a new TH1 within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the histogram
 * \param title Title of the histogram
 * \param xbins array of bin limits (contains also number of bins)
 * \param opt Additonal options (s for sumw2)
 */
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

/**
 * Create a new TH2 within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the histogram
 * \param title Title of the histogram
 * \param nbinsx number of bins in x-direction
 * \param xmin min. value of the range in x-direction
 * \param xmax max. value of the range in x-direction
 * \param nbinsy number of bins in y-direction
 * \param ymin min. value of the range in y-direction
 * \param ymax max. value of the range in y-direction
 */
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

/**
 * Create a new TH2 within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the histogram
 * \param title Title of the histogram
 * \param nbinsx number of bins in x-direction
 * \param xbins array of bin limits in x-direction
 * \param nbinsy number of bins in y-direction
 * \param ybins array of bin limits in y-direction
 */
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

/**
 * Create a new TH2 within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the histogram
 * \param title Title of the histogram
 * \param xbins array of bin limits in x-direction (contains also the number of bins)
 * \param ybins array of bin limits in y-direction (contains also the number of bins)
 */
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

/**
 * Create a new TH3 within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param nbinsx number of bins in x-direction
 * \param xmin min. value of the range in x-direction
 * \param xmax max. value of the range in x-direction
 * \param nbinsy number of bins in y-direction
 * \param ymin min. value of the range in y-direction
 * \param ymax max. value of the range in y-direction
 * \param nbinsz number of bins in z-direction
 * \param zmin min. value of the range in z-direction
 * \param zmax max. value of the range in z-direction
 */
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

/**
 * Create a new TH3 within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the histogram
 * \param title Title of the histogram
 * \param nbinsx number of bins in x-direction
 * \param xbins array of bin limits in x-direction
 * \param nbinsy number of bins in y-direction
 * \param ybins array of bin limits in y-direction
 * \param nbinsz number of bins in z-direction
 * \param zbins array of bin limits in z-direction
 */
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

/**
 * Create a new TH3 within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the histogram
 * \param title Title of the histogram
 * \param xbins array of bin limits in x-direction (contains also the number of bins)
 * \param ybins array of bin limits in y-direction (contains also the number of bins)
 * \param zbins array of bin limits in z-direction (contains also the number of bins)
 */
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

/**
 * Create a new THnSparse within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the histogram
 * \param title Title of the histogram
 * \param ndim Number of dimensions
 * \param nbins Number of bins per dimension
 * \param min min. value of the range for each dimension
 * \param max max. value of the range for each dimension
 */
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

/**
 * Create a new THnSparse within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the histogram
 * \param title Title of the histogram
 * \param ndim Number of dimensions
 * \param axes Array of pointers to TAxis for containing the axis definition for each dimension
 */
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

/**
 * Create a new TProfile within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the profile histogram
 * \param title Title of the profile histogram
 * \param nbinsX Number of bins in x-direction
 * \param xmin min. value in x-direction
 * \param xmax max. value in x-direction
 * \param opt Further options
 */
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

/**
 * Create a new TProfile within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the profile histogram
 * \param title Title of the profile histogram
 * \param nbinsX Number of bins in x-direction
 * \param xmbins binning in x-direction
 * \param opt Further options
 */
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

/**
 * Create a new TProfile within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the profile histogram
 * \param title Title of the profile histogram
 * \param nbinsX Number of bins in x-direction
 * \param xmbins binning in x-direction
 * \param opt Further options
 */
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

/**
 * Set a new group into the container into the parent group
 *
 * \param o the object ot be included
 */
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

/**
 * Fill a 1D histogram within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the histogram
 * \param x x-coordinate
 * \param weight optional weight of the entry (default 1)
 */
void THistManager::FillTH1(const char *name, double x, double weight) {
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
	hist->Fill(x, weight);
}

void THistManager::FillTH1(const char *name, const char *label, double weight) {
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
  hist->Fill(label, weight);
}

/**
 * Fill a 2D histogram within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the histogram
 * \param x x-coordinate
 * \param y y-coordinate
 * \param weight optional weight of the entry (default 1)
 */
void THistManager::FillTH2(const char *name, double x, double y, double weight) {
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
	hist->Fill(x, y, weight);
}

/**
 * Fill a 2D histogram within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the histogram
 * \param point coordinates of the data
 * \param weight optional weight of the entry (default 1)
 */
void THistManager::FillTH2(const char *name, double *point, double weight) {
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
	hist->Fill(point[0], point[1], weight);
}

/**
 * Fill a 3D histogram within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the histogram
 * \param x x-coordinate
 * \param y y-coordinate
 * \param z z-coordinate
 * \param weight optional weight of the entry (default 1)
 */
void THistManager::FillTH3(const char* name, double x, double y, double z, double weight) {
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
	hist->Fill(x, y, z, weight);
}

/**
 * Fill a 3D histogram within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the histogram
 * \param point 3D-coordinate (x,y,z) of the point to be filled
 * \param weight optional weight of the entry (default 1)
 */
void THistManager::FillTH3(const char* name, const double* point, double weight) {
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
	hist->Fill(point[0], point[1], point[2], weight);
}

/**
 * Fill a  nD histogram within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the histogram
 * \param x coordinates of the data
 * \param weight optional weight of the entry (default 1)
 */
void THistManager::FillTHnSparse(const char *name, const double *x, double weight) {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent){
		Fatal("THistManager::FillTH3", "Parent group %s does not exist", dirname.Data());
		return;
	}
	THnSparseD *hist = dynamic_cast<THnSparseD *>(parent->FindObject(hname.Data()));
	if(!hist){
		Fatal("THistManager::FillTH3", "Histogram %s not found in parent group %s", hname.Data(), dirname.Data());
		return;
	}
	hist->Fill(x, weight);
}

/**
 * Fill a profile histogram within the container. The histogram name also contains the parent group(s) according to the common
 * group notation.
 *
 * \param name Name of the profile histogram
 * \param x x-coordinate
 * \param y y-coordinate
 * \param weight optional weight of the entry (default 1)
 */
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

/**
 * Find an object inside the container. The object can also be within a
 * histogram group. For this the name has to follow the common notation
 *
 * \param name Name of the object to find inside the container
 * \return pointer to the object (NULL if not found)
 */
TObject *THistManager::FindObject(const char *name) const {
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent) return NULL;
	return parent->FindObject(hname);
}

/**
 * Find and object inside the container. The object name is expected to contain the
 * full path of the histogram object, including parent groups
 *
 * \param obj the object to find
 * \return pointer to the object (NULL if not found)
 */
TObject* THistManager::FindObject(const TObject* obj) const {
	TString dirname(basename(obj->GetName())), hname(histname(obj->GetName()));
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent) return NULL;
	return parent->FindObject(hname);
}

/**
 * Find histogram group. Name is using common notation
 *
 * \param dirname Path of the group (treat empty path as top node
 * \return TList of objects (NULL if group does not exist)
 */
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

/**
 * Tokenizes a string. Results are stored inside the vector listoftokens
 *
 * \param name string to be tokenised
 * \param delim delimiter string
 * \param listoftokens list of tokens (C++ strings)
 */
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

/**
 * Helper function extracting the basename from a given histogram path.
 *
 * \param path histogram path
 * \return basename extracted
 */
TString THistManager::basename(const char *path) const {
	TString s(path);
	int index = s.Last('/');
	if(index < 0) return "";  // no directory structure
	return TString(s(0, index)).Data();
}

/**
 * Helper function extracting the histogram name from a given histogram path.
 *
 * \param path histogram path
 * \return basename extracted
 */
TString THistManager::histname(const char *path) const {
	TString s(path);
	int index = s.Last('/');
	if(index < 0) return path;    // no directory structure
	return TString(s(index+1, s.Length() - (index+1))).Data();
}

/**
 * Create forward iterator starting at the beginning of the
 * container
 * \return Forward iterator at the beginning of the container
 */
THistManager::iterator THistManager::begin() const {
  return iterator(this, 0, iterator::kTHMIforward);
}

/**
 * Create backward iterator starting at the end of the
 * container
 * \return Backward iterator after the end of the container
 */
THistManager::iterator THistManager::rbegin() const {
  return iterator(this, fHistos->GetEntries()-1, iterator::kTHMIbackward);
}

/**
 * Create forward iterator starting behind the end of the container.
 * Used to terminate the iteration.
 * \return Forward iterator behind the end of the histogram manager
 */
THistManager::iterator THistManager::end() const {
  return iterator(this, fHistos->GetEntries(), iterator::kTHMIforward);
}

/**
 * Create backward iterator starting before the beginning of the container.
 * Used to terminate the iteration.
 * \return Backward iterator starting before the beginning of the histogram manager
 */
THistManager::iterator THistManager::rend() const {
  return iterator(this, -1, iterator::kTHMIbackward);
}

//////////////////////////////////////////////////////////
///                                                    ///
/// Implementation of THistManager::iterator           ///
///                                                    ///
//////////////////////////////////////////////////////////

/**
 * Constructor. Initializing the histogram manager to be iterated over, the starting
 * position for the iteration, and the direction of the iteration.
 * \param[in] hmgr Histogram manager to be iterated over, containing the data
 * \param[in] currentpos Starting position for the iteration
 * \param[in] dir Direction of the iteration
 */
THistManager::iterator::iterator(const THistManager * hmgr, Int_t currentpos, THMIDirection_t dir):
    fkArray(hmgr),
    fCurrentPos(),
    fNext(),
    fDirection(dir)
{}

/**
 * Copy constructor. Initializing all
 * values from the reference.
 * \param[in] ref Reference for the copy
 */
THistManager::iterator::iterator(const THistManager::iterator &ref):
    fkArray(ref.fkArray),
    fCurrentPos(ref.fCurrentPos),
    fNext(ref.fNext),
    fDirection(ref.fDirection)
{}

/**
 * Assignment operator. Initializing all
 * values from the reference.
 * \param[in] ref
 * \return iterator after assignment.
 */
THistManager::iterator &THistManager::iterator::operator=(const THistManager::iterator &ref){
  if(this != &ref){
    fkArray = ref.fkArray;
    fCurrentPos = ref.fCurrentPos;
    fNext = ref.fNext;
    fDirection = ref.fDirection;
  }
  return *this;
}

/**
 * Comparison operator for unequalness. The comparison is
 * performed based on the current position of the iterator.
 * \param other iterator to compare to
 * \return True if the iterators are not equal
 */
bool THistManager::iterator::operator!=(const THistManager::iterator &other) const{
  return fCurrentPos == other.fCurrentPos;
}

/**
 * Prefix increment operator. Incrementing / decrementing
 * the current position of the iterator based on the
 * direction.
 * \return State after incrementation
 */
THistManager::iterator &THistManager::iterator::operator++(){
  if(fDirection == kTHMIforward)
    fCurrentPos++;
  else
    fCurrentPos--;
  return *this;
}

/**
 * Postfix increment operator. Incrementing / decrementing
 * the current position of the iterator based on the
 * direction, but returning the state before iteration
 * \return State before incrementation
 */
THistManager::iterator THistManager::iterator::operator++(int){
  iterator tmp(*this);
  operator++();
  return tmp;
}

/**
 * Prefix decrement operator. Decrementing / incrementing
 * the current position of the iterator based on the
 * direction.
 * \return State after decrementation
 */
THistManager::iterator &THistManager::iterator::operator--(){
  if(fDirection == kTHMIforward)
    fCurrentPos--;
  else
    fCurrentPos++;
  return *this;
};

/**
 * Postfix decrement operator. Decrementing / incrementing
 * the current position of the iterator based on the
 * direction, but returning the state before iteration
 * \return State before decrementation
 */
THistManager::iterator THistManager::iterator::operator--(int){
  iterator tmp(*this);
  operator--();
  return tmp;
};

/**
 * Dereferncing operator. Providing access to the object
 * connected to the current state of the iterator. Objects
 * might be histograms or THashLists in case of wheter groups
 * are defined or not.
 * @return Object connected ot the current state of the iterator
 */
TObject *THistManager::iterator::operator*() const{
  if(fCurrentPos >=0 && fCurrentPos < fkArray->GetListOfHistograms()->GetEntries())
    return fkArray->GetListOfHistograms()->At(fCurrentPos);
  return NULL;
}
