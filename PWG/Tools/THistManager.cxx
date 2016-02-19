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

/*
 * Container class for histogram objects. Currenly can handle
 *   TH1
 *   TH2
 *   TH3
 *   THnSparse
 * Histograms can be stored in groups. For this the parent group is 
 * included inside the histogram name, i.e. /base/inheriting/histogram.
 * In case just the histogram name is given, it is assumed that the
 * histogram is stored at the top level.
 *
 *   Author: Markus Fasel
 */

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
#include <TString.h>

#include "THistManager.h"

ClassImp(THistManager)

//______________________________________________________________________________
THistManager::THistManager():
		TNamed(),
		fHistos(NULL),
		fIsOwner(true)
{
	/*
	 * Default constructor, only initialising pointers with 0
	 */
}

//______________________________________________________________________________
THistManager::THistManager(const char *name):
		TNamed(name, Form("Histogram container %s", name)),
		fHistos(NULL),
		fIsOwner(true)
{
	/*
	 * Main constructor, creating also a list for the histograms
	 *
	 * @param name: Name of the object (list named accordingly)
	 */
	fHistos = new THashList();
	fHistos->SetName(Form("histos%s", name));
	fHistos->SetOwner();
}

//______________________________________________________________________________
THistManager::~THistManager(){
	/*
	 * Destructor, deletes the list of histograms if it is the owner
	 */
	if(fHistos && fIsOwner) delete fHistos;
}

//______________________________________________________________________________
THashList* THistManager::CreateHistoGroup(const char *groupname, const char *parent) {
	/*
	 * Create a new group of histograms within a parent group. Groups are represented as list. The default parent is
	 * always the top list. List name structure accouding to unix paths (i.e. top list /, hirarchies separated by /).
	 *
	 * @param groupname: Name of the new group
	 * @param parent (@default "/"): Name of the parent group
	 * @throw HistManagerException
	 */
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

//______________________________________________________________________________
TH1* THistManager::CreateTH1(const char *name, const char *title, int nbins, double xmin, double xmax){
	/*
	 * Create a new TH1 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param name: Name of the histogram
	 * @param title: Title of the histogram
	 * @param nbins: number of bins
	 * @param xmin: min. value of the range
	 * @param xmax: max. value of the range
	 * Raises fatals in case the parent group does not exist or the object is attempted to be duplicated within the group
	 */
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
	parent->Add(h);
	return h;
}

//______________________________________________________________________________
TH1* THistManager::CreateTH1(const char *name, const char *title, int nbins, const double *xbins){
	/*
	 * Create a new TH1 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param name: Name of the histogram
	 * @param title: Title of the histogram
	 * @param nbins: number of bins
	 * @param xbins: array of bin limits
	 * Raises fatals in case the parent group does not exist or the object is attempted to be duplicated within the group
	 */
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
	parent->Add(h);
	return h;
}

//______________________________________________________________________________
TH1* THistManager::CreateTH1(const char *name, const char *title, const TArrayD &xbins){
	/*
	 * Create a new TH1 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param name: Name of the histogram
	 * @param title: Title of the histogram
	 * @param xbins: array of bin limits (contains also number of bins)
	 * Raises fatals in case the parent group does not exist or the object is attempted to be duplicated within the group
	 */
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
	parent->Add(h);
	return h;
}

//______________________________________________________________________________
TH2* THistManager::CreateTH2(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax){
	/*
	 * Create a new TH2 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param name: Name of the histogram
	 * @param title: Title of the histogram
	 * @param nbinsx: number of bins in x-direction
	 * @param xmin: min. value of the range in x-direction
	 * @param xmax: max. value of the range in x-direction
	 * @param nbinsy: number of bins in y-direction
	 * @param ymin: min. value of the range in y-direction
	 * @param ymax: max. value of the range in y-direction
	 * Raises fatals in case the parent group does not exist or the object is attempted to be duplicated within the group
	 */
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
	parent->Add(h);
	return h;
}

//______________________________________________________________________________
TH2* THistManager::CreateTH2(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy, const double *ybins){
	/*
	 * Create a new TH2 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param name: Name of the histogram
	 * @param title: Title of the histogram
	 * @param nbinsx: number of bins in x-direction
	 * @param xbins: array of bin limits in x-direction
	 * @param nbinsy: number of bins in y-direction
	 * @param ybins: array of bin limits in y-direction
	 * Raises fatals in case the parent group does not exist or the object is attempted to be duplicated within the group
	 */
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
	parent->Add(h);
	return h;
}

//______________________________________________________________________________
TH2* THistManager::CreateTH2(const char *name, const char *title, const TArrayD &xbins, const TArrayD &ybins){
	/*
	 * Create a new TH2 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param name: Name of the histogram
	 * @param title: Title of the histogram
	 * @param xbins: array of bin limits in x-direction (contains also the number of bins)
	 * @param ybins: array of bin limits in y-direction (contains also the number of bins)
	 * Raises fatals in case the parent group does not exist or the object is attempted to be duplicated within the group
	 */
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
	parent->Add(h);
	return h;
}

//______________________________________________________________________________
TH3* THistManager::CreateTH3(const char* name, const char* title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, int nbinsz, double zmin, double zmax) {
	/*
	 * Create a new TH3 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param nbinsx: number of bins in x-direction
	 * @param xmin: min. value of the range in x-direction
	 * @param xmax: max. value of the range in x-direction
	 * @param nbinsy: number of bins in y-direction
	 * @param ymin: min. value of the range in y-direction
	 * @param ymax: max. value of the range in y-direction
	 * @param nbinsz: number of bins in z-direction
	 * @param zmin: min. value of the range in z-direction
	 * @param zmax: max. value of the range in z-direction
	 * Raises fatals in case the parent group does not exist or the object is attempted to be duplicated within the group
	 */
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
	parent->Add(h);
	return h;
}

//______________________________________________________________________________
TH3* THistManager::CreateTH3(const char* name, const char* title, int nbinsx, const double* xbins, int nbinsy, const double* ybins, int nbinsz, const double* zbins) {
	/*
	 * Create a new TH3 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param name: Name of the histogram
	 * @param title: Title of the histogram
	 * @param nbinsx: number of bins in x-direction
	 * @param xbins: array of bin limits in x-direction
	 * @param nbinsy: number of bins in y-direction
	 * @param ybins: array of bin limits in y-direction
	 * @param nbinsz: number of bins in z-direction
	 * @param zbins: array of bin limits in z-direction
	 * Raises fatals in case the parent group does not exist or the object is attempted to be duplicated within the group
	 */
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
	parent->Add(h);
	return h;
}

//______________________________________________________________________________
TH3* THistManager::CreateTH3(const char* name, const char* title, const TArrayD& xbins, const TArrayD& ybins, const TArrayD& zbins) {
	/*
	 * Create a new TH3 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param name: Name of the histogram
	 * @param title: Title of the histogram
	 * @param xbins: array of bin limits in x-direction (contains also the number of bins)
	 * @param ybins: array of bin limits in y-direction (contains also the number of bins)
	 * @param zbins: array of bin limits in z-direction (contains also the number of bins)
	 * Raises fatals in case the parent group does not exist or the object is attempted to be duplicated within the group
	 */
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
	parent->Add(h);
	return h;
}

//______________________________________________________________________________
THnSparse* THistManager::CreateTHnSparse(const char *name, const char *title, int ndim, const int *nbins, const double *min, const double *max) {
	/*
	 * Create a new THnSparse within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param name: Name of the histogram
	 * @param title: Title of the histogram
	 * @param ndim: Number of dimensions
	 * @param nbins: Number of bins per dimension
	 * @param min: min. value of the range for each dimension
	 * @param max: max. value of the range for each dimension
	 * Raises fatals in case the parent group does not exist or the object is attempted to be duplicated within the group
	 */
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
	parent->Add(h);
	return h;
}

//______________________________________________________________________________
THnSparse* THistManager::CreateTHnSparse(const char *name, const char *title, int ndim, const TAxis **axes) {
	/*
	 * Create a new THnSparse within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param name: Name of the histogram
	 * @param title: Title of the histogram
	 * @param ndim: Number of dimensions
	 * @param axes: Array of pointers to TAxis for containing the axis definition for each dimension
	 * Raises fatals in case the parent group does not exist or the object is attempted to be duplicated within the group
	 */
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
	parent->Add(hsparse);
	return hsparse;
}

//______________________________________________________________________________
void THistManager::SetObject(TObject * const o, const char *group) {
	/*
	 * Set a new group into the container into the parent group
	 *
	 * @param o: the object ot be included
	 * Raises fatals in case the parent group is not found, the object is attempted to be duplicated in the group, or the type of
	 * the object is not a histogram type
	 */
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

//______________________________________________________________________________
void THistManager::FillTH1(const char *name, double x, double weight) {
	/*
	 * Fill a 1D histogram within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param name: Name of the histogram
	 * @param x: x-coordinate
	 * @param weight (@default 1): optional weight of the entry
	 * Raises fatals in case the parent group is not found or the histogram is not found in the parent group
	 */
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

//______________________________________________________________________________
void THistManager::FillTH1(const char *name, const char *label, double weight) {
  /*
   * Fill a 1D histogram within the container. The histogram name also contains the parent group(s) according to the common
   * group notation.
   *
   * @param name: Name of the histogram
   * @param label: bin label
   * @param weight (@default 1): optional weight of the entry
   * Raises fatals in case the parent group is not found or the histogram is not found in the parent group
   */
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

//______________________________________________________________________________
void THistManager::FillTH2(const char *name, double x, double y, double weight) {
	/*
	 * Fill a 2D histogram within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param name: Name of the histogram
	 * @param x: x-coordinate
	 * @param y: y-coordinate
	 * @param weight (@default 1): optional weight of the entry
	 * Raises fatals in case the parent group is not found or the histogram is not found in the parent group
	 */
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

//______________________________________________________________________________
void THistManager::FillTH2(const char *name, double *point, double weight) {
	/*
	 * Fill a 2D histogram within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param name: Name of the histogram
	 * @param point: coordinates of the data
	 * @param weight (@default 1): optional weight of the entry
	 * Raises fatals in case the parent group is not found or the histogram is not found in the parent group
	 */
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

//______________________________________________________________________________
void THistManager::FillTH3(const char* name, double x, double y, double z, double weight) {
	/*
	 * Fill a 3D histogram within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param name: Name of the histogram
	 * @param x: x-coordinate
	 * @param y: y-coordinate
	 * @param z: z-coordinate
	 * @param weight (@default 1): optional weight of the entry
	 * Raises fatals in case the parent group is not found or the histogram is not found in the parent group
	 */
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

//______________________________________________________________________________
void THistManager::FillTH3(const char* name, const double* point, double weight) {
	/*
	 * Fill a 3D histogram within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param name: Name of the histogram
	 * @param point: 3D-coordinate of the point
	 * @param weight (@default 1): optional weight of the entry
	 * Raises fatals in case the parent group is not found or the histogram is not found in the parent group
	 */
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


//______________________________________________________________________________
void THistManager::FillTHnSparse(const char *name, const double *x, double weight) {
	/*
	 * Fill a  nD histogram within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 *
	 * @param name: Name of the histogram
	 * @param x: coordinates of the data
	 * @param weight (@default 1): optional weight of the entry
	 * Raises fatals in case the parent group is not found or the histogram is not found in the parent group
	 */
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

//______________________________________________________________________________
TObject *THistManager::FindObject(const char *name) const {
	/*
	 * Find an object inside the container. The object can also be within a
	 * histogram group. For this the name has to follow the common notation
	 *
	 * @param name: Name of the object to find inside the container
	 * @return: pointer to the object (NULL if not found)
	 */
	TString dirname(basename(name)), hname(histname(name));
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent) return NULL;
	return parent->FindObject(hname);
}

//______________________________________________________________________________
TObject* THistManager::FindObject(const TObject* obj) const {
	/*
	 * Find and object inside the container. The object name is expected to contain the
	 * full path of the histogram object, including parent groups
	 *
	 * @param obj: the object to find
	 * @return: pointer to the object (NULL if not found)
	 */
	TString dirname(basename(obj->GetName())), hname(histname(obj->GetName()));
	THashList *parent(FindGroup(dirname.Data()));
	if(!parent) return NULL;
	return parent->FindObject(hname);
}

//______________________________________________________________________________
THashList *THistManager::FindGroup(const char *dirname) const {
	/*
	 * Find histogram group. Name is using common notation
	 *
	 * @param dirname: Path of the group (treat empty path as top node
	 * @return: TList of objects (NULL if group does not exist)
	 */
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

//______________________________________________________________________________
void THistManager::TokenizeFilename(const char *name, const char *delim, std::vector<std::string> &listoftokens) const {
	/*
	 * Tokenizes a string. Results are stored inside the vector listoftokens
	 *
	 * @ param name: string to be tokenised
	 * @ param delim: delimiter string
	 * @ param listoftokens: list of tokens (C++ strings)
	 */
	TString s(name);
	TObjArray *arr = s.Tokenize(delim);
	TObjString *ostr(NULL);
	TIter toks(arr);
	while((ostr = dynamic_cast<TObjString *>(toks()))){
		listoftokens.push_back(std::string(ostr->String().Data()));
	}
	delete arr;
}

//______________________________________________________________________________
TString THistManager::basename(const char *path) const {
	/*
	 * Helper function extracting the basename from a given histogram path.
	 *
	 * @param path: histogram path
	 * @return: basename extracted
	 */
	TString s(path);
	int index = s.Last('/');
	if(index < 0) return "";  // no directory structure
	return TString(s(0, index)).Data();
}

//______________________________________________________________________________
TString THistManager::histname(const char *path) const {
	/*
	 * Helper function extracting the histogram name from a given histogram path.
	 *
	 * @param path: histogram path
	 * @return: basename extracted
	 */
	TString s(path);
	int index = s.Last('/');
	if(index < 0) return path;    // no directory structure
	return TString(s(index+1, s.Length() - (index+1))).Data();
}
