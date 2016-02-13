#ifndef THISTMANAGER_H
#define THISTMANAGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include <TNamed.h>

class TArrayD;
class TAxis;
class TList;
class THashList;
class TH1;
class TH2;
class TH3;
class THnSparse;

class THistManager : public TNamed {
public:
	THistManager();
	THistManager(const char *name);
	~THistManager();
	void ReleaseOwner() { fIsOwner = kFALSE; };

	THashList* CreateHistoGroup(const char *groupname, const char *parent = "/");

	TH1* CreateTH1(const char *name, const char *title, int nbins, double xmin, double xmax);
	TH1* CreateTH1(const char *name, const char *title, int nbins, const double *xbins);
	TH1* CreateTH1(const char *name, const char *title, const TArrayD &xbins);
	TH2* CreateTH2(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax);
	TH2* CreateTH2(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy, const double *ybins);
	TH2* CreateTH2(const char *name, const char *title, const TArrayD &xbins, const TArrayD &ybins);
	TH3* CreateTH3(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, int nbinsz, double zmin, double zmax);
	TH3* CreateTH3(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy, const double *ybins, int nbinsz, const double *zbins);
	TH3* CreateTH3(const char *name, const char *title, const TArrayD &xbins, const TArrayD &ybins, const TArrayD &zbins);
	THnSparse* CreateTHnSparse(const char *name, const char *title, int ndim, const int *nbins, const double *min, const double *max);
	THnSparse* CreateTHnSparse(const char *name, const char *title, int ndim, const TAxis **axes);
	void SetObject(TObject * const o, const char *group = "/");
	void FillTH1(const char *hname, double x, double weight = 1.);
	void FillTH2(const char *hname, double x, double y, double weight = 1.);
	void FillTH2(const char *hname, double *point, double weight = 1.);
	void FillTH3(const char *hname, double x, double y, double z, double weight = 1.);
	void FillTH3(const char *hname, const double *point, double weight = 1.);
	void FillTHnSparse(const char *name, const double *x, double weight = 1.);

	THashList *GetListOfHistograms() { return fHistos; }
	TObject *FindObject(const char *name) const;
	virtual TObject *FindObject(const TObject *obj) const;

private:
	THistManager(const THistManager &);
	THistManager &operator=(const THistManager &);
	THashList *FindGroup(const char *dirname) const;
	void TokenizeFilename(const char *name, const char *delim, std::vector<std::string> &listoftokens) const;
	TString basename(const char *path) const;
	TString histname(const char *path) const;

	THashList *fHistos;                   // List of histograms
	bool fIsOwner;                        // Set the ownership

	ClassDef(THistManager, 1);  // Container for histograms
};

#endif
