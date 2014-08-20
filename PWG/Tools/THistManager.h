#ifndef THISTMANAGER_H
#define THISTMANAGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include <cstring>
#include <exception>
#include <string>
#include <sstream>
#include <TNamed.h>

class TArrayD;
class TAxis;
class TList;
class THashList;

class HistManagerException: public std::exception {
	/*
	 * Error handling class for the histogram container
	 */
public:
	enum ExceptionType_t {
		kHistNotFoundException = 0,
		kTypeException = 1,
		kHistDuplicationException = 2,
		kGroupException = 3
	};

	HistManagerException(const char *histname, const char *hgroup, ExceptionType_t etype):
		fHistname(),
		fGroup(),
		fExceptionType(etype)
	{
		if(histname) fHistname = histname;
		if(hgroup) fGroup = hgroup;
	}
	virtual ~HistManagerException() throw() {}

	virtual const char *what() const throw() {
		std::stringstream msgbuilder;
		switch(fExceptionType) {
		case kHistNotFoundException:
			msgbuilder << "Histogram " << fHistname << " not found in";
			if(strlen(fGroup.c_str())) msgbuilder << " group " << fGroup;
			else msgbuilder << " the list of histograms.";
			break;
		case kTypeException:
			msgbuilder << "Object " << fHistname << " is of wrong type.";
			break;
		case kHistDuplicationException:
			msgbuilder << "Histogram " << fHistname << " already exists in";
			if(strlen(fGroup.c_str())) msgbuilder << " group " << fGroup;
			else msgbuilder << " the list of histograms.";
			break;
		case kGroupException:
			msgbuilder << "Group " << fGroup << " not found.";
			break;
		};
		return msgbuilder.str().c_str();
	}

	const char * GetErrorHistogramName() const { return fHistname.c_str(); }
	ExceptionType_t GetExceptionType() const { return fExceptionType; }

private:
	std::string           fHistname;            // Name of the histogram producing the exception
	std::string           fGroup;               // Group of objects producing the exception
	ExceptionType_t       fExceptionType;       // type of the exception

};

class THistManager : public TNamed {
public:
	THistManager();
	THistManager(const char *name);
	~THistManager();
	void ReleaseOwner() { fIsOwner = kFALSE; };

	void CreateHistoGroup(const char *groupname, const char *parent = "/") throw(HistManagerException);

	void CreateTH1(const char *name, const char *title, int nbins, double xmin, double xmax) throw(HistManagerException);
	void CreateTH1(const char *name, const char *title, int nbins, const double *xbins) throw(HistManagerException);
	void CreateTH1(const char *name, const char *title, const TArrayD &xbins) throw(HistManagerException);
	void CreateTH2(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax) throw(HistManagerException);
	void CreateTH2(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy, const double *ybins) throw(HistManagerException);
	void CreateTH2(const char *name, const char *title, const TArrayD &xbins, const TArrayD &ybins) throw(HistManagerException);
	void CreateTH3(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, int nbinsz, double zmin, double zmax) throw (HistManagerException);
	void CreateTH3(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy, const double *ybins, int nbinsz, const double *zbins) throw (HistManagerException);
	void CreateTH3(const char *name, const char *title, const TArrayD &xbins, const TArrayD &ybins, const TArrayD &zbins) throw(HistManagerException);
	void CreateTHnSparse(const char *name, const char *title, int ndim, const int *nbins, const double *min, const double *max) throw(HistManagerException);
	void CreateTHnSparse(const char *name, const char *title, int ndim, const TAxis **axes) throw(HistManagerException);
	void SetObject(TObject * const o, const char *group = "/") throw(HistManagerException);
	void FillTH1(const char *hname, double x, double weight = 1.) throw(HistManagerException);
	void FillTH2(const char *hname, double x, double y, double weight = 1.) throw(HistManagerException);
	void FillTH2(const char *hname, double *point, double weight = 1.) throw(HistManagerException);
	void FillTH3(const char *hname, double x, double y, double z, double weight = 1.) throw(HistManagerException);
	void FillTH3(const char *hname, const double *point, double weight = 1.) throw(HistManagerException);
	void FillTHnSparse(const char *name, const double *x, double weight = 1.) throw(HistManagerException);

	THashList *GetListOfHistograms() { return fHistos; }
	TObject *FindObject(const char *name) const;

private:
	THistManager(const THistManager &);
	THistManager &operator=(const THistManager &);
	THashList *FindGroup(const char *dirname) const;
	void TokenizeFilename(const char *name, const char *delim, std::vector<std::string> &listoftokens) const;
	const char *basename(const char *path) const;
	const char *histname(const char *path) const;

	THashList *fHistos;                   // List of histograms
	bool fIsOwner;                        // Set the ownership

	ClassDef(THistManager, 1);  // Container for histograms
};

#endif
