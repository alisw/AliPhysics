#ifndef ALIEMCALHISTOCONTAINER_H
#define ALIEMCALHISTOCONTAINER_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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

/**
 * \namespace EMCalTriggerPtAnalysis
 * \brief Analysis of high-\f$ p_{t} \f$ tracks in triggered events
 *
 * This namespace contains classes for the analysis of high-\f$ p_{t} \f$ tracks in
 * triggered events.
 */
namespace EMCalTriggerPtAnalysis {

class HistoContainerContentException : public std::exception {
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

  HistoContainerContentException(const char *histname, const char *hgroup, ExceptionType_t etype):
    fHistname(),
    fGroup(),
    fErrorMessage(),
    fExceptionType(etype)
  {
    if(histname) fHistname = histname;
    if(hgroup) fGroup = hgroup;

    CreateErrorMessage();
  }
  virtual ~HistoContainerContentException() throw() {}

  virtual const char *what() const throw() {
    return fErrorMessage.c_str();
  }

  const char * GetErrorHistogramName() const { return fHistname.c_str(); }
  ExceptionType_t GetExceptionType() const { return fExceptionType; }

private:
  void CreateErrorMessage(){
    /*
     * Create error message with the histogram name, the histogram group, and the error type
     */
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
    fErrorMessage = msgbuilder.str();
  }

  std::string           fHistname;            // Name of the histogram producing the exception
  std::string           fGroup;               // Group of objects producing the exception
  std::string			  fErrorMessage;		// container for the error message produced in the what function
  ExceptionType_t       fExceptionType;       // type of the exception

};

class AliEMCalHistoContainer : public TNamed {
public:
  AliEMCalHistoContainer();
  AliEMCalHistoContainer(const char *name);
  ~AliEMCalHistoContainer();
  void ReleaseOwner() { fIsOwner = kFALSE; };

  void CreateHistoGroup(const char *groupname, const char *parent = "/") throw(HistoContainerContentException);

  void CreateTH1(const char *name, const char *title, int nbins, double xmin, double xmax, Option_t *opt = "") throw(HistoContainerContentException);
  void CreateTH1(const char *name, const char *title, int nbins, const double *xbins, Option_t *opt = "") throw(HistoContainerContentException);
  void CreateTH1(const char *name, const char *title, const TArrayD &xbins, Option_t *opt = "") throw(HistoContainerContentException);
  void CreateTH2(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, Option_t *opt = "") throw(HistoContainerContentException);
  void CreateTH2(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy, const double *ybins, Option_t *opt = "") throw(HistoContainerContentException);
  void CreateTH2(const char *name, const char *title, const TArrayD &xbins, const TArrayD &ybins, Option_t *opt = "") throw(HistoContainerContentException);
  void CreateTH3(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, int nbinsz, double zmin, double zmax, Option_t *opt = "") throw (HistoContainerContentException);
  void CreateTH3(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy, const double *ybins, int nbinsz, const double *zbins, Option_t *opt = "") throw (HistoContainerContentException);
  void CreateTH3(const char *name, const char *title, const TArrayD &xbins, const TArrayD &ybins, const TArrayD &zbins, Option_t *opt = "") throw(HistoContainerContentException);
  void CreateTHnSparse(const char *name, const char *title, int ndim, const int *nbins, const double *min, const double *max, Option_t *opt = "") throw(HistoContainerContentException);
  void CreateTHnSparse(const char *name, const char *title, int ndim, const TAxis **axes, Option_t *opt = "") throw(HistoContainerContentException);
  void SetObject(TObject * const o, const char *group = "/") throw(HistoContainerContentException);
  void FillTH1(const char *hname, double x, double weight = 1.) throw(HistoContainerContentException);
  void FillTH2(const char *hname, double x, double y, double weight = 1.) throw(HistoContainerContentException);
  void FillTH2(const char *hname, double *point, double weight = 1.) throw(HistoContainerContentException);
  void FillTH3(const char *hname, double x, double y, double z, double weight = 1.) throw(HistoContainerContentException);
  void FillTH3(const char *hname, const double *point, double weight = 1.) throw(HistoContainerContentException);
  void FillTHnSparse(const char *name, const double *x, double weight = 1.) throw(HistoContainerContentException);

  THashList *GetListOfHistograms() { return fHistos; }
  virtual TObject *FindObject(const char *name) const;
  virtual TObject *FindObject(const TObject *obj) const;

private:
  AliEMCalHistoContainer(const AliEMCalHistoContainer &);
  AliEMCalHistoContainer &operator=(const AliEMCalHistoContainer &);
  THashList *FindGroup(const char *dirname) const;
  void TokenizeFilename(const char *name, const char *delim, std::vector<std::string> &listoftokens) const;
  const char *basename(const char *path) const;
  const char *histname(const char *path) const;

  THashList *fHistos;                   // List of histograms
  bool fIsOwner;                        // Set the ownership

  ClassDef(AliEMCalHistoContainer, 1);  // Container for histograms
};

}
#endif
