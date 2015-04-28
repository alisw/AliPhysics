/**
 * \file AliEMCalHistoContainer.h
 * \brief Declarartion of class AliEMCalHistoContainer
 *
 * In this file the class AliEMCalHistoContainer and the corresponding exception class,
 * HistoContainerContentException, are declared. The histogram container is a tool storing
 * histograms of different types and providing easy access via their names.
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Labortatory
 * \date Aug 5, 2014
 */
#ifndef ALIEMCALHISTOCONTAINER_H
#define ALIEMCALHISTOCONTAINER_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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

/**
 * \class HistoContainerContentException
 * \brief Exception thrown by the histogram container in case of problems
 *
 * HistoContainerContentException is thrown by the histogram container in case
 * of problems appearing when accessing or filling histograms. Problems handled by
 * this class are:
 * - Required histogram not defined
 * - Histogram types not matching
 * - Histogram would be duplicated
 * - Group not existing or duplicated when creating
 *
 * Error handling class for the histogram container
 */
class HistoContainerContentException : public std::exception {
public:
  /**
   * \enum ExceptionType_t
   * \brief Definition of exception types thrown by the histogram container
   *
   * This enumeration defines possible exeption types handled by the exception class
   */
  enum ExceptionType_t {
    kHistNotFoundException = 0,   //!< Histogram with name not found in the container
    kTypeException = 1,           //!< Histogram type mismatch
    kHistDuplicationException = 2,//!< Histogram with name duplicated
    kGroupException = 3           //!< Group error (not existing or duplicated)
  };

  /**
   * Constuctor, defining the exception. Called when a HistoContainerContentException is thrown
   * \param histname Name of the histogram throwing the exception
   * \param hgroup Group throwing the exception
   * \param etype
   */
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

  /**
   * Destructor
   */
  virtual ~HistoContainerContentException() throw() {}

  /**
   * Get error message associated with the histogram
   * \return
   */
  virtual const char *what() const throw() {
    return fErrorMessage.c_str();
  }

  /**
   * Get the name of the histogram raising the exception
   * \return Name of the histogram
   */
  const char * GetErrorHistogramName() const { return fHistname.c_str(); }
  /**
   * Get the type of the exception
   * \return Type of the exception
   */
  ExceptionType_t GetExceptionType() const { return fExceptionType; }

private:
  /**
   * Create error message with the histogram name, the histogram group, and the error type
   */
  void CreateErrorMessage(){
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

  std::string           fHistname;            ///< Name of the histogram producing the exception
  std::string           fGroup;               ///< Group of objects producing the exception
  std::string			      fErrorMessage;		    ///< container for the error message produced in the what function
  ExceptionType_t       fExceptionType;       ///< type of the exception

};

/**
 * \class AliEMCalHistoContainer
 * \brief Container class for histograms for the high-\f$ p_{t} \f$ charged particle analysis
 *
 * Container class for histogram objects. Currenly can handle
 *   TH1
 *   TH2
 *   TH3
 *   THnSparse
 * Histograms can be stored in groups. For this the parent group is
 * included inside the histogram name, i.e. /base/inheriting/histogram.
 * In case just the histogram name is given, it is assumed that the
 * histogram is stored at the top level.
 */
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
  void CreateTProfile(const char *name, const char *title, int nbinsX, double xmin, double xmax, Option_t *opt = "") throw(HistoContainerContentException);
  void CreateTProfile(const char *name, const char *title, int nbinsX, const double *xbins, Option_t *opt = "") throw(HistoContainerContentException);
  void CreateTProfile(const char *name, const char *title, const TArrayD &xbins, Option_t *opt = "") throw(HistoContainerContentException);
  void SetObject(TObject * const o, const char *group = "/") throw(HistoContainerContentException);
  void FillTH1(const char *hname, double x, double weight = 1.) throw(HistoContainerContentException);
  void FillTH2(const char *hname, double x, double y, double weight = 1.) throw(HistoContainerContentException);
  void FillTH2(const char *hname, double *point, double weight = 1.) throw(HistoContainerContentException);
  void FillTH3(const char *hname, double x, double y, double z, double weight = 1.) throw(HistoContainerContentException);
  void FillTH3(const char *hname, const double *point, double weight = 1.) throw(HistoContainerContentException);
  void FillTHnSparse(const char *name, const double *x, double weight = 1.) throw(HistoContainerContentException);
  void FillProfile(const char *name, double x, double y, double weight = 1.) throw(HistoContainerContentException);

  /**
   * Get the list of histograms
   * \return The list of histograms
   */
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

  THashList *fHistos;                   ///< List of histograms
  bool fIsOwner;                        ///< Set the ownership

  /// \cond CLASSIMP
  ClassDef(AliEMCalHistoContainer, 1);  // Container for histograms
  /// \endcond
};

}
#endif
