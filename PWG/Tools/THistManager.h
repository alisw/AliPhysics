#ifndef THISTMANAGER_H
#define THISTMANAGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <THashList.h>
#include <TIterator.h>
#include <TNamed.h>
#include <iterator>

class TArrayD;
class TAxis;
class TList;
class TH1;
class TH2;
class TH3;
class THnSparse;
class TProfile;

/**
 * \class THistManager
 * \brief Container class for histograms for the high-\f$ p_{t} \f$ charged particle analysis
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 *
 * Container class for histogram objects. Currently can handle
 *   TH1
 *   TH2
 *   TH3
 *   THnSparse
 *   TProfile
 * Histograms can be stored in groups. For this the parent group is
 * included inside the histogram name, i.e. /base/inheriting/histogram.
 * In case just the histogram name is given, it is assumed that the
 * histogram is stored at the top level.
 */
class THistManager : public TNamed {
public:
  /**
   * \class iterator
   * \brief stl-iterator for the histogram manager
   *
   * stl-type iterator for the histogram manager. Iterating
   * over primary content of the histogram manager. In case
   * histograms are organized in groups, the primary content
   * will be histogram groups, with the data structure
   * THashList. The iterator is implemented as bidirectional
   * iterator, providing forward and backward iteration.
   */
  class iterator :  public std::iterator<std::bidirectional_iterator_tag,
                                                 TObject, std::ptrdiff_t,
                                                 TObject **, TObject *&>{
  public:
    /**
     * \enum THMIDirection_t
     * \brief Direction for the iteration
     *
     * Switch for the firection of the iteration.
     * Currently the iterator is implemented as
     * bi-directional iterator, providing only
     * forward and backward iteration.
     */
    enum THMIDirection_t {
      kTHMIforward = 0,//!< Forward iteration
      kTHMIbackward = 1//!< Backward iteration
    };
    iterator(const THistManager *hmgr, Int_t currentpos, THMIDirection_t dir = kTHMIforward);
    iterator(const iterator &ref);
    ~iterator() { }

    iterator              &operator=(const iterator &rhs);
    Bool_t                operator!=(const iterator &aIter) const;
    iterator              &operator++();
    iterator              operator++(int);
    iterator              &operator--();
    iterator              operator--(int);
    TObject               *operator*() const;

  private:
    const THistManager          *fkArray;
    Int_t                       fCurrentPos;
    Int_t                       fNext;
    THMIDirection_t             fDirection;

    iterator();
  };

	THistManager();
	THistManager(const char *name);
	~THistManager();
	void ReleaseOwner() { fIsOwner = kFALSE; };

	THashList* CreateHistoGroup(const char *groupname, const char *parent = "/");

	TH1* CreateTH1(const char *name, const char *title, int nbins, double xmin, double xmax, Option_t *opt = "");
	TH1* CreateTH1(const char *name, const char *title, int nbins, const double *xbins, Option_t *opt = "");
	TH1* CreateTH1(const char *name, const char *title, const TArrayD &xbins, Option_t *opt = "");
	TH2* CreateTH2(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, Option_t *opt = "");
	TH2* CreateTH2(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy, const double *ybins, Option_t *opt = "");
	TH2* CreateTH2(const char *name, const char *title, const TArrayD &xbins, const TArrayD &ybins, Option_t *opt = "");
	TH3* CreateTH3(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, int nbinsz, double zmin, double zmax, Option_t *opt = "");
	TH3* CreateTH3(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy, const double *ybins, int nbinsz, const double *zbins, Option_t *opt = "");
	TH3* CreateTH3(const char *name, const char *title, const TArrayD &xbins, const TArrayD &ybins, const TArrayD &zbins, Option_t *opt = "");
	THnSparse* CreateTHnSparse(const char *name, const char *title, int ndim, const int *nbins, const double *min, const double *max, Option_t *opt = "");
	THnSparse* CreateTHnSparse(const char *name, const char *title, int ndim, const TAxis **axes, Option_t *opt = "");
  void CreateTProfile(const char *name, const char *title, int nbinsX, double xmin, double xmax, Option_t *opt = "");
  void CreateTProfile(const char *name, const char *title, int nbinsX, const double *xbins, Option_t *opt = "");
  void CreateTProfile(const char *name, const char *title, const TArrayD &xbins, Option_t *opt = "");
	void SetObject(TObject * const o, const char *group = "/");
	void FillTH1(const char *hname, double x, double weight = 1.);
	void FillTH1(const char *name, const char *label, double weight = 1.);
	void FillTH2(const char *hname, double x, double y, double weight = 1.);
	void FillTH2(const char *hname, double *point, double weight = 1.);
	void FillTH3(const char *hname, double x, double y, double z, double weight = 1.);
	void FillTH3(const char *hname, const double *point, double weight = 1.);
	void FillTHnSparse(const char *name, const double *x, double weight = 1.);
  void FillProfile(const char *name, double x, double y, double weight = 1.);

  iterator begin() const;
  iterator end() const;
  iterator rbegin() const;
  iterator rend() const;

  /**
   * Get the list of histograms
   * \return The list of histograms
   */
	THashList *GetListOfHistograms() const { return fHistos; }
	TObject *FindObject(const char *name) const;
	virtual TObject *FindObject(const TObject *obj) const;

private:
	THistManager(const THistManager &);
	THistManager &operator=(const THistManager &);
	THashList *FindGroup(const char *dirname) const;
	void TokenizeFilename(const char *name, const char *delim, std::vector<std::string> &listoftokens) const;
	TString basename(const char *path) const;
	TString histname(const char *path) const;

	THashList *fHistos;                   ///< List of histograms
	bool fIsOwner;                        ///< Set the ownership

  /// \cond CLASSIMP
	ClassDef(THistManager, 1);  // Container for histograms
  /// \endcond
};
#endif
