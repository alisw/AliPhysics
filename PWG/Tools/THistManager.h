#ifndef THISTMANAGER_H
#define THISTMANAGER_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <THashList.h>
#include <TIterator.h>
#include <TNamed.h>
#include <iterator>

class TArrayD;
class TAxis;
class TBinning;
class TList;
class TH1;
class TH2;
class TH3;
class THnSparse;
class TProfile;

/**
 * @defgroup Histmanager
 * @brief Histogram manager and components needed to make it work.
 */

/**
 * @class THistManager
 * @brief Container class for histograms for the high-\f$ p_{t} \f$ charged particle analysis
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @ingroup Histmanager
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
   * @class iterator
   * @brief stl-iterator for the histogram manager
   * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
   * @ingroup Histmanager
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
     * @enum THMIDirection_t
     * @brief Direction for the iteration
     *
     * Switch for the firection of the iteration.
     * Currently the iterator is implemented as
     * bi-directional iterator, providing only
     * forward and backward iteration.
     */
    enum THMIDirection_t {
      kTHMIforward = 0,   //!< Forward iteration
      kTHMIbackward = 1   //!< Backward iteration
    };

    /**
     * Constructor. Initializing the histogram manager to be iterated over, the starting
     * position for the iteration, and the direction of the iteration.
     * @param[in] hmgr Histogram manager to be iterated over, containing the data
     * @param[in] currentpos Starting position for the iteration
     * @param[in] dir Direction of the iteration
     */
    iterator(const THistManager *hmgr, Int_t currentpos, THMIDirection_t dir = kTHMIforward);

    /**
     * Copy constructor. Initializing all
     * values from the reference.
     * @param[in] ref Reference for the copy
     */
    iterator(const iterator &ref);

    /**
     * Destructor
     */
    ~iterator() { }

    /**
     * Assignment operator. Initializing all
     *   values from the reference.
     * @param[in] ref
     * @return iterator after assignment.
     */
    iterator              &operator=(const iterator &rhs);

    /**
     * Comparison operator for unequalness. The comparison is
     * performed based on the current position of the iterator.
     * @param[in] other iterator to compare to
     * @return True if the iterators are not equal
     */
    Bool_t                operator!=(const iterator &aIter) const;

    /**
     * Prefix increment operator. Incrementing / decrementing
     * the current position of the iterator based on the
     * direction.
     * @return State after incrementation
     */
    iterator              &operator++();

    /**
     * Postfix increment operator. Incrementing / decrementing
     * the current position of the iterator based on the
     * direction, but returning the state before iteration
     * @return State before incrementation
     */
    iterator              operator++(int);

    /**
     * Prefix decrement operator. Decrementing / incrementing
     * the current position of the iterator based on the
     * direction.
     * @return State after decrementation
     */
    iterator              &operator--();

    /**
     * Postfix decrement operator. Decrementing / incrementing
     * the current position of the iterator based on the
     * direction, but returning the state before iteration
     * @return State before decrementation
     */
    iterator              operator--(int);

    /**
     * Dereferncing operator. Providing access to the object
     * connected to the current state of the iterator. Objects
     * might be histograms or THashLists in case of wheter groups
     * are defined or not.
     * @return Object connected ot the current state of the iterator
     */
    TObject               *operator*() const;

  private:
    const THistManager          *fkArray;
    Int_t                       fCurrentPos;
    Int_t                       fNext;
    THMIDirection_t             fDirection;

    iterator();
  };

  /**
   * Default constructor, only initialising pointers with 0
   */
	THistManager();

	/**
	 * Main constructor, creating also a list for the histograms
	 * @param name Name of the object (list named accordingly)
	 */
	THistManager(const char *name);

	/**
	 * Destructor, deletes the list of histograms if it is the owner
	 */
	~THistManager();

	void ReleaseOwner() { fIsOwner = kFALSE; };

	/**
	 * Create a new group of histograms within a parent group. Groups are represented as list. The default parent is
	 * always the top list. List name structure accouding to unix paths (i.e. top list /, hirarchies separated by /).
	 * @param groupname Name of the new group
	 * @param parent (default "/") Name of the parent group
	 * @throw HistoContainerContentException
	 */
	THashList* CreateHistoGroup(const char *groupname, const char *parent = "/");

	/**
	 * Create a new TH1 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param name Name of the histogram
	 * @param title Title of the histogram
	 * @param nbins number of bins
	 * @param xmin min. value of the range
	 * @param xmax max. value of the range
	 * @param opt Additonal options (s for sumw2)
	 */
	TH1* CreateTH1(const char *name, const char *title, int nbins, double xmin, double xmax, Option_t *opt = "");

	/**
	 * Create a new TH1 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] nbins number of bins
	 * @param[in] xbins array of bin limits
	 * @param[in] opt Additonal options (s for sumw2)
	 */
	TH1* CreateTH1(const char *name, const char *title, int nbins, const double *xbins, Option_t *opt = "");

	/**
	 * Create a new TH1 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] xbins array of bin limits (contains also number of bins)
	 * @param[in] opt Additonal options (s for sumw2)
	 */
	TH1* CreateTH1(const char *name, const char *title, const TArrayD &xbins, Option_t *opt = "");

	/**
	 * Create a new TH1 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] xbins User Binning
	 * @param[in] opt Additonal options (s for sumw2)
	 */
	TH1* CreateTH1(const char *name, const char *title, const TBinning &binning, Option_t *opt);

	/**
	 * Create a new TH2 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] nbinsx number of bins in x-direction
	 * @param[in] xmin min. value of the range in x-direction
	 * @param[in] xmax max. value of the range in x-direction
	 * @param[in] nbinsy number of bins in y-direction
	 * @param[in] ymin min. value of the range in y-direction
	 * @param[in] ymax max. value of the range in y-direction
	 */
	TH2* CreateTH2(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, Option_t *opt = "");

	/**
	 * Create a new TH2 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] nbinsx number of bins in x-direction
	 * @param[in] xmin min. value of the range in x-direction
	 * @param[in] xmax max. value of the range in x-direction
	 * @param[in] nbinsy number of bins in y-direction
	 * @param[in] ymin min. value of the range in y-direction
	 * @param[in] ymax max. value of the range in y-direction
	 */
	TH2* CreateTH2(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy, const double *ybins, Option_t *opt = "");

	/**
	 * Create a new TH2 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] xbins array of bin limits in x-direction (contains also the number of bins)
	 * @param[in] ybins array of bin limits in y-direction (contains also the number of bins)
	 */
	TH2* CreateTH2(const char *name, const char *title, const TArrayD &xbins, const TArrayD &ybins, Option_t *opt = "");

	/**
	 * Create a new TH2 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] User binning in x-direction
	 * @param[in] User binning in y-direction
	 */
	TH2* CreateTH2(const char *name, const char *title, const TBinning &xbins, const TBinning &ybins, Option_t *opt = "");

	/**
	 * Create a new TH2 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] nbinsx number of bins in x-direction
	 * @param[in] xmin min. value of the range in x-direction
	 * @param[in] xmax max. value of the range in x-direction
	 * @param[in] nbinsy number of bins in y-direction
	 * @param[in] ymin min. value of the range in y-direction
	 * @param[in] ymax max. value of the range in y-direction
	 */
	TH3* CreateTH3(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, int nbinsz, double zmin, double zmax, Option_t *opt = "");

	/**
	 * Create a new TH3 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] nbinsx number of bins in x-direction
	 * @param[in] xbins array of bin limits in x-direction
	 * @param[in] nbinsy number of bins in y-direction
	 * @param[in] ybins array of bin limits in y-direction
	 * @param[in] nbinsz number of bins in z-direction
	 * @param[in] zbins array of bin limits in z-direction
	 */
	TH3* CreateTH3(const char *name, const char *title, int nbinsx, const double *xbins, int nbinsy, const double *ybins, int nbinsz, const double *zbins, Option_t *opt = "");

	/**
	 * Create a new TH3 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] xbins array of bin limits in x-direction (contains also the number of bins)
	 * @param[in] ybins array of bin limits in y-direction (contains also the number of bins)
	 * @param[in] zbins array of bin limits in z-direction (contains also the number of bins)
	 */
	TH3* CreateTH3(const char *name, const char *title, const TArrayD &xbins, const TArrayD &ybins, const TArrayD &zbins, Option_t *opt = "");

	/**
	 * Create a new TH3 within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] User binning in x-direction
	 * @param[in] User binning in y-direction
	 * @param[in] User binning in z-direction
	 */
	TH3* CreateTH3(const char *name, const char *title, const TBinning &xbins, const TBinning &ybins, const TBinning &zbins, Option_t *opt = "");

	/**
	 * Create a new THnSparse within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] ndim Number of dimensions
	 * @param[in] nbins Number of bins per dimension
	 * @param[in] min min. value of the range for each dimension
	 * @param[in] max max. value of the range for each dimension
	 */
	THnSparse* CreateTHnSparse(const char *name, const char *title, int ndim, const int *nbins, const double *min, const double *max, Option_t *opt = "");

	/**
	 * Create a new THnSparse within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] ndim Number of dimensions
	 * @param[in] axes Array of pointers to TAxis for containing the axis definition for each dimension
	 */
	THnSparse* CreateTHnSparse(const char *name, const char *title, int ndim, const TAxis **axes, Option_t *opt = "");

	/**
	 * Create a new TProfile within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the profile histogram
	 * @param[in] title Title of the profile histogram
	 * @param[in] nbinsX Number of bins in x-direction
	 * @param[in] xmin min. value in x-direction
	 * @param[in] xmax max. value in x-direction
	 * @param[in] opt Further options
	 */
  void CreateTProfile(const char *name, const char *title, int nbinsX, double xmin, double xmax, Option_t *opt = "");

  /**
   * Create a new TProfile within the container. The histogram name also contains the parent group(s) according to the common
   * group notation.
   * @param[in] name Name of the profile histogram
   * @param[in] title Title of the profile histogram
   * @param[in] nbinsX Number of bins in x-direction
   * @param[in] xbins binning in x-direction
   * @param[in] opt Further options
   */
  void CreateTProfile(const char *name, const char *title, int nbinsX, const double *xbins, Option_t *opt = "");

  /**
   * Create a new TProfile within the container. The histogram name also contains the parent group(s) according to the common
   * group notation.
   * @param[in] name Name of the profile histogram
   * @param[in] title Title of the profile histogram
   * @param[in] xbins binning in x-direction
   * @param[in] opt Further options
   */
  void CreateTProfile(const char *name, const char *title, const TArrayD &xbins, Option_t *opt = "");

  /**
   * Create a new TProfile within the container. The histogram name also contains the parent group(s) according to the common
   * group notation.
   * @param[in] name Name of the profile histogram
   * @param[in] title Title of the profile histogram
   * @param[in] xbins User binning
   * @param[in] opt Further options
   */
  void CreateTProfile(const char *name, const char *title, const TBinning &xbins, Option_t *opt = "");

  /**
   * Set a new group into the container into the parent group
   * @param[in] o the object to be included
   */
	void SetObject(TObject * const o, const char *group = "/");

	/**
	 * Fill a 1D histogram within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] x x-coordinate
	 * @param[in] weight optional weight of the entry (default 1)
	 */
	void FillTH1(const char *hname, double x, double weight = 1., Option_t *opt = "");

	void FillTH1(const char *name, const char *label, double weight = 1., Option_t *opt = "");

	/**
	 * Fill a 2D histogram within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] x x-coordinate
	 * @param[in] y y-coordinate
	 * @param[in] weight optional weight of the entry (default 1)
	 */
	void FillTH2(const char *hname, double x, double y, double weight = 1., Option_t *opt = "");

	/**
	 * Fill a 2D histogram within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] point coordinates of the data
	 * @param[in] weight optional weight of the entry (default 1)
	 */
	void FillTH2(const char *hname, double *point, double weight = 1., Option_t *opt = "");

	/**
	 * Fill a 3D histogram within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] x x-coordinate
	 * @param[in] y y-coordinate
	 * @param[in] z z-coordinate
	 * @param[in] weight optional weight of the entry (default 1)
	 */
	void FillTH3(const char *hname, double x, double y, double z, double weight = 1., Option_t *opt = "");

	/**
	 * Fill a 3D histogram within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] point 3D-coordinate (x,y,z) of the point to be filled
	 * @param[in] weight optional weight of the entry (default 1)
	 */
	void FillTH3(const char *hname, const double *point, double weight = 1., Option_t *opt = "");

	/**
	 * Fill a  nD histogram within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] x coordinates of the data
	 * @param[in] weight optional weight of the entry (default 1)
	 */
	void FillTHnSparse(const char *name, const double *x, double weight = 1., Option_t *opt = "");

	/**
	 * Fill a profile histogram within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the profile histogram
	 * @param[in] x x-coordinate
	 * @param[in] y y-coordinate
	 * @param[in] weight optional weight of the entry (default 1)
	 */
  void FillProfile(const char *name, double x, double y, double weight = 1.);

  /**
   * Create forward iterator starting at the beginning of the
   * container
   * @return Forward iterator at the beginning of the container
   */
  inline iterator begin() const;

  /**
   * Create forward iterator starting at the end of the
   * container
   * @return Forward iterator after the end of the container
   */
  inline iterator end() const;

  /**
   * Create backward iterator starting behind the end of the container.
   * Used to terminate the iteration.
   * @return Backward iterator behind the end of the histogram manager
   */
  inline iterator rbegin() const;

  /**
   * Create backward iterator starting before the beginning of the container.
   * Used to terminate the iteration.
   * @return Backward iterator starting before the beginning of the histogram manager
   */
  inline iterator rend() const;

  /**
   * Get the list of histograms
   * @return The list of histograms
   */
	THashList *GetListOfHistograms() const { return fHistos; }

	/**
	 * Find an object inside the container. The object can also be within a
	 * histogram group. For this the name has to follow the common notation
	 * @param[in] name Name of the object to find inside the container
	 * @return pointer to the object (NULL if not found)
	 */
	TObject *FindObject(const char *name) const;

	/**
	 * Find and object inside the container. The object name is expected to contain the
	 * full path of the histogram object, including parent groups
	 * @param[in] obj the object to find
	 * @return pointer to the object (NULL if not found)
	 */
	virtual TObject *FindObject(const TObject *obj) const;

private:
	THistManager(const THistManager &);
	THistManager &operator=(const THistManager &);


	/**
	 * Find histogram group. Name is using common notation
	 * @param[in] dirname Path of the group (treat empty path as top node
	 * @return TList of objects (NULL if group does not exist)
	 */
	THashList *FindGroup(const char *dirname) const;

	/**
	 * Tokenizes a string. Results are stored inside the vector listoftokens
	 * @param[in] name string to be tokenised
	 * @param[in] delim delimiter string
	 * @param[in] listoftokens list of tokens (C++ strings)
	 */
	void TokenizeFilename(const char *name, const char *delim, std::vector<std::string> &listoftokens) const;

	/**
	 * Helper function extracting the basename from a given histogram path.
	 * @param[in] path histogram path
	 * @return basename extracted
	 */
	TString basename(const char *path) const;

	/**
	 * Helper function extracting the histogram name from a given histogram path.
	 * @param[in] path histogram path
	 * @return basename extracted
	 */
	TString histname(const char *path) const;

	THashList *fHistos;                   ///< List of histograms
	bool fIsOwner;                        ///< Set the ownership

  /// \cond CLASSIMP
	ClassDef(THistManager, 1);  // Container for histograms
  /// \endcond
};

THistManager::iterator THistManager::begin() const {
  return iterator(this, 0, iterator::kTHMIforward);
}

THistManager::iterator THistManager::rbegin() const {
  return iterator(this, fHistos->GetEntries()-1, iterator::kTHMIbackward);
}

THistManager::iterator THistManager::end() const {
  return iterator(this, fHistos->GetEntries(), iterator::kTHMIforward);
}

THistManager::iterator THistManager::rend() const {
  return iterator(this, -1, iterator::kTHMIbackward);
}
#endif
