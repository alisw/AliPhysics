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
 * @brief Container class for histograms
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @ingroup Histmanager
 *
 * # A container class for histogram handling
 *
 * Usual analyses contain a set of histograms to be handled in analysis tasks.
 * They are usually created in a UserCreateOutputObjects method, handled as
 * data members, and filled in certain methods of the UserExec event loop.
 *
 * The THistManager provides functionality handling histograms for user
 * analyses in a simple, transparent and dynamic way. Histograms are created
 * and automatically added to the histogram manager in Create methods. Histograms
 * are handled via their names. Several Fill methods are available in order to
 * fill a histogram with a certain name.
 *
 * Currently the histogram manager can handle the following types:
 * - TH1
 * - TH2
 * - TH3
 * - THnSparse
 * - TProfile
 * For histograms with multiple data types always the double precision
 * version is used.
 *
 * #Structurizing histogram output
 *
 * Histograms can be stored in groups. For this the parent group is
 * included inside the histogram name, i.e. /base/inheriting/histogram.
 * In case just the histogram name is given, it is assumed that the
 * histogram is stored at the top level.
 *
 * # Creating Histograms
 *
 * Creating histograms is done using the Create method for the various
 * histogram types. Create methods need a name and a binning definition.
 * It is recommended to use different @ref{TBinning} classes to initialize
 * the binning. Once created the histogram is added to the histmanager
 * and can be filled with the corresponding Fill method.
 *
 * In the following example we create a 1-dimensional histogram:
 *
 * ~~~{.cxx}
 * THistManager mgr(testhists);
 * mgr.CreateTH1("hPt", "pt-distribution", TLinearBinning(100, 0., 100.));
 * ~~~
 *
 * # Filling histograms
 *
 * Once histograms are created, they can be filled with the corresponding
 * fill method of a histogram type. For some histograms different Fill
 * methods are provided to cover various use cases.
 *
 * In the following example we fill the histogram created in the Create section
 * with random values of an exponential distribution.
 *
 * ~~~{.cxx}
 * for(auto en : ROOT::TSeqI(0, 10000) {
 *   double pt = gRandom->Exp(-1);
 *   mgr.FillTH1("hPt", pt);
 * }
 * ~~~
 *
 * ## Optional automatic correction of the bin width
 *
 * Correction for the bin width can be automatically handled by the histogram
 * manager when filling the histogram. For this purpose the Fill methods provide
 * an argument for options. Automatic correction for the bin width is done when
 * specifying the argument *W*, followed by the direction. Adding multiple directions
 * the weight is calculated for all directions at the same time.
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
     * @brief Constructor.
     *
     * Initializing the histogram manager to be iterated over, the starting
     * position for the iteration, and the direction of the iteration.
     * @param[in] hmgr Histogram manager to be iterated over, containing the data
     * @param[in] currentpos Starting position for the iteration
     * @param[in] dir Direction of the iteration
     */
    iterator(const THistManager *hmgr, Int_t currentpos, THMIDirection_t dir = kTHMIforward);

    /**
     * @brief Copy constructor.
     *
     * Initializing all values from the reference.
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
     * @brief Comparison operator for unequalness.
     *
     * The comparison is performed based on the current
     * position of the iterator.
     * @param[in] other iterator to compare to
     * @return True if the iterators are not equal
     */
    Bool_t                operator!=(const iterator &aIter) const;

    /**
     * @brief Prefix increment operator.
     *
     * Incrementing / decrementing the current position of the
     * iterator based on the direction.
     * @return State after incrementation
     */
    iterator              &operator++();

    /**
     * @brief Postfix increment operator.
     *
     * Incrementing / decrementing the current position of the
     * iterator based on the direction, but returning the state
     * before iteration
     * @return State before incrementation
     */
    iterator              operator++(int);

    /**
     * @brief Prefix decrement operator.
     *
     * Decrementing / incrementing the current position
     * of the iterator based on the direction.
     * @return State after decrementation
     */
    iterator              &operator--();

    /**
     * @brief Postfix decrement operator.
     *
     * Decrementing / incrementing the current position of
     * the iterator based on the direction, but returning
     * the state before iteration
     * @return State before decrementation
     */
    iterator              operator--(int);

    /**
     * @brief Dereferncing operator.
     *
     * Providing access to the object connected to the current
     * state of the iterator. Objects might be histograms or
     * THashLists in case of wheter groups are defined or not.
     * @return Object connected ot the current state of the iterator
     */
    TObject               *operator*() const;

  private:
    const THistManager          *fkArray;             ///< Underlying histmanager to iterate over
    Int_t                       fCurrentPos;          ///< Current position of the iterator in the histmanager
    Int_t                       fNext;                ///< Next position in the histmanager
    THMIDirection_t             fDirection;           ///< Direction of the iterator

    iterator();
  };

  /**
   * @brief Default constructor.
   *
   * Only initialising pointers with 0
   */
	THistManager();

	/**
	 * @brief Main constructor.
	 *
	 * Creating also a list for the histograms
	 * @param name Name of the object (list named accordingly)
	 */
	THistManager(const char *name);

	/**
	 * @brief Destructor.
	 *
	 * Deletes the list of histograms if it is the owner
	 */
	~THistManager();

	void ReleaseOwner() { fIsOwner = kFALSE; };

	/**
	 * @brief Create a new group of histograms within a parent group.
	 *
	 * Groups are represented as list. The default parent is always
	 * the top list. List name structure accouding to unix paths
	 * (i.e. top list /, hirarchies separated by /).
	 * @param groupname Name of the new group
	 * @param parent (default "/") Name of the parent group
	 * @throw HistoContainerContentException
	 */
	THashList* CreateHistoGroup(const char *groupname);

	/**
	 * @brief Create a new TH1 within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
	 * @param name Name of the histogram
	 * @param title Title of the histogram
	 * @param nbins number of bins
	 * @param xmin min. value of the range
	 * @param xmax max. value of the range
	 * @param opt Additonal options (s for sumw2)
	 */
	TH1* CreateTH1(const char *name, const char *title, int nbins, double xmin, double xmax, Option_t *opt = "");

	/**
	 * @brief Create a new TH1 within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] nbins number of bins
	 * @param[in] xbins array of bin limits
	 * @param[in] opt Additonal options (s for sumw2)
	 */
	TH1* CreateTH1(const char *name, const char *title, int nbins, const double *xbins, Option_t *opt = "");

	/**
	 * @brief Create a new TH1 within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] xbins array of bin limits (contains also number of bins)
	 * @param[in] opt Additonal options (s for sumw2)
	 */
	TH1* CreateTH1(const char *name, const char *title, const TArrayD &xbins, Option_t *opt = "");

	/**
	 * @brief Create a new TH1 within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] xbins User Binning
	 * @param[in] opt Additonal options (s for sumw2)
	 */
	TH1* CreateTH1(const char *name, const char *title, const TBinning &binning, Option_t *opt = "");

	/**
	 * @brief Create a new TH2 within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
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
	 * @brief Create a new TH2 within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
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
	 * @brief Create a new TH2 within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] xbins array of bin limits in x-direction (contains also the number of bins)
	 * @param[in] ybins array of bin limits in y-direction (contains also the number of bins)
	 */
	TH2* CreateTH2(const char *name, const char *title, const TArrayD &xbins, const TArrayD &ybins, Option_t *opt = "");

	/**
	 * @brief Create a new TH2 within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] User binning in x-direction
	 * @param[in] User binning in y-direction
	 */
	TH2* CreateTH2(const char *name, const char *title, const TBinning &xbins, const TBinning &ybins, Option_t *opt = "");

	/**
	 * @brief Create a new TH2 within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
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
	 * @brief Create a new TH3 within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
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
	 * @brief Create a new TH3 within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] xbins array of bin limits in x-direction (contains also the number of bins)
	 * @param[in] ybins array of bin limits in y-direction (contains also the number of bins)
	 * @param[in] zbins array of bin limits in z-direction (contains also the number of bins)
	 */
	TH3* CreateTH3(const char *name, const char *title, const TArrayD &xbins, const TArrayD &ybins, const TArrayD &zbins, Option_t *opt = "");

	/**
	 * @brief Create a new TH3 within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] User binning in x-direction
	 * @param[in] User binning in y-direction
	 * @param[in] User binning in z-direction
	 */
	TH3* CreateTH3(const char *name, const char *title, const TBinning &xbins, const TBinning &ybins, const TBinning &zbins, Option_t *opt = "");

	/**
	 * @brief Create a new THnSparse within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] ndim Number of dimensions
	 * @param[in] nbins Number of bins per dimension
	 * @param[in] min min. value of the range for each dimension
	 * @param[in] max max. value of the range for each dimension
	 */
	THnSparse* CreateTHnSparse(const char *name, const char *title, int ndim, const int *nbins, const double *min, const double *max, Option_t *opt = "");

	/**
	 * @brief Create a new THnSparse within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] title Title of the histogram
	 * @param[in] ndim Number of dimensions
	 * @param[in] axes Array of pointers to TAxis for containing the axis definition for each dimension
	 */
	THnSparse* CreateTHnSparse(const char *name, const char *title, int ndim, const TAxis **axes, Option_t *opt = "");

  /**
   * @brief Create a new THnSparse within the container.
   *
   * The histogram name also contains the parent group(s)
   * according to the common group notation.
   * @param[in] name Name of the histogram
   * @param[in] title Title of the histogram
   * @param[in] ndim Number of dimensions
   * @param[in] axes Array of pointers to TAxis for containing the axis definition for each dimension
   */
  THnSparse* CreateTHnSparse(const char *name, const char *title, int ndim, const TBinning **axes, Option_t *opt = "");


	/**
	 * @brief Create a new TProfile within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
	 * @param[in] name Name of the profile histogram
	 * @param[in] title Title of the profile histogram
	 * @param[in] nbinsX Number of bins in x-direction
	 * @param[in] xmin min. value in x-direction
	 * @param[in] xmax max. value in x-direction
	 * @param[in] opt Further options
	 */
  void CreateTProfile(const char *name, const char *title, int nbinsX, double xmin, double xmax, Option_t *opt = "");

  /**
   * @brief Create a new TProfile within the container.
   *
   * The histogram name also contains the parent group(s)
   * according to the common group notation.
   * @param[in] name Name of the profile histogram
   * @param[in] title Title of the profile histogram
   * @param[in] nbinsX Number of bins in x-direction
   * @param[in] xbins binning in x-direction
   * @param[in] opt Further options
   */
  void CreateTProfile(const char *name, const char *title, int nbinsX, const double *xbins, Option_t *opt = "");

  /**
   * @brief Create a new TProfile within the container.
   *
   * The histogram name also contains the parent group(s)
   * according to the common group notation.
   * @param[in] name Name of the profile histogram
   * @param[in] title Title of the profile histogram
   * @param[in] xbins binning in x-direction
   * @param[in] opt Further options
   */
  void CreateTProfile(const char *name, const char *title, const TArrayD &xbins, Option_t *opt = "");

  /**
   * @brief Create a new TProfile within the container.
   *
   * The histogram name also contains the parent group(s)
   * according to the common group notation.
   * @param[in] name Name of the profile histogram
   * @param[in] title Title of the profile histogram
   * @param[in] xbins User binning
   * @param[in] opt Further options
   */
  void CreateTProfile(const char *name, const char *title, const TBinning &xbins, Option_t *opt = "");

  /**
   * @brief Set a new group into the container into the parent group
   * @param[in] o the object to be included
   */
	void SetObject(TObject * const o, const char *group = "/");

	/**
	 * @brief Fill a 1D histogram within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] x x-coordinate
	 * @param[in] weight optional weight of the entry (default 1)
	 * @param[in] option Optional filling arguments
	 */
	void FillTH1(const char *hname, double x, double weight = 1., Option_t *opt = "");

  /**
   * @brief Fill a 1D histogram within the container.
   *
   * Instead of an x-value a bin label is used instead.
   * The histogram name also contains the parent group(s)
   * according to the common group notation.
   * @param[in] name Name of the histogram
   * @param[in] label Label of the bin to fill
   * @param[in] weight optional weight of the entry (default 1)
   * @param[in] option Optional filling arguments
   */
	void FillTH1(const char *name, const char *label, double weight = 1., Option_t *opt = "");

	/**
	 * @brief Fill a 2D histogram within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] x x-coordinate
	 * @param[in] y y-coordinate
	 * @param[in] weight optional weight of the entry (default 1)
   * @param[in] option Optional filling arguments
	 */
	void FillTH2(const char *hname, double x, double y, double weight = 1., Option_t *opt = "");

  /**
   * @brief Fill a 2D histogram within the container.
   *
   * The histogram name also contains the parent group(s)
   * according to the common group notation.
   * @param[in] name Name of the histogram
   * @param[in] labelX x-coordinate
   * @param[in] labelY y-coordinate
   * @param[in] weight optional weight of the entry (default 1)
   * @param[in] option Optional filling arguments
   */
	void FillTH2(const char *name, const char *labelX, const char *labelY, double weight = 1., Option_t *opt = "");

	/**
	 * @brief Fill a 2D histogram within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] point coordinates of the data
	 * @param[in] weight optional weight of the entry (default 1)
   * @param[in] option Optional filling arguments
	 */
	void FillTH2(const char *hname, double *point, double weight = 1., Option_t *opt = "");

	/**
	 * @brief Fill a 3D histogram within the container.
	 *
	 * The histogram name also contains the parent group(s)
	 * according to the common group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] x x-coordinate
	 * @param[in] y y-coordinate
	 * @param[in] z z-coordinate
	 * @param[in] weight optional weight of the entry (default 1)
   * @param[in] option Optional filling arguments
	 */
	void FillTH3(const char *hname, double x, double y, double z, double weight = 1., Option_t *opt = "");

	/**
	 * Fill a 3D histogram within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] point 3D-coordinate (x,y,z) of the point to be filled
	 * @param[in] weight optional weight of the entry (default 1)
   * @param[in] option Optional filling arguments
	 */
	void FillTH3(const char *hname, const double *point, double weight = 1., Option_t *opt = "");

	/**
	 * Fill a  nD histogram within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the histogram
	 * @param[in] x coordinates of the data
	 * @param[in] weight optional weight of the entry (default 1)
   * @param[in] option Optional filling arguments
	 */
	void FillTHnSparse(const char *name, const double *x, double weight = 1., Option_t *opt = "");

	/**
	 * Fill a profile histogram within the container. The histogram name also contains the parent group(s) according to the common
	 * group notation.
	 * @param[in] name Name of the profile histogram
	 * @param[in] x x-coordinate
	 * @param[in] y y-coordinate
	 * @param[in] weight optional weight of the entry (default 1)
   * @param[in] option Optional filling arguments
	 */
  void FillProfile(const char *name, double x, double y, double weight = 1.);

  /**
   * @brief Create forward iterator starting at the beginning of the
   * container
   * @return Forward iterator at the beginning of the container
   */
  inline iterator begin() const;

  /**
   * @brief Create forward iterator starting at the end of the
   * container
   * @return Forward iterator after the end of the container
   */
  inline iterator end() const;

  /**
   * @brief Create backward iterator starting behind the end of the container.
   *
   * Used to terminate the iteration.
   * @return Backward iterator behind the end of the histogram manager
   */
  inline iterator rbegin() const;

  /**
   * @brief Create backward iterator starting before the beginning of the container.
   *
   * Used to terminate the iteration.
   * @return Backward iterator starting before the beginning of the histogram manager
   */
  inline iterator rend() const;

  /**
   * @brief Get the list of histograms.
   * @return The list of histograms
   */
	THashList *GetListOfHistograms() const { return fHistos; }

	/**
	 * @brief Find an object inside the container.
	 *
	 * The object can also be within a histogram group.
	 * For this the name has to follow the common notation.
	 * @param[in] name Name of the object to find inside the container
	 * @return pointer to the object (NULL if not found)
	 */
	TObject *FindObject(const char *name) const;

	/**
	 * @brief Find and object inside the container.
	 *
	 * The object name is expected to contain the full
	 * path of the histogram object, including parent
	 * groups.
	 * @param[in] obj the object to find
	 * @return pointer to the object (NULL if not found)
	 */
	virtual TObject *FindObject(const TObject *obj) const;

private:
	THistManager(const THistManager &);
	THistManager &operator=(const THistManager &);


	/**
	 * @brief Find histogram group.
	 *
	 * Name is using common notation
	 * @param[in] dirname Path of the group (treat empty path as top node
	 * @return TList of objects (NULL if group does not exist)
	 */
	THashList *FindGroup(const char *dirname) const;

	/**
	 * @brief Extracting the basename from a given histogram path.
	 * @param[in] path histogram path
	 * @return basename extracted
	 */
	TString basename(const TString &path) const;

	/**
	 * @brief Extracting the histogram name from a given histogram path.
	 * @param[in] path histogram path
	 * @return basename extracted
	 */
	TString histname(const TString &path) const;

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

/**
 * @namespace TestTHistManager
 * @brief Collection of simple test for the THistManager
 * @ingroup Histmanager
 */
namespace TestTHistManager {

/**
 * @class THistManagerTestSuite
 * @brief Collection of tests for the THistManager
 * @ingroup Histmanager
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Nov 23, 2016
 *
 * Test suite for histogram manager. Currently implemented tests:
 * - Build simple histograms
 * - Build histrogram in groups
 * - Simple fill
 * - Fill histograms in groups
 */
class THistManagerTestSuite {
public:

  /**
   * Constructor
   */
  THistManagerTestSuite() {}

  /**
   * Destructor
   */
  virtual ~THistManagerTestSuite() {}

  /**
   * Purpose of the test: Check whether histmanager builds histogram correctly (no grouping)
   *
   * Create 1 histogram of each type
   * - TH1
   * - TH2
   * - TH3
   * - THnSparse
   * - TProfile
   * Find the histogram according to its name in the list of histograms
   *
   * Test passed: All histograms need to be found in the list of histograms, and type must match
   * @return 0 if test is passed, 1 if it failed
   */
  int TestBuildSimpleHistograms();

  /**
   * Purpose of the test: Check whether histmanager assings histograms correctly into groups
   * Relies on: TestBuildSimpleHistograms
   *
   * Creating 3 groups with 2 histograms
   * - Group 1 has only 1D histograms
   * - Group 2 has only 2D histograms
   * - Group 3 has only 3D histograms
   * Names of the histograms are the same for all groups.
   *
   * In addition: 1 group with a subgroup and a TProfile.
   *
   * Test passed:
   * - All Groups need to be found, and all histograms within the group with the matching type
   * - Subgroup need to be found in the subgroup test
   * @return 0 if test is passed, 1 if it failed
   */
  int TestBuildGroupedHistograms();

  /**
   * Purpose of the test: Check whether histograms are found correctly by the automatic procedure, and
   * whether fill is correctly propagated
   * Relies on: TestBuildSimpleHistograms
   *
   * Creating histograms of all types
   * - TH1
   * - TH2
   * - TH3
   * - THnSparse
   * - TProfile
   * with 1 bin per dimension and each with 100 times the same value.
   *
   * Test passed:
   * - All histograms need to have in its 1 bin the bin content 100
   * @return 0 if test is passed, 1 if it failed
   */
  int TestFillSimpleHistograms();

  /**
   * Purpose of the test: Test access of histograms in groups, check whether histograms are filled properly
   * Relies on: TestBuildSimpleHistograms, TestBuildGroupedHistograms, TestFillSimpleHistograms
   *
   * Fill 2 test histograms in 2 groups
   * - Group1: TH1
   * - Group2: TH2
   * each 100 times for bin 1. In addition Fill TProfile in Group3 with Subgroup1 100 times with weight1
   *
   * Test passed:
   * - All Histograms are found in the groups and subgroups
   * - All Histograms have the expected value (100 for histograms, 1 for profile)
   * @return 0 if test is passed, 1 if it failed
   */
  int TestFillGroupedHistograms();
};

/**
 * Runs all tests for THistManager. See @ref THistManagerTestSuite
 * for details.
 * @return O if all tests are passed, 1 if tests fail
 */
int TestRunAll();

/**
 * Run the test for building histograms. See @ref THistManagerTestSuite
 * for details.
 * @return 0 if test is passed, 1 if failed
 */
int TestRunBuildSimple();

/**
 * Run the test for building histograms in groups. See @ref THistManagerTestSuite
 * for details.
 * @return 0 if test is passed, 1 if failed
 */
int TestRunBuildGrouped();

/**
 * Run the test for filling histograms. See @ref THistManagerTestSuite
 * for details.
 * @return 0 if test is passed, 1 if failed
 */
int TestRunFillSimple();

/**
 * Run the test for building histograms. See @ref THistManagerTestSuite
 * for details.
 * @return 0 if test is passed, 1 if failed
 */
int TestRunFillGrouped();

}
#endif
