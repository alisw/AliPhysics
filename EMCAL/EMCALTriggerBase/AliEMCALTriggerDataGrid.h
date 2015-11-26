#ifndef ALIEMCALTRIGGERDATAGRID_H
#define ALIEMCALTRIGGERDATAGRID_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <exception>
#include <sstream>
#include <string>

/**
 * \class AliEmcalTriggerDataGrid
 * \brief Container for ADC / Amplitudes from the EMCAL triggers
 *
 * Dynamical-size container for ADC values from the FASTOR
 */
template<typename T>
class AliEMCALTriggerDataGrid : public TObject {
public:
  /**
   * \class UninitException
   * \brief Error handling for uninitialized grid
   */
  class UninitException : public std::exception{
  public:
    /**
     * Constructor
     */
    UninitException():
      std::exception()
    {}
    /**
     * Destructor, nothing to do
     */
    virtual ~UninitException() throw() {}

    /**
     * Access error message connected to the exception
     * @return Error message
     */
    virtual const char *what() const throw() { return "Trigger channel map not initialized"; }
  };

  /**
   * \class OutOfBoundsException
   * \brief Exception class handling access to non-existing container element
   */
  class OutOfBoundsException : public std::exception{
  public:
    /**
     * Definition of directions
     */
    enum Direction_t {
      kColDir     = 0,///< Column direction (eta)
      kRowDir     = 1,///< Row direction (phi)
      kUndef      = 2 ///< Not defined
    };
    /**
     * Dumny constructor
     */
    OutOfBoundsException():
      std::exception(),
      fMessage(),
      fDir(kUndef),
      fSize(0),
      fIndex()
    {
    }
    /**
     * Regular constructor, to be called when exception is thrown
     * @param dir Direction (col or row)
     * @param index Index for which exception is thrown
     * @param size Size of the grid in direction
     */
    OutOfBoundsException(Direction_t dir, int index, int size):
      std::exception(),
      fMessage(""),
      fDir(dir),
      fSize(size),
      fIndex(index)
    {
      std::stringstream errormessage;
      errormessage << "Out-of-bounds access in " << fDir << "Direction: Element " << fIndex << ", Size " << fSize;
      fMessage = errormessage.str();
    }
    /**
     * Destructor
     */
    virtual ~OutOfBoundsException() throw() {}

    /**
     * Get error message
     * @return error message (created in constructor)
     */
    const char *what() const throw() { return fMessage.c_str(); }

    /**
     * Get the size of the grid in direction for which exception is thrown
     * @return Size of the grid in direction
     */
    int GetSize() const { return fSize; }

    /**
     * Get index for which exception is thrown
     * @return Index for which exception is thrown
     */
    int GetIndex() const { return fIndex;}

    /**
     * Get the direction for which exception is thrown.
     * @return Direction for which exception is thrown.
     */
    Direction_t GetDirection() const { return fDir; }

  private:
    std::string         fMessage; ///< Error message, accessible via "what"
    Direction_t         fDir;     ///< Direction
    int                 fSize;    ///< size of the container in direction
    int                 fIndex;   ///< Index requested
  };

  /**
   * Dummy constructor, does not allocate anything
   */
  AliEMCALTriggerDataGrid();

  /**
   * Constructror
   * Allocates also for storage for the ADC values / amplitudes
   * @param cols Number of cols
   * @param rows Number of rows
   */
  AliEMCALTriggerDataGrid(Int_t cols, Int_t rows);

  /**
   * Copy constructor
   * New channel map will get its own storage. The content
   * of the ref storage will be copied into this storage
   * @param ref Reference for the copy
   */
  AliEMCALTriggerDataGrid(const AliEMCALTriggerDataGrid<T> &ref);

  /**
   * Assignment operator
   * New channel map will get its own storage. The content
   * of the ref storage will be copied into this storage
   * @param ref Reference for the copy
   * @return This channel map
   */
  AliEMCALTriggerDataGrid<T> &operator=(const AliEMCALTriggerDataGrid<T> &ref);

  /**
   * Constant acces operator at position (col, row)
   * @param col Column
   * @param row Row
   * @return Constant reference to entry at that position (can not modify the entry)
   */
  const T &operator()(Int_t col, Int_t row) const;

  /**
   * Access operator at position (col, row)
   * @param col Column
   * @param row Row
   * @return Reference to entry at that position (can modify the entry)
   */
  T &operator()(Int_t col, Int_t row);

  /**
   * Destructor
   */
  virtual ~AliEMCALTriggerDataGrid();

  /**
   * Set the ADC values stored in the 2D map again to 0
   */
  void Reset();


  Bool_t IsAllocated() const { return fValues != NULL; }

  void Allocate(Int_t ncols, Int_t nrows);

  /**
   * Set ADC value for position (col, row). Checks for boundary.
   * @param col Column of the position
   * @param row Row of the position
   * @param ADC The value to set
   * @throw OutOfBoundsException in case the index in any direction is out of bounds
   * @throw UninitException in case the grid is not initialized (allocated)
   */
  void SetADC(Int_t col, Int_t row, const T &adc);

  /**
   * Get ADC value at position (col, row). Checks for boundary.
   * @param col
   * @param row
   * @return
   * @throw OutOfBoundsException in case the index in any direction is out of bounds
   * @throw UninitException in case the grid is not initialized (allocated)
   */
  const T &GetADC(Int_t col, Int_t row) const;

  /**
   * Get the number of columns in the map
   * @return The number of colums
   */
  Int_t GetNumberOfCols() const { return fNCols; }

  /**
   * Get the number of rows in the map
   * @return
   */
  Int_t GetNumberOfRows() const { return fNRows; }

protected:
  /**
   * Get grid index in the ADC value list
   * @param col Column of the grid
   * @param row Row of the grid
   * @return Grid index
   */
  Int_t GetIndex(Int_t col, Int_t row) const;

  Int_t                     fNCols;           ///< Number of columns
  Int_t                     fNRows;           ///< Number of rows
  T                         *fValues;         ///< Array of Trigger ADC values

  /// \cond
  ClassDef(AliEMCALTriggerDataGrid, 1);
  /// \endcond
};


#endif
