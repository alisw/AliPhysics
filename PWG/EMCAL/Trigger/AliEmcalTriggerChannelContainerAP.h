/**
 * @file AliEmcalTriggerChannelContainerAP.h
 * @since
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
#ifndef AliEmcalTriggerChannelContainerAPAP_H
#define AliEmcalTriggerChannelContainerAPAP_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <TSortedList.h>

/**
 * @struct AliEmcalTriggerChannelContainerAP
 * @brief Structure for position of trigger channels
 *
 * This structure is a container for trigger channels in col-row space with
 * a given mask. Channels can only be added to the container, or it can be
 * checked whether the channel is listed in the container.
 */
class AliEmcalTriggerChannelContainerAP : public TObject {
public:

  /**
   * @struct AliEmcalTriggerChannelPositionAP
   * @brief 2D position of a trigger channel on the EMCAL surface
   *
   * This class represents the position of a trigger channel in a 2D coordinate
   * system consisting of column and row.
   */
  class AliEmcalTriggerChannelPositionAP : public TObject{
  public:

    /**
     * Dummy (I/O) constructor, not to be used
     */
    AliEmcalTriggerChannelPositionAP(): fCol(-1), fRow(-1) { }

    /**
     * Main constuctor, setting the position in column and row
     * @param col Column of the trigger channel
     * @param row Row of the trigger channel
     */
    AliEmcalTriggerChannelPositionAP(int col, int row): fCol(col), fRow(row) { }

    /**
     * Destructor, nothing to do
     */
    virtual ~AliEmcalTriggerChannelPositionAP() {}

    /**
     * Get the column of the channel
     * @return  The column of the channel
     */
    int GetCol() const { return fCol; }

    /**
     * Get the row of the channel
     * @return The row of the channel
     */
    int GetRow() const { return fRow; }

    /**
     * Set the colummn of the channel
     * \param col The column of the channel
     */
    void SetCol(int col) { fCol = col; }

    /**
     * Set the row of the channel
     * @param row The row of the channel
     */
    void SetRow(int row) { fRow = row; }

    /**
     * Check if the object is equal to object ref. Object can only be equal if ref is of
     * the same type (AliEmcalTrigger channel position). If this is the case, col and row
     * of the two objects have to match.
     * @param ref The object to check
     * @return True if objects are equal, false otherwise
     */
    virtual Bool_t IsEqual(const TObject *ref);

    /**
     * Compare objects. If objects differ, return always greater (+1). Otherwise compare col and
     * row of the object. Col has priority with respect to row.
     * @param ref The object ot comparte to
     * @return 0 if objects are equal, -1 if this object is smaller, +1 if this object is larger.
     */
    virtual Int_t Compare(const TObject *ref);

  private:
    Int_t                    fCol;            ///< Column of the trigger channel
    Int_t                    fRow;            ///< Row of the trigger channel

    /// \cond CLASSIMP
    ClassDef(AliEmcalTriggerChannelPositionAP, 1);
    /// \endcond
  };

  /**
   * Constructor
   */
  AliEmcalTriggerChannelContainerAP(): fChannels() {}

  /**
   * Destructor, cleans up the container
   */
  virtual ~AliEmcalTriggerChannelContainerAP(){}

  /**
   * Add a new channel with the postion in column and row to the container, In case the channel
   * is already listed in the trigger channel container we don't add it again.
   * @param col Column of the channel
   * @param row Row of the channel
   */
  void AddChannel(int col, int row);

  /**
   * Check whether channel with the position (col, row) is listed in the trigger channel container
   * @param col Column of the channel
   * @param row Row of the channel
   * @return True if the channel is listed, false otherwise
   */
  bool HasChannel(int col, int row);

private:
  TSortedList             fChannels;      ///< Container for listed channels

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerChannelContainerAP, 1);
  /// \endcond
};

#endif
