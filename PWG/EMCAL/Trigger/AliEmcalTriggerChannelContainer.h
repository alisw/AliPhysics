#ifndef ALIEMCALTRIGGERCHANNELPOSITION_H
#define ALIEMCALTRIGGERCHANNELPOSITION_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <TSortedList.h>

/**
 * \struct AliEmcalTriggerChannelContainer
 * \brief Structure for position of trigger channels
 *
 * This structure is a container for trigger channels in col-row space with
 * a given mask. Channels can only be added to the container, or it can be
 * checked whether the channel is listed in the container.
 */
class AliEmcalTriggerChannelContainer : public TObject {
public:
  /**
   * \struct AliEmcalTriggerChannelPosition
   * \brief 2D position of a trigger channel on the EMCAL surface
   *
   * This class represents the position of a trigger channel in a 2D coordinate
   * system consisting of column and row.
   */
  class AliEmcalTriggerChannelPosition : public TObject{
  public:
    /**
     * Dummy (I/O) constructor, not to be used
     */
    AliEmcalTriggerChannelPosition(): fCol(-1), fRow(-1) { }
    /**
     * Main constuctor, setting the position in column and row
     * \param col Column of the trigger channel
     * \param row Row of the trigger channel
     */
    AliEmcalTriggerChannelPosition(int col, int row): fCol(col), fRow(row) { }
    /**
     * Destructor, nothing to do
     */
    virtual ~AliEmcalTriggerChannelPosition() {}

    /**
     * Get the column of the channel
     * \return  The column of the channel
     */
    int GetCol() const { return fCol; }
    /**
     * Get the row of the channel
     * \return The row of the channel
     */
    int GetRow() const { return fRow; }

    /**
     * Set the colummn of the channel
     * \param col The column of the channel
     */
    void SetCol(int col) { fCol = col; }
    /**
     * Set the row of the channel
     * \param row The row of the channel
     */
    void SetRow(int row) { fRow = row; }

    virtual Bool_t IsEqual(const TObject *ref);
    virtual Int_t Compare(const TObject *ref);
  private:
    Int_t                    fCol;            ///< Column of the trigger channel
    Int_t                    fRow;            ///< Row of the trigger channel

    /// \cond CLASSIMP
    ClassDef(AliEmcalTriggerChannelPosition, 1);
    /// \endcond
  };
  /**
   * Constructor
   */
  AliEmcalTriggerChannelContainer(): fChannels() {}

  /**
   * Destructor, cleans up the container
   */
  virtual ~AliEmcalTriggerChannelContainer(){}

  void AddChannel(int col, int row);
  bool HasChannel(int col, int row);

private:
  TSortedList             fChannels;      ///< Container for listed channels

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerChannelContainer, 1);
  /// \endcond
};

#endif
