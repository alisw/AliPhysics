/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#ifndef ALIEMCALTRIGGERPARTBADCHANNELCONTAINER_H
#define ALIEMCALTRIGGERPARTBADCHANNELCONTAINER_H

#include <vector>
#include <TObject.h>

namespace PWG {

namespace EMCAL {

namespace TriggerPart {

/**
 * @struct AliEmcalTriggerPartBadChannelContainer
 * @brief Structure for position of trigger channels
 *
 * This structure is a container for trigger channels in col-row space with
 * a given mask. Channels can only be added to the container, or it can be
 * checked whether the channel is listed in the container.
 */
class AliEmcalTriggerPartBadChannelContainer : public TObject {
public:

  /**
   * @struct TriggerChannelPosition
   * @brief 2D position of a trigger channel on the EMCAL surface
   *
   * This class represents the position of a trigger channel in a 2D coordinate
   * system consisting of column and row.
   */
  class TriggerChannelPosition {
  public:

    /**
     * Dummy (I/O) constructor, not to be used
     */
    TriggerChannelPosition(): fCol(-1), fRow(-1) { }

    /**
     * Main constuctor, setting the position in column and row
     * @param col Column of the trigger channel
     * @param row Row of the trigger channel
     */
    TriggerChannelPosition(int col, int row): fCol(col), fRow(row) { }

    /**
     * Destructor, nothing to do
     */
    virtual ~TriggerChannelPosition() {}

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
    bool operator==(const TriggerChannelPosition &ref) const;

    /**
     * Compare objects. If objects differ, return always greater (+1). Otherwise compare col and
     * row of the object. Col has priority with respect to row.
     * @param ref The object ot comparte to
     * @return 0 if objects are equal, -1 if this object is smaller, +1 if this object is larger.
     */
    bool operator<(const TriggerChannelPosition &other) const;

  private:
    int                     fCol;            ///< Column of the trigger channel
    int                   fRow;            ///< Row of the trigger channel
  };

  /**
   * Constructor
   */
  AliEmcalTriggerPartBadChannelContainer(): TObject(), fChannels() {}

  /**
   * Destructor, cleans up the container
   */
  virtual ~AliEmcalTriggerPartBadChannelContainer() {}

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

  /**
   * Get Channels
   */
  std::vector<PWG::EMCAL::TriggerPart::AliEmcalTriggerPartBadChannelContainer::TriggerChannelPosition> GetChannels() const
  {
    return fChannels;
  }

private:
  std::vector<PWG::EMCAL::TriggerPart::AliEmcalTriggerPartBadChannelContainer::TriggerChannelPosition>             fChannels;      ///< Container for listed channels

  ClassDef(AliEmcalTriggerPartBadChannelContainer, 1);
};

}
}
}

#endif
