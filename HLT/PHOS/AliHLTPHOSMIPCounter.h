
 /**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTPHOSMIPCOUNTER_H
#define ALIHLTPHOSMIPCOUNTER_H

#include "AliHLTPHOSBase.h"

class AliHLTPHOSDigitContainerDataStruct;
class TH2I;

/**
 @author Øystein Djuvsland <oystein.djuvsland@gmail.com>
*/

class AliHLTPHOSMIPCounter : public AliHLTPHOSBase
  {
    public:
      /**
       * Default constructor
       */
      AliHLTPHOSMIPCounter();
      
      /**
       * Default destructor
       */
      ~AliHLTPHOSMIPCounter();

      /**
       * This functions counts the number of approved MIP hits in one digit container
       * @param digitContainer ContainerStruct for digits
       * @return Number of approved MIP hits
       */
      Int_t CountMIPs(AliHLTPHOSDigitContainerDataStruct* digitContainer);

  
      /**
       * Sets the upper bound of the amplitude accepted for a MIP
       * @param val Upper bound value
       */
      void SetUpperBound ( const Float_t& val ) { fUpperBound = val; }
      
      /**
       * Gets the upper bound of the amplitude accepted for a MIP
       * @return Upper bound value
       */
      Float_t GetUpperBound() const { return fUpperBound; }
      
      /**
       * Sets the lower bound of the amplitude accepted for a MIP
       * @param val lower bound value
       */
      void SetLowerBound ( const Float_t& val ) { fLowerBound = val; }
      
      /**
       * Gets the lower bound of the amplitude accepted for a MIP
       * @return Lower bound value
       */
      Float_t GetLowerBound() const { return fLowerBound; }
      
      /**
       * Sets the upper bound for the start time of the signal, i.e. where the
       * baseline shouldn't be zero anymore
       * @param val Upper start time
       */
      void SetUpperStartTime ( const Float_t& val ) { fUpperStartTime = val; }
      
      /**
       * Sets the upper bound for the start time of the signal, i.e. where the 
       * baseline shouldn't be zero anymore
       * @return Upper start time
       */
      Float_t GetUpperStartTime() const { return fUpperStartTime; }
      
      /**
       * Sets the lower bound for the start time of the signal, i.e. where the 
       * baseline shouldn't be zero anymore
       * @param val Lower start time
      */
      void SetLowerStartTime ( const Float_t& val ) { fLowerStartTime = val; }
      
      /**
       * Gets the lower bound for the start time of the signal, i.e. where the 
       * baseline shouldn't be zero anymore
       * @return Lower start time
       */
      Float_t GetLowerStartTime() const { return fLowerStartTime; }
      
      /**
       * Sets the threshold for what is zero when checking the baseline
       * @param val Zero threshold
       */
      void SetZeroThreshold ( const Float_t& val ) { fZeroThreshold = val; } 
      
      /**
       * Gets the threshold for what is zero when checking the baseline
       * @return Zero threshold
       */
      Float_t GetZeroThreshold() const { return fZeroThreshold; }
      
      /**
       * Get the total count of MIPs which the MIP counter has seen
       * @return Total count of MIPs
       */
      Int_t GetMIPCountTotal() const { return fMIPCountTotal; }
      
      /**
       * Get the total count of MIPs in the last event
       * @return Total count of MIPs in last event
       */
      Int_t GetMIPCountEvent() const { return fMIPCountEvent; }
      
      /**
       * Sets a pointer to the TH2I histogram of channels with MIP hits
       * @param hist Pointer to the histogram
       */
      void SetChannelHistogram(TH2I *hist) { fChannelHistPtr = hist; }
      
      /**
       * Get the pointer to the TH2I histogram of channels with MIP hits
       * @return Pointer to the histogram
       */
      TH2I* GetChannelHistogram() { return fChannelHistPtr; }
      

    private:
      Int_t fMIPCountEvent;
      Int_t fMIPCountTotal;
      Float_t fMIPRate;
      Float_t fLowerBound;
      Float_t fUpperBound;
      Float_t fUpperStartTime;
      Float_t fLowerStartTime;
      Float_t fZeroThreshold;
      TH2I *fChannelHistPtr;
  };

#endif
