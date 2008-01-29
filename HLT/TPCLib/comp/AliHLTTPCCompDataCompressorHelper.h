// XEmacs -*-C++-*-
// $Id$

/**************************************************************************
 * TPCCompModelAnalysisright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Anders Vestbo <mailto:vestbo@fi.uib.no>                       *
 *          for The ALICE Off-line Project.                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCCompDataCompressorHelper.h
    @author Anders Vestbo
    @date   30-11-2006
    @brief  The Data Compressor helper for the Vestbo-compression for the TPC
*/

#ifndef AliHLTTPCComp_DataCompressorHelper
#define AliHLTTPCComp_DataCompressorHelper

#include "AliHLTTPCRootTypes.h"

/** @class   AliHLTTPCCompDataCompressorHelper
    @author Anders Vestbo
    @date   30-11-2006
    @brief  Class to implement compression functions that will be used in Vestbo-compression  
*/
class AliHLTTPCCompDataCompressorHelper {
  
 public:

  /** function to set numbers of bits for compression model 
   * @param pad    integer to define pad
   * @param time   integer to define time bin
   * @param charge integer to define total charge of the cluster
   * @param shape  integer to define shape of the cluster
   */
  static void SetBitNumbers(Int_t pad,Int_t time,Int_t charge,Int_t shape);

  /** function to define transverse resolutions
   * @param res1  float to define first resolution
   * @param res2  float to define second resolution
   * @param res3  float to define third resolution
   * @param width float to define width (set to 0.005 if nothing is specified)
   */
  static void SetTransverseResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width=0.005);

  /** function to define longitudinal resolutions
   * @param res1  float to define first resolution
   * @param res2  float to define second resolution
   * @param res3  float to define third resolution
   * @param width float to define width (set to 0.005 if nothing is specified)
   */
  static void SetLongitudinalResolutions(Float_t res1,Float_t res2,Float_t res3,Float_t width=0.005);

  /** function to set number of bits for remaining clusters
   * @param pad   integer to define pad
   * @param time  integer to define time bin
   * @param shape integer to define shape of the cluster
   */ 
  static void SetRemainingBitNumbers(Int_t pad,Int_t time,Int_t shape);

  /** function to get number of pad bits
   * @return fgNumPadBits integer for number of pad bits
   */
  static Int_t GetNPadBits() {return fgNumPadBits;}

  /** function to get number of time bits
   * @return fgNumTimeBits integer for number of time bin bits
   */
  static Int_t GetNTimeBits() {return fgNumTimeBits;}

  /** function to get number of charge bits (total charge)
   * @return fgNumChargeBits integer for number of charge bits
   */
  static Int_t GetNChargeBits() {return fgNumChargeBits;}

  /** function to get number of shape bits
   * @return fgNumShapeBits integer for number of shape bits
   */
  static Int_t GetNShapeBits() {return fgNumShapeBits;}

  /** function to get detector xy width step
   * @return fgXYWidthStep float for the x-y-width-step
   */
  static Float_t GetXYWidthStep() {return fgXYWidthStep;}

  /** function to get detector z width step
   * @return fgZWidthStep float for the z-width-step
   */
  static Float_t GetZWidthStep() {return fgZWidthStep;}

  /** function to get the total charge of the cluster
   * @return fgClusterCharge integer for the total charge of the cluster
   */
  static Int_t GetClusterCharge() {return fgClusterCharge;}

  /** function to get detector xy residual step
   * @param row integer to define the pad row 
   * @return fgXYResidualStepX float for the x-y-residual-step (X stands for 1 to 3, depending on the step)
   */
  static Float_t GetXYResidualStep(Int_t row);

  /** function to get detector z residual step
   * @param row integer to define the pad row 
   * @return fgZResidualStepX float for the z-residual-step (X stands for 1 to 3, depending on the step)
   */
  static Float_t GetZResidualStep(Int_t row);

  /** function to get the number of bits for pads from remaining clusters
   * @return fgNumPadBitsRemaining integer for the number of bits for pads from remaining clusters
   */
  static Int_t GetNPadBitsRemaining() {return fgNumPadBitsRemaining;}

  /** function to get the number of bits for time bins from remaining clusters
   * @return fgNumTimeBitsRemaining integer for the number of bits for time bins from remaining clusters
   */
  static Int_t GetNTimeBitsRemaining() {return fgNumTimeBitsRemaining;}

  /** function to get the number of bits for shapes of remaining clusters
   * @return fgNumShapeBitsRemaining integer for the number of bits for shapes of remaining clusters
   */
  static Int_t GetNShapeBitsRemaining() {return fgNumShapeBitsRemaining;}

  /** function to get pad precision factor
   * @return float that defines the number of bits (-> precision)
   */
  static Float_t GetPadPrecisionFactor();

  /** function to get time precision factor
   * @return float that defines the number of bits (-> precision)
   */
  static Float_t GetTimePrecisionFactor();

  //taken from TMath:

  /** function to convert double to int
   * @param x double value to be converted
   * @return integer created from double value
   */
  static Int_t Nint(Double_t x); 

  /** function to calculate the absolute of an integer
   * @param d integer to be taken absolute value from
   * @return positive integer
   */
  static Int_t Abs(Int_t d) { return (d > 0) ? d : -d; }

  /** function to calculate the absolute of a double
   * @param d double to be taken absolute value from
   * @return positive double
   */
  static Double_t Abs(Double_t d) { return (d > 0) ? d : -d; }

 private:

  /** standard constructor */
  AliHLTTPCCompDataCompressorHelper();

  /** standard destructor */
  virtual ~AliHLTTPCCompDataCompressorHelper() {};

  /** number of pad bits */
  static Int_t fgNumPadBits; // Number of pad bits
  /** number of time bits */
  static Int_t fgNumTimeBits; // Number of time bits
  /** number of charge bits */
  static Int_t fgNumChargeBits; // Number of charge bits
  /** number of shape bits */
  static Int_t fgNumShapeBits; // Number of shape bits
  /** number of remaining pad bits */
  static Int_t fgNumPadBitsRemaining; // Number of remaining pad bits
  /** number of remaining time bits */
  static Int_t fgNumTimeBitsRemaining; // Number of remaining time bits
  /** number of remaining shape bits */
  static Int_t fgNumShapeBitsRemaining; // Number of remaining shape bits

  /** xy residual at step 1 */
  static Float_t fgXYResidualStep1; // XY residual at step 1
  /** xy residual at step 2 */
  static Float_t fgXYResidualStep2; // XY residual at step 2
  /** xy residual at step 3 */
  static Float_t fgXYResidualStep3; // XY residual at step 3
  /** z residual at step 1 */
  static Float_t fgZResidualStep1; // Z residual at step 1
  /** z residual at step 2 */
  static Float_t fgZResidualStep2; // Z residual at step 2
  /** z residual at step 3 */
  static Float_t fgZResidualStep3; // Z residual at step 3
  /** width of xy step */
  static Float_t fgXYWidthStep; // Width of XY step
  /** width of z step */
  static Float_t fgZWidthStep;  // Width of Z step
  /** cluster charge */
  static Int_t fgClusterCharge; // Cluster charge


  ClassDef(AliHLTTPCCompDataCompressorHelper,1) 

};

#endif
