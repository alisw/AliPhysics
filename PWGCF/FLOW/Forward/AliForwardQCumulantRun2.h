//
// Helper class to calculate Q cumulant in forward & central regions
//
#ifndef AliForwardQCumulantRun2_H
#define AliForwardQCumulantRun2_H
/**
 * @file AliForwardQCumulantRun2.h
 * @author Freja Thoresen <freja.thoresen@cern.ch>
 *
 * @brief
 *
 * @ingroup pwgcf_forward_flow
 */
#include <TObject.h>
#include <TH2D.h>
#include <TH3D.h>
#include "TString.h"
#include "TList.h"
#include "TMath.h"
#include "TRandom.h"
#include <THn.h>
#include "TString.h"
#include "AliForwardSettings.h"
#include "TComplex.h"
/**
 * Class to handle cumulant calculations.
 */
class AliForwardQCumulantRun2
{
public:
  /*
  * Constructor
  */
  AliForwardQCumulantRun2();

  /**
   * Destructor
   */
  virtual ~AliForwardQCumulantRun2(){}


  AliForwardSettings fSettings;
  /**
   * Do cumulants calculations for current event with
   * centrality cent
   *
   * @param cent Event centrality
   */
  void CumulantsAccumulate(TH2D& dNdetadphi, TList* outputList, double cent,double vertexpos, TString detType);

  void saveEvent(TList* outputList, double cent, double vertexpos,UInt_t r);

  /**
   * Constants
   */
  enum {
    ktpcOnly = 128,        // TPC only tracks
    kphiAcceptanceBin = 21 // phi acceptance bin in the FMD histogram (dNdetadphi)
  };


  /**
   * Get the bin number of <<cos(nphi)>>
   *
   * @param n moment
   *
   * @return bin number
   */
  Int_t GetBinNumberCos(Int_t n = 0) const;


  /**
   * Get the bin number of <<sin(nphi)>>
   *
   * @param n moment
   *
   * @return bin number
   */
  Int_t GetBinNumberSin(Int_t n = 0) const;

  /**
   * Reset histograms
   */
  void reset();

  /**
   * Do 4' calculation
   *
   * @param dQnReA Real part of the Q vector
   * @param dQnImA Imaginary part of the Q vector
   * @param qnRe Real part of the q vector (vector made from POIs which are also RPs)
   * @param qnIm Imaginary part of the q vector (vector made from POIs which are also RPs)
   * @param pnRe Real part of the p vector (vector made from POIs)
   * @param pnIm Imaginary part of the p vector (vector made from POIs)
   * @param q2nRe Real part of the q vector in the second mode
   * @param q2nIm Imaginary part of the q vector in the second mode
   * @param dQ2nReA Real part of the Q vector in the second mode
   * @param dQ2nImA Imaginary part of the Q vector in the second mode
   * @param mq Multiplicity of the POIs which are also RPs
   * @param multA Weights of the RPs
   *
   * @return four prime
   */
  Double_t calcFourPrime(double dQnReA, double dQnImA, double qnRe,double qnIm, double pnRe,double pnIm, double q2nRe, double q2nIm, double dQ2nReA, double dQ2nImA, double mq, double multA) const {
    Double_t temp = (TMath::Power(dQnReA,2.)+TMath::Power(dQnImA,2.))*(pnRe*dQnReA+pnIm*dQnImA)
    - q2nRe*(TMath::Power(dQnReA,2.)-TMath::Power(dQnImA,2.))
    - 2.*q2nIm*dQnReA*dQnImA
    - pnRe*(dQnReA*dQ2nReA+dQnImA*dQ2nImA)
    + pnIm*(dQnImA*dQ2nReA-dQnReA*dQ2nImA)
    - 2.*multA*(pnRe*dQnReA+pnIm*dQnImA)
    - 2.*(TMath::Power(dQnReA,2.)+TMath::Power(dQnImA,2.))*mq
    + 6.*(qnRe*dQnReA+qnIm*dQnImA)
    + 1.*(q2nRe*dQ2nReA+q2nIm*dQ2nImA)
    + 2.*(pnRe*dQnReA+pnIm*dQnImA)
    + 2.*mq*multA
    - 6.*mq;
    return temp;
  };


  /**
   * Calculate <cos(psi1 + phi2 + phi3)>
   *
   * @param pnRe Real part of the p vector (vector made from POIs)
   * @param dQnImA Imaginary part of the Q vector
   * @param dQnReA Real part of the Q vector
   * @param multA Weights of the RPs
   * @param q2nRe Real part of the q vector in the second mode
   * @param q2nIm Imaginary part of the q vector in the second mode
   * @param mq Multiplicity of the POIs which are also RPs
   * @param qnRe Real part of the q vector (vector made from POIs which are also RPs)
   *
   * @return <cos(psi1 + phi2 + phi3)>
   */
  Double_t calcCosPsi1Phi2Phi3p(double pnRe, double dQnImA, double dQnReA, double multA, double q2nRe, double q2nIm, double mq, double qnRe) const {
    Double_t temp;
    temp = pnRe*(TMath::Power(dQnImA,2.)+TMath::Power(dQnReA,2.)-multA)
    - 1.*(q2nRe*dQnReA+q2nIm*dQnImA)
    - mq*dQnReA+2.*qnRe;
    return temp;
  };

  /**
   * Calculate <sin(psi1 + phi2 + phi3)>
   *
   * @param pnIm Imaginary part of the p vector (vector made from POIs)
   * @param dQnImA Imaginary part of the Q vector
   * @param dQnReA Real part of the Q vector
   * @param multA Weights of the RPs
   * @param mq Multiplicity of the POIs which are also RPs
   * @param qnIm Imaginary part of the q vector (vector made from POIs which are also RPs)
   * @param q2nRe Real part of the q vector in the second mode
   * @param q2nIm Imaginary part of the q vector in the second mode
   *
   * @return <cos(psi1 + phi2 + phi3)>
   */
  Double_t calcSinPsi1Phi2Phi3p(double pnIm, double dQnImA, double dQnReA, double multA,double  mq, double qnIm, double q2nIm, double q2nRe) const {
    Double_t temp;
    temp = pnIm*(TMath::Power(dQnImA,2.)+TMath::Power(dQnReA,2.)-multA)
    - 1.*(q2nIm*dQnReA-q2nRe*dQnImA)
    - mq*dQnImA+2.*qnIm;
    return temp;
  };


  /**
   * Calculate <cos(psi1 - phi2 - phi3)>
   *
   * @param pnRe Real part of the p vector (vector made from POIs)
   * @param dQnReA Real part of the Q vector
   * @param dQnImA Imaginary part of the Q vector
   * @param pnIm Imaginary part of the p vector (vector made from POIs)
   * @param dQ2nReA Real part of the Q vector in the second mode
   * @param dQ2nImA Imaginary part of the Q vector in the second mode
   * @param mq Multiplicity of the POIs which are also RPs
   * @param qnRe Real part of the q vector (vector made from POIs which are also RPs)
   *
   * @return <cos(psi1 - phi2 - phi3)>
   */
  Double_t calcCosPsi1Phi2Phi3m(double pnRe, double dQnReA, double dQnImA, double pnIm,double  dQ2nReA,double  dQ2nImA,double mq, double qnRe) const {
    Double_t temp;
    temp = pnRe*(TMath::Power(dQnReA,2.)-TMath::Power(dQnImA,2.))+2.*pnIm*dQnReA*dQnImA
    - 1.*(pnRe*dQ2nReA+pnIm*dQ2nImA)
    - 2.*mq*dQnReA+2.*qnRe;
    return temp;
  };


  /**
   * Calculate <sin(psi1 - phi2 - phi3)>
   *
   * @param pnIm Imaginary part of the p vector (vector made from POIs)
   * @param dQnReA Real part of the Q vector
   * @param dQnImA Imaginary part of the Q vector
   * @param pnRe Real part of the p vector (vector made from POIs)
   * @param dQ2nReA Real part of the Q vector in the second mode
   * @param dQ2nImA Imaginary part of the Q vector in the second mode
   * @param mq Multiplicity of the POIs which are also RPs
   * @param qnIm Imaginary part of the q vector (vector made from POIs which are also RPs)
   *
   * @return <sin(psi1 - phi2 - phi3)>
   */
  Double_t calcSinPsi1Phi2Phi3m(double pnIm, double dQnReA, double dQnImA, double pnRe, double dQ2nReA,double dQ2nImA, double mq, double qnIm) const {
    Double_t temp;
    temp = pnIm*(TMath::Power(dQnReA,2.)-TMath::Power(dQnImA,2.))-2.*pnRe*dQnReA*dQnImA
    - 1.*(pnIm*dQ2nReA-pnRe*dQ2nImA)
    + 2.*mq*dQnImA-2.*qnIm;
    return temp;
  };


  /**
   * Calculate <cos(phi1 + phi2 + phi3)>
   *
   * @param dQnReA Real part of the Q vector
   * @param dQnImA Imaginary part of the Q vector
   * @param dQ2nReA Real part of the Q vector in the second mode
   * @param dQ2nImA Imaginary part of the Q vector in the second mode
   * @param multA Weights of the RPs
   *
   * @return <cos(phi1 + phi2 + phi3)>
   */
  Double_t calcCosPhi1Phi2Phi3m(double dQnReA, double dQnImA, double dQ2nReA, double dQ2nImA, double multA) const {
   return (dQnReA*(TMath::Power(dQnReA,2)+TMath::Power(dQnImA,2))
    -dQnReA*dQ2nReA-dQnImA*dQ2nImA-2.*(multA-1)*dQnReA);
  };


  /**
   * Calculate <sin(phi1 + phi2 + phi3)>
   *
   * @param dQnReA Real part of the Q vector
   * @param dQnImA Imaginary part of the Q vector
   * @param dQ2nReA Real part of the Q vector in the second mode
   * @param dQ2nImA Imaginary part of the Q vector in the second mode
   * @param multA Weights of the RPs
   *
   * @return <sin(phi1 + phi2 + phi3)>
   */
  Double_t calcSinPhi1Phi2Phi3m(double dQnReA, double dQnImA,double  dQ2nReA, double dQ2nImA,double  multA) const {
    return (-dQnImA*(TMath::Power(dQnReA,2)+TMath::Power(dQnImA,2))
      + dQnReA*dQ2nImA - dQnImA*dQ2nReA + 2.*(multA-1)*dQnImA);
  };



  Int_t fMaxMoment; // maximum mode of vn
  Int_t nHBins;
  Int_t nRefBins; // bins in accummulating reference particles in eta
  Int_t nEtaBins; // bins in accummulating differential particles in eta


  TH2D fCumuRef;     // Accumulated reference particles
  TH2D fCumuDiff;    // Accumulated differential particles
  bool useEvent;


  ClassDef(AliForwardQCumulantRun2, 1); // object for eta dependent cumulants ananlysis
};

#endif
