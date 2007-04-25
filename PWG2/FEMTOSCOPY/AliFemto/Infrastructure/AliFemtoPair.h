/***************************************************************************
 *
 * $Id: AliFemtoPair.h,v 1.17
 *
 * Author: Brian Laziuk, Yale University
 *         slightly modified by Mike Lisa
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *    the Pair object is passed to the PairCuts for verification, and
 *    then to the AddRealPair and AddMixedPair methods of the
 *    Correlation Functions
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.19  2003/01/14 09:44:00  renault
 * corrections on average separation calculation for tracks which doesn't cross
 * all 45 padrows.
 *
 * Revision 1.18  2002/11/19 23:33:18  renault
 * Enable average separation calculation for all combinaisons of
 * V0 daughters and tracks
 *
 * Revision 1.16  2002/02/28 14:18:36  rcwells
 * Added emissionAngle function to AliFemtoPair
 *
 * Revision 1.15  2001/12/14 23:11:30  fretiere
 * Add class HitMergingCut. Add class fabricesPairCut = HitMerginCut + pair purity cuts. Add TpcLocalTransform function which convert to local tpc coord (not pretty). Modify AliFemtoTrack, AliFemtoParticle, AliFemtoHiddenInfo, AliFemtoPair to handle the hit information and cope with my code
 *
 * Revision 1.14  2001/04/03 21:04:36  kisiel
 *
 *
 *   Changes needed to make the Theoretical code
 *   work. The main code is the ThCorrFctn directory.
 *   The most visible change is the addition of the
 *   HiddenInfo to AliFemtoPair.
 *
 * Revision 1.13  2001/03/28 22:35:23  flierl
 * changes and bugfixes in qYKP*
 * add pairrapidity
 *
 * Revision 1.12  2001/01/22 22:56:40  laue
 * Yano-Koonin-Podgoretskii Parametrisation added
 *
 * Revision 1.11  2000/10/26 16:09:16  lisa
 * Added OpeningAngle PairCut class and method to AliFemtoPair
 *
 * Revision 1.10  2000/10/05 23:09:05  lisa
 * Added kT-dependent radii to mixed-event simulator AND implemented AverageSeparation Cut and CorrFctn
 *
 * Revision 1.9  2000/07/17 20:03:17  lisa
 * Implemented tools for addressing and assessing trackmerging
 *
 * Revision 1.8  2000/04/04 16:13:09  lisa
 * AliFemtoPair:quality() now returns normalized value (and so is double) and add a CorrFctn which looks at quality()
 *
 * Revision 1.7  2000/04/03 22:09:12  rcwells
 * Add member function ... quality().
 *
 * Revision 1.6  2000/02/13 21:13:34  lisa
 * changed ambiguous AliFemtoPair::fourMomentum() to fourMomentumSum() and fourMomentumDiff() and fixed related bug in QvecCorrFctn
 *
 * Revision 1.5  2000/01/25 17:35:17  laue
 * I. In order to run the stand alone version of the AliFemtoMaker the following
 * changes have been done:
 * a) all ClassDefs and ClassImps have been put into #ifdef __ROOT__ statements
 * b) unnecessary includes of StMaker.h have been removed
 * c) the subdirectory AliFemtoMaker/doc/Make has been created including everything
 * needed for the stand alone version
 *
 * II. To reduce the amount of compiler warning
 * a) some variables have been type casted
 * b) some destructors have been declared as virtual
 *
 * Revision 1.4  1999/07/29 16:16:34  lisa
 * Selemons upgrade of AliFemtoPair class
 *
 * Revision 1.3  1999/07/22 18:49:10  lisa
 * Implement idea of Fabrice to not create and delete AliFemtoPair all the time
 *
 * Revision 1.2  1999/07/06 22:33:22  lisa
 * Adjusted all to work in pro and new - dev itself is broken
 *
 * Revision 1.1.1.1  1999/06/29 16:02:57  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/

#ifndef ST_HBT_PAIR_HH
#define ST_HBT_PAIR_HH

#include <utility>

#include "Infrastructure/AliFemtoParticle.h"
#include "Infrastructure/AliFemtoTypes.h"

class AliFemtoPair {
public:
  AliFemtoPair();
  AliFemtoPair(AliFemtoParticle*, AliFemtoParticle*);
  

  ~AliFemtoPair();
  //AliFemtoPair(const AliFemtoPair&);
  //AliFemtoPair& operator=(const AliFemtoPair&);

  // track Gets:
  AliFemtoParticle* track1() const;
  AliFemtoParticle* track2() const;
  // track Sets:
  void SetTrack1(const AliFemtoParticle* trkPtr);
  void SetTrack2(const AliFemtoParticle* trkPtr);

  AliFemtoLorentzVector fourMomentumDiff() const;
  AliFemtoLorentzVector fourMomentumSum() const;
  double qInv() const;
  double kT()   const;
  double mInv() const;
  // pair rapidity
  double rap() const;
  double emissionAngle() const;

  // Bertsch-Pratt momentum components in Pair Frame - written by Bekele/Humanic
  double qSidePf() const;
  double qOutPf() const;
  double qLongPf() const;
   
  // Bertsch-Pratt momentum components in Local CMS (longitudinally comoving) frame
  // - written by Bekele/Humanic
  double qSideCMS() const;
  double qOutCMS() const;
  double qLongCMS() const;

  double dKSide() const;
  double dKOut() const;
  double dKLong() const;

  // Bertsch-Pratt momentum components in a longitudinally boosted frame
  // the argument is the beta of the longitudinal boost (default is 0.0, meaning lab frame)
  // - written by Bekele/Humanic
  double qSideBf(double beta=0.0) const;
  double qOutBf(double beta=0.0) const;
  double qLongBf(double beta=0.0) const;

  // Yano-Koonin-Podgoretskii Parametrisation 
  // source rest frame (usually lab frame)
  void qYKPCMS(double& qP, double& qT, double& q0) const ;
  // longitudinal comoving frame
  void qYKPLCMS(double& qP, double& qT, double& q0) const ;
  // pair rest frame
  void qYKPPF(double& qP, double& qT, double& q0) const ;


  double quality() const;

  // the following two methods calculate the "nominal" separation of the tracks 
  // at the inner field cage (EntranceSeparation) and when they exit the TPC,
  // which may be at the outer field cage, or at the endcaps.
  // "nominal" means that the tracks are assumed to start at (0,0,0).  Making this
  // assumption is important for the Event Mixing-- it is not a mistake. - MALisa
  double NominalTpcExitSeparation() const;
  double NominalTpcEntranceSeparation() const;
  double NominalTpcAverageSeparation() const;
  // adapted calculation of Entrance/Exit/Average Tpc separation to V0 daughters
  double TpcExitSeparationTrackV0Pos() const;
  double TpcEntranceSeparationTrackV0Pos() const;
  double TpcAverageSeparationTrackV0Pos() const; 

  double TpcExitSeparationTrackV0Neg() const;
  double TpcEntranceSeparationTrackV0Neg() const;
  double TpcAverageSeparationTrackV0Neg() const; 

  double TpcExitSeparationV0PosV0Pos() const;
  double TpcEntranceSeparationV0PosV0Pos() const;
  double TpcAverageSeparationV0PosV0Pos() const; 

  double TpcExitSeparationV0PosV0Neg() const;
  double TpcEntranceSeparationV0PosV0Neg() const;
  double TpcAverageSeparationV0PosV0Neg() const; 
 
  double TpcExitSeparationV0NegV0Pos() const;
  double TpcEntranceSeparationV0NegV0Pos() const;
  double TpcAverageSeparationV0NegV0Pos() const; 
  
  double TpcExitSeparationV0NegV0Neg() const;
  double TpcEntranceSeparationV0NegV0Neg() const;
  double TpcAverageSeparationV0NegV0Neg() const; 

  double pInv() const;
  double KStar() const;
  double KStarFlipped() const;
  double CVK() const;
  double CVKFlipped() const;
  double qInvFlippedXY() const;

  double OpeningAngle() const;

  // Fabrice Private <<<
  double KStarSide() const;
  double KStarOut() const;
  double KStarLong() const;

  float PionPairProbability() const;
  float ElectronPairProbability() const;
  float KaonPairProbability() const;
  float ProtonPairProbability() const;
  float KaonPionPairProbability() const;

  double dcaInsideTpc() const;
  double quality2() const;

 /* double KStarGlobal() const;
  double CVKGlobal() const;
  double KStarSideGlobal() const;
  double KStarOutGlobal() const;
  double KStarLongGlobal() const;*/

  void setMergingPar(double aMaxDuInner, double aMaxDzInner,
		     double aMaxDuOuter, double aMaxDzOuter);
  void setDefaultHalfFieldMergingPar();
  void setDefaultFullFieldMergingPar();
  double getFracOfMergedRow() const;
  double getClosestRowAtDCA() const;
  double getWeightedAvSep() const;
  // >>>
  double getFracOfMergedRowTrkV0Pos() const;
  double getClosestRowAtDCATrkV0Pos() const;

  double getFracOfMergedRowTrkV0Neg() const;
  double getClosestRowAtDCATrkV0Neg() const;

  double getFracOfMergedRowV0PosV0Neg() const;
  double getFracOfMergedRowV0NegV0Pos() const;
  double getFracOfMergedRowV0PosV0Pos() const;
  double getFracOfMergedRowV0NegV0Neg() const;

private:
  AliFemtoParticle* fTrack1;
  AliFemtoParticle* fTrack2;

  mutable short fNonIdParNotCalculated;
  mutable double fDKSide;
  mutable double fDKOut;
  mutable double fDKLong;
  mutable double fCVK;
  mutable double kStarCalc;
  void calcNonIdPar() const;

  mutable short fNonIdParNotCalculatedGlobal;
 /* mutable double fDKSideGlobal;
  mutable double fDKOutGlobal;
  mutable double fDKLongGlobal;
  mutable double kStarCalcGlobal;
  mutable double fCVKGlobal;*/
  //void calcNonIdParGlobal() const;

  mutable short fMergingParNotCalculated;
  mutable double fWeightedAvSep;
  mutable double fFracOfMergedRow;
  mutable double fClosestRowAtDCA;

  mutable short fMergingParNotCalculatedTrkV0Pos;
  mutable double fFracOfMergedRowTrkV0Pos;
  mutable double fClosestRowAtDCATrkV0Pos;

  mutable short fMergingParNotCalculatedTrkV0Neg;
  mutable double fFracOfMergedRowTrkV0Neg;
  mutable double fClosestRowAtDCATrkV0Neg;

  mutable short fMergingParNotCalculatedV0PosV0Neg;
  mutable double fFracOfMergedRowV0PosV0Neg;
  mutable double fClosestRowAtDCAV0PosV0Neg;

  mutable short fMergingParNotCalculatedV0NegV0Pos;
  mutable double fFracOfMergedRowV0NegV0Pos;
  mutable double fClosestRowAtDCAV0NegV0Pos;

  mutable short fMergingParNotCalculatedV0PosV0Pos;
  mutable double fFracOfMergedRowV0PosV0Pos;
  mutable double fClosestRowAtDCAV0PosV0Pos;

  mutable short fMergingParNotCalculatedV0NegV0Neg;
  mutable double fFracOfMergedRowV0NegV0Neg;
  mutable double fClosestRowAtDCAV0NegV0Neg;

  static double fMaxDuInner;
  static double fMaxDzInner;
  static double fMaxDuOuter;
  static double fMaxDzOuter;
  void calcMergingPar() const;

  void CalcMergingParFctn(short* tmpMergingParNotCalculatedFctn,
			  float* tmpZ1,float* tmpU1,
			  float* tmpZ2,float* tmpU2,
			  int *tmpSect1,int *tmpSect2,
			  double* tmpFracOfMergedRow,
			  double* tmpClosestRowAtDCA
			  ) const;

  void resetParCalculated();
};

inline void AliFemtoPair::resetParCalculated(){
  fNonIdParNotCalculated=1;
  fNonIdParNotCalculatedGlobal=1;
  fMergingParNotCalculated=1;
  fMergingParNotCalculatedTrkV0Pos=1;
  fMergingParNotCalculatedTrkV0Neg=1;
  fMergingParNotCalculatedV0PosV0Pos=1;
  fMergingParNotCalculatedV0NegV0Pos=1;
  fMergingParNotCalculatedV0PosV0Neg=1;
  fMergingParNotCalculatedV0NegV0Neg=1;
}

inline void AliFemtoPair::SetTrack1(const AliFemtoParticle* trkPtr){
  fTrack1=(AliFemtoParticle*)trkPtr;
  resetParCalculated();
}
inline void AliFemtoPair::SetTrack2(const AliFemtoParticle* trkPtr){
  fTrack2=(AliFemtoParticle*)trkPtr;
  resetParCalculated();
}

inline AliFemtoParticle* AliFemtoPair::track1() const {return fTrack1;}
inline AliFemtoParticle* AliFemtoPair::track2() const {return fTrack2;}

inline double AliFemtoPair::dKSide() const{
  if(fNonIdParNotCalculated) calcNonIdPar();
  return fDKSide;
}
inline double AliFemtoPair::dKOut() const{
  if(fNonIdParNotCalculated) calcNonIdPar();
  return fDKOut;
}
inline double AliFemtoPair::dKLong() const{
  if(fNonIdParNotCalculated) calcNonIdPar();
  return fDKLong;
}
inline double AliFemtoPair::KStar() const{
  if(fNonIdParNotCalculated) calcNonIdPar();
  return kStarCalc;
}
inline double AliFemtoPair::qInv() const {
  AliFemtoLorentzVector tDiff = (fTrack1->FourMomentum()-fTrack2->FourMomentum());
  return ( -1.* tDiff.m());
}

// Fabrice private <<<
inline double AliFemtoPair::KStarSide() const{
  if(fNonIdParNotCalculated) calcNonIdPar();
  return fDKSide;//mKStarSide;
}
inline double AliFemtoPair::KStarOut() const{
  if(fNonIdParNotCalculated) calcNonIdPar();
  return fDKOut;//mKStarOut;
}
inline double AliFemtoPair::KStarLong() const{
  if(fNonIdParNotCalculated) calcNonIdPar();
  return fDKLong;//mKStarLong;
}
inline double AliFemtoPair::CVK() const{
  if(fNonIdParNotCalculated) calcNonIdPar();
  return fCVK;
}

/*inline double AliFemtoPair::KStarGlobal() const{
  if(fNonIdParNotCalculatedGlobal) calcNonIdParGlobal();
  return kStarCalcGlobal;
}
inline double AliFemtoPair::KStarSideGlobal() const{
  if(fNonIdParNotCalculatedGlobal) calcNonIdParGlobal();
  return fDKSideGlobal;//mKStarSide;
}
inline double AliFemtoPair::KStarOutGlobal() const{
  if(fNonIdParNotCalculatedGlobal) calcNonIdParGlobal();
  return fDKOutGlobal;//mKStarOut;
}
inline double AliFemtoPair::KStarLongGlobal() const{
  if(fNonIdParNotCalculatedGlobal) calcNonIdParGlobal();
  return fDKLongGlobal;//mKStarLong;
}
inline double AliFemtoPair::CVKGlobal() const{
  if(fNonIdParNotCalculatedGlobal) calcNonIdParGlobal();
  return fCVKGlobal;
}*/


inline float AliFemtoPair::PionPairProbability() const{
  return (fTrack1->Track()->PidProbPion()) * 
         (fTrack2->Track()->PidProbPion());
}
inline float AliFemtoPair::ElectronPairProbability() const{
  return (fTrack1->Track()->PidProbElectron()) * 
         (fTrack2->Track()->PidProbElectron());
}
inline float AliFemtoPair::KaonPairProbability() const{
  return (fTrack1->Track()->PidProbKaon()) * 
         (fTrack2->Track()->PidProbKaon());
}
inline float AliFemtoPair::ProtonPairProbability() const{
  return (fTrack1->Track()->PidProbProton()) * 
         (fTrack2->Track()->PidProbProton());
}
inline float AliFemtoPair::KaonPionPairProbability() const{
  return (fTrack1->Track()->PidProbKaon()) * 
         (fTrack2->Track()->PidProbPion());
}

inline double AliFemtoPair::getFracOfMergedRow() const{
  if(fMergingParNotCalculated) calcMergingPar();
  return fFracOfMergedRow;
}
inline double AliFemtoPair::getClosestRowAtDCA() const { 
  if(fMergingParNotCalculated) calcMergingPar();
  return fClosestRowAtDCA;
}
inline double AliFemtoPair::getWeightedAvSep() const {
  if(fMergingParNotCalculated) calcMergingPar();
  return fWeightedAvSep;
}


inline double AliFemtoPair::getFracOfMergedRowTrkV0Pos() const{
  if(fMergingParNotCalculatedTrkV0Pos)
    CalcMergingParFctn(&fMergingParNotCalculatedTrkV0Pos,
		       &(fTrack1->fZ[0]),&(fTrack1->fU[0]),
		       &(fTrack2->fZ[0]),&(fTrack2->fU[0]),
		       &(fTrack1->fSect[0]),&(fTrack2->fSect[0]),
		       &(fFracOfMergedRowTrkV0Pos),&(fClosestRowAtDCATrkV0Pos)
		       );
  return fFracOfMergedRowTrkV0Pos;
}
inline double AliFemtoPair::getClosestRowAtDCATrkV0Pos() const{
  if(fMergingParNotCalculatedTrkV0Pos)
    CalcMergingParFctn(&fMergingParNotCalculatedTrkV0Pos,
		       &(fTrack1->fZ[0]),&(fTrack1->fU[0]),
		       &(fTrack2->fZ[0]),&(fTrack2->fU[0]),
		       &(fTrack1->fSect[0]),&(fTrack2->fSect[0]),
		       &fFracOfMergedRowTrkV0Pos,&fClosestRowAtDCATrkV0Pos
		       );
  return fClosestRowAtDCATrkV0Pos;
}
inline double AliFemtoPair::getFracOfMergedRowTrkV0Neg() const{
  if(fMergingParNotCalculatedTrkV0Neg)
    CalcMergingParFctn(&fMergingParNotCalculatedTrkV0Neg,
		       &(fTrack1->fZ[0]),&(fTrack1->fU[0]),
		       &(fTrack2->fV0NegZ[0]),&(fTrack2->fV0NegU[0]),
		       &(fTrack1->fSect[0]),&(fTrack2->fV0NegSect[0]),
		       &(fFracOfMergedRowTrkV0Neg),&(fClosestRowAtDCATrkV0Neg)
		       );
  return fFracOfMergedRowTrkV0Neg;
}
inline double AliFemtoPair::getClosestRowAtDCATrkV0Neg() const{
  if(fMergingParNotCalculatedTrkV0Neg)
    CalcMergingParFctn(&fMergingParNotCalculatedTrkV0Neg,
		       &(fTrack1->fZ[0]),&(fTrack1->fU[0]),
		       &(fTrack2->fV0NegZ[0]),&(fTrack2->fV0NegU[0]),
		       &(fTrack1->fSect[0]),&(fTrack2->fV0NegSect[0]),
		       &fFracOfMergedRowTrkV0Neg,&fClosestRowAtDCATrkV0Neg
		       );
  return fClosestRowAtDCATrkV0Neg;
}
inline double AliFemtoPair::getFracOfMergedRowV0PosV0Neg() const{
  if(fMergingParNotCalculatedV0PosV0Neg)
    CalcMergingParFctn(&fMergingParNotCalculatedV0PosV0Neg,
		       &(fTrack1->fZ[0]),&(fTrack1->fU[0]),
		       &(fTrack2->fV0NegZ[0]),&(fTrack2->fV0NegU[0]),
		       &(fTrack1->fSect[0]),&(fTrack2->fV0NegSect[0]),
		       &(fFracOfMergedRowV0PosV0Neg),
		       &(fClosestRowAtDCAV0PosV0Neg)
		       );
  return fFracOfMergedRowV0PosV0Neg;
}
inline double AliFemtoPair::getFracOfMergedRowV0NegV0Pos() const{
  if(fMergingParNotCalculatedV0NegV0Pos)
    CalcMergingParFctn(&fMergingParNotCalculatedV0NegV0Pos,
		       &(fTrack1->fV0NegZ[0]),&(fTrack1->fV0NegU[0]),
		       &(fTrack2->fZ[0]),&(fTrack2->fU[0]),
		       &(fTrack1->fV0NegSect[0]),
		       &(fTrack2->fSect[0]),
		       &(fFracOfMergedRowV0NegV0Pos),
		       &(fClosestRowAtDCAV0NegV0Pos)
		       );
  return fFracOfMergedRowV0NegV0Pos;
}
inline double AliFemtoPair::getFracOfMergedRowV0PosV0Pos() const{
  if(fMergingParNotCalculatedV0PosV0Pos)
    CalcMergingParFctn(&fMergingParNotCalculatedV0PosV0Pos,
		       &(fTrack1->fZ[0]),&(fTrack1->fU[0]),
		       &(fTrack2->fZ[0]),&(fTrack2->fU[0]),
		       &(fTrack1->fSect[0]),
		       &(fTrack2->fSect[0]),
		       &(fFracOfMergedRowV0PosV0Pos),
		       &(fClosestRowAtDCAV0PosV0Pos)
		       );
  return fFracOfMergedRowV0PosV0Pos;
}
inline double AliFemtoPair::getFracOfMergedRowV0NegV0Neg() const{
  if(fMergingParNotCalculatedV0NegV0Neg)
    CalcMergingParFctn(&fMergingParNotCalculatedV0NegV0Neg,
		       &(fTrack1->fV0NegZ[0]),&(fTrack1->fV0NegU[0]),
		       &(fTrack2->fV0NegZ[0]),&(fTrack2->fV0NegU[0]),
		       &(fTrack1->fV0NegSect[0]),
		       &(fTrack2->fV0NegSect[0]),
		       &(fFracOfMergedRowV0NegV0Neg),
		       &(fClosestRowAtDCAV0NegV0Neg)
		       );
  return fFracOfMergedRowV0NegV0Neg;
}

#endif
