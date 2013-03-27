///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  AliFemtoEvent: hold the information specific to the event and a      //
//  track list                                                           //
//  AliFemtoEvent is the "transient microDST"  Objects of this class are //
//   generated from the input data by a Reader, and then presented to	 //
//   the Cuts of the various active Analyses.				 //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOEVENT_H
#define ALIFEMTOEVENT_H

#include "AliFemtoTypes.h"
#include "AliFemtoTrackCollection.h"
#include "AliFemtoV0Collection.h"
#include "AliFemtoXiCollection.h"
#include "AliFemtoKinkCollection.h"

class AliFemtoTrackCut;
class AliFemtoV0Cut;
class AliFemtoXiCut;
class AliFemtoKinkCut;
class AliEventplane;

#ifdef __ROOT__
// the following encapsulation by malisa 21apr2006
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
class StMuDst;
#endif
#endif

class AliFemtoEvent{
public:
  AliFemtoEvent();
#ifdef __ROOT__
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
//
#endif
#endif
  
  AliFemtoEvent(const AliFemtoEvent& ev, AliFemtoTrackCut* tCut=0, AliFemtoV0Cut* vCut=0,  AliFemtoXiCut* xCut=0, AliFemtoKinkCut* kCut=0); // copy constructor with track and v0 cuts
  AliFemtoEvent(const AliFemtoEvent& ev); // copy constructor
  ~AliFemtoEvent();
  AliFemtoEvent& operator=(const AliFemtoEvent& aEvent);

  unsigned short EventNumber() const;
  int RunNumber() const;
  unsigned short NumberOfTracks() const;
  AliFemtoThreeVector PrimVertPos() const;
  const double* PrimVertCov() const;
  AliFemtoV0Collection* V0Collection() const;
  AliFemtoXiCollection* XiCollection() const;
  AliFemtoKinkCollection* KinkCollection() const;
  AliFemtoTrackCollection* TrackCollection() const;
  double MagneticField() const;
  bool IsCollisionCandidate() const;

  //functions for alice variables
  float ZDCN1Energy() const;      
  float ZDCP1Energy() const;      
  float ZDCN2Energy() const;      
  float ZDCP2Energy() const;      
  float ZDCEMEnergy() const;    
  unsigned int ZDCParticipants() const; 

  unsigned long int     TriggerMask() const;     
  unsigned char      TriggerCluster() const;  

  float ReactionPlaneAngle() const;
  AliEventplane* EP() const;
  
  void SetEventNumber(const unsigned short& s);
  void SetRunNumber(const int& i);
  void SetNumberOfTracks(const unsigned short& s);
  void SetNormalizedMult(const int& i);
  void SetMultiplicityEstimateITSTPC(const unsigned short &s);
  void SetMultiplicityEstimateTracklets(const unsigned short &s);
  void SetMultiplicityEstimateITSPure(const unsigned short &s);


  void SetCentralityV0(const float &c);
  void SetCentralityV0A(const float &c);
  void SetCentralityV0C(const float &c);
  void SetCentralityZNA(const float &c);
  void SetCentralityZNC(const float &c);
  void SetCentralityCL1(const float &c);
  void SetCentralityCL0(const float &c);
  void SetCentralityTKL(const float &c);
  void SetCentralityFMD(const float &c);
  void SetCentralityTrk(const float &c);
  void SetCentralityCND(const float &c);
  void SetCentralityNPA(const float &c);
  void SetCentralitySPD1(const float &c);


  void SetSPDMult(const int& i);
  void SetPrimVertPos(const AliFemtoThreeVector& v);
  void SetPrimVertCov(const double* v);
  void SetMagneticField(const double& x);
  void SetIsCollisionCandidate(const bool& is);

   //functions for alice variables
  void SetZDCN1Energy(const float& x);      
  void SetZDCP1Energy(const float& x);      
  void SetZDCN2Energy(const float& x);      
  void SetZDCP2Energy(const float& x);      
  void SetZDCEMEnergy(const float& x);    
  void SetZDCParticipants(const unsigned int& i);
  
  void SetTriggerMask(const unsigned long int& i);     
  void SetTriggerCluster(const unsigned char& c); 

  void SetReactionPlaneAngle(const float& a);
  void SetEP(AliEventplane* ep);
  
  int UncorrectedNumberOfNegativePrimaries() const;
  int UncorrectedNumberOfPrimaries() const;
  int SPDMultiplicity() const;
  int NumberOfV0s() const;

  unsigned short MultiplicityEstimateITSTPC() const;
  unsigned short MultiplicityEstimateTracklets() const;
  unsigned short MultiplicityEstimateITSPure() const;


  float CentralityV0() const; 
  float CentralityV0A() const;
  float CentralityV0C() const;
  float CentralityZNA() const;
  float CentralityZNC() const;
  float CentralityCL1() const;
  float CentralityCL0() const;
  float CentralityTKL() const;
  float CentralityFMD() const;
  float CentralityTrk() const;
  float CentralityCND() const;
  float CentralityNPA() const;
  float CentralitySPD1() const;

private:
  unsigned short fEventNumber;           // Event number in file
  unsigned short fRunNumber;             // run number the event belong to
  unsigned short fNumberOfTracks;        // total number of TPC tracks
  int   fNormalizedMult;                 // normalized multiplicity
  int   fSPDMult;                        // Multiplicity of SPD tracklets
  unsigned short fEstimateITSTPC;        // Official multiplicity estimate ITS+TPC
  unsigned short fEstimateTracklets;     // Official multiplicity estimate Tracklets
  unsigned short fEstimateITSPure;       // Official multiplicity estimate ITS SA


  float fCentralityV0; 
  float fCentralityV0A;
  float fCentralityV0C;
  float fCentralityZNA;
  float fCentralityZNC;
  float fCentralityCL1;
  float fCentralityCL0;
  float fCentralityTKL;
  float fCentralityFMD;
  float fCentralityTrk;
  float fCentralityCND;
  float fCentralityNPA;
  float fCentralitySPD1;

  double fMagneticField;                 // magnetic field in Z direction
  bool fIsCollisionCandidate;            // is collision candidate
  
  AliFemtoThreeVector fPrimVertPos;      // primary vertex position
  double fPrimVertCov[6];                // primary vertex covariances
  AliFemtoTrackCollection* fTrackCollection; // collection of tracks
  AliFemtoV0Collection* fV0Collection;   // collection of V0s
  AliFemtoXiCollection* fXiCollection;   // collection of Xis
  AliFemtoKinkCollection* fKinkCollection; // collection of kinks

  //for alice changed by Marek Chojnacki
  float      fZDCN1Energy;      // reconstructed energy in the neutron ZDC
  float      fZDCP1Energy;      // reconstructed energy in the proton ZDC
  float      fZDCN2Energy;      // reconstructed energy in the neutron ZDC
  float      fZDCP2Energy;      // reconstructed energy in the proton ZDC
  float      fZDCEMEnergy;     // reconstructed energy in the electromagnetic ZDC
  unsigned int        fZDCParticipants; // number of participants estimated by the ZDC
  
  unsigned long int     fTriggerMask;     // Trigger Type (mask)
  unsigned char      fTriggerCluster;  // Trigger cluster (mask)

  float      fReactionPlaneAngle; // reconstructed reaction plane angle
  AliEventplane*  fEP; // pointer to full event plane information
};



#endif 
