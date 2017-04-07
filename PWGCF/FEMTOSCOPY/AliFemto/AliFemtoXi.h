///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoXi: special type of particle desling with the specifics       //
// of the Xi type of particle                                            //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOXI_H
#define ALIFEMTOXI_H

#include "AliFemtoVector.h" //same as in AliFemtoTrack.h
#include "AliFemtoV0.h"


class AliFemtoXi : public AliFemtoV0 {
public:
  AliFemtoXi();
  AliFemtoXi(const AliFemtoV0* aV0);
  virtual ~AliFemtoXi(){/* no-op */}

  AliFemtoXi(const AliFemtoXi& aXi);
  AliFemtoXi& operator=(const AliFemtoXi& aXi);

  void UpdateXi();
  float DecayLengthXi() const;            // 3-d decay distance
  AliFemtoThreeVector DecayVertexXi() const; // Coordinates of decay vertex
  float DecayVertexXiX() const;           // Coordinates of decay vertex
  float DecayVertexXiY() const;           // Coordinates of decay vertex
  float DecayVertexXiZ() const;           // Coordinates of decay vertex
  float DcaXiDaughters() const;           // DCA of xi daughters at decay vertex
  float DcaXiToPrimVertex() const;        // DCA of xi to primary vertex
  float DcaBacToPrimVertex() const;       // DCA of bachelor  xi daughter to pri vertex
  AliFemtoThreeVector MomBac() const;        // Momentum components of bac. daughter
  float MomBacX() const;                  // Momentum components of bac. daughter
  float MomBacY() const;                  // Momentum components of bac. daughter
  float MomBacZ() const;                  // Momentum components of bac. daughter

  int   TpcHitsBac() const;               // Number of TPC hits on bac. daughter
  unsigned long TrackTopologyMapBac(unsigned int i) const;

  AliFemtoThreeVector MomXi() const ;        // Momentum components of Xi
  float MomXiX() const ;                  // Momentum components of Xi
  float MomXiY() const ;                  // Momentum components of Xi
  float MomXiZ() const ;                  // Momentum components of Xi
  float AlphaXi() const ;                 // Armenteros-Podolanski variable
  float PtArmXi() const ;                 // Armenteros-Podolanski variable
  float EXi() const ;                     // Energy assuming xi hypothesis		     
  float EOmega() const ;                  // Energy assuming omega hypothesis		     
  float EBacKaon() const ;                // Energy of bac. daughter assuming kaon	     
  float EBacPion() const ;                // Energy of bac. daughter assuming pion	     
  float MassXi() const ;                  // Mass assuming Xi hypothesis		     
  float MassOmega() const ;               // Mass assuming Omega hypothesis		     
  float RapXi() const ;                   // Rapidity assuming (anti) xi		     
  float RapOmega() const ;                // Rapidity assuming (anti) omega		     
  float CTauXi() const ;                  // Lifetime (ctau) const assuming (anti) xi	     
  float CTauOmega() const ;               // Lifetime (ctau) const assuming (anti) omega     
  float PtXi() const ;                    // Transverse momentum			     
  float PtotXi() const ;                  // Total momentum				     
  float PtBac() const ;                   // Transverse momentum of bac. daughter	     
  float PtotBac() const ;                 // Total momentum of bac. daughter		     
  float DedxBac() const;                  // dedx of Bac track				     
  unsigned short   IdBac() const;         // Id of bac. track				     
  unsigned short   KeyBac() const;        // Id of bac. track                                

  void SetdecayLengthXi(const float x);  
  void SetdecayVertexXi(const AliFemtoThreeVector v);  
  void SetdecayVertexXiX(const float x);
  void SetdecayVertexXiY(const float x);
  void SetdecayVertexXiZ(const float x);
  void SetdcaXiDaughters(const float x); 
  void SetdcaXiToPrimVertex(const float x);  
  void SetdcaBacToPrimVertex(const float x); 
  void SetmomBac(const AliFemtoThreeVector v);  
  void SetmomBacX(const float x);  
  void SetmomBacY(const float x);  
  void SetmomBacZ(const float x);  

  void SettpcHitsBac(const int& word);      

  void SetTrackTopologyMapBac(unsigned int word, const unsigned long& m);

  void SetmomXi( AliFemtoThreeVector v);
  void SetmomXiX( float x);
  void SetmomXiY( float x);
  void SetmomXiZ( float x);
  void SetalphaXi( float x);       
  void SetptArmXi( float x);       
  void SeteXi( float x);     
  void SeteOmega( float x);    
  void SeteBacPion( float x);  
  void SeteBacKaon( float x);    
  void SetmassXi( float x);  
  void SetmassOmega( float x);
  void SetrapXi( float x);    
  void SetrapOmega( float x);   
  void SetcTauXi( float x);   
  void SetcTauOmega( float x);  
  void SetptXi( float x);         
  void SetptotXi( float x);       
  void SetptBac( float x);        
  void SetptotBac( float x);      
  void SetidBac(const unsigned short& s);
  void SetdedxBac(float x);
  void SetkeyBac(const unsigned short& s);



  //!!!!!!!!!!! new below 1.09.2015
  void SetCosPointingAngleXi(double x);
  void SetCosPointingAngleV0toXi(double x);
  void SetPhiXi(double x);
  void SetEtaXi(double x);

  void SetTPCNclsBac(int x);
  void SetNdofBac(int x);
  void SetStatusBac(unsigned long x);
  void SetEtaBac(float x);

  void SetBacNSigmaTPCK(float x);
  void SetBacNSigmaTPCPi(float x);
  void SetBacNSigmaTPCP(float x);

  void SetBacNSigmaTOFK(float x);
  void SetBacNSigmaTOFPi(float x);
  void SetBacNSigmaTOFP(float x);

  void SetTPCMomentumBac(double x);

  void SetTOFProtonTimeBac(double x);
  void SetTOFPionTimeBac(double x);
  void SetTOFKaonTimeBac(double x);

  void SetChargeXi(int x);



  double CosPointingAngleXi() const;
  double CosPointingAngleV0toXi() const;
  double EtaXi() const;
  double PhiXi() const;

  int TPCNclsBac() const;
  int NdofBac() const;
  unsigned long StatusBac() const;
  float EtaBac() const;

  double GetTPCMomentumBac() const;

  float BacNSigmaTPCK() const;
  float BacNSigmaTPCPi() const;
  float BacNSigmaTPCP() const;

  float BacNSigmaTOFK() const;
  float BacNSigmaTOFPi() const;
  float BacNSigmaTOFP() const;

  int ChargeXi() const;


protected: 
  int   fCharge;                           // Charge                
  float fDecayLengthXi;                    // Decay Length of th Xi
  AliFemtoThreeVector fDecayVertexXi;      // Xi decay vertex location
  float fDcaXiDaughters;                   // Dca of Xi daughters
  float fDcaXiToPrimVertex;                // Dca of Xi to primary vertex
  float fDcaBachelorToPrimVertex;          // Dca of bachelor to primary vertex
  AliFemtoThreeVector fMomBachelor;        // Momentum of the bachelor

  unsigned int   fTopologyMapBachelor[2];  // TPC topology map for the bachelor
  unsigned short fKeyBachelor;             // Unique key for the bachelor

  int   fTpcHitsBac;                       // Number of TPC hits for bachelor

  float fChi2Xi;                           // Fit quality for Xi
  float fClXi;				   // Confidence level for Xi
  float fChi2Bachelor;			   // Fit quality for bachelor
  float fClBachelor;			   // Confidence level for bachelor

  float fDedxBachelor;                     // dEdx for bachelor
  unsigned short fNufDedxBachelor;         // dEdx for bachelor

  // the following variables are not in the persistent version and can be calculated via UpdateXi();
  AliFemtoThreeVector fMomXi; // Momentum of the Xi
  float fAlphaXi;             // Armenteros-Podolanski variable
  float fPtArmXi;	      // Armenteros-Podolanski variable

  float fEXi;                 // Energy assuming xi hypothesis		     
  float fEOmega;	      // Energy assuming omega hypothesis		     
  float fEBacPion;	      // Energy of bac. daughter assuming kaon	     
  float fEBacKaon;	      // Energy of bac. daughter assuming pion	     
  float fMassXi;	      // Mass assuming Xi hypothesis		     
  float fMassOmega;	      // Mass assuming Omega hypothesis		     
  float fRapXi;		      // Rapidity assuming (anti) xi		     
  float fRapOmega;	      // Rapidity assuming (anti) omega		     
  float fCTauXi;	      // Lifetime (ctau) const assuming (anti) xi	     
  float fCTauOmega;	      // Lifetime (ctau) const assuming (anti) omega     
  float fPtXi;		      // Transverse momentum			     
  float fPtotXi;	      // Total momentum				     
  float fPtBac;		      // Transverse momentum of bac. daughter	     
  float fPtotBac;	      // Total momentum of bac. daughter		     

  unsigned short   fKeyBac;   // Key of bac. track


  //!!!!!!!!! new below 1.09.2015
  double fCosPointingAngleXi; // Cosine of Xi pointing angle
  double fCosPointingAngleV0toXi; //Cosing of V0 pointing angle to Xi decay vertex 
    //AliFemtoV0 has fCosPointingAngle, but this is the pointing angle to the primary vertex.
    //What is more important here is the V0 pointing angle to the Xi decay vertex
  double fPhiXi;  //Phi angle of Xi
  double fEtaXi;

  //bachelor track properties
  int fTPCNclsBac; //Bachelor number of TPC clusters
  int   fNdofBac;
  unsigned long fStatusBac; //Status (TPC refit, ITS refit, ...) of bachelor track
  float fEtaBac; //Eta of bachelor track
  int fIdBac; //id of bachelor track

  float fBacNSigmaTPCK;
  float fBacNSigmaTPCPi;
  float fBacNSigmaTPCP;

  float fBacNSigmaTOFK;
  float fBacNSigmaTOFPi;
  float fBacNSigmaTOFP;

  double fTPCMomentumBac;

  double fTOFProtonTimeBac;
  double fTOFPionTimeBac;
  double fTOFKaonTimeBac;
  
};

inline float AliFemtoXi::DecayLengthXi() const { return fDecayLengthXi; }
inline AliFemtoThreeVector AliFemtoXi::DecayVertexXi() const { return fDecayVertexXi; } 
inline float AliFemtoXi::DecayVertexXiX() const { return fDecayVertexXi.x(); } 
inline float AliFemtoXi::DecayVertexXiY() const { return fDecayVertexXi.y(); } 
inline float AliFemtoXi::DecayVertexXiZ() const { return fDecayVertexXi.z(); } 
inline float AliFemtoXi::DcaXiDaughters() const { return fDcaXiDaughters; }
inline float AliFemtoXi::DcaXiToPrimVertex() const { return fDcaXiToPrimVertex; }
inline float AliFemtoXi::DcaBacToPrimVertex() const { return fDcaBachelorToPrimVertex; }
inline AliFemtoThreeVector AliFemtoXi::MomBac() const { return fMomBachelor; }
inline float AliFemtoXi::MomBacX() const { return fMomBachelor.x(); }
inline float AliFemtoXi::MomBacY() const { return fMomBachelor.y(); }
inline float AliFemtoXi::MomBacZ() const { return fMomBachelor.z(); }
inline AliFemtoThreeVector AliFemtoXi::MomXi() const { return fMomXi; }
inline float AliFemtoXi::MomXiX() const { return fMomXi.x(); }
inline float AliFemtoXi::MomXiY() const { return fMomXi.y(); }
inline float AliFemtoXi::MomXiZ() const { return fMomXi.z(); }
inline float AliFemtoXi::AlphaXi() const { return fAlphaXi; }
inline float AliFemtoXi::PtArmXi() const {return fPtArmXi;}
inline float AliFemtoXi::EXi() const {return fEXi;}
inline float AliFemtoXi::EOmega() const {return fEOmega;}
inline float AliFemtoXi::EBacPion() const {return fEBacPion;}
inline float AliFemtoXi::EBacKaon() const {return fEBacKaon;}
inline float AliFemtoXi::MassXi() const {return fMassXi;}
inline float AliFemtoXi::MassOmega() const {return fMassOmega;}
inline float AliFemtoXi::RapXi() const {return fRapXi;}
inline float AliFemtoXi::RapOmega() const {return fRapOmega;}
inline float AliFemtoXi::CTauXi() const {return fCTauXi;}
inline float AliFemtoXi::CTauOmega() const {return fCTauOmega;}
inline float AliFemtoXi::PtXi() const {return fPtXi;}
inline float AliFemtoXi::PtotXi() const {return fPtotXi;}
inline float AliFemtoXi::PtBac() const {return fPtBac;}
inline float AliFemtoXi::PtotBac() const {return fPtotBac;}
inline int   AliFemtoXi::TpcHitsBac() const
             { return fTpcHitsBac; }
inline float AliFemtoXi::DedxBac() const {return fDedxBachelor;}

inline unsigned long    AliFemtoXi::TrackTopologyMapBac(unsigned int word) const { return fTopologyMapBachelor[word]; }
inline unsigned short   AliFemtoXi::IdBac() const { return fKeyBac; }
inline unsigned short   AliFemtoXi::KeyBac() const { return fKeyBac; }

inline void AliFemtoXi::SetdecayLengthXi(const float x){ fDecayLengthXi= x;}   
inline void AliFemtoXi::SetdecayVertexXiX(const float x){ fDecayVertexXi.SetX(x);}
inline void AliFemtoXi::SetdecayVertexXiY(const float x){ fDecayVertexXi.SetY(x);}
inline void AliFemtoXi::SetdecayVertexXiZ(const float x){ fDecayVertexXi.SetZ(x);}
inline void AliFemtoXi::SetdecayVertexXi(const AliFemtoThreeVector v){ fDecayVertexXi = v; }
inline void AliFemtoXi::SetdcaXiDaughters(const float x){fDcaXiDaughters= x;} 
inline void AliFemtoXi::SetdcaXiToPrimVertex(const float x){fDcaXiToPrimVertex= x;}   
inline void AliFemtoXi::SetdcaBacToPrimVertex(const float x){ fDcaBachelorToPrimVertex = x;} 
inline void AliFemtoXi::SetmomBac(const AliFemtoThreeVector v){fMomBachelor = v; }
inline void AliFemtoXi::SetmomBacX(const float x){fMomBachelor.SetX(x);}
inline void AliFemtoXi::SetmomBacY(const float x){fMomBachelor.SetY(x);}
inline void AliFemtoXi::SetmomBacZ(const float x){fMomBachelor.SetZ(x);}
inline void AliFemtoXi::SetTrackTopologyMapBac(unsigned int word, const unsigned long& m){fTopologyMapBachelor[word]=m;} 
inline void AliFemtoXi::SetmomXi(AliFemtoThreeVector v){fMomXi= v; }
inline void AliFemtoXi::SetmomXiX(const float x){fMomXi.SetX(x);}
inline void AliFemtoXi::SetmomXiY(const float x){fMomXi.SetY(x);}
inline void AliFemtoXi::SetmomXiZ(const float x){fMomXi.SetZ(x);}

inline void AliFemtoXi::SetalphaXi( float x){fAlphaXi= x;}
inline void AliFemtoXi::SetptArmXi( float x){fPtArmXi = x;}
inline void AliFemtoXi::SeteXi( float x){fEXi= x;}       
inline void AliFemtoXi::SeteOmega( float x){fEOmega= x;}
inline void AliFemtoXi::SeteBacPion( float x){fEBacPion= x;}
inline void AliFemtoXi::SeteBacKaon( float x){fEBacKaon= x;}
inline void AliFemtoXi::SetmassXi( float x){fMassXi = x;} 
inline void AliFemtoXi::SetmassOmega( float x){fMassOmega= x;}  
inline void AliFemtoXi::SetrapXi( float x){fRapXi= x;}
inline void AliFemtoXi::SetrapOmega( float x){fRapOmega = x;}   
inline void AliFemtoXi::SetcTauXi( float x){fCTauXi = x;}   
inline void AliFemtoXi::SetcTauOmega( float x){fCTauOmega = x;}   
inline void AliFemtoXi::SetptXi( float x){fPtXi = x;}          
inline void AliFemtoXi::SetptotXi( float x){fPtotXi = x;}
inline void AliFemtoXi::SetptBac( float x){fPtBac = x;}
inline void AliFemtoXi::SetptotBac( float x){fPtotBac = x;}    
inline void AliFemtoXi::SetidBac(const unsigned short& s){ fKeyBac= s;}
inline void AliFemtoXi::SetkeyBac(const unsigned short& s){ fKeyBac= s;}
inline void AliFemtoXi::SettpcHitsBac(const int& i){fTpcHitsBac=i;} 
inline void AliFemtoXi::SetdedxBac(float x){fDedxBachelor=x;}

//!!!!!!!!! new below 1.09.2015
inline void AliFemtoXi::SetCosPointingAngleXi(double x) {fCosPointingAngleXi = x;}
inline void AliFemtoXi::SetCosPointingAngleV0toXi(double x) {fCosPointingAngleV0toXi = x;}
inline void AliFemtoXi::SetPhiXi(double x) {fPhiXi = x;}
inline void AliFemtoXi::SetEtaXi(double x) {fEtaXi = x;}

inline void AliFemtoXi::SetTPCNclsBac(int x) {fTPCNclsBac=x;}
inline void AliFemtoXi::SetNdofBac(int x) {fNdofBac=x;}
inline void AliFemtoXi::SetStatusBac(unsigned long x) {fStatusBac=x;}
inline void AliFemtoXi::SetEtaBac(float x) {fEtaBac=x;}

inline void AliFemtoXi::SetBacNSigmaTPCK(float x) {fBacNSigmaTPCK = x;}
inline void AliFemtoXi::SetBacNSigmaTPCPi(float x) {fBacNSigmaTPCPi = x;}
inline void AliFemtoXi::SetBacNSigmaTPCP(float x) {fBacNSigmaTPCP = x;}

inline void AliFemtoXi::SetBacNSigmaTOFK(float x) {fBacNSigmaTOFK = x;}
inline void AliFemtoXi::SetBacNSigmaTOFPi(float x) {fBacNSigmaTOFPi = x;}
inline void AliFemtoXi::SetBacNSigmaTOFP(float x) {fBacNSigmaTOFP = x;}

inline void AliFemtoXi::SetTPCMomentumBac(double x) {fTPCMomentumBac = x;}

inline void AliFemtoXi::SetTOFProtonTimeBac(double x) {fTOFProtonTimeBac = x;}
inline void AliFemtoXi::SetTOFPionTimeBac(double x) {fTOFPionTimeBac = x;}
inline void AliFemtoXi::SetTOFKaonTimeBac(double x) {fTOFKaonTimeBac = x;}

inline void AliFemtoXi::SetChargeXi(int x) {fCharge = x;}

inline double AliFemtoXi::CosPointingAngleXi() const {return fCosPointingAngleXi;}
inline double AliFemtoXi::CosPointingAngleV0toXi() const {return fCosPointingAngleV0toXi;}
inline double AliFemtoXi::EtaXi() const {return fEtaXi;}
inline double AliFemtoXi::PhiXi() const {return fPhiXi;}

inline int AliFemtoXi::TPCNclsBac() const {return fTPCNclsBac;}
inline int AliFemtoXi::NdofBac() const {return fNdofBac;}
inline unsigned long AliFemtoXi::StatusBac() const {return fStatusBac;}
inline float AliFemtoXi::EtaBac() const {return fEtaBac;}

inline double AliFemtoXi::GetTPCMomentumBac() const {return fTPCMomentumBac;}

inline float AliFemtoXi::BacNSigmaTPCK() const {return fBacNSigmaTPCK;}
inline float AliFemtoXi::BacNSigmaTPCPi() const {return fBacNSigmaTPCPi;}
inline float AliFemtoXi::BacNSigmaTPCP() const {return fBacNSigmaTPCP;}

inline float AliFemtoXi::BacNSigmaTOFK() const {return fBacNSigmaTOFK;}
inline float AliFemtoXi::BacNSigmaTOFPi() const {return fBacNSigmaTOFPi;}
inline float AliFemtoXi::BacNSigmaTOFP() const {return fBacNSigmaTOFP;}

inline int AliFemtoXi::ChargeXi() const {return fCharge;}
 

#endif


/***********************************************************************
 *
 * $Log$
 * Revision 1.1.2.1  2007/10/05 09:38:17  akisiel
 * Fix stray colons
 *
 * Revision 1.1  2007/05/16 10:22:12  akisiel
 * Making the directory structure of AliFemto flat. All files go into one common directory
 *
 * Revision 1.2  2007/05/03 09:42:29  akisiel
 * Fixing Effective C++ warnings
 *
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.3  2003/09/02 17:58:33  perev
 * gcc 3.2 updates + WarnOff
 *
 * Revision 1.2  2001/12/05 15:10:33  laue
 * Boris' updates (mainly access functions)
 *
 *
 ***********************************************************************/
















