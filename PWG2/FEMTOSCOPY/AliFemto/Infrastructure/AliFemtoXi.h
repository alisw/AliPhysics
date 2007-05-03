/***********************************************************************
 *
 * $Id: AliFemtoXi.h,v 1.0 1999/09/07
 *
 * Authors: Frank Laue
 *
 ***********************************************************************
 *
 * Description: Xi class with members copied from StXihbtXi.h
 *
 ***********************************************************************/
#ifndef AliFemtoXi_hh
#define AliFemtoXi_hh

#include "Infrastructure/AliFemtoVector.h" //same as in AliFemtoTrack.h
#include "Infrastructure/AliFemtoV0.h"

#ifdef __ROOT__
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
#include "StStrangeMuDstMaker/StXiMuDst.h"
#endif
#endif

class AliFemtoXi : public AliFemtoV0 {
public:
  AliFemtoXi();
#ifdef __ROOT__
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
  AliFemtoXi(StXiMuDst&); // from strangeness Xi micro dst structure
#endif
#endif
  virtual ~AliFemtoXi(){/* no-op */}

  void UpdateXi();
  float decayLengthXi() const;            // 3-d decay distance
  AliFemtoThreeVector decayVertexXi() const; // Coordinates of decay vertex
  float decayVertexXiX() const;           // Coordinates of decay vertex
  float decayVertexXiY() const;           // Coordinates of decay vertex
  float decayVertexXiZ() const;           // Coordinates of decay vertex
  float dcaXiDaughters() const;           // DCA of xi daughters at decay vertex
  float dcaXiToPrimVertex() const;        // DCA of xi to primary vertex
  float dcaBacToPrimVertex() const;       // DCA of bachelor  xi daughter to pri vertex
  AliFemtoThreeVector momBac() const;        // Momentum components of bac. daughter
  float momBacX() const;                  // Momentum components of bac. daughter
  float momBacY() const;                  // Momentum components of bac. daughter
  float momBacZ() const;                  // Momentum components of bac. daughter

  int   tpcHitsBac() const;               // Number of TPC hits on bac. daughter
  unsigned long trackTopologyMapBac(unsigned int) const;

  AliFemtoThreeVector momXi() const ;        // Momentum components of Xi
  float momXiX() const ;                  // Momentum components of Xi
  float momXiY() const ;                  // Momentum components of Xi
  float momXiZ() const ;                  // Momentum components of Xi
  float alphaXi() const ;                 // Armenteros-Podolanski variable
  float ptArmXi() const ;                 // Armenteros-Podolanski variable
  float eXi() const ;                     // Energy assuming xi hypothesis
  float eOmega() const ;                  // Energy assuming omega hypothesis
  float eBacKaon() const ;                // Energy of bac. daughter assuming kaon
  float eBacPion() const ;                // Energy of bac. daughter assuming pion
  float massXi() const ;                  // Mass assuming Xi hypothesis
  float massOmega() const ;               // Mass assuming Omega hypothesis
  float rapXi() const ;                   // Rapidity assuming (anti) xi
  float rapOmega() const ;                // Rapidity assuming (anti) omega
  float cTauXi() const ;                  // Lifetime (ctau) const assuming (anti) xi
  float cTauOmega() const ;               // Lifetime (ctau) const assuming (anti) omega
  float ptXi() const ;                    // Transverse momentum
  float ptotXi() const ;                  // Total momentum
  float ptBac() const ;                   // Transverse momentum of bac. daughter
  float ptotBac() const ;                 // Total momentum of bac. daughter
  float dedxBac() const;                  // dedx of Bac track
  unsigned short   idBac() const;         // Id of bac. track
  unsigned short   keyBac() const;        // Id of bac. track

  void SetdecayLengthXi(const float);  
  void SetdecayVertexXi(const AliFemtoThreeVector);  
  void SetdecayVertexXiX(const float);
  void SetdecayVertexXiY(const float);
  void SetdecayVertexXiZ(const float);
  void SetdcaXiDaughters(const float); 
  void SetdcaXiToPrimVertex(const float);  
  void SetdcaBacToPrimVertex(const float); 
  void SetmomBac(const AliFemtoThreeVector);  
  void SetmomBacX(const float);  
  void SetmomBacY(const float);  
  void SetmomBacZ(const float);  

  void SettpcHitsBac(const int&);      

  void SetTrackTopologyMapBac(unsigned int, const unsigned long&);

  void SetmomXi( AliFemtoThreeVector);
  void SetmomXiX( float);
  void SetmomXiY( float);
  void SetmomXiZ( float);
  void SetalphaXi( float);       
  void SetptArmXi( float);       
  void SeteXi( float);     
  void SeteOmega( float);    
  void SeteBacPion( float);  
  void SeteBacKaon( float);    
  void SetmassXi( float);  
  void SetmassOmega( float);
  void SetrapXi( float);    
  void SetrapOmega( float);   
  void SetcTauXi( float);   
  void SetcTauOmega( float);  
  void SetptXi( float);         
  void SetptotXi( float);       
  void SetptBac( float);        
  void SetptotBac( float);      
  void SetidBac(const unsigned short&);
  void SetdedxBac(float);
  void SetkeyBac(const unsigned short&);


protected:
  int   fCharge;                
  float fDecayLengthXi;
  AliFemtoThreeVector fDecayVertexXi;
  float fDcaXiDaughters; 
  float fDcaXiToPrimVertex;
  float fDcaBachelorToPrimVertex;
  AliFemtoThreeVector fMomBachelor;

  unsigned int   fTopologyMapBachelor[2];
  unsigned short fKeyBachelor;

  int   fTpcHitsBac;

  float fChi2Xi;
  float fClXi;
  float fChi2Bachelor;
  float fClBachelor;

  float fDedxBachelor;
  unsigned short fNufDedxBachelor;

  // the following variables are not in the persistent version and can be calculated via UpdateXi();
  AliFemtoThreeVector fMomXi;
  float fAlphaXi;
  float fPtArmXi;

  float fEXi;
  float fEOmega;
  float fEBacPion;
  float fEBacKaon;
  float fMassXi;
  float fMassOmega;
  float fRapXi;
  float fRapOmega;
  float fCTauXi;
  float fCTauOmega;
  float fPtXi;
  float fPtotXi;
  float fPtBac;
  float fPtotBac;

  unsigned short   fKeyBac;
};

inline float AliFemtoXi::decayLengthXi() const { return fDecayLengthXi; }
inline AliFemtoThreeVector AliFemtoXi::decayVertexXi() const { return fDecayVertexXi; } 
inline float AliFemtoXi::decayVertexXiX() const { return fDecayVertexXi.x(); } 
inline float AliFemtoXi::decayVertexXiY() const { return fDecayVertexXi.y(); } 
inline float AliFemtoXi::decayVertexXiZ() const { return fDecayVertexXi.z(); } 
inline float AliFemtoXi::dcaXiDaughters() const { return fDcaXiDaughters; }
inline float AliFemtoXi::dcaXiToPrimVertex() const { return fDcaXiToPrimVertex; }
inline float AliFemtoXi::dcaBacToPrimVertex() const { return fDcaBachelorToPrimVertex; }
inline AliFemtoThreeVector AliFemtoXi::momBac() const { return fMomBachelor; }
inline float AliFemtoXi::momBacX() const { return fMomBachelor.x(); }
inline float AliFemtoXi::momBacY() const { return fMomBachelor.y(); }
inline float AliFemtoXi::momBacZ() const { return fMomBachelor.z(); }
inline AliFemtoThreeVector AliFemtoXi::momXi() const { return fMomXi; }
inline float AliFemtoXi::momXiX() const { return fMomXi.x(); }
inline float AliFemtoXi::momXiY() const { return fMomXi.y(); }
inline float AliFemtoXi::momXiZ() const { return fMomXi.z(); }
inline float AliFemtoXi::alphaXi() const { return fAlphaXi; }
inline float AliFemtoXi::ptArmXi() const {return fPtArmXi;}
inline float AliFemtoXi::eXi() const {return fEXi;}
inline float AliFemtoXi::eOmega() const {return fEOmega;}
inline float AliFemtoXi::eBacPion() const {return fEBacPion;}
inline float AliFemtoXi::eBacKaon() const {return fEBacKaon;}
inline float AliFemtoXi::massXi() const {return fMassXi;}
inline float AliFemtoXi::massOmega() const {return fMassOmega;}
inline float AliFemtoXi::rapXi() const {return fRapXi;}
inline float AliFemtoXi::rapOmega() const {return fRapOmega;}
inline float AliFemtoXi::cTauXi() const {return fCTauXi;}
inline float AliFemtoXi::cTauOmega() const {return fCTauOmega;}
inline float AliFemtoXi::ptXi() const {return fPtXi;}
inline float AliFemtoXi::ptotXi() const {return fPtotXi;}
inline float AliFemtoXi::ptBac() const {return fPtBac;}
inline float AliFemtoXi::ptotBac() const {return fPtotBac;}
inline int   AliFemtoXi::tpcHitsBac() const
             { return fTpcHitsBac; }
inline float AliFemtoXi::dedxBac() const {return fDedxBachelor;}

inline unsigned long   AliFemtoXi::trackTopologyMapBac(unsigned int word) const { return fTopologyMapBachelor[word]; }
inline unsigned short   AliFemtoXi::idBac() const { return fKeyBac; }; 
inline unsigned short   AliFemtoXi::keyBac() const { return fKeyBac; }

inline void AliFemtoXi::SetdecayLengthXi(const float x){ fDecayLengthXi= x;}   
inline void AliFemtoXi::SetdecayVertexXiX(const float x){ fDecayVertexXi.setX(x);}
inline void AliFemtoXi::SetdecayVertexXiY(const float x){ fDecayVertexXi.setY(x);}
inline void AliFemtoXi::SetdecayVertexXiZ(const float x){ fDecayVertexXi.setZ(x);}
inline void AliFemtoXi::SetdecayVertexXi(const AliFemtoThreeVector v){ fDecayVertexXi = v; }
inline void AliFemtoXi::SetdcaXiDaughters(const float x){fDcaXiDaughters= x;} 
inline void AliFemtoXi::SetdcaXiToPrimVertex(const float x){fDcaXiToPrimVertex= x;}   
inline void AliFemtoXi::SetdcaBacToPrimVertex(const float x){ fDcaBachelorToPrimVertex = x;} 
inline void AliFemtoXi::SetmomBac(const AliFemtoThreeVector v){fMomBachelor = v; }
inline void AliFemtoXi::SetmomBacX(const float x){fMomBachelor.setX(x);}
inline void AliFemtoXi::SetmomBacY(const float x){fMomBachelor.setY(x);}
inline void AliFemtoXi::SetmomBacZ(const float x){fMomBachelor.setZ(x);}
inline void AliFemtoXi::SetTrackTopologyMapBac(unsigned int word, const unsigned long& m){fTopologyMapBachelor[word]=m;} 
inline void AliFemtoXi::SetmomXi(AliFemtoThreeVector v){fMomXi= v; }
inline void AliFemtoXi::SetmomXiX(const float x){fMomXi.setX(x);}
inline void AliFemtoXi::SetmomXiY(const float x){fMomXi.setY(x);}
inline void AliFemtoXi::SetmomXiZ(const float x){fMomXi.setZ(x);}

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

#endif


/***********************************************************************
 *
 * $Log$
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
















