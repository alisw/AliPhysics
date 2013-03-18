#ifndef ALIAODHMPIDRINGS_H
#define ALIAODHMPIDRINGS_H


//
// Class to handle the AOD tracks with good HMPID rings 
// Author: Levente Molnar
// levente.molnar@cern.ch , March 2012
// 



//___ROOT includes
#include <TMath.h>
//___AliRoot includes
#include "AliPID.h"


class AliAODHMPIDrings : public TObject {
  
 public:
  
  AliAODHMPIDrings();
  AliAODHMPIDrings(
                    Int_t trkId,
                    Int_t qn, 
                    Int_t cluIdx,
                    Double32_t  trkTheta,
                    Double32_t trkPhi,
                    Double32_t signal,
                    Double32_t occ,
                    Double32_t chi2,
                    Double32_t trkX,
                    Double32_t trkY,
                    Double32_t mipX,
                    Double32_t mipY,
                    Double32_t *pid,
                    Double32_t *p                  );      //              
      
      
  AliAODHMPIDrings(const AliAODHMPIDrings& hmpidAOD);//
  AliAODHMPIDrings &operator=(const AliAODHMPIDrings& hmpidAOD);//
  virtual ~AliAODHMPIDrings() {};
    
  //___ Getters
  Double32_t GetHmpTrkID()                const { return fHmpidAODtrkId; }
  
  Double32_t GetHmpMipCharge()            const { return fHmpidAODqn%1000000; }
  Double32_t GetHmpNumOfPhotonClusters()  const { return fHmpidAODqn/1000000;}
 
  Int_t      GetHmpChamber()              const { return fHmpidAODcluIdx/1000000; }
  
  Int_t      GetHmpCluIdx()               const { return fHmpidAODcluIdx; }
  
  Double32_t GetHmpTrackTheta()           const { return fHmpidAODtrkTheta;}
  Double32_t GetHmpTrackPhi()             const { return fHmpidAODtrkPhi;}
  
  Double32_t GetHmpSignal()               const { return fHmpidAODsignal;}
  Double32_t GetHmpOccupancy()            const { return fHmpidAODocc;}
  
  Double32_t GetHmpChi2()                 const { return fHmpidAODchi2;}

  Double32_t GetHmpTrackX()               const { return fHmpidAODtrkX;}
  Double32_t GetHmpTrackY()               const { return fHmpidAODtrkY;}

  Double32_t GetHmpMipX()                 const { return fHmpidAODmipX;}
  Double32_t GetHmpMipY()                 const { return fHmpidAODmipY;}

  Double32_t GetHmpDX()                   const { return fHmpidAODmipX - fHmpidAODtrkX;}
  Double32_t GetHmpDY()                   const { return fHmpidAODmipY - fHmpidAODtrkY;}
  Double32_t GetHmpDist()                 const { return TMath::Sqrt((fHmpidAODmipX - fHmpidAODtrkX)*(fHmpidAODmipX - fHmpidAODtrkX) + (fHmpidAODmipY - fHmpidAODtrkY)*(fHmpidAODmipY - fHmpidAODtrkY));}
  
  
  void GetHmpPidProbs(Double32_t *pid) const;   //defined in cxx
  void GetHmpMom(Double32_t *mom)      const;   //defined in cxx
  
  //___ Setters
  
  void SetHmpMipCharge(Int_t q)               { fHmpidAODqn = q; }
  void SetHmpCluIdx(Int_t ch,Int_t idx)       { fHmpidAODcluIdx=ch*1000000+idx;}
  
  void SetHmpNumOfPhotonClusters(Int_t nph)   { fHmpidAODqn = 1000000 * nph;}
  
  void SetHmpTrackTheta(Double32_t trkTheta)  { fHmpidAODtrkTheta = trkTheta;}
  void SetHmpTrackPhi(Double32_t trkPhi)      { fHmpidAODtrkPhi = trkPhi;}
  
  void SetHmpSignal(Double32_t thetaC)        { fHmpidAODsignal = thetaC;}
  void SetHmpOccupancy(Double32_t occ)        { fHmpidAODocc =  occ;}
  
  void SetHmpChi2(Double32_t chi2)            { fHmpidAODchi2 = chi2;}

  void SetHmpTrackX(Double32_t trkX)          { fHmpidAODtrkX = trkX;}
  void SetHmpTrackY(Double32_t trkY)          { fHmpidAODtrkY = trkY;}

  void SetHmpMipX(Double32_t mipX)            { fHmpidAODmipX = mipX;}
  void SetHmpMipY(Double32_t mipY)            { fHmpidAODmipY = mipY;}
 
  void SetHmpPidProbs(Double32_t *pid);       
  void SetHmpMom(Double32_t *mom);        
  
  
  // blablabla
  
  
  
 protected:
  
  Int_t       fHmpidAODtrkId;                      // Unique track id as in ESD
  Int_t       fHmpidAODqn;                         // 1000000*number of photon clusters + QDC
  Int_t       fHmpidAODcluIdx;                     // 1000000*chamber id + cluster idx of the assigned MIP cluster
  
  Double32_t  fHmpidAODtrkTheta;                   // [-2*pi,2*pi,16] theta of the track extrapolated to the HMPID, LORS
  Double32_t  fHmpidAODtrkPhi;                     // [-2*pi,2*pi,16] theta of the track extrapolated to the HMPID, LORS
  Double32_t  fHmpidAODsignal;                     // [0,0.9,8] HMPID signal (Theta ckov, rad)
  Double32_t  fHmpidAODocc;                        // [0,0,,8]  chamber occupancy where the track passed through: number of pads
  Double32_t  fHmpidAODchi2;                       // [0.,0.,8] chi2 in the HMPID  
  Double32_t  fHmpidAODtrkX;                       // [0.,0.,8] x of the track impact, LORS 
  Double32_t  fHmpidAODtrkY;                       // [0.,0.,8] y of the track impact, LORS 
  Double32_t  fHmpidAODmipX;                       // [0.,0.,8] x of the MIP in LORS
  Double32_t  fHmpidAODmipY;                       // [0.,0.,8] y of the MIP in LORS
  Double32_t  fHmpidAODpid[AliPID::kSPECIES];      // [0.,0.,8] "detector response probabilities" (for the PID)
  Double32_t  fHMPIDmom[3];                          // [0.,0.,8] track momentum at the HMPID ring reconstruction
  
  ClassDef(AliAODHMPIDrings,1)  
        
};
#endif
