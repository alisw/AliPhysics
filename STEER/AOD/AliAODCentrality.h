#ifndef AliAODCentrality_H
#define AliAODCentrality_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------------------------
//     AOD centrality class
//     Author: Alberica Toia, CERN, Alberica.Toia@cern.ch
//-------------------------------------------------------------------------

#include <TNamed.h>

class AliAODCentrality : public TNamed {
  
 public :
  
  AliAODCentrality();
  virtual ~AliAODCentrality();
  AliAODCentrality(const AliAODCentrality& cnt); 
  AliAODCentrality& operator=(const AliAODCentrality& cnt);
  
  void SetxVertex     	    (Double_t x) {fxVertex	 = x; }
  void SetyVertex 	    (Double_t x) {fyVertex	 = x; }
  void SetzVertex	    (Double_t x) {fzVertex	 = x; }
  void SetVertexer3d        (Bool_t   x) {fVertexer3d	 = x; }
  void SetbMC               (Double_t x) {fbMC = x; } 
  void SetNpartTargMC       (Int_t    x) {fNpartTargMC = x; }
  void SetNpartProjMC       (Int_t    x) {fNpartProjMC = x; }
  void SetNNColl            (Int_t    x) {fNNColl = x; }     
  void SetNNwColl           (Int_t    x) {fNNwColl = x; }    
  void SetNwNColl           (Int_t    x) {fNwNColl = x; }    
  void SetNwNwColl          (Int_t    x) {fNwNwColl = x; }   
  void SetNTracklets        (Int_t    x) {fNTracklets = x; } 
  void SetNSingleClusters   (Int_t    x) {fNSingleClusters = x; }
  void SetNClusters         (Int_t x0, Int_t x1, Int_t x2, Int_t x3, Int_t x4, Int_t x5) {
    fNClusters[0] = x0; 
    fNClusters[1] = x1; 
    fNClusters[2] = x2; 
    fNClusters[3] = x3; 
    fNClusters[4] = x4; 
    fNClusters[5] = x5; 
  }   
  void SetNChips            (Int_t x0, Int_t x1) {
    fNChips[0] = x0;       
    fNChips[1] = x1; 
  }      
  void SetbZDC              (Double_t x) {fbZDC = x; }           
  void SetNpartZDC          (Int_t    x) {fNpartZDC = x; }       
  void SetbZDCA             (Double_t x) {fbZDCA = x; }          
  void SetNpartZDCA         (Int_t    x) {fNpartZDCA = x; }      
  void SetbZDCC             (Double_t x) {fbZDCC = x; }          
  void SetNpartZDCC         (Int_t    x) {fNpartZDCC = x; }      
  void SetESDFlag  	  (UInt_t   x) {fESDFlag = x; } 
  void SetZNCEnergy 	  (Float_t  x) {fZNCEnergy = x; }
  void SetZPCEnergy	  (Float_t  x) {fZPCEnergy = x; }
  void SetZNAEnergy	  (Float_t  x) {fZNAEnergy = x; }
  void SetZPAEnergy	  (Float_t  x) {fZPAEnergy = x; }
  void SetZEM1Energy	  (Float_t  x) {fZEM1Energy = x; }
  void SetZEM2Energy	  (Float_t  x) {fZEM2Energy = x; }
  void SetZNCtower          (Float_t x0, Float_t x1, Float_t x2, Float_t x3, Float_t x4) {
    fZNCtower[0] = x0; fZNCtower[1] = x1; fZNCtower[2] = x2; fZNCtower[3] = x3; fZNCtower[4] = x4; 
  }   
  void SetZPCtower          (Float_t x0, Float_t x1, Float_t x2, Float_t x3, Float_t x4) {
    fZPCtower[0] = x0; fZPCtower[1] = x1; fZPCtower[2] = x2; fZPCtower[3] = x3; fZPCtower[4] = x4; 
  }   
  void SetZNAtower          (Float_t x0, Float_t x1, Float_t x2, Float_t x3, Float_t x4) {
    fZNAtower[0] = x0; fZNAtower[1] = x1; fZNAtower[2] = x2; fZNAtower[3] = x3; fZNAtower[4] = x4; 
  }   
  void SetZPAtower          (Float_t x0, Float_t x1, Float_t x2, Float_t x3, Float_t x4) {
    fZPAtower[0] = x0; fZPAtower[1] = x1; fZPAtower[2] = x2; fZPAtower[3] = x3; fZPAtower[4] = x4; 
  }   
  void SetCentrZNC          (Float_t x0, Float_t x1) {
    fCentrZNC[0] = x0; fCentrZNC[1] = x1; 
  }   
  void SetCentrZNA          (Float_t x0, Float_t x1) {
    fCentrZNA[0] = x0; fCentrZNA[1] = x1; 
  }   
  void SetNTracks           (Int_t    x) {fNTracks = x; }    
  void SetNPmdTracks        (Int_t    x) {fNPmdTracks = x; } 
  void SetMultV0A           (Double_t x) {fMultV0A = x; }    
  void SetMultV0C           (Double_t x) {fMultV0C = x; }    
  void SetMultFMDA          (Float_t  x) {fMultFMDA = x; }    
  void SetMultFMDC          (Float_t  x) {fMultFMDC = x; }
  
  Double_t   GetxVertex	          () const {return fxVertex	 ; }
  Double_t   GetyVertex	          () const {return fyVertex	 ; }
  Double_t   GetzVertex	          () const {return fzVertex	 ; }
  Bool_t     GetVertexer3d        () const {return fVertexer3d	 ; }
  Double_t   GetbMC               () const {return  fbMC 	  ;}               
  Int_t      GetNpartTargMC       () const {return  fNpartTargMC  ;}               
  Int_t      GetNpartProjMC       () const {return  fNpartProjMC  ;}               
  Int_t      GetNNColl            () const {return  fNNColl       ;}               
  Int_t      GetNNwColl           () const {return  fNNwColl      ;}               
  Int_t      GetNwNColl           () const {return  fNwNColl      ;}               
  Int_t      GetNwNwColl          () const {return  fNwNwColl     ;}               
  Int_t      GetNTracklets        () const {return  fNTracklets   ;}               
  Int_t      GetNSingleClusters   () const {return  fNSingleClusters;}               
  Int_t      GetNClusters0        () const {return  fNClusters[0];}               
  Int_t      GetNClusters1        () const {return  fNClusters[1];}               
  Int_t      GetNClusters2        () const {return  fNClusters[2];}               
  Int_t      GetNClusters3        () const {return  fNClusters[3];}               
  Int_t      GetNClusters4        () const {return  fNClusters[4];}               
  Int_t      GetNClusters5        () const {return  fNClusters[5];}              
  Int_t      GetNChip0            () const {return  fNChips[0];}               
  Int_t      GetNChip1            () const {return  fNChips[1];}               
  Double_t   GetbZDC              () const {return  fbZDC         ;}               
  Int_t      GetNpartZDC          () const {return  fNpartZDC     ;}               
  Double_t   GetbZDCA             () const {return  fbZDCA        ;}               
  Int_t      GetNpartZDCA         () const {return  fNpartZDCA    ;}               
  Double_t   GetbZDCC             () const {return  fbZDCC        ;}               
  Int_t      GetNpartZDCC         () const {return  fNpartZDCC    ;}               
  UInt_t     GetESDFlag  	  () const {return  fESDFlag      ;}               
  Float_t    GetZNCEnergy         () const {return  fZNCEnergy    ;}               
  Float_t    GetZPCEnergy	  () const {return  fZPCEnergy    ;}               
  Float_t    GetZNAEnergy	  () const {return  fZNAEnergy    ;}               
  Float_t    GetZPAEnergy	  () const {return  fZPAEnergy    ;}               
  Float_t    GetZEM1Energy        () const {return  fZEM1Energy   ;}               
  Float_t    GetZEM2Energy	  () const {return  fZEM2Energy   ;}               
  Float_t    GetZNCtower0         () const {return  fZNCtower[0]  ;}               
  Float_t    GetZNCtower1         () const {return  fZNCtower[1]  ;}               
  Float_t    GetZNCtower2         () const {return  fZNCtower[2]  ;}               
  Float_t    GetZNCtower3         () const {return  fZNCtower[3]  ;}               
  Float_t    GetZNCtower4         () const {return  fZNCtower[4]  ;}               
  Float_t    GetZPCtower0         () const {return  fZPCtower[0]  ;}               
  Float_t    GetZPCtower1         () const {return  fZPCtower[1]  ;}               
  Float_t    GetZPCtower2         () const {return  fZPCtower[2]  ;}               
  Float_t    GetZPCtower3         () const {return  fZPCtower[3]  ;}               
  Float_t    GetZPCtower4         () const {return  fZPCtower[4]  ;}               
  Float_t    GetZNAtower0         () const {return  fZNAtower[0]  ;}               
  Float_t    GetZNAtower1         () const {return  fZNAtower[1]  ;}               
  Float_t    GetZNAtower2         () const {return  fZNAtower[2]  ;}               
  Float_t    GetZNAtower3         () const {return  fZNAtower[3]  ;}               
  Float_t    GetZNAtower4         () const {return  fZNAtower[4]  ;}               
  Float_t    GetZPAtower0         () const {return  fZPAtower[0]  ;}               
  Float_t    GetZPAtower1         () const {return  fZPAtower[1]  ;}               
  Float_t    GetZPAtower2         () const {return  fZPAtower[2]  ;}               
  Float_t    GetZPAtower3         () const {return  fZPAtower[3]  ;}               
  Float_t    GetZPAtower4         () const {return  fZPAtower[4]  ;}               
  Float_t    GetCentrZNC0         () const {return  fCentrZNC[0]  ;}               
  Float_t    GetCentrZNC1         () const {return  fCentrZNC[1]  ;}               
  Float_t    GetCentrZNA0         () const {return  fCentrZNA[0]  ;}               
  Float_t    GetCentrZNA1         () const {return  fCentrZNA[1]  ;}               
  Int_t      GetNTracks           () const {return  fNTracks      ;}               
  Int_t      GetNPmdTracks        () const {return  fNPmdTracks   ;}               
  Double_t   GetMultV0A           () const {return  fMultV0A      ;}               
  Double_t   GetMultV0C           () const {return  fMultV0C      ;}               
  Float_t    GetMultFMDA          () const {return  fMultFMDA     ;}               
  Float_t    GetMultFMDC          () const {return  fMultFMDC     ;}               
  
  void         Print(Option_t* option = "") const;
  const char*  GetOutputFileName() const {return TNamed::GetName();}
  void         SetOutputFileName(const char* fname) {TNamed::SetName(fname);}

  
 private:
  Double_t fxVertex;		//  X vertex from ITS
  Double_t fyVertex;		//  Y vertex from ITS
  Double_t fzVertex;		//  Z vertex from ITS
  Bool_t   fVertexer3d;		//  Is vertex from 3d vertexer?
  //
  Double_t fbMC;                // impact parameter from MC
  Int_t    fNpartTargMC;        // no. of participants for target nucleus from MC
  Int_t    fNpartProjMC;        // no. of participants for project nucleus from MC
  Int_t    fNNColl;             // Number of N-N collisions
  Int_t    fNNwColl;            // Number of N-Nwounded collisions
  Int_t    fNwNColl;            // Number of Nwounded-N collisons
  Int_t    fNwNwColl;           // Number of Nwounded-Nwounded collisions
  //
  Int_t    fNTracklets;         //  no. tracklets
  Int_t    fNSingleClusters;    //  no. single clusters
  Int_t    fNClusters[6];       //  no. clusters on 6 ITS layers
  Int_t    fNChips[2];          //  no. chips on 2 SPD layers
  //
  Double_t fbZDC;               // impact parameter from ZDC reco
  Int_t    fNpartZDC;           // no. of participants from ZDC reco
  Double_t fbZDCA;              // impact parameter from ZDC reco side A
  Int_t    fNpartZDCA;          // no. of part. from ZDC reco side A
  Double_t fbZDCC;              // impact parameter from ZDC reco side C
  Int_t    fNpartZDCC;          // no. of part. from ZDC reco side C
  //
  UInt_t   fESDFlag;            //  ZDC ESD flags
  Float_t  fZNCEnergy;          //  ZNC Energy
  Float_t  fZPCEnergy;          //  ZPC Energy
  Float_t  fZNAEnergy;          //  ZNA Energy
  Float_t  fZPAEnergy;          //  ZPA Energy
  Float_t  fZEM1Energy;         //  ZEM1 Energy
  Float_t  fZEM2Energy;         //  ZEM2 Energy
  Float_t  fZNCtower[5];        //  ZNC 5 tower signals
  Float_t  fZPCtower[5];        //  ZPC 5 tower signals
  Float_t  fZNAtower[5];        //  ZNA 5 tower signals
  Float_t  fZPAtower[5];        //  ZPA 5 tower signals
  Float_t  fCentrZNC[2];        //  centroid over ZNC
  Float_t  fCentrZNA[2];        //  centroid over ZNA
  
  Int_t    fNTracks;            //  no. tracks
  Int_t    fNPmdTracks;         //  no. PMD tracks
  Double_t fMultV0A;            //  multiplicity from V0 reco side A
  Double_t fMultV0C;            //  multiplicity from V0 reco side C
  Float_t  fMultFMDA;      //  multiplicity from FMD on detector A
  Float_t  fMultFMDC;      //  multiplicity from FMD on detector C
  
  ClassDef(AliAODCentrality, 1);
};


#endif
