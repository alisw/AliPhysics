#ifndef ALIANALYSISTASKREASNUCLEXAOD_H
#define ALIANALYSISTASKREASNUCLEXAOD_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskReadNuclexAOD class
//-----------------------------------------------------------------

class TList;
class TH1F;
class TH2F;
class TH3F;
class TNtuple;
class AliAODPid;
class AliESDcascade;
//class AliCascadeVertexer; 

#include "TString.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskReadNuclexAOD : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskReadNuclexAOD();
  AliAnalysisTaskReadNuclexAOD(const char *name);
  virtual ~AliAnalysisTaskReadNuclexAOD();
  
  virtual void  UserCreateOutputObjects();
  virtual void  Init();
  virtual void LocalInit() {Init();}
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);
  
  void SetCollidingSystems(Short_t collidingSystems = 0)     {fCollidingSystems = collidingSystems;}
  void SetAnalysisType    (const char* analysisType = "ESD") {fAnalysisType = analysisType;}
  void SetDataType    (const char* dataType = "REAL") {fDataType = dataType;}

  //Double_t BetheBloch(Double_t bg,Double_t Charge,Bool_t optMC);
  Double_t BetheBloch(Double_t bg,Double_t Charge,Bool_t isPbPb);

  
 private:

  TString fAnalysisType;	     //! "ESD" or "AOD" analysis type	

  Short_t fCollidingSystems;	     //! 0 = pp collisions or 1 = AA collisions
  TString fDataType;		     //! "REAL" or "SIM" data type	
 

  TList	*fListHist;	             //! List of  histograms

  
  TH1F *fHistEventMultiplicity;
  TH2F *fHistTrackMultiplicity;
  TH2F *fHistTrackMultiplicityCent;
  TH2F *fHistTrackMultiplicitySemiCent;
  TH2F *fHistTrackMultiplicityMB;
  TH2F *fHistTrackMultiplicityPVCent;
  TH2F *fHistTrackMultiplicityPVSemiCent;
  TH2F *fHistTrackMultiplicityPVMB;
  TH2F *fhBB;
  TH2F *fhBBPions;
  TH2F *fhBBHe;
    
  TTree  *aodTree;
  TClonesArray *Nuclei;
  TClonesArray *SecondaryVertices;
  TClonesArray *DaughterTracks;

  Int_t nevent; //event of the reduced tree

  TTree *fNtuple1;                  //! Tree Pairs Pi/Proton "standard"
  
  Float_t trunNumber;
  Float_t tbunchcross;
  Float_t torbit;
  Float_t tperiod;
  Float_t teventtype;
  Float_t tTrackNumber;
  Float_t tpercentile;
  Float_t txPrimaryVertex;
  Float_t tyPrimaryVertex;
  Float_t tzPrimaryVertex;
  Float_t txSecondaryVertex;
  Float_t tySecondaryVertex;
  Float_t tzSecondaryVertex;
  Float_t tdcaTracks;
  Float_t tCosPointingAngle;
  Float_t tDCAV0toPrimaryVertex;
  Float_t tHeSign;
  Float_t tHepInTPC;
  Float_t tHeTPCsignal;
  Float_t tDcaHeToPrimVertex;
  Float_t tHeEta;
  Float_t tmomHex;
  Float_t tmomHey;
  Float_t tmomHez;
  Float_t tmomHeAtSVx;
  Float_t tmomHeAtSVy;
  Float_t tmomHeAtSVz;
  Float_t tHeTPCNcls;
  Float_t tHeimpactXY;
  Float_t tHeimpactZ;
  Float_t tHeITSClusterMap;
  Float_t tIsHeITSRefit;
  Float_t tPionSign;
  Float_t tPionpInTPC;
  Float_t tPionTPCsignal;
  Float_t tDcaPionToPrimVertex;
  Float_t tPionEta;
  Float_t tmomPionx;
  Float_t tmomPiony;
  Float_t tmomPionz;
  Float_t tmomNegPionAtSVx;
  Float_t tmomNegPionAtSVy;
  Float_t tmomNegPionAtSVz;
  Float_t tPionTPCNcls;
  Float_t tPionimpactXY;
  Float_t tPionimpactZ;
  Float_t tPionITSClusterMap;
  Float_t tIsPiITSRefit;
  Float_t txn;
  Float_t txp;
  Float_t tchi2He;
  Float_t tchi2Pi;

  
  TTree *fNtuple4;

  Float_t t3LHlrunNumber;        
  Float_t t3LHlBCNumber;         
  Float_t t3LHlOrbitNumber;      
  Float_t t3LHlPeriodNumber;     
  Float_t t3LHleventtype;        
  Float_t t3LHlpercentile;       
  Float_t t3LHlPx;               
  Float_t t3LHlPy;               
  Float_t t3LHlPz;               
  Float_t t3LHlEta;              
  Float_t t3LHlY;      
  Float_t t3LHlM;      
  Float_t t3LHlxPrimaryVertex;   
  Float_t t3LHlyPrimaryVertex;   
  Float_t t3LHlzPrimaryVertex;   
  Float_t t3LHlp0px;    
  Float_t t3LHlp0py;
  Float_t t3LHlp0pz;
  Float_t t3LHlp0ch;
  Float_t t3LHlp1px;
  Float_t t3LHlp1py;
  Float_t t3LHlp1pz;
  Float_t t3LHlp1ch;




 static const Int_t fgNrot;
 

  AliAnalysisTaskReadNuclexAOD(const AliAnalysisTaskReadNuclexAOD&);            // not implemented
  AliAnalysisTaskReadNuclexAOD& operator=(const AliAnalysisTaskReadNuclexAOD&); // not implemented
  
  ClassDef(AliAnalysisTaskReadNuclexAOD, 0);
};

#endif
