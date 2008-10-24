#ifndef AliAnaPartCorrBaseClass_H
#define AliAnaPartCorrBaseClass_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */

//_________________________________________________________________________
// Base class for analysis algorithms
//-- Author: Gustavo Conesa (INFN-LNF)

//ROOT
class TClonesArray ;
#include <TList.h>
#include <TObject.h>

//Analysis
class AliAODCaloCluster;
class AliAODCaloCells;
#include "AliAODParticleCorrelation.h"
class AliCaloTrackReader ;   
#include "AliCaloPID.h"
class AliFidutialCut ;
class AliIsolationCut ;
class AliNeutralMesonSelection ;
/* class AliStack ; */
/* class AliHeader ; */
/* class AliGenEventHeader ; */
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"

class AliAnaPartCorrBaseClass : public TObject {
       
  public: 
       
       AliAnaPartCorrBaseClass() ; // default ctor
       AliAnaPartCorrBaseClass(const AliAnaPartCorrBaseClass & g) ; // cpy ctor
       AliAnaPartCorrBaseClass & operator = (const AliAnaPartCorrBaseClass & g) ;//cpy assignment
       virtual ~AliAnaPartCorrBaseClass() ; //virtual dtor
              
       virtual void AddAODCaloCluster(AliAODCaloCluster calo) ;
       virtual void AddAODParticleCorrelation(AliAODParticleCorrelation pc) ;

       virtual void ConnectAODCaloClusters();
       virtual void ConnectAODPHOSCells();
       virtual void ConnectAODEMCALCells();

       virtual TList * GetCreateOutputObjects(){ return (new TList) ;}
       
       virtual void Print(const Option_t * ) const {;}
       
       virtual void MakeAnalysisFillAOD()  {;}
       
       virtual void MakeAnalysisFillHistograms() {;}
  
       virtual Int_t GetDebug() const  { return fDebug ; }
       virtual void SetDebug(Int_t d)   { fDebug = d ; }

       virtual AliCaloTrackReader * GetReader() const {return fReader ; }
       virtual void SetReader(AliCaloTrackReader * reader) { fReader = reader ; }

       virtual TClonesArray* GetAODBranch() const {return fAODBranch ;}
       virtual TClonesArray* GetAODCaloClusters() const {return fAODCaloClusters ;}
       virtual AliAODCaloCells* GetAODCaloCells() const {return fAODCaloCells ;}

       virtual TClonesArray* GetAODCTS() const ;
       virtual TClonesArray* GetAODEMCAL() const ;
       virtual TClonesArray* GetAODPHOS() const ;
       
       virtual TNamed * GetEMCALCells() const ;
       virtual TNamed * GetPHOSCells() const ;

       virtual AliStack * GetMCStack() const ;
       virtual AliHeader* GetMCHeader() const ;
       virtual AliGenEventHeader* GetMCGenEventHeader() const ;

       virtual void SetAODBranch(TClonesArray * tca) { fAODBranch = tca ; }

       virtual AliCaloPID * GetCaloPID() const {return  fCaloPID ;}
       virtual void SetCaloPID(AliCaloPID * pid) { fCaloPID = pid ;}

       virtual AliFidutialCut * GetFidutialCut() const {return  fFidCut ;}
       virtual void SetFidutialCut(AliFidutialCut * fc) { fFidCut = fc ;}

       virtual AliIsolationCut * GetIsolationCut() const {return  fIC ;}
       virtual void SetIsolationCut(AliIsolationCut * fc) { fIC = fc ;}
       
       virtual AliNeutralMesonSelection * GetNeutralMesonSelection() const {return  fNMS ;}
       virtual void SetNeutralMesonSelection(AliNeutralMesonSelection * nms) { fNMS = nms ;}

       virtual Bool_t     IsDataMC() const {return fDataMC ; }
       virtual void SwitchOnDataMC()    {fDataMC = kTRUE ; }
       virtual void SwitchOffDataMC()    {fDataMC = kFALSE ; }
       
       virtual Bool_t IsFidutialCutOn() {return fCheckFidCut ; }
       virtual void SwitchOnFidutialCut() { fCheckFidCut = kTRUE;}
       virtual void SwitchOffFidutialCut() { fCheckFidCut = kFALSE;}
       
       virtual Bool_t IsCaloPIDOn() {return fCheckCaloPID ; }
       virtual void SwitchOnCaloPID() { fCheckCaloPID = kTRUE;}
       virtual void SwitchOffCaloPID() { fCheckCaloPID = kFALSE;}
       
       virtual Bool_t IsCaloPIDRecalculationOn() {return fRecalculateCaloPID ; }
       virtual void SwitchOnCaloPIDRecalculation() { fRecalculateCaloPID  = kTRUE;}
       virtual void SwitchOffCaloPIDRecalculation() { fRecalculateCaloPID  = kFALSE;}
       
       virtual Float_t    GetMaxPt()         const {return fMaxPt ; }
       virtual Float_t    GetMinPt()         const {return fMinPt ; }
       virtual void SetMaxPt(Float_t pt)              {fMaxPt = pt ; }
       virtual void SetMinPt(Float_t pt)              {fMinPt = pt ; }
       void SetPtCutRange(Double_t ptmin, Double_t ptmax)
       {  fMaxPt=ptmax;   fMinPt=ptmin;}
       
       virtual void InitParameters() ;
    
 private:    

       Bool_t fDataMC ; //Flag to access MC data when using ESD or AOD     
       Int_t fDebug ; // Debug level
       Bool_t fCheckFidCut ; // Do analysis for clusters in defined region         
       Bool_t fCheckCaloPID ; // Do analysis for calorimeters
       Bool_t fRecalculateCaloPID ; //Recalculate PID or use PID weights in calorimeters
       Float_t fMinPt ; //Maximum pt of (trigger) particles in the analysis
       Float_t fMaxPt ; //Minimum pt of (trigger) particles in the analysis
 
       AliCaloTrackReader * fReader; //Acces to ESD/AOD/MC data

       TClonesArray* fAODBranch ;        //! selected particles branch
       TClonesArray* fAODCaloClusters ;     //! selected PHOS/EMCAL CaloClusters
       AliAODCaloCells * fAODCaloCells ;     //! selected PHOS/EMCAL CaloCells
       AliCaloPID * fCaloPID; // PID calculation
       AliFidutialCut * fFidCut; //Acceptance cuts
       AliIsolationCut * fIC; // Isolation cut 
       AliNeutralMesonSelection * fNMS; // Neutral Meson Selection

       ClassDef(AliAnaPartCorrBaseClass,1)
 } ;


#endif //AliAnaPartCorrBaseClass_H



