PWGLF/SPECTRA/ChargedHadrons/dNdPt/Framework/AliAnalysisTaskUEStudy.hPWGLF/SPECTRA/ChargedHadrons/dNdPt/Framework/AlidNdPtTools.h/// \class AliAnalysisTaskUEStudy
/// \brief Task to study efficiency and contamination in MC
///
/// Fills THnSparse to study pt vs. nch in underlying event
///
/// \author Michael Linus Knichel <michael.linus.knichel@cern.ch>, CERN
/// \date Mar 8, 2019

#ifndef AliAnalysisTaskUEStudy_H
#define AliAnalysisTaskUEStudy_H

#include "AliAnalysisTaskMKBase.h"

class AliESDtrackCuts;
class AliVEvent;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class AliStack;
class AliHeader;
class AliGenEventHeader;
class AliESDtrack;
class AliMCParticle;

class AliAnalysisTaskUEStudy : public AliAnalysisTaskMKBase
{
    public:
                                AliAnalysisTaskUEStudy();
                                AliAnalysisTaskUEStudy(const char *name);
        virtual                 ~AliAnalysisTaskUEStudy();
        
        virtual void            AddOutput();                     //called at the beginning
        virtual Bool_t          IsEventSelected();               //called for each event
        virtual void            AnaEvent();                      //called once for every event                
        virtual void            AnaTrack(Int_t flag = 0);        //called once for every track
        virtual void            AnaParticleMC(Int_t flag = 0);   //called once for every mc particle
        
        
        static AliAnalysisTaskUEStudy* AddTaskUEStudy(const char* name = "TaskUEStudy", const char* outfile = 0);

    protected:    
        Double_t                fPtMax;                     //!<! largest Pt track in the event
        Double_t                fPtMaxPhi;                  //!<! Phi of the largest pt
        Double_t                fMCPtMax;                   //!<! largest Pt particle in the mc event
        Double_t                fMCPtMaxPhi;                //!<! Phi of the largest pt particle
        Int_t                   fNTracksTowards;            //!<! ntracks in towards region
        Int_t                   fNTracksTransverse;         //!<! ntracks in transvers region
        Int_t                   fNTracksAway;               //!<! ntracks in away region
        Int_t                   fMCNChTowards;              //!<! MC nCH in towards region
        Int_t                   fMCNChTransverse;           //!<! MC nCH in transvers region
        Int_t                   fMCNChAway;                 //!<! MC nCH in away region                   
        THnSparseF*             fHistUETracks;              //-> underlying event histogram in data
        THnSparseF*             fHistUE;                    //-> underlying event histogram in data
        THnSparseF*             fHistUETracksMC;            //-> underlying event histogram in mc
        THnSparseF*             fHistUEMC;                  //-> underlying event histogram in mc
        THnSparseF*             fHistUEPhiRes;              //-> underlying event resolution phi
        

        
    private:
        AliAnalysisTaskUEStudy(const AliAnalysisTaskUEStudy&); // not implemented
        AliAnalysisTaskUEStudy& operator=(const AliAnalysisTaskUEStudy&); // not implemented
        
    /// \cond CLASSIMP    
        ClassDef(AliAnalysisTaskUEStudy, 2);
    /// \endcond        
};

#endif
