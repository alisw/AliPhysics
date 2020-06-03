/// \class AliAnalysisTaskPtResStudy
/// \brief Analsysis task to study pT resolution correction
///
/// to study pt Resolution in MC and Data
///
/// \author Michael Linus Knichel <michael.linus.knichel@cern.ch>, CERN
/// \date Apr 8, 2019

#ifndef AliAnalysisTaskPtResStudy_H
#define AliAnalysisTaskPtResStudy_H

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

class AliAnalysisTaskPtResStudy : public AliAnalysisTaskMKBase
{
    public:
                                AliAnalysisTaskPtResStudy();
                                AliAnalysisTaskPtResStudy(const char *name);
        virtual                 ~AliAnalysisTaskPtResStudy();

        virtual void            AddOutput();                     //called at the beginning
        virtual Bool_t          IsEventSelected();               //called for each event
        virtual void            AnaEvent();                      //called once for every selected event        
        virtual void            AnaTrack(Int_t flag = 0);        //called once for every track
        
        static AliAnalysisTaskPtResStudy* AddTaskPtResStudy(const char* name = "TaskPtResStudy", const char* outfile = 0);

    protected:    
        THnSparseF*             fHistPtResCov;        //-> pt resolution from covariance matrix
        THnSparseF*             fHistPtResCovHighPt;  //-> pt resolution from covariance matrix smaller binning for highest pt
        THnSparseF*             fHistPtResMC;         //-> pt resolution from mc
        THnSparseF*             fHistPtRes;           //-> pt resolution, combination of covariance and track fit
        
    private:
        AliAnalysisTaskPtResStudy(const AliAnalysisTaskPtResStudy&); // not implemented
        AliAnalysisTaskPtResStudy& operator=(const AliAnalysisTaskPtResStudy&); // not implemented
        
    /// \cond CLASSIMP    
        ClassDef(AliAnalysisTaskPtResStudy, 2);
    /// \endcond        
};

#endif
