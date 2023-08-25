#ifndef AliAnalysisTaskPHOSCalibSelection_H
#define AliAnalysisTaskPHOSCalibSelection_H

#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliVCaloCells.h"

class AliAnalysisTaskPHOSCalibSelection : public AliAnalysisTaskSE
{
    public: 

        AliAnalysisTaskPHOSCalibSelection(); 
        AliAnalysisTaskPHOSCalibSelection(const char * name);
        virtual ~AliAnalysisTaskPHOSCalibSelection();

        void    UserCreateOutputObjects();

        virtual void   UserExec(Option_t *option);

        void    ProcessCells();
        
        void     ResetBuffer();


        void    SetCellMinimumEnergy(float emin)    {fCellEmin = emin;}
        void    SetSaveFullTree()                  {fSaveFullTree = true;}
        void    SetSaveCells()                     {fSaveCells = true; }

    private:

        float   fCellEmin;
        AliVEvent*          fInputEvent;
        AliVCaloCells       *fPHOSCells;
        bool                fSaveCells;
        bool                fSaveFullTree;
        TList             * fOutputContainer;

     

    protected:

        std::vector<int>      fCellIDVector;
        std::vector<float>    fCellTimeVector;
        std::vector<float>    fCellEnergyVector;

        TTree*        fCellTree;

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskPHOSCalibSelection,2) ;
    /// \endcond





};
#endif //AliAnalysisTaskPHOSCalibSelection_H
