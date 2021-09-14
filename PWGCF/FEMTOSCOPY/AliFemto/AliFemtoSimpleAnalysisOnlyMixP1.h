#ifndef AliFemtoSimpleAnalysisOnlyMixP1_hh
#define AliFemtoSimpleAnalysisOnlyMixP1_hh


//////////////////////////////////////////////////////////////////////////////////
// dongfang.wang@cern.ch                                                        //
// AliFemtoSimpleAnalysisOnlyMixP1: Only keep update/save first particle        //
// The key is in line 200(if current event has first particle, but not have     //
// second particle, then save it)                                               //
// and line 260(Only allow first particle in mixing pool making pair            //
// with second particle in current event                                        //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoSimpleAnalysis.h"
class AliFemtoSimpleAnalysisOnlyMixP1 : public AliFemtoSimpleAnalysis{

    public:
        AliFemtoSimpleAnalysisOnlyMixP1();
        AliFemtoSimpleAnalysisOnlyMixP1(const AliFemtoSimpleAnalysisOnlyMixP1 &OriAnalysis);
        virtual ~AliFemtoSimpleAnalysisOnlyMixP1();
        AliFemtoSimpleAnalysisOnlyMixP1& operator =(const AliFemtoSimpleAnalysisOnlyMixP1 &OriAnalysis);

        
        virtual void ProcessEvent(const AliFemtoEvent* hbtEvent);

        // wdf for test
        virtual bool Pass(float inputcut);
        void Test(float InputNumber);
    private:
        float  number;

};


#endif
