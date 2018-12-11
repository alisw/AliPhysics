#ifndef ALIANALYSISTASKBSEMBEDDING_H
#define ALIANALYSISTASKBSEMBEDDING_H
#include "AliAnalysisTaskSE.h"
using namespace std;
                                                                               
class AliAnalysisTaskBSEmbedding : public AliAnalysisTaskSE{
  public:
    AliAnalysisTaskBSEmbedding();                                          
    AliAnalysisTaskBSEmbedding(const char *name, const char *option);
    ~AliAnalysisTaskBSEmbedding();
    AliAnalysisTaskBSEmbedding(const AliAnalysisTaskBSEmbedding& ap);             
    AliAnalysisTaskBSEmbedding& operator =(const AliAnalysisTaskBSEmbedding& ap); 

    
    //AliAnalysisTaskMyf0f2(const AliAnalysisTaskMyf0f2& ap);           
    //AliAnalysisTaskMyf0f2& operator=(const AliAnalysisTaskMyf0f2& ap);
    //virtual void UserCreateOutputObjects();                     
    virtual void UserExec(Option_t* option);
    static AliAnalysisTaskBSEmbedding * AddTaskBSEmbedding();

  
  private:
    TString fOption ="";
  ClassDef(AliAnalysisTaskBSEmbedding, 10)


};
#endif

