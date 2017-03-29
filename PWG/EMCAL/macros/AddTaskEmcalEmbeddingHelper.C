/**
 * \file AddTaskEmcalEmbeddingHelper.C
 * \brief AddTask for AliAnalysisTaskEmcalEmbeddingHelper
 *
 * The main implementation of the AddTask is in the class. This is simply a wrapper, primarily intended
 * for use on the LEGO train.
 *
 * \author Raymond Ehlers <raymond.ehlers@cern.ch>, Yale University
 * \date Aug 25, 2016
 */
AliAnalysisTaskEmcalEmbeddingHelper * AddTaskEmcalEmbeddingHelper()
{  
  return AliAnalysisTaskEmcalEmbeddingHelper::AddTaskEmcalEmbeddingHelper();
}
