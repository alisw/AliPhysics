/// \class AliJetEmbeddingSelRhoTask
/// \brief Track embedding into an event with rho in a given range
///
/// \ingroup PWGJEUSER
/// The class inherits from AliJetEmbeddingTask and takes care of the event
/// selection. Otherwise all is in AliJetEmbeddingTask 

///
/// \author Chiara Bianchin

#ifndef ALIJETEMBEDDINGSELRHOTASK_H
#define ALIJETEMBEDDINGSELRHOTASK_H

// $Id$

#include "AliJetEmbeddingTask.h"

class AliJetEmbeddingSelRhoTask : public AliJetEmbeddingTask{
	
public:
	
	AliJetEmbeddingSelRhoTask();
	AliJetEmbeddingSelRhoTask(const char *name);
	virtual ~AliJetEmbeddingSelRhoTask()                         {}
	
	void           UserCreateOutputObjects();
	
	void     SetRhoRange(Double_t min, Double_t max)             { fRhoMin = min; fRhoMax = max; }
	Double_t GetRhoMin() const                                   { return fRhoMin ;}
	Double_t GetRhoMax() const                                   { return fRhoMax ;}
	
	void     SetRhoName(TString rhoname)                         { fRhoName = rhoname; }
	
	
protected:
	void     Run();
	void     FillHistograms();
	void     Terminate(Option_t *option="");
	
private:

	Double_t fRhoMin;                           ///< Minimum Rho accepted
	Double_t fRhoMax;                           ///< Maximum Rho accepted
    TString  fRhoName;                          ///< Name of rho to be read
    
    // histograms
    TH1F     *fhQARhoEventRejection;            //!<! Events accepted and rejected
    TH1F     *fhQARho;                          //!<! Rho distribution of the accepted events
	/// \cond CLASSIMP
	ClassDef(AliJetEmbeddingSelRhoTask, 1) /// Jet embedding task with rho event selection
	/// \endcond
};

#endif
