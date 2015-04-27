#include "AliMultiplictyLoaderTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "TParameter.h"
#include "AliESDtrackCuts.h"
#include "TChain.h"
#include "AliCentrality.h"
#include "AliPPVsMultUtils.h"
#include "AliESDtrackCuts.h"
using namespace std;

ClassImp(AliMultiplictyLoaderTask)

AliMultiplictyLoaderTask::AliMultiplictyLoaderTask(const char *name ):
AliAnalysisTaskSE(name),
fESD(0),
fAliPPVsMultUtils(0),
fCentEstimator(""),
fUseAliPPVsMultUtils(0),
fcentvalue(0),
fncharged05value(0),
fncharged08value(0),
fFirstEvent(kFALSE),
fDonotusetrackelts(kFALSE)		
{
	 DefineInput(0, TChain::Class());
	//DefineOutput(1, TList::Class());
}
//___________________________________________________________________________
AliMultiplictyLoaderTask:: ~AliMultiplictyLoaderTask()
{
	if(fAliPPVsMultUtils)
		delete fAliPPVsMultUtils;
}
//____________________________________________________________________
void AliMultiplictyLoaderTask::UserCreateOutputObjects()
{

	fAliPPVsMultUtils=new AliPPVsMultUtils();
	TString par1name("cent");
	TString par2name("Ncheta0dot5");
	TString par3name("Ncheta0dot8");
	par1name+=fCentEstimator;	
	fcentvalue=new TParameter<Double_t>(par1name.Data(),-10.0);
	fncharged05value=new TParameter<Int_t>(par2name.Data(),-10);
	fncharged08value=new TParameter<Int_t>(par3name.Data(),-10);


}
//_______________________________________________________________________
void AliMultiplictyLoaderTask::UserExec(Option_t *)
{
	fESD = dynamic_cast<AliESDEvent*> (InputEvent());
	if (!fESD) 
	{
		Printf("ERROR: fESD not available");
		return;
	}
	Double_t cent=-1.0;
	Int_t ncharged05=-10;
	Int_t ncharged08=-10;

	if(fUseAliPPVsMultUtils)
	{
        	cent=fAliPPVsMultUtils->GetMultiplicityPercentile(fESD,fCentEstimator.Data());
		ncharged08=AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD);
		if(ncharged08<1)
		{
			ncharged05=ncharged08;	
		}
		else
		{
			AliESDtrackCuts::MultEstTrackType estType = fESD->GetPrimaryVertexTracks()->GetStatus() ? AliESDtrackCuts::kTrackletsITSTPC : AliESDtrackCuts::kTracklets;
			if(fDonotusetrackelts)
			{
				estType = AliESDtrackCuts::kTrackletsITSTPC;
			}
			ncharged05=AliESDtrackCuts::GetReferenceMultiplicity(fESD,estType,0.5);
		}
				

	}	
	else
	{
		cent=fESD->GetCentrality()->GetCentralityPercentile(fCentEstimator.Data());		
		AliESDtrackCuts::MultEstTrackType estType = fESD->GetPrimaryVertexTracks()->GetStatus() ? AliESDtrackCuts::kTrackletsITSTPC : AliESDtrackCuts::kTracklets;
		if(fDonotusetrackelts)
		{
			estType = AliESDtrackCuts::kTrackletsITSTPC;
		}
		ncharged05=AliESDtrackCuts::GetReferenceMultiplicity(fESD,estType,0.5);
		ncharged08=AliESDtrackCuts::GetReferenceMultiplicity(fESD,estType,0.8);	
	}

	fcentvalue->SetVal(cent);
	fncharged05value->SetVal(ncharged05);
	fncharged08value->SetVal(ncharged08);
	if(!fFirstEvent)
	{
		fESD->AddObject(fcentvalue);	
		fESD->AddObject(fncharged05value);
		fESD->AddObject(fncharged08value);
		fFirstEvent=kTRUE;
	}	

} 
