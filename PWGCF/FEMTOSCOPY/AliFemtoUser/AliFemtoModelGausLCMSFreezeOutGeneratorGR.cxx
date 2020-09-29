////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelGausLCMSFreezeOutGeneratorGR - freeze-out                     ///
/// coordinates generator, generating a 3D gaussian ellipsoid in LCMS        ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoModelGausLCMSFreezeOutGeneratorGR, 1);
  /// \endcond
#endif

#include "math.h"
#include "AliFemtoModelGausLCMSFreezeOutGeneratorGR.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelGlobalHiddenInfo.h"
#include "AliFemtoLorentzVector.h"
#include <TMath.h>

//_______________________
AliFemtoModelGausLCMSFreezeOutGeneratorGR::AliFemtoModelGausLCMSFreezeOutGeneratorGR() :
fSizeOut(0), fSizeSide(0), fSizeLong(0),
fSelectPrimary(false)
{
	// Default constructor
	fRandom = new TRandom2();
}

//_______________________
AliFemtoModelGausLCMSFreezeOutGeneratorGR::AliFemtoModelGausLCMSFreezeOutGeneratorGR(const AliFemtoModelGausLCMSFreezeOutGeneratorGR &aModel):
AliFemtoModelFreezeOutGenerator(aModel),
fSizeOut(0), fSizeSide(0), fSizeLong(0),
fSelectPrimary(false)
{
	// Copy constructor
	fRandom = new TRandom2();
	SetSizeOut(aModel.GetSizeOut());
	SetSizeSide(aModel.GetSizeSide());
	SetSizeLong(aModel.GetSizeLong());
}
//_______________________
AliFemtoModelGausLCMSFreezeOutGeneratorGR::~AliFemtoModelGausLCMSFreezeOutGeneratorGR()
{
}
//_______________________
AliFemtoModelGausLCMSFreezeOutGeneratorGR& AliFemtoModelGausLCMSFreezeOutGeneratorGR::operator=(const AliFemtoModelGausLCMSFreezeOutGeneratorGR &aModel)
{
	if (this != &aModel) {
		fRandom = new TRandom2();
		SetSizeOut(aModel.GetSizeOut());
		SetSizeSide(aModel.GetSizeSide());
		SetSizeLong(aModel.GetSizeLong());
	}
	
	return *this;
}
//_______________________
void AliFemtoModelGausLCMSFreezeOutGeneratorGR::GenerateFreezeOut(AliFemtoPair *aPair)
{
	// Generate two particle emission points with respect
	// to their pair momentum
	// The source is the 3D Gaussian ellipsoid in the LCMS frame
	//AliFemtoModelHiddenInfo *inf1 = (AliFemtoModelHiddenInfo *) aPair->Track1()->HiddenInfo();
	//AliFemtoModelHiddenInfo *inf2 = (AliFemtoModelHiddenInfo *) aPair->Track2()->HiddenInfo();
	// AliFemtoModelGlobalHiddenInfo
	
	//ml
	AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
	AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();
	
	if ((!inf1) || (!inf2)) { cout << "Hidden info not created! "  << endl; exit(kFALSE); }
	
	//std::cout<<" we are in Freeze-out Generator inf1 inf2  "<<inf1<<"  "<<inf2<<std::endl;
	//std::cout<<" inf1 GetMass "<<((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetPDGPid()<<std::endl;
	//std::cout<<" true mom    " <<((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetTrueMomentum()->x()<<std::endl;
	
	if (fSelectPrimary) 
	{
		const AliFemtoTrack *infg1 = dynamic_cast<const AliFemtoTrack *>(aPair->Track1()->Track());
		const AliFemtoTrack *infg2 = dynamic_cast<const AliFemtoTrack *>(aPair->Track2()->Track());
		
		if ((infg1) && (infg2)) {
			// assume the emission point is in [cm] and try to judge if
			// both particles are primary
			// Double_t dist1 = infg1->GetGlobalEmissionPoint()->Perp();
			//  Double_t dist2 = infg2->GetGlobalEmissionPoint()->Perp();
			
			Double_t dist1 = ((AliFemtoModelGlobalHiddenInfo*)infg1->GetHiddenInfo())->GetGlobalEmissionPoint()->Perp();
			Double_t dist2 = ((AliFemtoModelGlobalHiddenInfo*)infg2->GetHiddenInfo())->GetGlobalEmissionPoint()->Perp();
			
			
			if ((dist1 > 0.05) && (dist2 > 0.05)) {
				// At least one particle is non primary
				if (!(((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetEmissionPoint())) 
				{
					AliFemtoLorentzVector tPos(-1000,1000,-500,0);
					inf1->SetEmissionPoint(&tPos);
					((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->SetEmissionPoint(&tPos);
					
				}
				else
				{
					inf1->SetEmissionPoint(-1000,1000,-500,0);
					((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->SetEmissionPoint(-1000,1000,-500,0);
				}
				if (!(((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetEmissionPoint())) {
					AliFemtoLorentzVector tPos(fRandom->Gaus(0,1000.0),fRandom->Gaus(0,1000),fRandom->Gaus(0,1000),0.0);
					inf2->SetEmissionPoint(&tPos);
					((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->SetEmissionPoint(&tPos);
				}
				else
				{
					inf2->SetEmissionPoint(fRandom->Gaus(0,1000), fRandom->Gaus(0,1000), fRandom->Gaus(0,1000), 0.0);
					((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->SetEmissionPoint(fRandom->Gaus(0,1000), fRandom->Gaus(0,1000), fRandom->Gaus(0,1000), 0.0);
					
				}
				
				return;
			}
		}
	}
	
	
	Double_t tPx = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetTrueMomentum()->x()  + ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetTrueMomentum()->x();
	Double_t tPy = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetTrueMomentum()->y()  + ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetTrueMomentum()->y();
	Double_t tPz = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetTrueMomentum()->z()  + ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetTrueMomentum()->z();
	
	
	// Double_t tPy = inf1->GetTrueMomentum()->y() + inf2->GetTrueMomentum()->y();
	// Double_t tPz = inf1->GetTrueMomentum()->z() + inf2->GetTrueMomentum()->z();
	
	//std::cout<<" tPx tPy tPz"<<tPx<<" "<<tPy<<" "<<tPz<<std::endl;
	if (!(tPx==0 && tPy==0 && tPz==0 )) {
		
		Double_t tM1 = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetMass();
		Double_t tM2 = ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetMass();
		Double_t tE1 = sqrt(tM1*tM1 + ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetTrueMomentum()->Mag2());
		Double_t tE2 = sqrt(tM2*tM2 + ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetTrueMomentum()->Mag2());
		Double_t tEs = tE1 + tE2;
		
		//std::cout<<" tM1 tM2 tE1 tE2"<<tM1<<" "<<tM2<<" "<<tE2<<std::endl;
		
		
		Double_t tPt = sqrt(tPx*tPx + tPy*tPy);
		
		Double_t tRout = fRandom->Gaus(0.0, fSizeOut);
		Double_t tRside = fRandom->Gaus(0.0, fSizeSide);
		Double_t tRlong = fRandom->Gaus(0.0, fSizeLong);
		
		Double_t tXout = (tPx * tRout + tPy * tRside)/tPt;
		Double_t tXside = (tPy * tRout - tPx * tRside)/tPt;
		Double_t tBetaz = tPz/tEs;
		Double_t tGammaz = 1.0/TMath::Sqrt(1-tBetaz*tBetaz);
		
		Double_t tXlong = tGammaz * (tRlong + tBetaz * 0);
		Double_t tXtime = tGammaz * (0 + tBetaz * tRlong);
		
		//std::cout<<" tXout tXside before hidden infor "<<tXout<<" "<<tXside<<std::endl;
		
		inf1->SetEmissionPoint(0,0,0,0);
		((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->SetEmissionPoint(0,0,0,0);
		
		inf2->SetEmissionPoint(tXout, tXside, tXlong, tXtime);
		((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->SetEmissionPoint(tXout, tXside, tXlong, tXtime);
		
		//std::cout<<" after all tXout tXside "<<tXout<<" "<<tXside<<std::endl;
		
	}
}

//_______________________
void AliFemtoModelGausLCMSFreezeOutGeneratorGR::SetSizeOut(Double_t aSizeOut)
{
	fSizeOut = aSizeOut;
}
//_______________________
void AliFemtoModelGausLCMSFreezeOutGeneratorGR::SetSizeSide(Double_t aSizeSide)
{
	fSizeSide = aSizeSide;
}
//_______________________
void AliFemtoModelGausLCMSFreezeOutGeneratorGR::SetSizeLong(Double_t aSizeLong)
{
	fSizeLong = aSizeLong;
}
//_______________________
Double_t AliFemtoModelGausLCMSFreezeOutGeneratorGR::GetSizeOut() const
{
	return fSizeOut;
}
//_______________________
Double_t AliFemtoModelGausLCMSFreezeOutGeneratorGR::GetSizeSide() const
{
	return fSizeSide;
}
//_______________________
Double_t AliFemtoModelGausLCMSFreezeOutGeneratorGR::GetSizeLong() const
{
	return fSizeLong;
}
//_______________________
AliFemtoModelFreezeOutGenerator* AliFemtoModelGausLCMSFreezeOutGeneratorGR::Clone() const
{
	return GetGenerator();
}
//_______________________
AliFemtoModelFreezeOutGenerator* AliFemtoModelGausLCMSFreezeOutGeneratorGR::GetGenerator() const
{
	AliFemtoModelFreezeOutGenerator* tModel = new AliFemtoModelGausLCMSFreezeOutGeneratorGR(*this);
	return tModel;
}
//_______________________
void AliFemtoModelGausLCMSFreezeOutGeneratorGR::SetSelectPrimaryFromHidden(bool aUse)
{
	fSelectPrimary = aUse;
}
Bool_t AliFemtoModelGausLCMSFreezeOutGeneratorGR::GetSelectPrimaryFromHidden()
{
	return fSelectPrimary;
}