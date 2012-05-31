#include "../AliHLTTriggerDetectorGeomRectangle.h"
#include "../AliHLTTriggerDecisionParameters.h"

void CreateDetectorGeomsPhos()
{

  bool mod[5] = {false, false, true, true, true};

  TObjArray *geomConf = new TObjArray();

  //  geomConf->SetName("GeomConfig");

  
  // Make the geometry objects
  for(int i = 0; i < 5; i++)
    {
      if(mod[i])
	{
	  AliHLTTriggerDetectorGeom *tmpDet= new AliHLTTriggerDetectorGeomRectangle();

	  char tmpName[128];
	  sprintf(tmpName, "PHOS_%d", i);
	  TString name = tmpName;

	  tmpDet->SetDetName(name);
	  tmpDet->SetEtaMin(-0.12);
	  tmpDet->SetEtaMax(0.12);
	  float phiModule = 230 + i*20;
	  float phiMin = (phiModule - 10)*TMath::DegToRad();
	  float phiMax = (phiModule + 10)*TMath::DegToRad();
	  tmpDet->SetPhiMin(phiMin);
	  tmpDet->SetPhiMax(phiMax);

	  Double_t nX = TMath::Cos(TMath::DegToRad()*phiModule);
	  Double_t nY = TMath::Sin(TMath::DegToRad()*phiModule);
	  Double_t nZ = 0;
	  Double_t normVector[3] = {nX, nY, nZ};

	  Double_t pX = TMath::Cos(TMath::DegToRad()*phiModule)*460;
	  Double_t pY = TMath::Sin(TMath::DegToRad()*phiModule)*460;
	  Double_t pZ = 0;
	  Double_t initialPoint[3] = {pX, pY, pZ};

	  tmpDet->SetNormVector(normVector);
	  tmpDet->SetInitialPoint(initialPoint);
	  tmpDet->PrintDetectorGeom(cout);
	  geomConf->AddLast(tmpDet);
	}
    }
  
  //  Make the decision object
  
  AliHLTTriggerDecisionParameters *pars = new AliHLTTriggerDecisionParameters();
  pars->SetTriggerName(TString("PHOSgeomtrigger"));
  pars->SetReadoutListParameter(AliHLTReadoutList::kPHOS);
  pars->SetDescription(TString("Track in PHOS"));

  geomConf->AddLast(pars);

  TFile *fin = TFile::Open("PHOSgeomtrigger.root", "RECREATE");

  geomConf->Write("GeomConf", TObject::kSingleKey);

}
