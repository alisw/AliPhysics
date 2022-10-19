#include "AliGammaConvEventMixing.h"

//________________________________________________________________________________
EventMixPoolMesonJets::EventMixPoolMesonJets()
{
  mixingPool.resize(vecJetPClasses.size() - 1);
  for (auto& i : mixingPool) {
    i.clear();
    i.resize(0);
  }
}

//________________________________________________________________________________
EventMixPoolMesonJets::EventMixPoolMesonJets(std::vector<float> vec)
{
  vecJetPClasses = vec;
  mixingPool.resize(vec.size() - 1);
  for (auto& i : mixingPool) {
    i.clear();
    i.resize(0);
  }
}

//________________________________________________________________________________
void EventMixPoolMesonJets::SetJetPtClasses(std::vector<float> vec)
{
  vecJetPClasses = vec;
  mixingPool.resize(vec.size() - 1);
}

//________________________________________________________________________________
int EventMixPoolMesonJets::getJetPIndex(float jetP)
{
  int index = -1;
  for (unsigned int i = 0; i < vecJetPClasses.size() - 1; ++i) {
    if (vecJetPClasses[i] < jetP && jetP < vecJetPClasses[i + 1]) {
      index = i;
      break;
    }
  }
  return index;
}

//________________________________________________________________________________
void EventMixPoolMesonJets::AddEvent(std::shared_ptr<EventWithJetAxis> ev, float jetP)
{
  int index = getJetPIndex(jetP);
  if (index < 0) {
    printf("ERROR: index out of range\n");
    return;
  }
  if (mixingPool[index].size() >= poolDepth) {
    // if (mixingPool[index][0])
    //   delete mixingPool[index][0];
    mixingPool[index].erase(mixingPool[index].begin());
  }
  mixingPool.at(index).push_back(ev);
}

//________________________________________________________________________________
unsigned int EventMixPoolMesonJets::GetNBGEvents(float jetP)
{
  int index = getJetPIndex(jetP);
  return mixingPool[index].size();
}

//________________________________________________________________________________
unsigned int EventMixPoolMesonJets::GetNGammasInEvt(float jetP, int evt, bool isCaloPhoton)
{
  int index = getJetPIndex(jetP);
  if (isCaloPhoton)
    return mixingPool[index][evt]->caloPhotons.size();
  return mixingPool[index][evt]->convPhotons.size();
}

//________________________________________________________________________________
std::vector<std::shared_ptr<AliAODConversionPhoton>> EventMixPoolMesonJets::getPhotonsRotated(unsigned int nEvt, float jetP, TVector3 jetAxis, bool isCaloPhoton)
{
  int index = getJetPIndex(jetP);
  if (index < 0) {
    printf("ERROR: index out of range\n");
    std::vector<std::shared_ptr<AliAODConversionPhoton>> tmpVec(1);
    tmpVec[0] = nullptr;
    return tmpVec;
  }

  if (nEvt >= mixingPool[index].size()) {
    printf("index for mixing pool out of range");
    std::vector<std::shared_ptr<AliAODConversionPhoton>> tmpVec(1);
    tmpVec[0] = nullptr;
    return tmpVec;
  }
  // calculate the shift
  double jetThetaCurEv = jetAxis.Theta();
  double jetPhiCurEv = jetAxis.Phi();

  double jetThetaMixEv = mixingPool[index][nEvt]->jetAxis.Theta();
  double jetPhiMixEv = mixingPool[index][nEvt]->jetAxis.Phi();

  double diffTheta = jetThetaCurEv - jetThetaMixEv;
  double diffPhi = jetPhiCurEv - jetPhiMixEv;

  unsigned int nGammas = mixingPool[index][nEvt]->getNPhotons(isCaloPhoton);
  std::vector<std::shared_ptr<AliAODConversionPhoton>> vecRotatedGammas(nGammas);
  TLorentzVector LVGammaRot;
  TVector3 tmpVec;
  for (unsigned int i = 0; i < nGammas; ++i) {
    // E has to stay the same

    double energy = mixingPool[index][nEvt]->getPhotonE(i, isCaloPhoton);
    double theta = mixingPool[index][nEvt]->getPhotonTheta(i, isCaloPhoton);
    double phi = mixingPool[index][nEvt]->getPhotonPhi(i, isCaloPhoton);

    theta += diffTheta;
    phi += diffPhi;

    tmpVec.SetMagThetaPhi(energy, theta, phi);
    LVGammaRot.SetPtEtaPhiM(tmpVec.Pt(), tmpVec.Eta(), tmpVec.Phi(), 0.);
    vecRotatedGammas[i] = std::make_shared<AliAODConversionPhoton>(&LVGammaRot);
  }
  return vecRotatedGammas;
}