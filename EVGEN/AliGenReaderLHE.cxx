#include <TVirtualMC.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TParticle.h>

#include "AliLog.h"
#include "AliGenReaderLHE.h"
#include "AliRun.h"
#include "AliStack.h"

ClassImp(AliGenReaderLHE);

AliGenReaderLHE::~AliGenReaderLHE() {
  if (fXMLDoc) {
    fXMLEngine.FreeDoc(fXMLDoc);
    fXMLDoc = NULL;
  }
}

void AliGenReaderLHE::Init()
{
  if (fXMLDoc)
    fXMLEngine.FreeDoc(fXMLDoc);

  fXMLDoc = fXMLEngine.ParseFile(fFileName);
  if (!fXMLDoc)
    AliFatalF("fXMLEngine.ParseFile('%s') failed", fFileName);

  XMLNodePointer_t xmlMainNode = fXMLEngine.DocGetRootElement(fXMLDoc);

  if (TString(fXMLEngine.GetNodeName(xmlMainNode)) != "LesHouchesEvents")
    AliFatal("main node != 'LesHouchesEvents'");

  XMLAttrPointer_t attr = fXMLEngine.GetFirstAttr(xmlMainNode);
  if (!attr)
    AliFatal("cannot determine version of Les Houches file format");

  if (TString(fXMLEngine.GetAttrName(attr))  != "version" ||
      TString(fXMLEngine.GetAttrValue(attr)) != "1.0")
    AliFatal("version != 1.0");

  Bool_t foundEvent = kFALSE;
  fXMLChild = fXMLEngine.GetChild(xmlMainNode);
  while (fXMLChild) {
    foundEvent = (TString(fXMLEngine.GetNodeName(fXMLChild)) == "event");
    if (foundEvent)
      break;
    fXMLChild  = fXMLEngine.GetNext(fXMLChild);
  }
  if (!fXMLChild || !foundEvent)
    AliFatal("no events found");
}

Int_t AliGenReaderLHE::NextEvent()
{
  if (!fXMLChild) {
    AliError("no new event found");
    return 0;
  }

  const char* content = fXMLEngine.GetNodeContent(fXMLChild);
  if (!content) {
    AliError("event is empty");
    return 0;
  }
  fStrStream.str(content);

  Int_t nup=0, idprup=0;
  Double_t xwgtup=0, scalup=0, aqedup=0, aqcdup=0;
  fStrStream >> nup >> idprup >> xwgtup >> scalup >> aqedup >> aqcdup;
  if (!fStrStream)
    AliFatal("malformed event record (1st line)");

  // todo: xwgtup, scalup, aqedup, aqcdup -> event header
  fPosTracksBegin = fStrStream.tellg();

  // advance to the next event
  fXMLChild  = fXMLEngine.GetNext(fXMLChild);

  return nup;
}

TParticle* AliGenReaderLHE::NextParticle()
{
  if (!fStrStream) {
    AliError("missing particle");
    return NULL;
  }
  Int_t idup=0, istup=0, mothup[2]={0,0}, icolup[2]={0,0};
  Double_t p[4] = { 0,0,0,0 }, mass=0, vtimup=0, spinup=0;
  fStrStream >> idup >> istup >> mothup[0] >> mothup[1] >> icolup[0] >> icolup[1]
	     >> p[0] >> p[1] >> p[2] >> p[3] >> mass >> vtimup >> spinup;
  if (!fStrStream) {
    AliError("malformed particle record");
    return NULL;
  }

  const TLorentzVector vp(p);
  fParticle.SetPdgCode(idup);
  fParticle.SetStatusCode(istup);
  fParticle.SetFirstMother(mothup[0] - fConvIndicesFortranToC);
  fParticle.SetLastMother (mothup[1] - fConvIndicesFortranToC);
  fParticle.SetMomentum(vp);
  fParticle.SetProductionVertex(mass ? vp*(vtimup/mass) : TLorentzVector(0,0,0,0)); // TODO: check unit
  fParticle.SetBit(kTransportBit, istup == 1);

  return &fParticle;
}

void AliGenReaderLHE::RewindEvent()
{
  fStrStream.seekg(fPosTracksBegin);
}
