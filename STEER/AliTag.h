#ifndef ALITAG
#define ALITAG

#include <stdlib.h>
#include <iostream.h>

#include "TObject.h"
#include "TClonesArray.h"

//______________________________________________________________________________
class AliEventTag : public TObject
{

	private:
		Int_t    fAliceEventId;                    //The event id
		Int_t    fGUID;                            //The unique identifier of the file
		Int_t    fNumberOfParticipants;            //Number of participants
		Float_t  fImpactParameter;                 //The impact parameter
		Float_t  fPrimaryVertexX;                  //Primary vertex - X coordinate
		Float_t  fPrimaryVertexY;                  //Primary vertex - Y coordinate
		Float_t  fPrimaryVertexZ;                  //Primary vertex - Z coordinate
		Int_t    fTriggerInfo;                     //Information from trigger
		Float_t  fZDCNeutronEnergy;                //ZDC info - neutron
		Float_t  fZDCProtonEnergy;                 //ZDC info - proton
		Float_t  fZDCEMEnergy;                     //ZDC info - em
		Float_t  fT0VertexZ;                       //T0 info
		Int_t    fNumberOfTracks;                  //Multiplicity
		Int_t    fNumberOfPositiveTracks;          //Multiplicity of positive tracks
		Int_t    fNumberOfNegativeTracks;          //Multiplicity of negative tracks
		Int_t    fNumberOfNeutralTracks;           //Multiplicity of neutral tracks
		Int_t    fNumberOfV0s;                     //Number of V0s
		Int_t    fNumberOfCascades;                //Number of cascades
		Int_t    fNumberOfKinks;                   //Number of kinks
		Int_t    fNumberOfPMDTracks;               //PMD tracks
		Int_t    fNumberOfPHOSClusters;            //PHOS clusters
		Int_t    fNumberOfEMCALClusters;           //EMCAL clusters
		Int_t    fNumberOfFMDTracks;               //FMD tracks
		Int_t    fNumberOfJetCandidates;           //Jet candidates
		Int_t    fNumberOfHardPhotonsCandidates;   //Hard photons candidates
		Int_t    fNumberOfElectrons;               //Number of electrons
		Int_t    fNumberOfMuons;                   //Number of muons
		Int_t    fNumberOfPions;                   //Number of pions
		Int_t    fNumberOfKaons;                   //Number of kaons
		Int_t    fNumberOfProtons;                 //Number of protons
		Int_t    fNumberOfLambdas;                 //Number of lambdas
		Int_t    fNumberOfJPsiCandidates;          //JPsi candidates
		Int_t    fNumberOfPsiPrimeCandidates;      //Psi prime candidates
		Int_t    fNumberOfUpsilonCandidates;       //Upsilon candidates
		Int_t    fNumberOfUpsilonPrimeCandidates;  //Upsilon prime candidates
		Int_t    fNumberOfUpsilonDoublePrimeCandidates;
		Int_t    fNumberOfCharmParticleCandidates;
		Int_t    fNumberOfBeautyParticleCandidates;
		Float_t  fK0PeakPosition;
		Float_t  fK0PeakWidth;
		Float_t  fTotalP;
		Float_t  fMeanPt;
		Float_t  fMaxPt;
		Float_t  fFlowV1;
		Float_t  fFlowV2;

		virtual void CopyTag(AliEventTag *EvTag);

	public:
		AliEventTag();
		AliEventTag(AliEventTag *t);
		virtual ~AliEventTag();

		void          SetEventId(Int_t Pid) {fAliceEventId = Pid;}
		void          SetGUID(Int_t Pid) {fGUID = Pid;}

		void          SetNumOfParticipants(Int_t P) {fNumberOfParticipants = P;}
		void          SetImpactParameter(Float_t Pimpact) {fImpactParameter = Pimpact;}

		void          SetVertexX(Float_t Pvx) {fPrimaryVertexX = Pvx;}
		void          SetVertexY(Float_t Pvy) {fPrimaryVertexY = Pvy;}
		void          SetVertexZ(Float_t Pvz) {fPrimaryVertexZ = Pvz;}

		void          SetTrigger(Int_t Ptr) {fTriggerInfo = Ptr;}

		void          SetZDCNeutronEnergy(Float_t Pen) {fZDCNeutronEnergy = Pen;}
		void          SetZDCProtonEnergy(Float_t Pen) {fZDCProtonEnergy = Pen;}
		void          SetZDCEMEnergy(Float_t Pen) {fZDCEMEnergy = Pen;}

		void          SetT0VertexZ(Float_t Pvz) {fT0VertexZ = Pvz;}

		void          SetNumOfTracks(Int_t Ptr) {fNumberOfTracks = Ptr;}
		void          SetNumOfPosTracks(Int_t Ptr) {fNumberOfPositiveTracks = Ptr;}
		void          SetNumOfNegTracks(Int_t Ptr) {fNumberOfNegativeTracks = Ptr;}
		void          SetNumOfNeutrTracks(Int_t Ptr) {fNumberOfNeutralTracks = Ptr;}

		void          SetNumOfV0s(Int_t Ptr) {fNumberOfV0s = Ptr;}
		void          SetNumOfCascades(Int_t Ptr) {fNumberOfCascades = Ptr;}
		void          SetNumOfKinks(Int_t Ptr) {fNumberOfKinks = Ptr;}

		void          SetNumOfPMDTracks(Int_t Ptr) {fNumberOfPMDTracks = Ptr;}
		void          SetNumOfPHOSClusters(Int_t Ptr) {fNumberOfPHOSClusters = Ptr;}
		void          SetNumOfEMCALClusters(Int_t Ptr) {fNumberOfEMCALClusters = Ptr;}
		void          SetNumOfFMDTracks(Int_t Ptr) {fNumberOfFMDTracks = Ptr;}

		void          SetNumOfJetCandidates(Int_t Ptr) {fNumberOfJetCandidates = Ptr;}
		void          SetNumOfHardPhotonsCandidates(Int_t Ptr) {fNumberOfHardPhotonsCandidates = Ptr;}
		void          SetNumOfJPsiCandidates(Int_t Ptr) {fNumberOfJPsiCandidates = Ptr;}
		void          SetNumOfPsiPrimeCandidates(Int_t Ptr) {fNumberOfPsiPrimeCandidates = Ptr;}
		void          SetNumOfUpsilonCandidates(Int_t Ptr) {fNumberOfUpsilonCandidates = Ptr;}
		void          SetNumOfUpsilonPrimeCandidates(Int_t Ptr) {fNumberOfUpsilonPrimeCandidates = Ptr;}
		void          SetNumOfUpsilonDoublePrimeCandidates(Int_t Ptr) {fNumberOfUpsilonDoublePrimeCandidates = Ptr;}
		void          SetNumOfCharmCandidates(Int_t Ptr) {fNumberOfCharmParticleCandidates = Ptr;}
		void          SetNumOfBeautyCandidates(Int_t Ptr) {fNumberOfBeautyParticleCandidates = Ptr;}

		void          SetNumOfElectrons(Int_t Ptr) {fNumberOfElectrons = Ptr;}
		void          SetNumOfMuons(Int_t Ptr) {fNumberOfMuons = Ptr;}
		void          SetNumOfPions(Int_t Ptr) {fNumberOfPions = Ptr;}
		void          SetNumOfKaons(Int_t Ptr) {fNumberOfKaons = Ptr;}
		void          SetNumOfProtons(Int_t Ptr) {fNumberOfProtons = Ptr;}
		void          SetNumOfLambdas(Int_t Ptr) {fNumberOfLambdas = Ptr;}

		void          SetK0Peak(Float_t Ppeak) {fK0PeakPosition = Ppeak;}
		void          SetK0Width(Float_t Pw) {fK0PeakWidth = Pw;}

		void          SetTotalMomentum(Float_t P) {fTotalP = P;}
		void          SetMeanPt(Float_t Pt) {fMeanPt = Pt;}
		void          SetMaxPt(Float_t Pt) {fMaxPt = Pt;}

		void          SetFlowV1(Float_t Pv1) {fFlowV1 = Pv1;}
		void          SetFlowV2(Float_t Pv2) {fFlowV2 = Pv2;}




		Int_t         GetEventId() {return fAliceEventId;}
		Int_t         GetGUID() {return fGUID;}

		Int_t         GetNumOfParticipants() {return fNumberOfParticipants;}
		Float_t       GetImpactParameter() {return fImpactParameter;}

		Float_t       GetVertexX() {return fPrimaryVertexX;}
		Float_t       GetVertexY() {return fPrimaryVertexY;}
		Float_t       GetVertexZ() {return fPrimaryVertexZ;}

		Int_t         GetTrigger() {return fTriggerInfo;}

		Float_t       GetZDCNeutronEnergy() {return fZDCNeutronEnergy;}
		Float_t       GetZDCProtonEnergy() {return fZDCProtonEnergy;}
		Float_t       GetZDCEMEnergy() {return fZDCEMEnergy;}

		Float_t       GetT0VertexZ() {return fT0VertexZ;}

		Int_t         GetNumOfTracks() {return fNumberOfTracks;}
		Int_t         GetNumOfPosTracks() {return fNumberOfPositiveTracks;}
		Int_t         GetNumOfNegTracks() {return fNumberOfNegativeTracks;}
		Int_t         GetNumOfNeutrTracks() {return fNumberOfNeutralTracks;}

		Int_t         GetNumOfV0s() {return fNumberOfV0s;}
		Int_t         GetNumOfCascades() {return fNumberOfCascades;}
		Int_t         GetNumOfKinks() {return fNumberOfKinks;}

		Int_t         GetNumOfPMDTracks() {return fNumberOfPMDTracks;}
		Int_t         GetNumOfFMDTracks() {return fNumberOfFMDTracks;}
		Int_t         GetNumOfPHOSClusters() {return fNumberOfPHOSClusters;}
		Int_t         GetNumOfEMCALClusters() {return fNumberOfEMCALClusters;}

		Int_t         GetNumOfJetCandidates() {return fNumberOfJetCandidates;}
		Int_t         GetNumOfHardPhotonsCandidates() {return fNumberOfHardPhotonsCandidates;}
		Int_t         GetNumOfJPsiCandidates() {return fNumberOfJPsiCandidates;}
		Int_t         GetNumOfPsiPrimeCandidates() {return fNumberOfPsiPrimeCandidates;}
		Int_t         GetNumOfUpsilonCandidates() {return fNumberOfUpsilonCandidates;}
		Int_t         GetNumOfUpsilonPrimeCandidates() {return fNumberOfUpsilonPrimeCandidates;}
		Int_t         GetNumOfUpsilonDoublePrimeCandidates() {return fNumberOfUpsilonDoublePrimeCandidates;}
		Int_t         GetNumOfCharmCandidates() {return fNumberOfCharmParticleCandidates;}
		Int_t         GetNumOfBeautyCandidates() {return fNumberOfBeautyParticleCandidates;}

		Int_t         GetNumOfElectrons() {return fNumberOfElectrons;}
		Int_t         GetNumOfMuons() {return fNumberOfMuons;}
		Int_t         GetNumOfPions() {return fNumberOfPions;}
		Int_t         GetNumOfKaons() {return fNumberOfKaons;}
		Int_t         GetNumOfProtons() {return fNumberOfProtons;}
		Int_t         GetNumOfLambdas() {return fNumberOfLambdas;}

		Float_t       GetK0Peak() {return fK0PeakPosition;}
		Float_t       GetK0Width() {return fK0PeakWidth;}

		Float_t       GetTotalMomentum() {return fTotalP;}
		Float_t       GetMeanPt() {return fMeanPt;}
		Float_t       GetMaxPt() {return fMaxPt;}

		Float_t       GetFlowV1() {return fFlowV1;}
		Float_t       GetFlowV2() {return fFlowV2;}

	ClassDef(AliEventTag,2)  //(ClassName, ClassVersion)
};
//______________________________________________________________________________


//______________________________________________________________________________
class AliLHCTag : public TObject
{

	private:
		Char_t   fLHCState[50];                 //LHC run conditions
		Float_t  fLHCLuminosity;                //the value of the luminosity

	public:
		AliLHCTag();
		virtual ~AliLHCTag();

		void          SetLHCState(char *type) {strcpy(fLHCState,type);}
		void          SetLuminosity(Float_t lumin) {fLHCLuminosity = lumin;}
		void          SetLHCTag(Float_t lumin, char *type) {fLHCLuminosity = lumin; strcpy(fLHCState,type); }

		char         *GetLHCState() {return fLHCState;}
		Float_t       GetLuminosity() {return fLHCLuminosity;}

	ClassDef(AliLHCTag,1)  //(ClassName, ClassVersion)
};
//______________________________________________________________________________


//______________________________________________________________________________
class AliDetectorTag : public TObject
{

	private:
		Bool_t   fITS;                 //ITS active = 1
		Bool_t   fTPC;                 //TPC active = 1
		Bool_t   fTRD;
		Bool_t   fTOF;
		Bool_t   fHMPID;
		Bool_t   fPHOS;
		Bool_t   fZDC;
		Bool_t   fMUON;
		Bool_t   fABSORBER;
		Bool_t   fPMD;
		Bool_t   fRICH;
		Bool_t   fEMCAL;
		Bool_t   fVZERO;
		Bool_t   fTZERO;

		virtual void CopyTag(AliDetectorTag *DetTag);

	public:
		AliDetectorTag();
		AliDetectorTag(AliDetectorTag *t);
		virtual ~AliDetectorTag();

		void          SetITS(Int_t n) {fITS = n;}
		void          SetTPC(Int_t n) {fTPC = n;}
		void          SetTRD(Int_t n) {fTRD = n;}
		void          SetTOF(Int_t n) {fTOF = n;}
		void          SetHMPID(Int_t n) {fHMPID = n;}
		void          SetPHOS(Int_t n) {fPHOS = n;}
		void          SetZDC(Int_t n) {fZDC = n;}
		void          SetMUON(Int_t n) {fMUON = n;}
		void          SetABSORBER(Int_t n) {fABSORBER = n;}
		void          SetPMD(Int_t n) {fPMD = n;}
		void          SetRICH(Int_t n) {fRICH = n;}
		void          SetEMCAL(Int_t n) {fEMCAL = n;}
		void          SetVZERO(Int_t n) {fVZERO = n;}
		void          SetTZERO(Int_t n) {fTZERO = n;}

		Bool_t        GetITS() {return fITS;}
		Bool_t        GetTPC() {return fTPC;}
		Bool_t        GetTRD() {return fTRD;}
		Bool_t        GetTOF() {return fTOF;}
		Bool_t        GetHMPID() {return fHMPID;}
		Bool_t        GetPHOS() {return fPHOS;}
		Bool_t        GetZDC() {return fZDC;}
		Bool_t        GetMUON() {return fMUON;}
		Bool_t        GetABSORBER() {return fABSORBER;}
		Bool_t        GetPMD() {return fPMD;}
		Bool_t        GetRICH() {return fRICH;}
		Bool_t        GetEMCAL() {return fEMCAL;}
		Bool_t        GetVZERO() {return fVZERO;}
		Bool_t        GetTZERO() {return fTZERO;}

	ClassDef(AliDetectorTag,1)  //(ClassName, ClassVersion)
};
//______________________________________________________________________________


//______________________________________________________________________________
class AliRunTag : public TObject
{

	private:
		Int_t    fAliceRunId;                   //the run id
		Float_t  fAliceMagneticField;           //value of the magnetic field
		Int_t    fAliceRunStartTime;            //run start date
		Int_t    fAliceRunStopTime;             //run stop date
		Int_t    fAliceReconstructionVersion;   //reco version
		Bool_t   fAliceRunQuality;              //validation script
		Float_t  fAliceBeamEnergy;              //beam energy cm
		Char_t   fAliceBeamType[5];             //run type (pp, AA, pA)
		Int_t    fAliceCalibrationVersion;      //calibration version

		Int_t    fNumEvents;                    //number of events per file
		Int_t    fNumDetectors;                 //number of detector configs per file
		TClonesArray  *fEventTag;               //array with all event tags
		TClonesArray  *fDetectorTag;            //array with all the detector tags

		AliLHCTag   fLHCTag;

		static TClonesArray *fgEvents;
		static TClonesArray *fgDetectors;

	public:
		AliRunTag();
		virtual ~AliRunTag();

		void          SetRunId(Int_t Pid) {fAliceRunId = Pid;}
		void          SetMagneticField(Float_t Pmag) {fAliceMagneticField = Pmag;}
		void          SetRunStartTime(Int_t Pt0) {fAliceRunStartTime = Pt0;}
		void          SetRunStopTime(Int_t Pt1) {fAliceRunStopTime = Pt1;}
		void          SetRecoVersion(Int_t Pn) {fAliceReconstructionVersion = Pn;}
		void          SetRunQuality(Int_t Pn) {fAliceRunQuality = Pn;}
		void          SetBeamEnergy(Float_t PE) {fAliceBeamEnergy = PE;}
		void          SetBeamType(char *Ptype) {strcpy(fAliceBeamType,Ptype);}
		void          SetCalibVersion(Int_t Pn) {fAliceCalibrationVersion = Pn;}

		void          SetNEvents(Int_t Pn) { fNumEvents = Pn; }

		void          SetLHCTag(Float_t Plumin, char *type);
		void          SetDetectorTag(AliDetectorTag *t);
		void          AddEventTag(AliEventTag *t);
		void          Clear();


		Int_t         GetRunId() {return fAliceRunId;}
		Float_t       GetMagneticField() {return fAliceMagneticField;}
		Int_t         GetRunStartTime() {return fAliceRunStartTime;}
		Int_t         GetRunStopTime() {return fAliceRunStopTime;}
		Int_t         GetRecoVersion() {return fAliceReconstructionVersion;}
		Int_t         GetRunQuality() {return fAliceRunQuality;}
		Float_t       GetBeamEnergy() {return fAliceBeamEnergy;}
		char         *GetBeamType() {return fAliceBeamType;}
		Int_t         GetCalibVersion() {return fAliceCalibrationVersion;}

		Int_t         GetNEvents() const {return fNumEvents;}

		AliLHCTag       *GetLHCTag() { return &fLHCTag; }

	ClassDef(AliRunTag,1)  //(ClassName, ClassVersion)
};
//______________________________________________________________________________

#endif
