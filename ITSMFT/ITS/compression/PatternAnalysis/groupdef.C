/*
In this macro we create groups of patterns with similar characteristics. First we load the pattern database with information about MC-hit.
We do not group patterns with a frequency below a treshold (frequencyTreshold).
We define 7 different ID. 2 for the sigma of the MChit-COG, x and z direction (sigmaXID and sigma ZID);
2 for the distance between COG and centre of the pixel, x and z (shiftXID and shiftZID);
2 for the distance between MChit and centre of the pixel, x and z direction (biasXID ad biasZID).
1 ID for the number of fired pixels.
We define the binning (number/width of bins) and assign the previous IDs.
We decide which the grouping method we want to use: sigma, shift and number of pixels (kShift) or sigma, bias and number of pixels (kBias)
Finally we assign group ID. Frequent patterns (above the treshold), form a one-pattern group. The tohers are in the same group if they have
the same IDs.
A fil.txt is print qith he inormation for each apttern, including the IDs and the pattID of the patterns in the same group.

*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TObjArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPavesText.h"
#include "TLatex.h"
#include "TBits.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "../ITS/UPGRADE/AliITSUClusterPix.h"
#include "../ITS/UPGRADE/AliITSURecoLayer.h"
#include "../ITS/UPGRADE/AliITSURecoDet.h"
#include "../ITS/UPGRADE/AliITSUHit.h"
#include "../ITS/UPGRADE/AliITSUGeomTGeo.h"
#include "AliITSsegmentation.h"
#include "AliGeomManager.h"

#endif

TObjArray histoArr;
TObjArray *pattDB=0; //it is an array with all the patterns in TBits format (forom the most to the least frequent)
TVectorF* pattFR=0; //frequency of the pattern in the database
TVectorF* xCentrePix=0; //coordinate of the centre of the pixel containing the COG for the down-left corner in fracion of the pitch
TVectorF* zCentrePix=0;
TVectorF* xCentreShift=0; //Difference between x coordinate fo COG and the centre of the pixel containing the COG, in fraction of the pitch
TVectorF* zCentreShift=0;
TVectorF* NPix=0;//Number of fired pixels
TVectorF* NRow=0;//Number of rows of the pattern
TVectorF* NCol=0;//Number of columns of the pattern
TVectorF* DeltaZmean=0; //mean value of the difference between MChit and COG z coordinate in micron derived from a gaussian fit
TVectorF* DeltaXmean=0; //mean value of the difference between MChit and COG x coordinate in micron derived from a gaussian fit
TVectorF* DeltaZsigma=0; //sigma of the difference between MChit and COG z coordinate in micron derived from a gaussian fit
TVectorF* DeltaXsigma=0; //sigma of the difference between MChit and COG x coordinate in micron derived from a gaussian fit
TVectorF* DeltaZmeanErr=0;
TVectorF* DeltaXmeanErr=0;
TVectorF* DeltaZsigmaErr=0;
TVectorF* DeltaXsigmaErr=0;

//Defining the ID to create the groups.
TArrayI sigmaXID;
TArrayI sigmaZID;
TArrayI shiftXID;
TArrayI shiftZID;
TArrayI biasXID;
TArrayI biasZID;
TArrayI NPixID;
TArrayI groupID;
Float_t totalGroupFreq=0;

void LoadDB(const char* namefile);

Int_t nPatterns=0;

//Defining the frequency treshold under which ti group patterns
Float_t frequencyTreshold = 0.001;
//Defining the bins
Int_t NumberofSigmaXbins=1000;
Int_t NumberofSigmaZbins=1000;
Int_t NumberofShiftXbins=20;
Int_t NumberofShiftZbins=20;
Int_t NumberofBiasXbins=2000;
Int_t NumberofBiasZbins=2000;
//Defining the boundaries of the bins concerning the number of pixel
const Int_t	limitsnumber = 4;
//Defining the width of the bins to store patterns with similar sigma/shift/bias (in micron)
Float_t sigmaXwidth=3;
Float_t sigmaZwidth=3;
Float_t biasXwidth=5;
Float_t biasZwidth=5;
Float_t shiftXwidth=1./NumberofShiftXbins; //(fraction of the pitch)
Float_t shiftZwidth=1./NumberofShiftZbins; //(fraction of the pitch)

Int_t nPixlimits[limitsnumber] = {4,10,25,50};

enum{kShift=0, kBias=1};
//Select kShift to group wrt sigma & COG-centreOfThePixel distance
//Select kBias to group wrt sigma & MChit-centreOfThePixel distance

Int_t groupingMethod = kShift;

Int_t tempgroupID=0;

void groupdef(){

	//Importing Data
	LoadDB("clusterTopology2.root");

	//Define moudule segmentation

	AliGeomManager::LoadGeometry("geometry.root");
	AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE);
	AliITSUClusterPix::SetGeom(gm);
	const AliITSsegmentation* segm = gm->GetSegmentation(0);
	Float_t pitchX = segm->Dpx(0)*10000; // for pitch in X in micron
	printf("pitchX: %f\n",pitchX);
	Float_t pitchZ = segm->Dpz(0)*10000; // for pitch in Z in micron
	printf("pitchX: %f\n",pitchX);
	delete gm;

	//Setting the number of patterns
	nPatterns = (pattDB->GetEntriesFast());

	sigmaXID.Set(nPatterns);
	sigmaZID.Set(nPatterns);

	//Assign -2 to frequent patterns, not to group, and -1 to rare clusters, -3 to patterns with sigma set to zero
	for(Int_t ID=0; ID<nPatterns; ID++){
		printf("Assign temporary sigma ID to pattern %d... ",ID);
 		if((*pattFR)[ID]>frequencyTreshold){
			sigmaXID[ID]=-2;//In order not to consider it within the groups
			sigmaZID[ID]=-2;//In order not to consider it within the groups
		}
		else{
			sigmaXID[ID]=-1;
			sigmaZID[ID]=-1;
			totalGroupFreq+=(*pattFR)[ID];
		}
		printf("done\n\n");
	}

	//Assign to similar patterns the same sigmaXID
	for(Int_t i=0; i<nPatterns; i++){
		printf("Assign sigmaXID to pattern %d... ",i);
		if(sigmaXID[i]==-1){
			if((*DeltaXsigma)[i]==0){
				sigmaXID[i]=-3; // In order not to cosider it within the groups
			}
			else{
				for(Int_t j=0; j<NumberofSigmaXbins; j++){
					if(j*sigmaXwidth < (*DeltaXsigma)[i] && (*DeltaXsigma)[i]<= (j+1)*sigmaXwidth){
						sigmaXID[i]=j+1;
						break;
					}
				}
			}
		}
		printf("done!!\n\n");
	}

	//Assign to similar patterns the same sigmaZID
	for(Int_t i=0; i<nPatterns; i++){
		printf("Assign sigmaZID to pattern %d... ",i);
		if(sigmaZID[i]==-1){
			if((*DeltaZsigma)[i]==0){
				sigmaZID[i]=-3; // In order not to cosider it within the groups
			}
			else{
				for(int j=0; j<NumberofSigmaZbins ; j++){
					if(j*sigmaZwidth < (*DeltaZsigma)[i] && (*DeltaZsigma)[i]<= (j+1)*sigmaZwidth){
						sigmaZID[i]=j+1;
						break;
					}
				}
			}
		}
		printf("done!!\n\n");
	}

	//assigning shiftID

	shiftXID.Set(nPatterns);
	shiftZID.Set(nPatterns);

	for(Int_t i=0; i<nPatterns; i++){
		printf("Assign shiftXID to pattern %d... ",i);
		
		for(int j=-NumberofShiftXbins/2; j<NumberofShiftXbins/2; j++){
			
			if(j*shiftXwidth < (*xCentreShift)[i] && (*xCentreShift)[i]<= (j+1)*shiftXwidth){
				shiftXID[i]=j+1;
				printf("done!!\n\n");
				break;
			}	
		}
	}

	for(Int_t i=0; i<nPatterns; i++){
		printf("Assign shiftZID to pattern %d... ",i);
		
		for(int j=-NumberofShiftZbins/2; j<NumberofShiftZbins/2; j++){
			if(j*shiftZwidth < (*zCentreShift)[i] && (*zCentreShift)[i]<= (j+1)*shiftZwidth){
				shiftZID[i]=j+1;
				printf("done!!\n\n");
				break;
			}
		}	
	}
	
	//assigning BiasID

	biasXID.Set(nPatterns);
	biasZID.Set(nPatterns);

	//Setting all the bias ID to zero

	for(Int_t i=0; i<nPatterns; i++){
		biasXID[i]=0;
		biasZID[i]=0;
	}

	for(Int_t i=0; i<nPatterns; i++){
		printf("Assign biasXID to pattern %d... ",i);
		for(Int_t j=-NumberofBiasXbins/2; j<NumberofBiasXbins/2; j++){
			if(j*biasXwidth < ((*DeltaXmean)[i]+((*xCentreShift)[i]*pitchX))
				&& ((*DeltaXmean)[i]+((*xCentreShift)[i]*pitchX))<= (j+1)*biasXwidth){
			biasXID[i]=j+1;
			break;
			}
		}
		printf("done!!\n\n");
	}

	for(Int_t i=0; i<nPatterns; i++){
		biasXID[i]=0;
		biasZID[i]=0;
	}

	for(Int_t i=0; i<nPatterns; i++){
		printf("Assign biasZID to pattern %d... ",i);
		for(Int_t j=-NumberofBiasZbins/2; j<NumberofBiasZbins/2; j++){
			if(j*biasZwidth < ((*DeltaZmean)[i]+((*zCentreShift)[i]*pitchZ)) 
				&& ((*DeltaZmean)[i]+((*zCentreShift)[i]*pitchZ))<= (j+1)*biasZwidth){
			biasZID[i]=j+1;
			break;
			}
		}
		printf("done!!\n\n");
	}

	

	//Assigning NPixID

	NPixID.Set(nPatterns);

	for(Int_t i=0; i<nPatterns; i++){
		printf("Assigning NPixID to pattern %d...", i);
		if((*NPix)[i]<=nPixlimits[0]){
			NPixID[i]=0;
			printf("done!!\n\n");
		}
		else if(nPixlimits[0]<(*NPix)[i] && (*NPix)[i]<=nPixlimits[limitsnumber-1]){
			for(Int_t j=0; j<nPatterns; j++ ){
				if(nPixlimits[j]<(*NPix)[i] && (*NPix)[i]<=nPixlimits[j+1]){
					NPixID[i]=j+1;
					printf("done!!\n\n");
					break;
				}
			}
		}
		else if((*NPix)[i]>nPixlimits[limitsnumber-1]){
			NPixID[i]=limitsnumber+1;
			printf("done!!\n\n");
		}
	}

	//Assigning groupID

	groupID.Set(nPatterns);

	//Assign -2 to frequent patterns, not to group, and -1 to rare clusters
	for(Int_t ID=0; ID<nPatterns; ID++){
		printf("Assign temporary group ID to pattern %d... ",ID);
		groupID[ID]=-1;
		printf("done\n\n");
	}
	Int_t k=0;
	while((*pattFR)[k]>frequencyTreshold){
		groupID[k]=tempgroupID;
		tempgroupID++;
		k++;
	}

	if(groupingMethod==kShift){
		for(Int_t i=0; i<nPatterns; i++){
			if(groupID[i]!=-1) continue;	
			groupID[i]=tempgroupID;
			printf("Assigning group ID %d... ",tempgroupID);
			for(Int_t j=i+1; j<nPatterns; j++){
				if(sigmaXID[j]==-3){
					groupID[j]=-1;
					continue;
				}
				else if(sigmaXID[j]==sigmaXID[i] && sigmaZID[j]==sigmaZID[i] && 
					shiftXID[j]==shiftXID[i] && shiftZID[j]==shiftZID[i] &&
					NPixID[i]==NPixID[j]) groupID[j]=tempgroupID;
			}
			printf("done!!\n\n");
			tempgroupID++;
		}
	}
	else if(groupingMethod==kBias){
		for(Int_t i=0; i<nPatterns; i++){
			if(groupID[i]!=-1) continue;
			groupID[i]=tempgroupID;
			printf("Assigning group ID %d... ",tempgroupID);
			for(Int_t j=i+1; j<nPatterns; j++){
				if(sigmaZID[j]==-3){
					groupID[j]=-1;
					continue;
				}
				else if(sigmaXID[j]==sigmaXID[i] && sigmaZID[j]==sigmaZID[i]
					&& biasXID[j]==biasXID[i] && biasZID[j]==biasZID[i]
					&& NPixID[i]==NPixID[j]) groupID[j]=tempgroupID;
			}
			printf("done!!\n\n");
			tempgroupID++;
		}
	}

	ofstream a("groupdef.txt");
		
	//setw(55) << "patterns in the group\n" << endl;

	a << Form("A NEGATIVE ID means that the pattern is not in a group.") << endl << endl <<
		Form("EXCEPTION: biasID CAN be negative") << endl <<
		"\n\n......................................................................................." << 
		"................................................................................................\n\n";

	for(int i=0; i<nPatterns; i++){

		printf("Writing info about pattern %d ...", i);

		a <<  setw(30) << Form("pattID: %d",i) <<  setw(30) << Form("freq: %f", (*pattFR)[i]) << setw(30) << Form("NPix: %d",Int_t((*NPix)[i]))<< setw(45) << Form("NRow: %d ",Int_t((*NRow)[i])) << setw(45) <<
		Form("NCol: %d", Int_t((*NCol)[i])) << endl << endl
		<< setw(45) << Form("DeltaXmean: %f (%f)",(*DeltaXmean)[i],(*DeltaXmeanErr)[i]) << setw(45) << 
		Form("DeltaZmean: %f (%f)",(*DeltaZmean)[i],(*DeltaZmeanErr)[i])<<
		setw(45) << Form("DeltaXsigma: %f (%f)",(*DeltaXsigma)[i],(*DeltaXsigmaErr)[i]) << setw(45) <<
		Form("DeltaZsigma: %f (%f)",(*DeltaZsigma)[i],(*DeltaZsigmaErr)[i]) << endl << endl << setw(30) << 
		Form("xShift: %f",(*xCentreShift)[i]) <<  setw(45) << Form("zShift: %f",(*zCentreShift)[i]) << setw(45) << Form("BiasX(MChit-centrePix): %f",(*DeltaXmean)[i]+((*xCentreShift)[i]*pitchX)) << setw(45) <<
		Form("BiasZ(MChit-centrePix): %f",(*DeltaZmean)[i]+((*zCentreShift)[i]*pitchZ)) << endl << endl <<
		setw(30) << Form("sigmaXID: %d",sigmaXID[i]) << setw(15) << Form("sigmaZID: %d",sigmaZID[i]) << setw(15) <<
		Form("shiftXID: %d",shiftXID[i]) << setw(15) << Form("shiftZID: %d",shiftZID[i]) << setw(15)<< Form("biasXID: %d",biasXID[i]) <<
		setw(15) << Form("biasZID: %d",biasZID[i]) << setw(15) << Form("NPixID: %d", NPixID[i]) << endl << endl << setw(30) << Form("groupID: %d", groupID[i])
		<< endl << endl << setw(30) << "patterns in this group:\n\n";

		Int_t matchcounter=0;

		for(int j=0; j<nPatterns;  j++){
			if(matchcounter!=0 && (matchcounter%20)==0){
				a << endl;
				matchcounter=0;
			}
			if(groupID[j]==groupID[i]) {
				a << j << ", ";
				matchcounter++;
			}
		}
		
		a << "\n\n......................................................................................." << 
		"................................................................................................\n\n";


		printf("done!\n\n");
	}

	a.close();

	printf("%d groups found!!!!\n\n",tempgroupID);

	printf("\n\nThe total frequency of the patterns in group is %f\n\n",totalGroupFreq);

	/*

	TVectorF groupDeltaX(tempGID+100);
	TVectorF groupShiftX(tempGID+100);
	TVectorF groupDeltaZ(tempGID+100);
	TVectorF groupShiftZ(tempGID+100);
	TVectorF patternNum(tempGID+100);
	groupDeltaX.Zero();
	groupShiftX.Zero();
	groupDeltaZ.Zero();
	groupShiftZ.Zero();
	patternNum.Zero();

	ofstream b("groupDB.txt");

	b<<setw(15)<<"groupID"<<setw(25)<<"patterns in the group"<<setw(15)<<
	"DeltaX"<<setw(15)<<"ShiftX"<<setw(15)<<"DeltaZ"<<setw(15)<<"ShiftZ\n"<<endl;

	TCanvas* c = new TCanvas("c","Patterns groups",900,600);
	TH1F* h = new TH1F("h","Patterns groups",tempGID,-0.5, tempGID-0.5);
	h->GetXaxis()->SetTitle("group ID");
	h->SetStats(0);

	for(Int_t gid=0; gid<tempGID; gid++){

		Int_t PatternNumberInGroup=0;
		Float_t freqSum=0.;//It is the sum of the frequencies of the patterns in the same group

		for(Int_t pattID=0; pattID<nPatterns; pattID++){
			if (groupID[pattID]==gid)
			{
				groupDeltaX[gid]+=(*DeltaXsigma)[pattID]*(*pattFR)[pattID];
				groupShiftX[gid]+=(*xCentreShift)[pattID]*(*pattFR)[pattID];
				groupDeltaZ[gid]+=(*DeltaZsigma)[pattID]*(*pattFR)[pattID];
				groupShiftZ[gid]+=(*zCentreShift)[pattID]*(*pattFR)[pattID];
				freqSum+=(*pattFR)[pattID];
				PatternNumberInGroup++;

				h->Fill(gid);
			}
		}

		groupDeltaX[gid]=groupDeltaX[gid]/freqSum;
		groupShiftX[gid]=groupShiftX[gid]/freqSum;
		groupDeltaZ[gid]=groupDeltaZ[gid]/freqSum;
		groupShiftZ[gid]=groupShiftZ[gid]/freqSum;
		patternNum[gid]=PatternNumberInGroup;

		b<<setw(15)<<gid<<setw(25)<<patternNum[gid]<<setw(15)<<
		groupDeltaX[gid]<<setw(15)<<groupShiftX[gid]<<setw(15)<<
		groupDeltaZ[gid]<<setw(15)<<groupShiftZ[gid]<<"\n"<<endl;
	}

	c->cd();
	h->Draw();

	TPavesText* info = new TPavesText(0.5,1,0.5,1,1,"nb");
	info->AddText(Form("Number of groups: %d",tempGID));
	info->AddText("#delta_{X}<5 & #delta_{Z}<5");
	info->AddText("#DeltaX<0.05 & #DeltaZ<0.05");
	c->cd();
	info->Draw();

	b.close();
	*/
}

void LoadDB(const char* fname){

  printf("\n\nLoading DB... ");
  // load database
  TFile* fl = TFile::Open(fname);
  if(!fl){printf("Could not find %s",fname); exit(1);}
  pattDB = (TObjArray*)fl->Get("TopDB");
  pattFR = (TVectorF*)fl->Get("TopFreq");
  xCentrePix =(TVectorF*)fl->Get("xCOG");
  zCentrePix =(TVectorF*)fl->Get("zCOG");
  xCentreShift =(TVectorF*)fl->Get("xShift");
  zCentreShift =(TVectorF*)fl->Get("zShift");
  NPix =(TVectorF*)fl->Get("NPix");
  NCol =(TVectorF*)fl->Get("NCol");
  NRow =(TVectorF*)fl->Get("NRow");
  DeltaXmean =(TVectorF*)fl->Get("DeltaXmean");
  DeltaZmean =(TVectorF*)fl->Get("DeltaZmean");
  DeltaXsigma =(TVectorF*)fl->Get("DeltaXsigma");
  DeltaZsigma =(TVectorF*)fl->Get("DeltaZsigma");
  DeltaXmeanErr =(TVectorF*)fl->Get("DeltaXmeanErr");
  DeltaZmeanErr =(TVectorF*)fl->Get("DeltaZmeanErr");
  DeltaZsigmaErr =(TVectorF*)fl->Get("DeltaXsigmaErr");
  DeltaXsigmaErr =(TVectorF*)fl->Get("DeltaZsigmaErr");
  printf("done!!\n\n");
}