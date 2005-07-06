//Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
//See cxx source for full Copyright notice
//
//The purpose of this class is to permorm the ITS tracking.
//
//      The constructor has the task to inizialize some private members.
// The method DoTracking is written to be called by a macro. It gets the event
// number, the minimum and maximum order number of TPC tracks that are to be
// tracked trough the ITS, and the file where the recpoints are registered.
//
// The method AliiTStracking is a recursive function that performs the
// tracking trough the ITS
//
// The method Intersection found the layer, ladder and detector where the
// intersection take place and caluclate the cohordinates of this
// intersection.  It returns an integer that is 0 if the intersection has
// been found successfully.
//
// The two mwthods Kalmanfilter and kalmanfiltervert operate the
// kalmanfilter without and with the verteximposition respectively.

#ifndef ALIITSTRACKERV1_H
#define ALIITSTRACKERV1_H

#include <TObject.h>

class TObjArray;
class TVector;
class TMatrix;
class AliITSTrackV1;
class AliITSRad;
class TStopwatch;

class AliITSTrackerV1 : public TObject {
 public:
    AliITSTrackerV1();
    AliITSTrackerV1(Int_t evnumber, Bool_t flag);
    AliITSTrackerV1(const AliITSTrackerV1 &cobj);
    ~AliITSTrackerV1();
    AliITSTrackerV1 &operator=(const AliITSTrackerV1 &obj);
//  void DelMatrix(Int_t NumOfModules);
    void DoTracking(Int_t evNumber, Int_t minTr, Int_t maxTr, TFile *file, Bool_t realmass=0);
    void RecursiveTracking(TList *trackITSlist);
    Int_t Intersection(AliITSTrackV1 &track, Int_t layer,Int_t &ladder,
		       Int_t &detector); 
    void KalmanFilter(AliITSTrackV1 *newtrack, TVector &cluster,
		      Double_t sigma[2]);
    void KalmanFilterVert(AliITSTrackV1 *newtrack, TVector &cluster,
			  Double_t sigma[2]);
    //void KalmanFilterVert(AliITSTrackV1 *newtrack, TVector &cluster,
    //                      Double_t sigma[2], Double_t chi2pred);

 private:
    //AliITS* fITS;              //! pointer to AliITS
    AliITSDetTypeRec *fDetTypeRec;//ITS obj. for reconstruction
    AliITSTrackV1 *fresult;    // result is a pointer to the final best track
    Double_t fPtref;           // transvers momentum obtained from TPC tracking
    Double_t fChi2max;         //  chi2 cut  
    //Double_t fepsphi;  //eps for definition window in phi 
   // Double_t fepsz;  //eps for definition window in z 
    TClonesArray  *frecPoints;    // pointer to RecPoints
    Int_t **fvettid;           // flag vector of used clusters
    Bool_t fflagvert;          // a flag to impose or not the vertex constraint
    AliITSRad *frl;            // pointer to get the radiation lenght matrix
///      To be put in a specific geometric class	 
    Int_t fNlad[6];            // Number of ladders for a given layer
    Int_t fNdet[6];            // Number of detector for a given layer
    Double_t fAvrad[6];        // Average radius for a given layer
    Double_t fDetx[6];// Semidimension of detectors along x axis for a given layer
    Double_t fDetz[6];// Semidimension of detectors along z axis for a given layer
    Double_t **fzmin;// Matrix of zmin for a given layer and a given detector
    Double_t **fzmax;// Matrix of zmax for a given layer and a given detector
    Double_t **fphimin;// Matrix of phimin for a given layer and a given ladder
    Double_t **fphimax;// Matrix of phimax for a given layer and a given ladder
    Double_t **fphidet; // azimuthal angle for a given layer and a given ladder
	 Int_t *fNRecPoints;        // pointer to the vector giving the number of recpoints for a given module
	 /*
	 Float_t **fRecCylR;        // pointer to the matrix giving the R cylindric cohordinate of a recpoint
	 Float_t **fRecCylPhi;      // pointer to the matrix giving the Phi cylindric cohordinate of a recpoint
	 Float_t **fRecCylZ;        // pointer to the matrix giving the Z cylindric cohordinate of a recpoint	
	 */
	 Double_t **fRecCylR;        // pointer to the matrix giving the R cylindric cohordinate of a recpoint
	 Double_t **fRecCylPhi;      // pointer to the matrix giving the Phi cylindric cohordinate of a recpoint
	 Double_t **fRecCylZ;        // pointer to the matrix giving the Z cylindric cohordinate of a recpoint		  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Double_t fFieldFactor;     // Magnetic field factor 
 
    ClassDef(AliITSTrackerV1,2)  
};

#endif
