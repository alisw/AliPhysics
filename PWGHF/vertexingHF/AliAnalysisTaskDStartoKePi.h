#ifndef ALIANALYSISTASKDSTARTOKEPI_H
#define ALIANALYSISTASKDSTARTOKEPI_H


class AliTRDdigitsManager;

#include "AliAnalysisTaskSE.h"

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

//#pragma link C++ class std::vector< std::vector< std::vector< std::vector< std::vector<TLorentzVector> > > > >+;


struct vectors_t
{
    std::vector< std::vector< std::vector< std::vector<Double_t> > > > vec_data;
};


//----------------------------------------------------------------------------------------
class AliDStarTrack : public TObject
{
private:
    // Track properties
    Int_t channel;
    TLorentzVector TLV_A; // TLorentzVector
    TLorentzVector TLV_B; // TLorentzVector
    TLorentzVector TLV_C; // TLorentzVector
    Float_t  dist_to_prim[2]; // decay length
    Float_t  dca_AB[2]; // dca between daughter tracks
    Float_t  dca_AB_to_prim[2]; // dca between mother particle and primary vertex
    Float_t  dca_A[2]; // dca of particle A to primary vertex
    Float_t  dca_B[2]; // dca of particle B to primary vertex
    Float_t  dca_C[2]; // dca of particle C to primary vertex
    Float_t  Pid_A[2]; // dE/dx, TOF
    Float_t  Pid_B[2]; // dE/dx, TOF
    Float_t  Pid_C[2]; // dE/dx, TOF
    UShort_t NITScls_A; // Number of ITS hits
    UShort_t NITScls_B; // Number of ITS hits
    UShort_t NITScls_C; // Number of ITS hits


public:
    AliDStarTrack() :
	channel(),TLV_A(),TLV_B(),TLV_C(),dist_to_prim(),dca_AB(),dca_AB_to_prim(),dca_A(),dca_B(),dca_C(),
	Pid_A(),Pid_B(),Pid_C(),NITScls_A(-1),NITScls_B(-1),NITScls_C(-1)
    {
    }
	~AliDStarTrack()
	{
	}

	// setters
	void set_channel(Int_t i)                 { channel = i;          }
        void set_TLV_part_A(TLorentzVector tlv)   { TLV_A   = tlv;        }
        void set_TLV_part_B(TLorentzVector tlv)   { TLV_B   = tlv;        }
        void set_TLV_part_C(TLorentzVector tlv)   { TLV_C   = tlv;        }
	void set_dist_to_prim(Float_t f1, Float_t f2)   { dist_to_prim[0]   = f1; dist_to_prim[1]   = f2; }
	void set_dca_AB(Float_t f1, Float_t f2)         { dca_AB[0]         = f1; dca_AB[1]         = f2; }
	void set_dca_AB_to_prim(Float_t f1, Float_t f2) { dca_AB_to_prim[0] = f1; dca_AB_to_prim[1] = f2; }
	void set_dca_A(Float_t f1, Float_t f2)          { dca_A[0]          = f1; dca_A[1]          = f2; }
	void set_dca_B(Float_t f1, Float_t f2)          { dca_B[0]          = f1; dca_B[1]          = f2; }
	void set_dca_C(Float_t f1, Float_t f2)          { dca_C[0]          = f1; dca_C[1]          = f2; }
	void set_Pid_A(Float_t f1, Float_t f2)          { Pid_A[0]          = f1; Pid_A[1]          = f2; }
	void set_Pid_B(Float_t f1, Float_t f2)          { Pid_B[0]          = f1; Pid_B[1]          = f2; }
        void set_Pid_C(Float_t f1, Float_t f2)          { Pid_C[0]          = f1; Pid_C[1]          = f2; }
	void set_NITScls_A(UShort_t s)            { NITScls_A = s;      }
	void set_NITScls_B(UShort_t s)            { NITScls_B = s;      }
	void set_NITScls_C(UShort_t s)            { NITScls_C = s;      }

	// getters
	Int_t get_channel() const                 { return channel;        }
	TLorentzVector get_TLV_part_A() const     { return TLV_A;          }
	TLorentzVector get_TLV_part_B() const     { return TLV_B;          }
	TLorentzVector get_TLV_part_C() const     { return TLV_C;          }
	Float_t get_dist_to_prim(Int_t i) const   { return dist_to_prim[i];   }
	Float_t get_dca_AB(Int_t i) const         { return dca_AB[i];         }
	Float_t get_dca_AB_to_prim(Int_t i) const { return dca_AB_to_prim[i]; }
	Float_t get_dca_A(Int_t i) const          { return dca_A[i];          }
	Float_t get_dca_B(Int_t i) const          { return dca_B[i];          }
	Float_t get_dca_C(Int_t i) const          { return dca_C[i];          }
	Float_t get_Pid_A(Int_t i) const          { return Pid_A[i];          }
	Float_t get_Pid_B(Int_t i) const          { return Pid_B[i];          }
	Float_t get_Pid_C(Int_t i) const          { return Pid_C[i];          }
	UShort_t get_NITScls_A() const            { return NITScls_A;      }
	UShort_t get_NITScls_B() const            { return NITScls_B;      }
	UShort_t get_NITScls_C() const            { return NITScls_C;      }

	Int_t     HasITShit_A_on_layer(Int_t ilayer) { return ((NITScls_A >> ilayer) & 1);}  // ITShit -> LOL
	Int_t     HasITShit_B_on_layer(Int_t ilayer) { return ((NITScls_B >> ilayer) & 1);}  // ITShit -> LOL
        Int_t     HasITShit_C_on_layer(Int_t ilayer) { return ((NITScls_C >> ilayer) & 1);}  // ITShit -> LOL

	ClassDef(AliDStarTrack,1);  // A simple track of a particle
};


class AliDStarEvent : public TObject
{
private:
    Float_t x; // Event vertex x
    Float_t y; // Event vertex y
    Float_t z; // Event vertex z
    Int_t   id; // Run id
    Int_t   N_tracks; // total number of tracks
    Int_t   N_TRD_tracklets; // total number of TRD tracklets
    Float_t   cent_class_ZNA; // ZDC neutral A
    Float_t   cent_class_ZNC; // ZDC neutral C
    Float_t   cent_class_V0A; // V0 A
    Float_t   cent_class_V0C; // V0 C
    Float_t   cent_class_V0M; // V0 average
    Float_t   cent_class_CL0; // clusters in layer 0
    Float_t   cent_class_CL1; // clusters in layer 1
    Float_t   cent_class_SPD; // SPD
    Float_t   cent_class_V0MEq; //
    Float_t   cent_class_V0AEq; //
    Float_t   cent_class_V0CEq; //


    Float_t BeamIntAA; // ZDC coincidence rate
    Float_t T0zVertex; // z-vertex position from VPD

    TString TriggerWord; // Trigger word

    UShort_t      fNumTracks; // number of tracks in event

    TClonesArray* fTracks;      //->

public:
    AliDStarEvent() :
	x(-1),y(-1),z(-1),id(-1),N_tracks(0),N_TRD_tracklets(0),
	cent_class_ZNA(0),cent_class_ZNC(0),cent_class_V0A(0),cent_class_V0C(0),cent_class_V0M(0),cent_class_CL0(0),cent_class_CL1(0),
	cent_class_SPD(0),cent_class_V0MEq(0),cent_class_V0AEq(0),cent_class_V0CEq(0),BeamIntAA(-1),T0zVertex(-1),TriggerWord(),fNumTracks(0)
    {
	fTracks      = new TClonesArray( "AliDStarTrack", 10 );
    }
	~AliDStarEvent()
	{
	    delete fTracks;
	    fTracks = NULL;
	}

	void       setx(Float_t r)                    { x = r;                         }
	Float_t    getx() const                       { return x;                      }

	void       sety(Float_t r)                    { y = r;                         }
	Float_t    gety() const                       { return y;                      }

	void       setz(Float_t r)                    { z = r;                         }
	Float_t    getz() const                       { return z;                      }

	void       setid(Int_t  r)                    { id = r;                        }
	Int_t      getid() const                      { return id;                     }

	void       setN_tracks(Int_t r)                 { N_tracks = r;                    }
	Int_t      getN_tracks() const                    { return N_tracks;                 }

	void       setN_TRD_tracklets(Int_t r)                 { N_TRD_tracklets = r;                    }
	Int_t      getN_TRD_tracklets() const                    { return N_TRD_tracklets;                 }

	void       setcent_class_ZNA(Float_t r)             { cent_class_ZNA = r;                }
	Float_t      getcent_class_ZNA() const              { return cent_class_ZNA;             }

	void       setcent_class_ZNC(Float_t r)             { cent_class_ZNC = r;                }
	Float_t      getcent_class_ZNC() const              { return cent_class_ZNC;             }

	void       setcent_class_V0A(Float_t r)             { cent_class_V0A = r;                }
	Float_t      getcent_class_V0A() const              { return cent_class_V0A;             }

	void       setcent_class_V0C(Float_t r)             { cent_class_V0C = r;                }
	Float_t      getcent_class_V0C() const              { return cent_class_V0C;             }

	void       setcent_class_V0M(Float_t r)             { cent_class_V0M = r;                }
	Float_t      getcent_class_V0M() const              { return cent_class_V0M;             }

	void       setcent_class_CL0(Float_t r)             { cent_class_CL0 = r;                }
	Float_t      getcent_class_CL0() const              { return cent_class_CL0;             }

	void       setcent_class_CL1(Float_t r)             { cent_class_CL1 = r;                }
	Float_t      getcent_class_CL1() const              { return cent_class_CL1;             }

	void       setcent_class_SPD(Float_t r)             { cent_class_SPD = r;                }
	Float_t      getcent_class_SPD() const              { return cent_class_SPD;             }

	void       setcent_class_V0MEq(Float_t r)             { cent_class_V0MEq = r;                }
	Float_t      getcent_class_V0MEq() const              { return cent_class_V0MEq;             }

	void       setcent_class_V0AEq(Float_t r)             { cent_class_V0AEq = r;                }
	Float_t      getcent_class_V0AEq() const              { return cent_class_V0AEq;             }

	void       setcent_class_V0CEq(Float_t r)             { cent_class_V0CEq = r;                }
	Float_t      getcent_class_V0CEq() const              { return cent_class_V0CEq;             }

	void       setBeamIntAA(Float_t r)                 { BeamIntAA = r;                      }
	Float_t    getBeamIntAA() const                    { return BeamIntAA;                   }

	void       setT0zVertex(Float_t r)            { T0zVertex = r;                     }
	Float_t    getT0zVertex() const               { return T0zVertex;                  }

	void       setTriggerWord(TString s)          { TriggerWord = s;}
	TString    getTriggerWord() const             { return TriggerWord; }

	AliDStarTrack* createTrack()
	{
	    if (fNumTracks == fTracks->GetSize())
		fTracks->Expand( fNumTracks + 10 );
	    if (fNumTracks >= 10000)
	    {
		Fatal( "AliDStarEvent::createTrack()", "ERROR: Too many tracks (>10000)!" );
		exit( 2 );
	    }

	    new((*fTracks)[fNumTracks++]) AliDStarTrack;
	    return (AliDStarTrack*)((*fTracks)[fNumTracks - 1]);
	}
	void clearTrackList()
	{
	    fNumTracks   = 0;
	    fTracks      ->Clear();
	}
	UShort_t getNumTracks() const
	{
	    return fNumTracks;
	}
	AliDStarTrack* getTrack(UShort_t i) const
	{
	    return i < fNumTracks ? (AliDStarTrack*)((*fTracks)[i]) : NULL;
}

ClassDef(AliDStarEvent,1);  // A simple event compiled of tracks
};
//----------------------------------------------------------------------------------------



class AliAnalysisTaskDStartoKePi : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskDStartoKePi()
	: AliAnalysisTaskSE(),
	fESD(0x0),
	fAOD(0x0),
	fDigitsInputFileName("TRD.FltDigits.root"), fDigitsInputFile(0),
	fDigitsOutputFileName(""), fDigitsOutputFile(0),
	fDigMan(0),fGeo(0),
	AS_Event(0),AS_Track(0),
	Tree_AS_Event(0),fEventNoInFile(-2), N_good_events(0), fDigitsLoadedFlag(kFALSE),
	h_delta_dca(0x0),h_invariant_mass_DStar(0x0),h_invariant_mass_DStar_nom(0x0),h_invariant_mass_D0(0x0),h_delta_invariant_mass(0x0),h_multiplicities(0x0),
	TP_pt_bins(0x0),TP_dca_vs_p(0x0),h2D_mass2_vs_p(0x0),TP_counts_pid(0x0),h_pT(0x0),TP_TPC_chi2_vs_pT(0x0),h_runnumber(0x0)
    {
	cout << "" << endl;
	cout << "***************************************************************************************" << endl;
	cout << "In AliAnalysisTaskDStartoKePi.h constructor" << endl;
	cout << "fDigitsInputFileName: " << fDigitsInputFileName << endl;
	cout << "***************************************************************************************" << endl;
	cout << "" << endl;
    }
	AliAnalysisTaskDStartoKePi(const char *name);
	//virtual ~AliAnalysisTaskDStartoKePi() {}

	virtual void   UserCreateOutputObjects();
	virtual Bool_t UserNotify();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);
        void SetIsAOD(Bool_t val){fIsAOD = val;};

	void FillHelix(AliESDtrack* track_in, Double_t magF_in);
	void FindDCAHelixPoint(TVector3 space_vec, AliHelix helixA, Double_t path_initA, Double_t path_initB,
						  Double_t &pathA, Double_t &dca_xy, Double_t &dca_z);
        Double_t calcDeterminant(TVector3& v1,TVector3& v2,TVector3& v3);
        TVector3 calculatePointOfClosestApproach(TVector3 &base1, TVector3 &dir1,TVector3 &base2, TVector3 &dir2);
        TVector3 calculateCrossPoint(TVector3 &base1, TVector3 &dir1,TVector3 &base2, TVector3 &dir2);
        TVector3 calcVertexAnalytical(TVector3 &base1, TVector3 &dir1,TVector3 &base2, TVector3 &dir2);
        Double_t calculateMinimumDistanceStraightToPoint(TVector3 &base, TVector3 &dir,TVector3 &point);
	Double_t calculateMinimumDistance(TVector3 &base1, TVector3 &dir1,TVector3 &base2, TVector3 &dir2);
        TVector3 calc_straight_to_point_dca_vector(TVector3 base1, TVector3 dir1,TVector3 point);
        std::vector< std::vector<TLorentzVector> > Rec_two_prong_primary_decay(std::vector< std::vector<vectors_t> > vec_track_info_in,
								       Int_t Pid_A_in, Int_t Pid_B_in, Int_t z_bin, Int_t centrality_bin, Int_t buffer_bin, Int_t buffer_bin_ME, Int_t SE_ME_mode,
									       Int_t N_cuts, std::vector< std::vector< std::vector<Int_t> > > &vec_track_id);
	std::vector< std::vector<TLorentzVector> > Rec_three_prong_primary_decay(std::vector< std::vector<vectors_t> > vec_track_info_in,
										 std::vector< std::vector<vectors_t> > vec_two_prong_info_in,
										 Int_t i_channel,Int_t N_cuts, Int_t Pid_C_in, Int_t z_bin, Int_t centrality_bin, Int_t buffer_bin, Int_t buffer_bin_ME, Int_t SE_ME_mode,
										 std::vector< std::vector<TLorentzVector> > &vec_TLV_AB
										);



	void SetDigitsInputFilename(TString x)
	{
	    fDigitsInputFileName=x;
	    cout << "" << endl;
	    cout << "***************************************************************************************" << endl;
	    cout << "fDigitsInputFileName: " << fDigitsInputFileName << endl;
	    cout << "***************************************************************************************" << endl;
	}
	void SetDigitsOutputFilename(TString x) {fDigitsOutputFileName=x;}

	AliHelix aliHelix;

    protected:

	Bool_t NextEvent(Bool_t preload=kFALSE);

	AliESDEvent     *fESD;
	AliAODEvent     *fAOD;
	TList           *fListOfHistos;       //! list of output histograms
        TTree           *fTree;               //! output tree
	AliPIDResponse  *fPIDResponse;        //! PID handling
	AliESDtrackCuts *EsdTrackCuts;        //!
        Bool_t           fIsAOD;

    private:

	TFile* OpenDigitsFile(TString inputfile, TString digfile, TString opt);

	TString fDigitsInputFileName;         //! Name of digits file for reading
	TFile*  fDigitsInputFile;             //! Digits file for reading
	TString fDigitsOutputFileName;        //! Name of digits file for writing
	TFile*  fDigitsOutputFile;            //! Digits file for writing

	AliTRDdigitsManager* fDigMan; //! digits manager
	AliTRDgeometry* fGeo; //! TRD geometry
	AliDStarEvent* AS_Event;
	AliDStarTrack*    AS_Track;
        TTree       *Tree_AS_Event;

	std::vector< std::vector<vectors_t> > vec_track_info; // this is needed since a 6D std::vector does not compile within an ali root class
	std::vector< std::vector<vectors_t> > vec_two_prong_info;
	std::vector< std::vector<Int_t> >     vec_buffer_counter;

        TVector3 TV3_prim_vertex;

        TH1D* h_delta_dca;
	std::vector<TH1D*> h_invariant_mass;
	std::vector< std::vector< std::vector< std::vector<TH1D*> > > > h_invariant_mass_DStar;
        std::vector< std::vector< std::vector< std::vector<TH1D*> > > > h_invariant_mass_DStar_nom;
	std::vector< std::vector< std::vector< std::vector<TH1D*> > > > h_invariant_mass_D0;
	std::vector< std::vector< std::vector< std::vector<TH1D*> > > > h_delta_invariant_mass;
	std::vector<TH1D*> h_multiplicities;
	TProfile* TP_pt_bins;


	TProfile* TP_dca_vs_p;
	std::vector< std::vector<TH2D*> > h2D_mass2_vs_p;
	std::vector<TProfile*> TP_counts_pid;
	std::vector< std::vector<TH1D*> > h_pT;
	std::vector<TProfile*> TP_TPC_chi2_vs_pT;
        TH1D* h_runnumber;
	std::vector< std::vector<Int_t> > ME_buffer_counter;
	//std::vector< std::vector< std::vector< std::vector< std::vector<Double_t> > > > > TLV_ME_tracks;
	//std::vector< std::vector< std::vector< std::vector< std::vector<TLorentzVector> > > > >  TLV_ME_tracks;
	//std::vector< std::vector< std::vector< std::vector< std::vector< std::vector<TLorentzVector> > > > > > TLV_ME_tracks;
	//std::vector< std::vector< std::vector< std::vector<TLorentzVector> > > >  TLV_ME_tracks;
        //std::vector< std::vector< std::vector< std::vector< std::vector<TLorentzVector> > > > > TLV_ME_tracks_D0;
	Int_t fEventNoInFile;
	Int_t N_good_events;
	Int_t fDigitsLoadedFlag;

	AliAnalysisTaskDStartoKePi(const AliAnalysisTaskDStartoKePi&); // not implemented
	AliAnalysisTaskDStartoKePi& operator=(const AliAnalysisTaskDStartoKePi&); // not implemented

	ClassDef(AliAnalysisTaskDStartoKePi, 1);
};

#endif
