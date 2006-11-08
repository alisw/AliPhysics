#ifndef ALITPCLASERTRACKS_H
#define ALITPCLASERTRACKS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
//              Container class for laser track positions                 //
////////////////////////////////////////////////////////////////////////////

class TString;
//class TPolyLine3D;
class TObjArray;
class TGraph;
class TVector3;

#include <TPolyLine3D.h>

////////////////////////////////////////////////////////////////////////
//              Class AliTPCLaserTracks
////////////////////////////////////////////////////////////////////////
class AliTPCLaserTracks : public TObject {

public:
    AliTPCLaserTracks();
    AliTPCLaserTracks(Int_t npoints);
    virtual ~AliTPCLaserTracks();
    AliTPCLaserTracks(const AliTPCLaserTracks &param); // copy constructor
    AliTPCLaserTracks &operator = (const AliTPCLaserTracks & param);

    Int_t GetId()    {return fId;     }
    Int_t GetSide()  {return fSide;   }
    Int_t GetRod()   {return fRod;    }
    Int_t GetBundle(){return fBundle; }
    Int_t GetBeam()  {return fBeam;   }

    Double_t GetX()    {return fX;    }
    Double_t GetY()    {return fY;    }
    Double_t GetZ()    {return fZ;    }
    Double_t GetPhi()  {return fPhi;  }
    Double_t GetTheta(){return fTheta;}

    TPolyLine3D *GetLine();
    Int_t SetPoint(Int_t point, Double_t x, Double_t y, Double_t z);



    void SetId    (Int_t id)    {fId     = id;    }
    void SetSide  (Int_t side)  {fSide   = side;  }
    void SetRod   (Int_t rod)   {fRod    = rod;   }
    void SetBundle(Int_t bundle){fBundle = bundle;}
    void SetBeam  (Int_t beam)  {fBeam   = beam;  }

    void SetX    (Double_t x)    {fX     = x;    }
    void SetY    (Double_t y)    {fY     = y;    }
    void SetZ    (Double_t z)    {fZ     = z;    }
    void SetPhi  (Double_t phi)  {fPhi   = phi;  }
    void SetTheta(Double_t theta){fTheta = theta;}

//    void SetLine(TPolyLine3D *l) {fLine  = l;    }

    void WriteTreeDesignData();

    Int_t FindMirror(Char_t *file, Double_t x, Double_t y, Double_t z, Double_t phi);

    TObjArray* GetLines(Char_t* file, Char_t *cuts="");

protected:
    Int_t fId;              //Laser beam id            (0-335)
    Int_t fSide;            //TPC side; 0:Shaft Side (A) -- 1:Muon Side (C)
    Int_t fRod;             //Laser Rod                (0-5)
    Int_t fBundle;          //Mirror bundle in the Rod (0-3)
    Int_t fBeam;            //Laser Beam in the bundle (0-6)

    Double_t fX;             //X-position of the mirror in the bundle
    Double_t fY;             //Y-position of the mirror in the bundle
    Double_t fZ;             //Z-position of the mirror in the bundle
    Double_t fTime;          //time-position of the mirror in the bundle
    Double_t fPhi;           //Phi direction of the laser beam
    Double_t fTheta;         //Theta direction of the laser beam

//    TPolyLine3D *fLine;      //Line of the track in the TPC
    Int_t fMaxSize;          //! Number of points in array
    Int_t fNpoints;          //Number of points

    Double_t  *fXarr;            //[fNpoints] array of a line (x,y,z)
    Double_t  *fYarr;            //[fNpoints] array of a line (x,y,z)
    Double_t  *fZarr;            //[fNpoints] array of a line (x,y,z)

    void InitPoints();

    Double_t FindBeamLength(TVector3 vS, TVector3 vP);

    ClassDef(AliTPCLaserTracks,1)
};



#endif

