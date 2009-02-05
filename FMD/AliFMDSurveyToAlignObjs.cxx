#include "AliFMDSurveyToAlignObjs.h"
#include "AliLog.h"
#include "AliSurveyPoint.h"
#include <TGraph2DErrors.h>
#include <TF2.h>
#include <TVector3.h>
#include <iostream>
#include <iomanip>
#include <TMath.h>
#include <TRotation.h>
#include <TGeoMatrix.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include "AliFMDGeometry.h"

//____________________________________________________________________
Double_t
AliFMDSurveyToAlignObjs::GetUnitFactor() const
{
  // Returns the conversion factor from the measured values to
  // centimeters. 
  if (!fSurveyObj) return 0;
  TString units(fSurveyObj->GetUnits());
  if      (units.CompareTo("mm", TString::kIgnoreCase) == 0) return .1;
  else if (units.CompareTo("cm", TString::kIgnoreCase) == 0) return 1.;
  else if (units.CompareTo("m",  TString::kIgnoreCase) == 0) return 100.;
  return 1;
}

//____________________________________________________________________
Bool_t
AliFMDSurveyToAlignObjs::GetPoint(const char* name, 
				  TVector3&   point, 
				  TVector3&   error) const
{
  // Get named point.   On return, point will contain the point
  // coordinates in centimeters, and error will contain the
  // meassurement errors in centimeters too.  If no point is found,
  // returns false, otherwise true. 
  if (!fSurveyPoints) return kFALSE;
  
  Double_t unit = GetUnitFactor();
  if (unit == 0) return kFALSE;
  
  TObject* obj  = fSurveyPoints->FindObject(name);
  if (!obj) return kFALSE;
  
  AliSurveyPoint* p = static_cast<AliSurveyPoint*>(obj);
  point.SetXYZ(unit * p->GetX(), 
	       unit * p->GetY(),
	       unit * p->GetZ());
  error.SetXYZ(unit * p->GetPrecisionX(),
	       unit * p->GetPrecisionY(),
	       unit * p->GetPrecisionZ());
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDSurveyToAlignObjs::CalculatePlane(const     TVector3& a, 
					const     TVector3& b,
					const     TVector3& c, 
					Double_t  depth,
					Double_t* trans,
					Double_t* rot) const
{
  // Vector a->b, b->c, and normal to plane defined by these two
  // vectors. 
  TVector3 ab(b-a), bc(c-b);
  
  // Normal vector to the plane of the fiducial marks obtained
  // as cross product of the two vectors on the plane d0^d1
  TVector3 nn(ab.Cross(bc));
  if (nn.Mag() < 1e-8) { 
    Info("CalculatePlane", "Normal vector is null vector");
    return kFALSE;
  }
  
  // We express the plane in Hessian normal form.
  //
  //   n x = -p,
  //
  // where n is the normalised normal vector given by 
  // 
  //   n_x = a / l,   n_y = b / l,   n_z = c / l,  p = d / l
  //
  // with l = sqrt(a^2+b^2+c^2) and a, b, c, and d are from the
  // normal plane equation 
  //
  //  ax + by + cz + d = 0
  // 
  // Normalize
  TVector3 n(nn.Unit());
  // Double_t p = - (n * a);
  
  // The center of the square with the fiducial marks as the
  // corners.  The mid-point of one diagonal - md.  Used to get the
  // center of the surveyd box. 
  TVector3 md(a + c);
  md *= 1/2.;

  // The center of the box. 
  TVector3 orig(md - depth * n);
  trans[0] = orig[0];
  trans[1] = orig[1];
  trans[2] = orig[2];
  
  // Normalize the spanning vectors 
  TVector3 uab(ab.Unit());
  TVector3 ubc(bc.Unit());
  
  for (size_t i = 0; i < 3; i++) { 
    rot[i * 3 + 0] = ubc[i];
    rot[i * 3 + 1] = uab[i];
    rot[i * 3 + 2] = n[i];
  }
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDSurveyToAlignObjs::FitPlane(const TObjArray& points, 
				  const TObjArray& errors,
				  Double_t         /* depth */,
				  Double_t*        trans,
				  Double_t*        rot) const
{
  Int_t nPoints = points.GetEntries();
  if (nPoints < 4) { 
    AliError(Form("Cannot fit a plane equation to less than 4 survey points, "
		  "got only %d", nPoints));
    return kFALSE;
  }
  
  TGraph2DErrors g;
  // Loop and fill graph 
  for (int i = 0; i < nPoints; i++) {
    TVector3* p = static_cast<TVector3*>(points.At(i));
    TVector3* e = static_cast<TVector3*>(errors.At(i));
  
    if (!p || !e) continue;
    
    g.SetPoint(i, p->X(), p->Y(), p->Z());
    g.SetPointError(i, e->X(), e->Y(), e->Z());
  }

  // Check that we have enough points
  if (g.GetN() < 4) { 
    AliError(Form("Only got %d survey points - no good for plane fit",
		  g.GetN()));
    return kFALSE;
  }

  // Next, declare fitting function and fit to graph. 
  // Fit to the plane equation: 
  // 
  //   ax + by + cz + d = 0
  //
  // or 
  // 
  //   z = - ax/c - by/c - d/c
  //
  TF2 f("plane", "-[0]*x-[1]*y-[2]", 
        g.GetXmin(), g.GetXmax(), g.GetYmin(), g.GetYmax());
  g.Fit(&f, "Q");
  
  // Now, extract the normal and offset
  TVector3 nv(f.GetParameter(0), f.GetParameter(1), 1);
  TVector3 n(nv.Unit());
  Double_t p = -f.GetParameter(2);
  
  // Create two vectors spanning the plane 
  TVector3 a(1, 0, f.Eval(1, 0)-p);
  TVector3 b(0, -1, f.Eval(0, -1)-p);
  TVector3 ua(a.Unit());
  TVector3 ub(b.Unit());
  // Double_t angAb = ua.Angle(ub);
  // PrintVector("ua: ", ua);
  // PrintVector("ub: ", ub);
  // std::cout << "Angle: " << angAb * 180 / TMath::Pi() << std::endl;
    
  for (size_t i = 0; i < 3; i++) { 
    rot[i * 3 + 0] = ua[i];
    rot[i * 3 + 1] = ub[i];
    rot[i * 3 + 2] = n[i];
  }
    
  // The intersection of the plane is given by (0, 0, -d/c)
  trans[0] = 0;
  trans[1] = 0;
  trans[2] = p;
  
  return kTRUE;
}


//____________________________________________________________________
Bool_t
AliFMDSurveyToAlignObjs::MakeDelta(const char*  path, 
				   Double_t*    rot, 
				   Double_t*    trans, 
				   TGeoHMatrix& delta) const
{
  // Make delta transformation
  if (!gGeoManager)           return kFALSE;
  if (!gGeoManager->cd(path)) return kFALSE;
  
  
  TGeoMatrix*  global = gGeoManager->GetCurrentMatrix();
#if 0
  PrintRotation(Form("%s rot:", global->GetName()),global->GetRotationMatrix());
  PrintVector(Form("%s trans:", global->GetName()),global->GetTranslation());
#endif

  return MakeDelta(global, rot, trans, delta);
}

//____________________________________________________________________
Bool_t
AliFMDSurveyToAlignObjs::MakeDelta(TGeoMatrix*  global,
				   Double_t*    rot, 
				   Double_t*    trans, 
				   TGeoHMatrix& delta) const
{
  // Make delta transformation
  TGeoHMatrix* geoM = new TGeoHMatrix;
  geoM->SetTranslation(trans);
  geoM->SetRotation(rot);

  delta = global->Inverse();
  delta.MultiplyLeft(geoM);
  
  return true;
}

//____________________________________________________________________
Bool_t
AliFMDSurveyToAlignObjs::GetFMD1Plane(Double_t* rot, Double_t* trans) const
{

  // The possile survey points 
  TVector3  icb, ict, ocb, oct, dummy;
  Int_t     missing = 0;
  if (!GetPoint("V0L_ICB", icb, dummy)) missing++;
  if (!GetPoint("V0L_ICT", ict, dummy)) missing++;
  if (!GetPoint("V0L_OCB", ocb, dummy)) missing++;
  if (!GetPoint("V0L_OCT", oct, dummy)) missing++;

  // Check that we have enough points
  if (missing > 1) { 
    AliWarning(Form("Only got %d survey points - no good for FMD1 plane",
		    4-missing));
    return kFALSE;
  }

#if 0
  const char* lidN = "FMD1_lid_mat0";
  TGeoMatrix* lidM = static_cast<TGeoMatrix*>(gGeoManager->GetListOfMatrices()
					      ->FindObject(lidN));
  if (!lidM) { 
    AliError(Form("Couldn't find FMD1 lid transformation %s", lidN));
    return kFALSE;
  }

  const Double_t* lidT = lidM->GetTranslation();
  Double_t        lidZ = lidT[2];
  Double_t        off  = lidZ-3.3;
#else
  Double_t        off  = 0;
#endif

  if (!CalculatePlane(ocb, icb, ict, off, trans, rot)) return kFALSE;
  // PrintRotation("FMD1 rotation:",  rot);
  // PrintVector("FMD1 translation:", trans);

  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDSurveyToAlignObjs::DoFMD1()
{
  // Do the FMD1 stuff
  Double_t rot[9], trans[3];
  if (!GetFMD1Plane(rot, trans)) return kFALSE;
  // const char* path = "/ALIC_1/F1MT_1/FMD1_lid_0";
  
  // TGeoHMatrix delta;
  TGeoTranslation global(0,0,324.670);
  if (!MakeDelta(&global, rot, trans, fFMD1Delta)) 
    return kFALSE;
  
  // PrintRotation("FMD1 delta rotation:",  fFMD1Delta.GetRotationMatrix());
  // PrintVector("FMD1 delta translation:", fFMD1Delta.GetTranslation());

  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDSurveyToAlignObjs::GetFMD2Plane(Double_t* rot, Double_t* trans) const
{

  // The possible survey points 
  const char*    names[] = { "FMD2_ITOP",  "FMD2_OTOP", 
			     "FMD2_IBOTM", "FMD2_OBOTM", 
			     "FMD2_IBOT",  "FMD2_OBOT", 
			     0 };
  const char**   name    = names;

  TObjArray points;
  TObjArray errors;
  
  // Loop and fill graph 
  while (*name) {
    TVector3 p, e;
    if (!GetPoint(*name++, p, e)) continue;
    
    points.Add(new TVector3(p));
    errors.Add(new TVector3(e));
  }
  if (points.GetEntries() < 4) { 
    AliWarning(Form("Only got %d survey points - no good for FMD2 plane",
		    points.GetEntries()));
    return kFALSE;
  }

  return FitPlane(points, errors, 0, trans, rot);
}

#define M(I,J) rot[(J-1) * 3 + (I-1)]
//____________________________________________________________________
Bool_t
AliFMDSurveyToAlignObjs::DoFMD2()
{
  // Do the FMD2 stuff
  Double_t rot[9], trans[3];
  if (!GetFMD2Plane(rot, trans)) return kFALSE;
  // PrintRotation("FMD2 rotation:",  rot);
  // PrintVector("FMD2 translation:", trans);

  for (int i = 0; i < 3; i++) { 
    for (int j = 0; j < 3; j++) { 
      rot[i*3+j] = (i == j ? 1 : 0);
    }
  }
  trans[0] = trans[1] = 0;
  trans[2] += 0.015;
  // PrintRotation("FMD2 rotation:",  rot);
  // PrintVector("FMD2 translation:", trans);
  
  // TGeoHMatrix delta;
  if (!MakeDelta("/ALIC_1/F2MT_2/FMD2_support_0/FMD2_back_cover_2", 
		 rot, trans, fFMD2Delta)) return kFALSE;
  
  // PrintRotation("FMD2 delta rotation:",  fFMD2Delta.GetRotationMatrix());
  // PrintVector("FMD2 delta translation:", fFMD2Delta.GetTranslation());

  return kTRUE;
}

//____________________________________________________________________
void
AliFMDSurveyToAlignObjs::Run()
{
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Init();
  geom->InitTransformations();
  
  DoFMD1();
  DoFMD2();
}

//____________________________________________________________________
Bool_t 
AliFMDSurveyToAlignObjs::CreateAlignObjs()
{
  TClonesArray& array = *fAlignObjArray;
  Int_t         n     = array.GetEntriesFast();

  if (!fFMD1Delta.IsIdentity()) { 
    new (array[n++]) AliAlignObjParams("FMD1/FMD1_T", 0, fFMD1Delta, kTRUE);
    new (array[n++]) AliAlignObjParams("FMD1/FMD1_B", 0, fFMD1Delta, kTRUE);
  }
  if (!fFMD2Delta.IsIdentity()) { 
    new (array[n++]) AliAlignObjParams("FMD2/FMD2_T", 0, fFMD2Delta, kTRUE);
    new (array[n++]) AliAlignObjParams("FMD2/FMD2_B", 0, fFMD2Delta, kTRUE);
  }
  // array.Print();
  
  return kTRUE;
}

//____________________________________________________________________
void 
AliFMDSurveyToAlignObjs::PrintVector(const char* text, const TVector3& v)
{
  // Print a vector
  Double_t va[] = { v.X(), v.Y(), v.Z() };
  PrintVector(text, va);
}
//____________________________________________________________________
void 
AliFMDSurveyToAlignObjs::PrintVector(const char* text, const Double_t* v)
{
  // Print a vector
  std::cout << text 
	    << std::setw(15) << v[0] 
	    << std::setw(15) << v[1]
	    << std::setw(15) << v[2]
	    << std::endl;
}


//____________________________________________________________________
void 
AliFMDSurveyToAlignObjs::PrintRotation(const char* text, const Double_t* rot)
{
  // Print a rotation matrix
  std::cout << text << std::endl;
  for (size_t i = 0; i < 3; i++) { 
    for (size_t j = 0; j < 3; j++) 
      std::cout << std::setw(15) << rot[i * 3 + j];
    std::cout << std::endl;
  }
}

//____________________________________________________________________
//
// EOF
//
