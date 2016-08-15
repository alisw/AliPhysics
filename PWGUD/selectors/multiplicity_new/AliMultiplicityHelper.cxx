#include "AliMultiplicityHelper.h"
#include <AliMultiplicity.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>
#include <AliMCEvent.h>
#include <AliESDEvent.h>
#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <TMath.h>
#include <TMathBase.h>
#include <AliLog.h>
#include <TRandom.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <AliGenCocktailEventHeader.h>
#include <iostream>
#include "AliAnalysisUtils.h"

ClassImp ( AliMultiplicityHelper )

AliMultiplicityHelper::AliMultiplicityHelper() {

}

AliMultiplicityHelper::~AliMultiplicityHelper() {

}

// MCProcessType AliMultiplicityHelper::GetProcessType ( const AliVEvent* mc ) {
// 	AliGenEventHeader* header = GetGenEventHeader ( mc );
// 	if ( !header )
// 		return kUndefined;
//
// 	AliGenPythiaEventHeader* pythiaGenHeader = 0;
// 	AliGenDPMjetEventHeader* dpmHeader = 0;
//
// 	pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*> ( header );
// 	dpmHeader = dynamic_cast<AliGenDPMjetEventHeader*> ( header );
//
// 	Bool_t isPythia8 = kFALSE;
// 	if ( pythiaGenHeader ) {
// 		isPythia8 = (pythiaGenHeader->ProcessType() > 100);
// 	}
//
// 	Int_t iSD1flag, iSD2flag, iDDflag;
// 	Int_t iMCProcessType = 0;
// 	if ( pythiaGenHeader ) {
// 		iMCProcessType = pythiaGenHeader->ProcessType();
// 		if ( !isPythia8 ) {
// 			iSD1flag = 92;
// 			iSD2flag = 93;
// 			iDDflag = 94;
// 		} else {
// 			iSD1flag = 103;
// 			iSD2flag = 104;
// 			iDDflag = 105;
// 		}
// 	} else
// 		if ( dpmHeader ) {
// 			iMCProcessType = dpmHeader->ProcessType();
// 			iSD1flag = 5;
// 			iSD2flag = 6;
// 			iDDflag = 7;
// 		} else
// 			return kUndefined;
//
// 	if ( iMCProcessType == iSD1flag || iMCProcessType == iSD2flag ) {
// 		AliInfoClass ( "Single-diffractive" );
// 		return kSingleDiffractive;
// 	} else
// 		if ( iMCProcessType == iDDflag ) {
// 			AliInfoClass ( "Double-diffractive" );
// 			return kDoubleDiffractive;
// 		} else {
// 			if ( (isPythia8 && ( iMCProcessType != 102 )) || (!isPythia8)  ) {
// 				AliInfoClass ( "INEL" );
// 				return kND;
// 			}
// 			AliInfoClass ( "Elastic?" );
// 			return kUndefined;
// 		}
// }
//_____________________________________________________________________________________________________________________________________________
const AliVertex* AliMultiplicityHelper::GetVertex ( VertexType vType, const AliVEvent* event ) {
	switch ( vType ) {
		case kVertexSPD:
			return dynamic_cast<const AliESDEvent*> ( event )->GetPrimaryVertexSPD();

		case kVertexTracks:
			return dynamic_cast<const AliESDEvent*> ( event )->GetPrimaryVertexTracks();

		case kVertexMC:
			return dynamic_cast<const AliVertex*> ( dynamic_cast<const AliMCEvent*> ( event )->GetPrimaryVertex() );

		case kPrimaryVertex:
			return dynamic_cast<const AliESDEvent*> ( event )->GetPrimaryVertex();
	}

	return 0;
}
//_____________________________________________________________________________________________________________________________________________
Bool_t AliMultiplicityHelper::IsQualityVertex ( const AliVertex* v, Bool_t& disp ) {
	if ( !v )
		return kFALSE;

	if ( !v->GetStatus() )
		return kFALSE;

	if ( v->IsFromVertexerZ() && ( v->GetDispersion() > 0.02 ) ) {
		disp = kTRUE;
		return kFALSE;
	}

	return kTRUE;
}
//_____________________________________________________________________________________________________________________________________________
Bool_t AliMultiplicityHelper::IsQualityVertex ( const AliVEvent* event ) {
	const AliESDEvent* esd = dynamic_cast<const AliESDEvent*> ( event );
	const AliVertex* vSPD = AliMultiplicityHelper::GetVertex ( kVertexSPD, esd );
	const AliVertex* vTracks = AliMultiplicityHelper::GetVertex ( kVertexTracks, esd );
	Bool_t disp = kFALSE;

	if ( IsQualityVertex ( vSPD, disp ) || IsQualityVertex ( vTracks, disp ) )
		return kTRUE;

	return kFALSE;
}
//_____________________________________________________________________________________________________________________________________________
// Bool_t AliMultiplicityHelper::IsSelected ( AliVEvent::EOfflineTriggerTypes trigger ) {
// 	AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*> ( AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler() );
// 	if ( !inputHandler ) {
// 		Printf ( "ERROR: Could not get input handler" );
// 		return kFALSE;
// 	}
// 	return inputHandler->IsEventSelected() & trigger;
// }
//_____________________________________________________________________________________________________________________________________________
Double_t AliMultiplicityHelper::GetTrackletChi2 ( const AliMultiplicity* m, Int_t iTracklet ) {

	if ( ( m->GetDPhiWindow2() < 1e-9 ) || ( m->GetDThetaWindow2() < 1e-9 ) )
		return -1;

	Double_t chi2phi = TMath::Abs ( m->GetDeltaPhi ( iTracklet ) ) - m->GetDPhiShift();
	chi2phi *= chi2phi;
	chi2phi /= m->GetDPhiWindow2();
	Double_t chi2theta = m->GetDeltaTheta ( iTracklet );
	chi2theta *= chi2theta;
	chi2theta /= m->GetDThetaWindow2();
	
	if ( !m->GetScaleDThetaBySin2T() ) {
		Double_t sin = TMath::Sin ( 0.5 * ( m->GetThetaSingleLr ( iTracklet, 0 ) + m->GetThetaSingleLr ( iTracklet, 1 ) ) );
		sin *= sin;
		chi2theta /= sin;
	}

	return chi2phi + chi2theta;
}
//_____________________________________________________________________________________________________________________________________________
// Bool_t AliMultiplicityHelper::CheckDiffraction ( const AliVEvent* mc ) {
// 	AliGenEventHeader* header = GetGenEventHeader ( mc );
//
// 	AliGenPythiaEventHeader* pythiaGenHeader = 0;
// 	AliGenDPMjetEventHeader* dpmHeader = 0;
//
// 	pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*> ( header );
// 	dpmHeader = dynamic_cast<AliGenDPMjetEventHeader*> ( header );
//
// 	if(pythiaGenHeader) return pythiaGenHeader->ch
//
// }
//_____________________________________________________________________________________________________________________________________________
AliGenEventHeader* AliMultiplicityHelper::GetGenEventHeader ( const AliVEvent* mc ) {
	const AliMCEvent* event =  dynamic_cast<const AliMCEvent*> ( mc );

	if ( !event )
		return 0;

	return event->GenEventHeader();
}

Int_t AliMultiplicityHelper::FindPrimaryMother ( AliStack* stack, Int_t label ) {
	Int_t nPrim  = stack->GetNprimary();

	while ( label >= nPrim ) {
		TParticle* particle = stack->Particle ( label );

		if ( !particle ) {
			AliWarningClassF ( "UNEXPECTED: particle with label %d not found in stack.", label );
			return -1;
		}

		// find mother
		if ( particle->GetMother ( 0 ) < 0 ) {
			AliWarningClassF ( "UNEXPECTED: Could not find mother of secondary particle %d.", label );
			return -1;
		}

		label = particle->GetMother ( 0 );
	}

	return label;
}

TParticle* AliMultiplicityHelper::GetMother ( AliStack* stack, Int_t label ) {
	Int_t mlabel = AliMultiplicityHelper::FindPrimaryMother ( stack, label );

	if ( mlabel < 0 || !stack->IsPhysicalPrimary ( mlabel ) )
		return 0;

	return stack->Particle ( mlabel );
}

Double_t AliMultiplicityHelper::GetDiffractiveMass ( AliStack* stack, Double_t energyCMS ) {
	Double_t M = -1.;
	Int_t np = stack->GetNprimary();

	Int_t iPart1 = -1;
	Int_t iPart2 = -1;

	Double_t y1 = 1e10;
	Double_t y2 = -1e10;

	const Int_t kNstable = 20;
	const Int_t pdgStable[20] = {
		22,             // Photon
		11,             // Electron
		12,             // Electron Neutrino
		13,             // Muon
		14,             // Muon Neutrino
		15,             // Tau
		16,             // Tau Neutrino
		211,            // Pion
		321,            // Kaon
		311,            // K0
		130,            // K0s
		310,            // K0l
		2212,           // Proton
		2112,           // Neutron
		3122,           // Lambda_0
		3112,           // Sigma Minus
		3222,           // Sigma Plus
		3312,           // Xsi Minus
		3322,           // Xsi0
		3334            // Omega
	};

	for ( Int_t i = 0; i < np; i++ ) {
		TParticle* part = stack->Particle ( i );

		Int_t statusCode = part->GetStatusCode();

		// Initial state particle
		if ( statusCode != 1 )
			continue;

		Int_t pdg = TMath::Abs ( part->GetPdgCode() );
		Bool_t isStable = kFALSE;

		for ( Int_t i1 = 0; i1 < kNstable; i1++ ) {
			if ( pdg == pdgStable[i1] ) {
				isStable = kTRUE;
				break;
			}
		}

		if ( !isStable )
			continue;

		Double_t y = part->Y();

		if ( y < y1 ) {
			y1 = y;
			iPart1 = i;
		}

		if ( y > y2 ) {
			y2 = y;
			iPart2 = i;
		}
	}

	if ( iPart1 < 0 || iPart2 < 0 )
		return M;

	y1 = TMath::Abs ( y1 );
	y2 = TMath::Abs ( y2 );

	TParticle*   part1 = stack->Particle ( iPart1 );
	TParticle*   part2 = stack->Particle ( iPart2 );
	Int_t pdg1 = part1->GetPdgCode();
	Int_t pdg2 = part2->GetPdgCode();


	Int_t iPart = -1;

	if ( pdg1 == 2212 && pdg2 == 2212 ) {
		if ( y1 > y2 )
			iPart = iPart1;
		else if ( y1 < y2 )
			iPart = iPart2;
		else {
			iPart = iPart1;

			if ( gRandom->Uniform ( 0., 1. ) > 0.5 )
				iPart = iPart2;
		}
	} else if ( pdg1 == 2212 )
		iPart = iPart1;
	else if ( pdg2 == 2212 )
		iPart = iPart2;

	if ( iPart > 0 ) {
		TParticle* part = stack->Particle ( iPart );
		Double_t E = part->Energy();
		Double_t P = part->P();
		M = TMath::Sqrt ( TMath::Abs ( ( energyCMS - E - P ) * ( energyCMS - E + P ) ) );
	}

	return M;
}

Bool_t AliMultiplicityHelper::CheckDiff ( AliStack* stack, Int_t& gen_type, Int_t& weighted_type, Double_t energyCMS, Double_t MMax ) {
	Double_t M = GetDiffractiveMass ( stack, energyCMS );

	AliInfoClassF ( "diff. mass = %.2f, cut = %.2f", M, MMax );

	if ( M > MMax )
		M = -1;

	Int_t proc0 = 2;

	if ( gen_type == 2 )
		proc0 = 1;

	if ( gen_type == 1 )
		proc0 = 0;

	Int_t proc = 2;

	if ( M > 0 )
		proc = 0;
	else if ( proc0 == 1 )
		proc = 1;

	if ( proc != 0 ) {
		if ( proc0 != 0 ) {
			weighted_type = gen_type;
		} else     {
			weighted_type = 0;
		}

		return kTRUE;
	}

	return kTRUE;
}

Bool_t AliMultiplicityHelper::CheckDiff ( AliStack* stack, Int_t& gen_type, Int_t& weighted_type, Double_t energyCMS, TH2D* hwSD, TH2D* hwNDDD ) {
	Double_t M = GetDiffractiveMass ( stack, energyCMS );

	Double_t Mmin, Mmax, wSD, wDD, wND;

	Int_t energy = 0;

	if ( TMath::Abs ( energyCMS - 900 ) < 1 ) energy = 1;
	else if ( TMath::Abs ( energyCMS - 2760 ) < 1 ) energy = 2;
	else if ( TMath::Abs ( energyCMS - 7000 ) < 1 ) energy = 3;
	else if ( TMath::Abs ( energyCMS - 8000 ) < 1 ) energy = 3;

	if ( !GetWeights ( M, Mmin, Mmax, wSD, wDD, wND, energy, hwSD, hwNDDD ) )
		return kFALSE;

	if ( M > -1 && M < Mmin )
		return kFALSE;

	if ( M > Mmax )
		M = -1;

	Int_t proc0 = 2;

	if ( gen_type == 2 )
		proc0 = 1;

	if ( gen_type == 1 )
		proc0 = 0;

	Int_t proc = 2;

	if ( M > 0 )
		proc = 0;
	else if ( proc0 == 1 )
		proc = 1;

	if ( proc == 1 && ( gRandom->Uniform ( 0., 1. ) > wDD ) )
		return kFALSE;

	if ( proc == 2 && ( gRandom->Uniform ( 0., 1. ) > wND ) )
		return kFALSE;


	if ( proc != 0 ) {
		if ( proc0 != 0 ) {
			weighted_type = gen_type;
		} else     {
			weighted_type = 0;
		}

		return kTRUE;
	}

	if ( wSD < 0 )
		AliErrorClass ( "wSD<0 ! \n" );

	if ( gRandom->Uniform ( 0., 1. ) > wSD )
		return kFALSE;

// 	if ( iPart == iPart1 ) {
// 		weighted_type = 1;
// 	} else if ( iPart == iPart2 ) {
// 		weighted_type = 1;
// 	}

	return kTRUE;
}

Bool_t AliMultiplicityHelper::GetWeights ( Double_t mass, Double_t& minMass, Double_t& maxMass, Double_t& wSD, Double_t& wDD, Double_t& wND, Int_t energyCMS, TH2D* hwSD, TH2D* hwNDDD ) {
	wDD = hwNDDD->GetBinContent ( energyCMS, 2 );
	wND = hwNDDD->GetBinContent ( energyCMS, 1 );
	wSD = -1;

	minMass = hwSD->GetXaxis()->GetBinLowEdge ( 1 );
	maxMass = hwSD->GetXaxis()->GetBinUpEdge ( hwSD->GetXaxis()->GetNbins() );

	if ( mass < minMass || mass > maxMass )
		return kTRUE;

	wSD = hwSD->GetBinContent ( hwSD->FindFixBin ( mass, energyCMS ) );

	return kTRUE;
}

Bool_t AliMultiplicityHelper::IsCocktail ( AliGenEventHeader* h ) {
	if ( dynamic_cast<AliGenCocktailEventHeader*> ( h ) ) return kTRUE;

	return kFALSE;
}

UInt_t AliMultiplicityHelper::Classify ( AliESDEvent* esd, const AliMCEvent* mc, Double_t consistencyCut, Double_t SPDPileUpThreshold ) {
	UInt_t mask = 0;

	AliAnalysisUtils* u = new AliAnalysisUtils();

	if ( !u->IsSPDClusterVsTrackletBG ( esd ) ) mask |= kNotClusterCutBG;

	delete u;

	//trigger
	AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*> ( AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler() );

	if ( !inputHandler ) {
		AliErrorClass ( "Could not get input handler." );
		return 0;
	}

	Bool_t selectedMBOR = inputHandler->IsEventSelected() & AliVEvent::kMB;
	Bool_t selectedV0AND = inputHandler->IsEventSelected() & AliVEvent::kINT7;
	Bool_t selectedHiMult = inputHandler->IsEventSelected() & AliVEvent::kHighMult;

	if ( selectedMBOR ) mask |= kMBOR;

	if ( selectedV0AND ) mask |= kV0AND;

	if ( selectedHiMult ) mask |= kHighMult;

	//SPD pileup
	if ( esd->IsPileupFromSPD ( 3, SPDPileUpThreshold ) ) {
		mask |= kSPDPileUp;
	} else {
		mask |= kSPDNoPileup;
	}

	//vertex
	if ( esd ) {
		Bool_t disp = kFALSE;
		const AliESDVertex* vSPD = ( AliESDVertex* ) ( AliMultiplicityHelper::GetVertex ( kVertexSPD, esd ) );
		const AliESDVertex* vTracks = ( AliESDVertex* ) ( AliMultiplicityHelper::GetVertex ( kVertexTracks, esd ) );
		Bool_t qTracks = AliMultiplicityHelper::IsQualityVertex ( vTracks, disp );
		Bool_t qSPD = AliMultiplicityHelper::IsQualityVertex ( vSPD, disp );


		Bool_t qC = kFALSE;

		if ( qSPD && qTracks ) {
			qC = ( TMath::Abs ( vSPD->GetZ() - vTracks->GetZ() ) < consistencyCut );
		}

		if ( qSPD ) {
			mask |= kVSPD;

			if ( vSPD->IsFromVertexer3D() ) mask |= kVSPD3D;

			if ( vSPD->IsFromVertexerZ() ) mask |= kVSPD1D;
		} else {
			mask |= kVNoSPD;

			if ( disp ) mask |= kVNoSPDDC;
		}

		if ( qTracks ) {
			mask |= kVGlobal;
		} else {
			mask |= kVNoGlobal;
		}

		if ( qC ) {
			mask |= kVConsistent;
		} else {
			mask |= kVInconsistent;
		}
	}

	//MC process type
	if ( mc ) {
		AliGenEventHeader* header = GetGenEventHeader ( mc );

		if ( !header ) {
			AliErrorClass ( "Could not get MC event header." );
			return 0;
		}

		AliGenPythiaEventHeader* pythiaGenHeader = 0;
		AliGenDPMjetEventHeader* dpmHeader = 0;

		pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*> ( header );
		dpmHeader = dynamic_cast<AliGenDPMjetEventHeader*> ( header );

		Bool_t isPythia8 = kFALSE;

		if ( pythiaGenHeader ) {
			isPythia8 = ( pythiaGenHeader->ProcessType() > 100 );
		}

		Int_t iSD1flag, iSD2flag, iDDflag;
		Int_t iMCProcessType = 0;

		if ( pythiaGenHeader ) {
			iMCProcessType = pythiaGenHeader->ProcessType();

			if ( !isPythia8 ) {
				iSD1flag = 92;
				iSD2flag = 93;
				iDDflag = 94;
			} else {
				iSD1flag = 103;
				iSD2flag = 104;
				iDDflag = 105;
			}
		} else if ( dpmHeader ) {
			iMCProcessType = dpmHeader->ProcessType();
			iSD1flag = 5;
			iSD2flag = 6;
			iDDflag = 7;
		} else {
			AliErrorClass ( "Unsupported generator type." );
			return 0;
		}

		if ( iMCProcessType == iSD1flag || iMCProcessType == iSD2flag ) {
			mask |= kSingleDiffractive;
			return mask;
		} else	if ( iMCProcessType == iDDflag ) {
			mask |= kDoubleDiffractive;
			return mask;
		} else {
// 			if (!isPythia8 ) {
// 				mask |= kNonDiffractive;
// 			}
// 			if ( iMCProcessType != 102 ) {
			mask |= kNonDiffractive;
			return mask;
// 			}
// 			return 0;
		}
	}

	return mask;
}

void AliMultiplicityHelper::Declassify ( UInt_t mask ) {
	TString s;

	for ( Int_t bit = 8; bit < 11; ++bit ) {
		s += Form ( "%d\t", ( Int_t ) ( ( BIT ( bit ) & mask ) != 0 ) );
	}

	AliInfoClassF ( "%s", s.Data() );
}
