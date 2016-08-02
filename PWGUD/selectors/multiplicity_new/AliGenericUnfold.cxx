#include "AliGenericUnfold.h"
#include <TString.h>
#include <TF1.h>
#include <TMath.h>
#include <TVirtualFitter.h>
#include <TFitter.h>
#include <TVectorD.h>
#include <TMatrixT.h>
#include <TVectorT.h>
#include <TMatrix.h>

#include <TLegend.h>
#include <TFile.h>
#include <TMinuitMinimizer.h>
#include <TObjString.h>
#include <TVectorDfwd.h>
#include <fstream>

void mFCN_wrapper ( Int_t&, Double_t*, Double_t& chi2, Double_t* params, Int_t ) {
	chi2 = AliGenericUnfold::GetInstance()->ChiSquared ( params );
}

ClassImp ( AliGenericUnfold )

AliGenericUnfold* AliGenericUnfold::fInstance = 0;

AliGenericUnfold::AliGenericUnfold ( Option_t* options ) :
	fOptions ( options ),
	fuseMinuit2 ( kFALSE ),
	fBayesIterationLimit ( 50 ),
	fBayesWeight ( 0.5 ),
	fChi2Weight ( 1 ),
	fChi2Weights ( 0 ),
	nX ( 1 ),
	nXstop ( 1 ),
	nY ( 1 ),
	bUseMeasuredDamp ( kFALSE ),
	fDampThreshold ( 10 ),
	regStart ( 0 ),
	mDampStart ( 0 ),
	mStrategy ( 2 ),
	mPrint ( 0 ),
	fMeasured ( 0 ),
	fMeasuredCopy ( 0 ),
	fMeasuredErr ( 0 ),
	fMeasuredCov ( 0 ),
	fMeasuredCovInv ( 0 ),
	fResponse ( 0 ),
	fUnfolded ( 0 ),
	fUnfoldedErr ( 0 ),
	fMeasuredIntegral ( 1.0 ),
	C ( 0 ),
	A ( 0 ),
	B ( 0 ),
	G ( 0 ),
	bias ( 0 ),
	Delta ( 0 ),
	deltaChi2 ( 0 ),
	totalChi2 ( 0 ),
	lastbin ( 1 ),
	M ( 0 ) {

}

AliGenericUnfold::~AliGenericUnfold() {

}

AliGenericUnfold* AliGenericUnfold::GetInstance() {
	if ( !fInstance )
		fInstance = new AliGenericUnfold();

	return fInstance;
}

Double_t AliGenericUnfold::ChiSquared ( Double_t* params ) {
	Double_t chi2total = 0;
	Double_t chi2proper = 0;
	Double_t chi2reg = 0;

	TVectorD paramsv ( nY, params );

	for ( Int_t i = 0; i < nY; ++i ) {
		paramsv[i] *= paramsv[i];
	}

	TVectorD temp = ( *fResponse ) * paramsv - ( *fMeasured );

	chi2proper = temp * ( ( *fMeasuredCovInv ) * temp );

	if ( chi2proper < 0 ) Printf ( "WTF!!?!?" );

	if ( fOptions.Contains ( optsd ) ) {
		for ( Int_t j = regStart + 1; j < nY - 1; j++ ) {
			Double_t sd = TMath::Power ( ( paramsv[j - 1] - 2 * paramsv[j]  + paramsv[j + 1] ) , 2 ) ;

			if ( fChi2Weights )
				sd *= ( *fChi2Weights ) [j];

			chi2reg += sd;
		}

		chi2reg *= nY;
	} else if ( fOptions.Contains ( optsdw ) ) {
		for ( Int_t j = regStart + 1; j < nY - 1; j++ ) {
			Double_t sd = TMath::Power ( ( paramsv[j - 1] - 2 * paramsv[j]  + paramsv[j + 1] )  , 2 ) ;
			Double_t w = paramsv[j];

			if ( bUseMeasuredDamp && ( j > mDampStart ) && ( j < nX - 1 ) ) {
				Double_t penalty = 0;
				Double_t sdm = TMath::Power ( ( *fMeasuredCopy ) [j - 1] - 2 * ( *fMeasuredCopy ) [j] + ( *fMeasuredCopy ) [j + 1] , 2 );
				Double_t wm = ( *fMeasuredCopy ) [j];

				if ( wm > 0 )
					penalty = sdm / wm;

				sd *= fDampThreshold / ( fDampThreshold + penalty * penalty ) ;
			}

			if ( w > 0 )
				sd /= w;

			if ( fChi2Weights )
				sd *= ( *fChi2Weights ) [j];

			if ( TMath::Finite ( sd ) )
				chi2reg += sd;
		}
	} else if ( fOptions.Contains ( optsd1w ) ) {
		for ( Int_t j = regStart + 1; j < nY - 1; j++ ) {
// 				std::cerr<<j<<"\t"<<paramsv[j-1]<<"\t"<<paramsv[j]<<"\t"<<paramsv[j+1]<<std::endl;
			Double_t sd = TMath::Power ( ( paramsv[j - 1] - 2 * paramsv[j]  + paramsv[j + 1] ) , 2 );
			Double_t weight = paramsv[j - 1] + 4 * paramsv[j] + paramsv[j + 1];

			if ( bUseMeasuredDamp && ( j > mDampStart ) && ( j < nX - 1 ) ) {
				Double_t penalty = 0;
				Double_t sdm = TMath::Power ( ( *fMeasuredCopy ) [j - 1] - 2 * ( *fMeasuredCopy ) [j] + ( *fMeasuredCopy ) [j + 1] , 2 );
				Double_t wm = ( *fMeasuredCopy ) [j - 1] + 4 * ( *fMeasuredCopy ) [j] + ( *fMeasuredCopy ) [j + 1];

				if ( wm > 0 )
					penalty = sdm / wm;

				sd *= fDampThreshold / ( fDampThreshold + penalty * penalty ) ;
			}

			if ( weight > 0 )
				chi2reg += sd / weight;
		}
	} else if ( fOptions.Contains ( optsd2w ) ) {
		for ( Int_t j = regStart + 1; j < nY - 1; j++ ) {
			Double_t sd = TMath::Power ( ( paramsv[j - 1] - 2 * paramsv[j]  + paramsv[j + 1] )  , 2 ) ;
			Double_t w = paramsv[j] * paramsv[j];

			if ( bUseMeasuredDamp && ( j > mDampStart ) && ( j < nX - 1 ) ) {
				Double_t penalty = 0;
				Double_t sdm = TMath::Power ( ( *fMeasuredCopy ) [j - 1] - 2 * ( *fMeasuredCopy ) [j] + ( *fMeasuredCopy ) [j + 1] , 2 );
				Double_t wm = ( *fMeasuredCopy ) [j];

				if ( wm > 0 )
					penalty = sdm / wm;

				sd *= fDampThreshold / ( fDampThreshold + penalty * penalty ) ;
			}

			if ( w > 0 )
				sd /= w;

			if ( fChi2Weights )
				sd *= ( *fChi2Weights ) [j];

			if ( TMath::Finite ( sd ) && ( w > 0 ) )
				chi2reg += sd;

// 				chi2reg *= nY;
		}
	} else if ( fOptions.Contains ( opttd ) ) {
		for ( Int_t j = regStart + 2; j < nY - 2; ++j ) {
			Double_t td = TMath::Power ( -1 / 2 * paramsv[j - 2] + paramsv[j - 1] + paramsv[j + 1] + 1 / 2 * paramsv[j + 2] , 2 );

// 				if ( j <= nDamp ) td /= fDampWeight;
			if ( paramsv[j] > 0 ) chi2reg += td / paramsv[j];
		}
	} //else

// 	if ( fOptions.Contains(optsdp) ) {
// 		for (Int_t j = 4; j < nY - 2; ++j ) {
// 			Double_t sd = TMath::Power( -1/12*paramsv[j-2] + 4/3*paramsv[j-1] - 5/2*paramsv[j] + 4/3*paramsv[j+1] - 1/12*paramsv[j+2] , 2);
// 			if ( j <= nDamp ) sd /= fDampWeight;
// 			chi2reg += sd;
// 		}
// 	}
	if ( !fChi2Weights )
		chi2reg *= fChi2Weight * nY;

	if ( fOptions.Contains ( optsd ) ) chi2reg /= fMeasuredIntegral;

	deltaChi2 = chi2reg;
	chi2total = chi2proper + chi2reg;
	totalChi2 = chi2total;
	return chi2total;

}

Option_t* AliGenericUnfold::GetOptions() {
	return fOptions.Data();
}

Double_t AliGenericUnfold::GetWeight ( Option_t* type ) {
	TString t = type;

	if ( t.Contains ( optbayes ) )
		return fBayesWeight;

	if ( t.Contains ( optchi2 ) )
		return fChi2Weight;

	LogMessage ( Form ( "WARNING: Wrong reg type: %s", type ) );
	return 0;
}

void AliGenericUnfold::SetRegStart ( Int_t start ) {
	if ( start >= 0 ) {
		regStart = start;
	} else
		LogMessage ( Form ( "WARNING: Starting index should be greater or equal 0: %d", start ) );
}

void AliGenericUnfold::SetmDampStart ( Int_t start ) {
	if ( start > 0 ) {
		mDampStart = start;
	} else
		LogMessage ( Form ( "WARNING: Starting index should be greater than 0: %d", start ) );
}


void AliGenericUnfold::SetOptions ( Option_t* options ) {
	fOptions = options;
}

void AliGenericUnfold::SetWeight ( Double_t weight, Option_t* type ) {
	TString t = type;

	if ( t.Contains ( optbayes ) ) {
		if ( ( 0. < weight ) && ( weight <= 1.0 ) )
			fBayesWeight = weight;
		else
			LogMessage ( Form ( "WARNING: Bayesian weight should be in (0:1] (%.2f)", weight ) ); //Printf ( "Bayesian weight should be in (0:1] (%.2f)", weight );

		return;
	}

	if ( t.Contains ( optchi2 ) ) {
		fChi2Weight = weight;
		return;
	}

	LogMessage ( Form ( "WARNING: Wrong reg type: %s", type ) );
}

void AliGenericUnfold::SetWeight ( TVectorD* weights ) {
	if ( weights && ( weights->GetNoElements() == nY ) ) {
		LogMessage ( "INFO: Using per-bin weights" );
		fChi2Weights = weights;
	} else {
		LogMessage ( "WARNING: per-bin weights vector not defined or has wrong size, not using" );
	}
}


void AliGenericUnfold::SetData ( TVectorD* m, TVectorD* me, TMatrixD* r, Int_t stop ) {
	if ( m ) {
		fMeasured = m;
		fMeasuredCopy = new TVectorD ( *fMeasured );
		nX = fMeasured->GetNoElements();
		fMeasuredIntegral = 0;

		for ( Int_t i = 0; i < nX; ++i ) {
			fMeasuredIntegral += ( *fMeasured ) [i];
		}
	}

	if ( r ) {
		fResponse = r;
		nY = fResponse->GetNcols();
	}

	if ( me ) {
		fMeasuredErr = me;
		fMeasuredCov = new TMatrixD ( nX, nX );
		fMeasuredCov->Zero();
		fMeasuredCovInv = new TMatrixD ( nX, nX );
		fMeasuredCovInv->Zero();
		Delta = new TMatrixD ( nY, nY );
		Delta->UnitMatrix();

		for ( Int_t i = 0; i < nX; i++ ) {
			Double_t e = ( *fMeasuredErr ) [i] > 0 ? ( *fMeasuredErr ) [i] : 1.0;
			( *fMeasuredCov ) [i][i] = e * e ;
			( *fMeasuredCovInv ) [i][i] = 1.0 / ( e * e );
		}

	}

	if ( !fMeasured )
		LogMessage ( "Measured not set" );

	if ( !fMeasuredErr )
		LogMessage ( "Measured error not set" );

	if ( !fResponse )
		LogMessage ( "Response not set" );

	if ( fMeasured && fMeasuredErr && fResponse ) {
		if ( fMeasuredCopy ) delete fMeasuredCopy;

		fMeasuredCopy = new TVectorD ( *fMeasured );

		if ( fUnfolded )
			delete fUnfolded;

		fUnfolded = new TVectorD ( nY );
		fUnfolded->Zero();

		if ( fUnfoldedErr )
			delete fUnfoldedErr;

		fUnfoldedErr = new TVectorD ( nY );
		fUnfoldedErr->Zero();

		if ( C )
			delete C;

		C = 0;

		if ( A )
			delete A;

		A = 0;

		if ( B )
			delete B;

		B = 0;

		if ( stop < 0 )
			nXstop = nX;
		else
			nXstop = stop;

		if ( M ) delete M;
	}
}

unsigned int AliGenericUnfold::GetBayesLimit() {
	return fBayesIterationLimit;
}

void AliGenericUnfold::SetBayesLimit ( unsigned int limit ) {
	fBayesIterationLimit = limit;
}

TVectorD* AliGenericUnfold::Unfold() {
	if ( !fMeasured || !fResponse ) {
		LogMessage ( "No data set" );
		return 0;
	}

	ResetUnfolded();

	if ( !fUnfolded ) {
		fUnfolded = new TVectorD ( nY );
		fUnfolded->Zero();
	}

	if ( !fUnfoldedErr ) {
		fUnfoldedErr = new TVectorD ( nY );
		fUnfoldedErr->Zero();
	}

	if ( fOptions.Contains ( optbayes ) ) {
		if ( !M ) M = new TMatrixD ( nY, nX );

		TMatrixD R1 ( nY, nX );

		for ( Int_t i = 0; i < nX ; i++ ) {
			( *fUnfolded ) [i] = ( *fMeasured ) [i];
		}

		TVectorD vTemp ( nX );
		TVectorD vPrev ( nY );
		Double_t diff = 0;

		Bool_t limit = ( GetBayesLimit() > 0 ) ? kTRUE : kFALSE;
		unsigned int nIterations = 0;

		do {
			diff = 0;
			nIterations++; //count iterations if we are limiting them
			vPrev = ( *fUnfolded );

			//smoothing
			for ( Int_t j = 1; j < nY - 1; j++ ) {
				( *fUnfolded ) [j] = ( 1 - fBayesWeight ) * ( *fUnfolded ) [j] + fBayesWeight * ( ( *fUnfolded ) [j - 1] + ( *fUnfolded ) [j] + ( *fUnfolded ) [j + 1] ) / 3;
			}

// 			Int_t nsmooth = 5;
// 			for ( Int_t j = nsmooth / 2; j < nY - nsmooth / 2; j++ ) {
// 				Double_t tmp = ( 1 - fBayesWeight ) * ( *fUnfolded ) [j];
// 				for (Int_t k = 0; k < nsmooth; ++k ) {
// 					tmp += (*fUnfolded)[k+j-2] / (Double_t)nsmooth;
// 				}
// 				(*fUnfolded)[j] = tmp;
// 			}
			( *fUnfolded ) [0] = ( 1 - fBayesWeight ) * ( ( *fUnfolded ) [0] + ( *fUnfolded ) [1] ) / 2;

			vTemp = ( *fResponse ) * ( *fUnfolded );

			for ( Int_t i = 0; i < nY; i++ ) {
				for ( Int_t j = 0; j < nX; j++ ) {
					R1[i][j] = ( *fResponse ) [j][i] * ( *fUnfolded ) [i] / vTemp[j];

					if ( !TMath::Finite ( R1[i][j] ) )
						R1[i][j] = 0;
				}
			}

			( *fUnfolded ) = R1 * ( *fMeasured );

			for ( Int_t i = 0; i < nXstop; i++ ) {
				diff += TMath::Power ( ( vPrev[i] - ( *fUnfolded ) [i] ) / ( *fMeasuredErr ) [i], 2 );
			}

			if ( limit && ( nIterations > GetBayesLimit() ) )
				break;
		} while ( diff > nY * 1e-6 );

		vTemp = ( *fResponse ) * ( *fUnfolded );

		for ( Int_t i = 0; i < nY; ++i ) {
			for ( Int_t j = 0; j < nX; ++j ) {
				( *M ) [i][j] = ( *fResponse ) [j][i] * ( *fUnfolded ) [i] / vTemp[j];

				if ( !TMath::Finite ( ( *M ) [i][j] ) )
					( *M ) [i][j] = 0;
			}
		}

		TVectorD temp = ( *fResponse ) * ( *fUnfolded ) - ( *fMeasured );
		totalChi2 = temp * ( ( *fMeasuredCovInv ) * temp );

		return fUnfolded;
	}

	if ( fOptions.Contains ( optchi2 ) ) {
		if ( fuseMinuit2 ) TVirtualFitter::SetDefaultFitter ( "Minuit2" );

		TVirtualFitter* minuit = TVirtualFitter::Fitter ( 0, nY );

		minuit->SetPrecision ( 1e-6 );
		minuit->SetFCN ( mFCN_wrapper ); //set ChiSq function

		Double_t arglist[10]; //aux arguments to Minuit
		Int_t status = 0;
		arglist[0] = mPrint;
		minuit->ExecuteCommand ( "SET PRINT", arglist, 1 );
		arglist[0] = mStrategy;
		minuit->ExecuteCommand ( "SET STRATEGY", arglist, 1 );

		TVectorD parameters ( nY );
		parameters.Zero();

		for ( Int_t i = 0; i < nX; ++i ) {
			parameters[i] = ( *fMeasured ) [i];
		}

		Double_t step = 0.1;

		for ( Int_t i = 0; i < nY; i++ ) {
			minuit->SetParameter ( i, Form ( "U%d", i ), parameters[i] == 0 ? 0.316228 : TMath::Sqrt ( parameters[i] ), parameters[i] == 0 ? step : step * TMath::Sqrt ( parameters[i] ), /*limit_down*/ 0, /*limit_up*/ 0 ); //set minuit parameters
		}

		arglist[0] = 1e8; //max iterations
		status = minuit->ExecuteCommand ( "MIGRAD", arglist, 1 ); //do the minimization
		LogMessage ( Form ( "Minuit status: %d", status ), kTRUE );
// 		if ( status < 0 ) return 0;

		for ( Int_t i = 0; i < nY; i++ ) { //get the result
			Double_t v = minuit->GetParameter ( i );
			Double_t ve = minuit->GetParError ( i );
			( *fUnfolded ) [i] = v * v;
			( *fUnfoldedErr ) [i]  = 2 * ve * v;
		}

		delete minuit;

		return fUnfolded;
	}

	LogMessage ( "Unfolding method not set in options" );
	return 0;
}

TVectorD* AliGenericUnfold::GetUnfoldErrorMinuit() {
	if ( fUnfoldedErr )
		return fUnfoldedErr;

	return 0;
}

TVectorD* AliGenericUnfold::Fold() {
	TVectorD* out = new TVectorD ( nX );
	( *out ) = ( *fResponse ) * ( *fUnfolded );
	return out;
}

void AliGenericUnfold::SetDampMeasured ( Double_t threshold ) {
	bUseMeasuredDamp = kTRUE;
	fDampThreshold = threshold;
}


TVectorD* AliGenericUnfold::GetUnfoldingBias() {
	TVectorD* unfoldingBias = new TVectorD ( nY );

	TVectorD temp = ( *fResponse ) * ( *fUnfolded ) - ( *fMeasured );

	if ( !C )
		CalculateCMatrix();

	( *unfoldingBias ) = ( *C ) * temp;

	return unfoldingBias;
}

TMatrixD* AliGenericUnfold::GetUnfoldingVariance() {
	TMatrixD* unfoldingVariance = new TMatrixD ( nY, nY );
	TMatrixD Ct ( nX, nY );

	if ( !C )
		CalculateCMatrix();

	Ct.Transpose ( *C );

	*unfoldingVariance = ( *C ) * ( *fMeasuredCov ) * Ct;

	return unfoldingVariance;
}

TMatrixD* AliGenericUnfold::GetUnfoldingBiasVariance() {
	TMatrixD* uvariance =  GetUnfoldingVariance();

	TMatrixD* biasVariance = new TMatrixD ( nY, nY );
	TMatrixD unit ( nY, nY );
	unit.UnitMatrix();
	TMatrixD CRmI = ( *C ) * ( *fResponse ) - unit;
	TMatrixD CRmIt ( nY, nY );
	CRmIt.Transpose ( CRmI );

	*biasVariance = CRmI * ( *uvariance ) * CRmIt;

	return biasVariance;
}

TMatrixD* AliGenericUnfold::GetFullVariance(Bool_t minus) {
	TMatrixD* fullvariance = new TMatrixD(nY, nY);
	
	TMatrixD* uvariance = GetUnfoldingVariance();
	TMatrixD CR = (*C) * (*fResponse);
	TMatrixD CRt(nY,nY); CRt.Transpose(CR);
	TMatrixD unit ( nY, nY );
	unit.UnitMatrix();
	if ( minus ) {	
		(*fullvariance) = (*uvariance) + CR * (*uvariance) - CR * (*uvariance) * CRt;
	} else {
		(*fullvariance) = ( 3.0 * unit + CR ) * (*uvariance) + (CR + 2.0 * unit) * (*uvariance) * CRt;
	}
	
	return fullvariance;
}


TMatrixD* AliGenericUnfold::CalculateCMatrix() {
	std::cout << "calculating C matrix..." << std::endl;
	TVectorD* result = 0;
	TVectorD* resultErr = 0;

	if ( fUnfolded ) {
		result = new TVectorD ( *fUnfolded );
		resultErr = new TVectorD ( *fUnfoldedErr );
	} else {
		result = Unfold();
		resultErr = GetUnfoldErrorMinuit();
	}

	C = new TMatrixD ( nY, nX );
	C->Zero();

	TMatrixD result_shift ( nY, 4 );
	TVectorD* tmp = 0;

	for ( Int_t i = 0; i < nXstop; ++i ) {
		Double_t backup_i = ( *fMeasured ) [i];

		Int_t nres = 0;

		for ( Float_t shift = -1; shift <= 1; shift += 0.5 ) {
			if ( shift == 0 )
				continue;

			( *fMeasured ) [i] = backup_i + shift * TMath::Sqrt ( backup_i );
			ResetUnfolded();

			tmp = Unfold();

			for ( Int_t k = 0; k < nY; ++k ) {
				result_shift[k][nres] = ( *tmp ) [k];
			}

			nres++;
		}

		for ( Int_t j = 0; j < nY; ++j ) {
			( *C ) [j][i] = 0.5 / TMath::Sqrt ( backup_i ) * ( 8. * ( ( result_shift[j][2] ) - ( result_shift[j][1] ) ) - ( ( result_shift[j][3] ) - ( result_shift[j][0] ) ) ) / 3;
		}

		( *fMeasured ) [i] = backup_i;
	}

	if ( result )
		fUnfolded = result;

	if ( resultErr )
		fUnfoldedErr = resultErr;

	return C;
}

TMatrixD* AliGenericUnfold::CalculateCMatrixDirect() {
	if ( C ) delete C;

	C = new TMatrixD ( nY, nX );

	if ( A ) delete A;

	CalculateA();

	if ( !B ) CalculateB();

	( *C ) = ( A->Invert() ) * ( *B );

	return C;
}

TMatrixD* AliGenericUnfold::CalculateA() {
	TVectorD muplus ( nY - 2 );
	TVectorD muminus ( nY - 2 );

	for ( Int_t i = 1; i < nY - 1; ++i ) {
		Double_t t = ( *fUnfolded ) [i - 1] + ( *fUnfolded ) [i + 1];
		muplus[i - 1] = t + 4 * ( *fUnfolded ) [i];
		muminus[i - 1] = t - 2 * ( *fUnfolded ) [i];
	}

	A = new TMatrixD ( nY, nY );
	A->Zero();

	for ( Int_t i = regStart; i < nY; ++i ) {
		for ( Int_t j = regStart; j < nY; ++j ) {
			Double_t t = 0;

			for ( Int_t m = 0; m < nY - 2; ++m ) {
				Double_t tmp = 0;

				if ( fOptions.Contains ( optsd ) ) {
					tmp = ( ( *Delta ) [j][m] - 2 * ( *Delta ) [j][m + 1] + ( *Delta ) [j][m + 2] ) * ( ( *Delta ) [i][m] - 2 * ( *Delta ) [i][m + 1] + ( *Delta ) [i][m + 2] );
					tmp /= nY;
				}

				if ( fOptions.Contains ( optsdw ) ) {
					tmp = ( ( ( *Delta ) [i][m] - 2 * ( *Delta ) [i][m + 1] + ( *Delta ) [i][m + 2] ) * ( ( *Delta ) [j][m] - 2 * ( *Delta ) [j][m + 1] + ( *Delta ) [j][m + 2] )
					        - muminus[m] * ( *Delta ) [j][m + 1] * ( ( *Delta ) [i][m] - 2 * ( *Delta ) [i][m + 1] + ( *Delta ) [i][m + 2] ) / ( ( *fUnfolded ) [m + 1] )
					        + ( *Delta ) [i][m + 1] * ( *Delta ) [j][m + 1] * ( muminus[m] * muminus[m] ) / ( ( *fUnfolded ) [m + 1] * ( *fUnfolded ) [m + 1] )
					      ) / ( *fUnfolded ) [m + 1];
				}

				if ( fOptions.Contains ( optsd1w ) ) {
					tmp = ( ( *Delta ) [j][m] - 2 * ( *Delta ) [j][m + 1] + ( *Delta ) [j][m + 2] ) * ( ( *Delta ) [i][m] - 2 * ( *Delta ) [i][m + 1] + ( *Delta ) [i][m + 2] ) / muplus[m]
					      + TMath::Power ( muminus[m] / muplus[m] , 2 ) / muplus[m] * ( ( *Delta ) [j][m] + 4 * ( *Delta ) [j][m + 1] + ( *Delta ) [j][m + 2] ) * ( ( *Delta ) [i][m] + 4 * ( *Delta ) [i][m + 1] + ( *Delta ) [i][m + 2] )
					      - muminus[m] / ( muplus[m] * muplus[m] ) * ( ( *Delta ) [j][m] + 4 * ( *Delta ) [j][m + 1] + ( *Delta ) [j][m + 2] ) * ( ( *Delta ) [i][m] - 2 * ( *Delta ) [i][m + 1] + ( *Delta ) [i][m + 2] )
					      - muminus[m] / ( muplus[m] * muplus[m] ) * ( ( *Delta ) [j][m] - 2 * ( *Delta ) [j][m + 1] + ( *Delta ) [j][m + 2] ) * ( ( *Delta ) [i][m] + 4 * ( *Delta ) [i][m + 1] + ( *Delta ) [i][m + 2] );
				}

				if ( fOptions.Contains ( optsd2w ) ) {
					tmp = ( 0.5 * ( ( *Delta ) [i][m] - 2 * ( *Delta ) [i][m + 1] + ( *Delta ) [i][m + 2] ) * ( ( *Delta ) [j][m] - 2 * ( *Delta ) [j][m + 1] + ( *Delta ) [j][m + 2] )
					        - muminus[m] / ( *fUnfolded ) [m + 1]
					        * ( ( ( *Delta ) [j][m] - 2 * ( *Delta ) [j][m + 1] + ( *Delta ) [j][m + 2] ) * ( *Delta ) [i][m + 1] + 3 * ( ( *Delta ) [i][m] - 2 * ( *Delta ) [i][m + 1] + ( *Delta ) [i][m + 2] ) * ( *Delta ) [j][m + 1] )
					        + 3 * ( muminus[m] * muminus[m] ) / ( ( *fUnfolded ) [m + 1] * ( *fUnfolded ) [m + 1] ) * ( *Delta ) [i][m + 1] * ( *Delta ) [j][m + 1]
					      ) / ( ( *fUnfolded ) [m + 1] * ( *fUnfolded ) [m + 1] );
// 					tmp *= nY;
				}

				if ( fChi2Weights )
					tmp *= ( *fChi2Weights ) [m + 1];

				if ( bUseMeasuredDamp && ( m > mDampStart ) ) {
					Double_t penalty = 0;
					Double_t sdm = TMath::Power ( ( *fMeasuredCopy ) [m] - 2 * ( *fMeasuredCopy ) [m + 1] + ( *fMeasuredCopy ) [m + 2] , 2 );
					Double_t wm = ( *fMeasuredCopy ) [m] + 4 * ( *fMeasuredCopy ) [m + 1] + ( *fMeasuredCopy ) [m + 2];

					if ( wm > 0 ) penalty = sdm / wm;

					tmp *= fDampThreshold / ( fDampThreshold + penalty * penalty ) ;
				}

				t += tmp;
			}

			t *= 2.0;

			if ( !fChi2Weights ) t *= fChi2Weight * nY;

			if ( fOptions.Contains ( optsd ) ) t /= fMeasuredIntegral;

			( *A ) [i][j] = t;
		}
	}

	TMatrixD Rt ( nY, nX );
	Rt.Transpose ( *fResponse );

	TMatrixD A2 = 2.0 * Rt * ( *fMeasuredCovInv ) * ( *fResponse );

	*A += A2;

	if ( ( *A ) [0][0] < 1e-16 ) ( *A ) [0][0] = 1e-15; //makes matrix non-singular in the case when regularization doesn't start from 0

	return A;
}

TMatrixD* AliGenericUnfold::CalculateB() {
	B = new TMatrixD ( nY, nX );
	B->Zero();

	for ( Int_t i = 0; i < nY; ++i ) {
		for ( Int_t j = 0 ; j < nX; ++j ) {
			( *B ) [i][j] = - 2.0 * ( *fResponse ) [j][i] * ( *fMeasuredCovInv ) [j][j];
		}
	}

	return B;
}


void AliGenericUnfold::ScanWeight ( Double_t wlow, Double_t whi, Int_t nsteps, char* filename, Bool_t draw, Bool_t blog ) {
	if ( nsteps <= 0 ) {
		LogMessage ( "Number of steps must be > 0" );
		return;
	}

	if ( wlow > whi ) {
		LogMessage ( "Upper weight limit must be higher than lower." );
		return;
	}

	Double_t wstep;

	if ( blog ) {
		wstep = ( TMath::Log10 ( whi / wlow ) ) / ( Double_t ) nsteps;
		std::cout << "a = " << wstep << std::endl;
	} else {
		wstep = ( whi - wlow ) / ( Double_t ) nsteps;
		std::cout << "d = " << wstep << std::endl;
	}

	TVectorD* result = 0;
	TVectorD* biasv = 0;
	TMatrixD* variance = 0;
	TMatrixD* bvariance = 0;

	Double_t tbias [nsteps];
	Double_t tbiasw [nsteps];
	Double_t tvar [nsteps];
	Double_t tvarw [nsteps];
	Double_t mse [nsteps];
	Double_t wmse [nsteps];
	Double_t dchi2[nsteps];
	Double_t tchi2[nsteps];
	Double_t chi2b[nsteps];

	for ( Int_t i = 0; i <= nsteps; ++i ) {

		tbias [i] = 0;
		tbiasw [i] = 0;
		tvar [i] = 0;
		tvarw [i] = 0;
		mse [i] = 0;
		wmse [i] = 0;
		chi2b[i] = 0;

		Double_t weight;

		if ( blog ) {
			weight = wlow * TMath::Power ( 10, wstep * i );
		} else {
			weight = wlow + i * wstep;
		}

		SetWeight ( weight, "Chi2" );

		result = Unfold();

		if ( ! result ) continue;

		dchi2[i] = deltaChi2;
		tchi2[i] = totalChi2;

		CalculateCMatrixDirect();
		biasv = GetUnfoldingBias();
		variance = GetUnfoldingVariance();
		bvariance = GetUnfoldingBiasVariance();

		for ( Int_t j = 0; j < nY; ++j ) {
			tbias [i] += ( *biasv ) [j] * ( *biasv ) [j];
			tvar [i] += ( *variance ) [j][j];
			mse [i] += ( *variance ) [j][j] + ( *biasv ) [j] * ( *biasv ) [j];

			tbiasw [i] += ( *biasv ) [j] * ( *biasv ) [j] / ( *result ) [j];
			tvarw [i] += ( *variance ) [j][j] / ( *result ) [j];
			wmse [i] += ( ( *variance ) [j][j] + ( *biasv ) [j] * ( *biasv ) [j] ) / ( *result ) [j];

			chi2b[i] += ( *biasv ) [j] * ( *biasv ) [j] / ( *bvariance ) [j][j];
		}

		tbias[i] /= nY;
		tbiasw[i] /= nY;
		tvar[i] /= nY;
		tvarw[i] /= nY;
		mse[i] /= nY;
		wmse[i] /= nY;

		delete biasv;
		delete variance;
		delete bvariance;
	}

// 	TFile* outroot = TFile::Open ( filename, "RECREATE" );
// 	Double_t xlow = blog ? 0.9 * wlow : wlow - 0.5 * wstep;
// 	Double_t xhi = blog ? 1.1 * whi : whi + 0.5 * wstep;
// 	TH1D hbias ( "hbias", "Total bias", nsteps, xlow, xhi );
// 	TH1D hwbias ( "hwbias", "Total weighted bias", nsteps, xlow, xhi );
// 	TH1D hvar ( "hvar", "Total variance", nsteps, xlow, xhi );
// 	TH1D hwvar ( "hwvar", "Total weighted variance", nsteps, xlow, xhi );
// 	TH1D hmse ( "hmse", "MSE", nsteps, xlow, xhi );
// 	TH1D hwmse ( "hwmse", "Weighted MSE", nsteps, xlow, xhi );
// 	TH1D hdchi2 ( "hdchi2", "#Delta #chi^{2}", nsteps, xlow, xhi );
// 	TH1D htchi2 ( "htchi2", "Total #chi^{2}", nsteps,  xlow, xhi );
// 	TH1D hchi2b ( "hchi2b", "#chi^{2}_{b}", nsteps, xlow, xhi );

// 	TH2D* hlcurve = new TH2D ( "hlcurve", "L-curve", 100, 1., 1.0e3, 100, 1., 1.0e3 );
// 	TGraph* glcurve  = new TGraph( nsteps  );
// 	glcurve->SetMarkerStyle(kFullDiamond);
// 	glcurve->SetMarkerColor(kRed);
//
// 	TGraph* gweight = new TGraph ( nsteps );
// 	gweight->SetMarkerStyle(kDiamond);
// 	gweight->SetMarkerColor(kGreen);
//
// 	for ( Int_t i = 0; i < nsteps; ++i ) {
// 		hbias.SetBinContent ( i + 1, tbias[i] );
// 		hwbias.SetBinContent ( i + 1, tbiasw[i] );
// 		hvar.SetBinContent ( i + 1, tvar[i] );
// 		hwvar.SetBinContent ( i + 1, tvarw[i] );
// 		hmse.SetBinContent ( i + 1, mse[i] );
// 		hwmse.SetBinContent ( i + 1, wmse[i] );
// 		hdchi2.SetBinContent ( i + 1, dchi2[i] );
// 		htchi2.SetBinContent ( i+1, tchi2[i] );
// 		hchi2b.SetBinContent ( i + 1, chi2b[i] );
// // 		hlcurve->Fill(tchi2[i],dchi2[i]/(wlow + i * wstep));
// 		Double_t x = tchi2[i] - dchi2[i];
// 		Double_t y = dchi2[i] / ( blog ? wlow * TMath::Power(10, wstep * i) : wlow + i * wstep );
// // 		std::cout << i << "\t" << tchi2[i] - dchi2[i] << "\t" << x << "\t" << y << std::endl;
// 		glcurve->SetPoint(i, x, y);
// 		gweight->SetPoint(i, x, blog ? wlow * TMath::Power(10, wstep * i) : wlow + i * wstep);
// 	}

// 	hbias.SetLineColor ( kRed );
// 	hbias.SetLineStyle ( kDashed );
// 	hbias.SetLineWidth ( 2 );
//
// 	hwbias.SetLineColor ( kRed );
// 	hwbias.SetLineStyle ( kDashed );
// 	hwbias.SetLineWidth ( 2 );
//
// 	hvar.SetLineColor ( kGreen );
// 	hvar.SetLineStyle ( kDashDotted );
// 	hvar.SetLineWidth ( 2 );
//
// 	hwvar.SetLineColor ( kGreen );
// 	hwvar.SetLineStyle ( kDashDotted );
// 	hwvar.SetLineWidth ( 2 );
//
// 	hmse.SetLineColor ( kBlue );
// 	hmse.SetLineStyle ( kSolid );
// 	hmse.SetLineWidth ( 2 );
//
// 	hwmse.SetLineColor ( kBlue );
// 	hwmse.SetLineStyle ( kSolid );
// 	hwmse.SetLineWidth ( 2 );
//
// 	hdchi2.SetLineColor ( kCyan );
// 	hdchi2.SetLineStyle ( kSolid );
// 	hdchi2.SetLineWidth ( 2 );
//
// 	htchi2.SetLineColor ( kAzure );
// 	htchi2.SetLineStyle ( kDashed );
// 	htchi2.SetLineWidth ( 2 );
//
// 	hchi2b.SetLineColor ( kMagenta );
// 	hchi2b.SetLineStyle ( kDashDotted );
// 	hchi2b.SetLineWidth ( 2 );
//
// 	hbias.Write();
// 	hwbias.Write();
// 	hvar.Write();
// 	hwvar.Write();
// 	hmse.Write();
// 	hwmse.Write();
// 	hdchi2.Write();
// 	htchi2.Write();
// 	hchi2b.Write();
// 	glcurve->Write();

// 	TObjString* regs = new TObjString(fOptions.Data());
// 	regs->Write();

// 	outroot->Close();

// 	if ( draw ) {
// 		TCanvas* display = new TCanvas ( "display", "display", 1610, 800 );
// 		display->Divide ( 2, 1 );
//
// 		display->cd ( 1 )->SetGrid ( 1, 1 );
// 		gPad->SetLogy();
// 		if ( blog ) gPad->SetLogx();
//
// 		hmse.DrawCopy ( "hist" );
// 		hvar.DrawCopy ( "hist same" );
// 		hbias.DrawCopy ( "hist same" );
// // 		hdchi2.DrawCopy ( "his same" );
// // 		hchi2b.DrawCopy ( "hist same" );
//
// 		display->cd ( 1 )->BuildLegend();
//
// 		display->cd ( 2 )->SetGrid ( 1, 1 );
// 		gPad->SetLogy();
// 		if ( blog ) gPad->SetLogx();
//
// 		hwmse.DrawCopy ( "hist" );
// 		hwvar.DrawCopy ( "hist same" );
// 		hwbias.DrawCopy ( "hist same" );
//
// 		display->cd ( 2 )->BuildLegend();
//
// 		TCanvas* c2 = new TCanvas();
// 		c2->cd();
// 		c2->SetLogx();
// 		c2->SetLogy();
// // 		glcurve->Print();
// 		glcurve->Draw("AP");
// 		gweight->Draw("P same");
// 	}


	std::ofstream outtxt ( TString::Format ( "%s.txt", filename ) );

	outtxt << "step\tdelta_chi2\ttotal bias\ttotal bias (w)\ttotal variance\ttotal variance (w)\tMSE\tMSE (w)\tbias chi2\ttotal chi2\ttotal chi2 - delta chi2"<<std::endl;
	for ( Int_t i = 0; i < nsteps; ++i ) {
		outtxt << ( blog ? ( wlow * TMath::Power ( 10, wstep * i ) ) : ( wlow + i * wstep ) ) << "\t" << dchi2[i] << "\t" << tbias[i] << "\t" << tbiasw[i] << "\t" << tvar[i] << "\t" << tvarw[i] << "\t" << mse[i] << "\t" << wmse[i] << "\t" << chi2b[i] << "\t" << tchi2[i] << "\t" << tchi2[i] - dchi2[i] << std::endl;
	}

	outtxt.close();
}

void AliGenericUnfold::ResetUnfolded () {
	if ( fUnfolded ) {
		delete fUnfolded;
		fUnfolded = 0;
	}
}

TH2D* AliGenericUnfold::GetC() {
	if ( !C )
		CalculateCMatrix();

	TH2D* hC = new TH2D ( "C", "", nX, -0.5, nX - 0.5, nY, -0.5, nY - 0.5 );

	for ( Int_t i = 0; i < nX; ++i ) {
		for ( Int_t j = 0; j < nY; ++j ) {
			hC->SetBinContent ( i + 1, j + 1, ( *C ) [j][i] );
		}
	}

	return hC;
}

TH2D* AliGenericUnfold::GetCR() {
	if ( !C )
		CalculateCMatrix();

	TMatrixD CR = ( *C ) * ( *fResponse );

	TH2D* hCR = new TH2D ( "CR", "", nX, -0.5, nX - 0.5, nX, -0.5, nX - 0.5 );

	for ( Int_t i = 0; i < nX; ++i ) {
		for ( Int_t j = 0; j < nX; ++j ) {
			hCR->SetBinContent ( i + 1, j + 1, CR[i][j] );
		}
	}

	return hCR;
}


void AliGenericUnfold::SetStrategy ( Int_t s ) {
	mStrategy = s;
}

void AliGenericUnfold::SetPrint ( Int_t p ) {
	mPrint = p;
}

void AliGenericUnfold::LogMessage ( const char* message, Bool_t bcerr ) {
	if ( bcerr ) {
		std::cerr << GetName() << ": " << message << std::endl;
	} else {
		std::cout << GetName() << ": " << message << std::endl;
	}
}

TVectorD* AliGenericUnfold::GetResiduals() {
	TVectorD* fld = Fold();
	TVectorD* residuals = new TVectorD ( nX );

	for ( Int_t i = 0; i < nX; ++i ) {
		( *residuals ) [i] = ( ( *fld ) [i] - ( *fMeasured ) [i] ) *  TMath::Sqrt ( ( *fMeasuredCovInv ) [i][i] );
	}

	delete fld;
	return residuals;
}

Double_t AliGenericUnfold::GetDeltaChi2() {
	return deltaChi2;
}

void AliGenericUnfold::SetUseMinuit2() {
	fuseMinuit2 = kTRUE;
}

Int_t AliGenericUnfold::GetLastBin() {
	return lastbin;
}

TMatrixD* AliGenericUnfold::GetBayesianVariance() {
	TMatrixD* V = new TMatrixD ( nY, nY );
	TMatrixD Mt ( nX, nY );
	Mt.Transpose ( *M );
	*V = ( *M ) * ( *fMeasuredCov ) * Mt;

	return V;
}

void AliGenericUnfold::ScanWeightBayes ( Int_t nsteps, char* filename, Bool_t draw ) {
	if ( nsteps <= 0 ) {
		LogMessage ( "Negative or 0 number of steps" );
		return;
	}

	Double_t wstep = 1.0 / ( Double_t ) nsteps;

	TVectorD* result = 0;
	TMatrixD* variance = 0;
	Double_t tvar [nsteps];
	Double_t tvarw [nsteps];
	Double_t tchi2[nsteps];

	for ( Int_t i = 0; i < nsteps; ++i ) {
		tvar [i] = 0;
		tvarw [i] = 0;

		Double_t weight = i * wstep;
		SetWeight ( weight, "Bayes" );

		result = Unfold();

		if ( ! result ) continue;

		tchi2[i] = totalChi2;

		variance = GetBayesianVariance();

		for ( Int_t j = 0; j < nY; ++j ) {
			tvar [i] += ( *variance ) [j][j];

			if ( ( *result ) [j] > 1e-10 ) tvarw [i] += ( *variance ) [j][j] / ( *result ) [j];
		}

		tvar[i] /= nY;
		tvarw[i] /= nY;
		delete variance;
	}

	/*TFile* outroot = */TFile::Open ( filename, "RECREATE" );
	TH1D hvar ( "hvar", "Total variance", nsteps,  - wstep / 2.0, 1.0 + wstep / 2.0 );
	TH1D hwvar ( "hwvar", "Total weighted variance", nsteps, 0 - wstep / 2.0, 1.0 + wstep / 2.0 );
	TH1D htchi2 ( "htchi2", "Total #chi^{2}", nsteps,   - wstep / 2.0, 1.0 + wstep / 2.0 );

	for ( Int_t i = 0; i < nsteps; ++i ) {
		hvar.SetBinContent ( i + 1, tvar[i] );
		hwvar.SetBinContent ( i + 1, tvarw[i] );
		htchi2.SetBinContent ( i + 1, tchi2[i] );
	}

	hvar.SetLineColor ( kGreen );
	hvar.SetLineStyle ( kDashDotted );
	hvar.SetLineWidth ( 2 );

	hwvar.SetLineColor ( kGreen );
	hwvar.SetLineStyle ( kDashDotted );
	hwvar.SetLineWidth ( 2 );

	htchi2.SetLineColor ( kAzure );
	htchi2.SetLineStyle ( kDashed );
	htchi2.SetLineWidth ( 2 );

	hvar.Write();
	hwvar.Write();
	htchi2.Write();

	TObjString* regs = new TObjString ( fOptions.Data() );
	regs->Write();

	if ( draw ) {
		TCanvas* display = new TCanvas ( "display", "display", 1610, 800 );
		display->Divide ( 2, 1 );

		display->cd ( 1 )->SetGrid ( 1, 1 );
		gPad->SetLogy();

		hvar.DrawCopy ( "hist" );
		htchi2.DrawCopy ( "hist same" );

		display->cd ( 1 )->BuildLegend();

		display->cd ( 2 )->SetGrid ( 1, 1 );
		gPad->SetLogy();

		hwvar.DrawCopy ( "hist" );

		display->cd ( 2 )->BuildLegend();
	}


	std::ofstream outtxt ( TString::Format ( "%s.txt", filename ) );

	for ( Int_t i = 0; i < nsteps; ++i ) {
		outtxt << i* wstep << "\t" << tvar[i] << "\t" << tvarw[i] << "\t" << tchi2[i] << std::endl;
	}

	outtxt.close();
}

Double_t AliGenericUnfold::GetDeltaBeta() {
	if ( !A ) CalculateA();

	CalculateG();
	TMatrixD uvar = GetUnfoldingVariance()->Invert();
	TMatrixD Rt ( nY, nX ); Rt.Transpose ( *fResponse );

	TVectorD* vr = new TVectorD ( nY );

	if ( fOptions.Contains ( optsd2w ) ) {
		*vr = ( -1.0 ) * ( A->Invert() ) * ( ( *G ) * ( *fUnfolded ) );
	}

	TVectorD temp = ( Rt * ( *fMeasuredCovInv ) ) * ( ( *fResponse ) * ( *fUnfolded )  - ( *fMeasured ) );

	Double_t db =

	    ( Double_t ) nY / (
	        ( *vr ) * temp
	    );
// 	delete uvar;
	delete vr;
	return db;
}

TMatrixD* AliGenericUnfold::CalculateG() {
	if ( G ) delete G;

	G = new TMatrixD ( nY, nY );
	G->Zero();

	for ( Int_t i = 0; i < nY; ++i ) {
		for ( Int_t j = 0; j < nY; ++j ) {
// 			std::cout<<i<<"\t"<<j<<std::endl;
			if ( i == j && i > 1 && i < nY - 2 ) ( *G ) [i][j] = 6 / ( ( *fUnfolded ) [i] * ( *fUnfolded ) [i] );

			if ( i == j && ( i == 1 || i == nY - 2 ) ) ( *G ) [i][j] = 5 / ( ( *fUnfolded ) [i] * ( *fUnfolded ) [i] );

			if ( i == j && ( i == 0 || i == nY - 1 ) ) ( *G ) [i][j] = 1 / ( ( *fUnfolded ) [i] * ( *fUnfolded ) [i] );

			if ( ( i == 1 && j == 0 ) || ( i == nY - 1 && j == nY - 2 ) ) ( *G ) [i][j] = - 2 / ( ( *fUnfolded ) [j] * ( *fUnfolded ) [j] );

			if ( ( i == 0 && j == 1 ) || ( i == nY - 2 && j == nY - 1 ) ) ( *G ) [i][j] = - 2 / ( ( *fUnfolded ) [i] * ( *fUnfolded ) [i] );

			if ( ( i >= 2 && i <= nY - 1 ) && ( i - j ) == 2 ) ( *G ) [i][j] = 1 / ( ( *fUnfolded ) [i - 1] * ( *fUnfolded ) [i - 1] );

			if ( ( j >= 2 && j <= nY - 1 ) && ( j - i ) == 2 ) ( *G ) [i][j] = 1 / ( ( *fUnfolded ) [j - 1] * ( *fUnfolded ) [j - 1] );

// 			if ( TMath::IsNaN ( (*G)[i][j] ) ) (*G)[i][j] = 0;
			( *G ) [i][j] *= nY;
		}
	}

	return G;
}

