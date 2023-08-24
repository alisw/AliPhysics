/*
 * Maintainer: Mingrui Zhao
 */
#include "CorrelationCalculator.h"

ClassImp(CorrelationCalculator);
//_____________________________________________________________________
TComplex CorrelationCalculator::Q(int n, int p)
{

	if(n>=0) return Qvector[n][p];
	else return TComplex::Conjugate(Qvector[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap10M(int n, int p)
{

	if(n>=0) return Qvector10M[n][p];
	else return TComplex::Conjugate(Qvector10M[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap10P(int n, int p)
{

	if(n>=0) return Qvector10P[n][p];
	else return TComplex::Conjugate(Qvector10P[-n][p]);
}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap12M(int n, int p)
{

	if(n>=0) return Qvector12M[n][p];
	else return TComplex::Conjugate(Qvector12M[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap12P(int n, int p)
{

	if(n>=0) return Qvector12P[n][p];
	else return TComplex::Conjugate(Qvector12P[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap0M(int n, int p)
{

	if(n>=0) return Qvector0M[n][p];
	else return TComplex::Conjugate(Qvector0M[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap0P(int n, int p)
{

	if(n>=0) return Qvector0P[n][p];
	else return TComplex::Conjugate(Qvector0P[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap2M(int n, int p)
{

	if(n>=0) return Qvector2M[n][p];
	else return TComplex::Conjugate(Qvector2M[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap2P(int n, int p)
{

	if(n>=0) return Qvector2P[n][p];
	else return TComplex::Conjugate(Qvector2P[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap4M(int n, int p)
{

	if(n>=0) return Qvector4M[n][p];
	else return TComplex::Conjugate(Qvector4M[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap4P(int n, int p)
{

	if(n>=0) return Qvector4P[n][p];
	else return TComplex::Conjugate(Qvector4P[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap6M(int n, int p)
{

	if(n>=0) return Qvector6M[n][p];
	else return TComplex::Conjugate(Qvector6M[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap6P(int n, int p)
{

	if(n>=0) return Qvector6P[n][p];
	else return TComplex::Conjugate(Qvector6P[-n][p]);

}

//____________________________________________________________________
TComplex CorrelationCalculator::QGap8M(int n, int p)
{

	if(n>=0) return Qvector8M[n][p];
	else return TComplex::Conjugate(Qvector8M[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap8P(int n, int p)
{

	if(n>=0) return Qvector8P[n][p];
	else return TComplex::Conjugate(Qvector8P[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::pGap10M(int n, int p)
{

	if(n>=0) return pvectorM[n][p];
	else return TComplex::Conjugate(pvectorM[n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::pGap10P(int n, int p)
{

	if(n>=0) return pvectorP[n][p];
	else return TComplex::Conjugate(pvectorP[n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap14M(int n, int p)
{

	if(n>=0) return Qvector14M[n][p];
	else return TComplex::Conjugate(Qvector14M[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap14P(int n, int p)
{

	if(n>=0) return Qvector14P[n][p];
	else return TComplex::Conjugate(Qvector14P[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap16M(int n, int p)
{

	if(n>=0) return Qvector16M[n][p];
	else return TComplex::Conjugate(Qvector16M[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap16P(int n, int p)
{

	if(n>=0) return Qvector16P[n][p];
	else return TComplex::Conjugate(Qvector16P[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap18M(int n, int p)
{

	if(n>=0) return Qvector18M[n][p];
	else return TComplex::Conjugate(Qvector18M[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QGap18P(int n, int p)
{

	if(n>=0) return Qvector18P[n][p];
	else return TComplex::Conjugate(Qvector18P[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QsubLeft(int n, int p)
{

	if(n>=0) return QvectorSubLeft[n][p];
	else return TComplex::Conjugate(QvectorSubLeft[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QsubRight(int n, int p)
{

	if(n>=0) return QvectorSubRight[n][p];
	else return TComplex::Conjugate(QvectorSubRight[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QsubMiddle(int n, int p)
{

	if(n>=0) return QvectorSubMiddle[n][p];
	else return TComplex::Conjugate(QvectorSubMiddle[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QsubGap2Left(int n, int p)
{

	if(n>=0) return QvectorSubGap2Left[n][p];
	else return TComplex::Conjugate(QvectorSubGap2Left[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QsubGap2Right(int n, int p)
{

	if(n>=0) return QvectorSubGap2Right[n][p];
	else return TComplex::Conjugate(QvectorSubGap2Right[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::QsubGap2Middle(int n, int p)
{

	if(n>=0) return QvectorSubGap2Middle[n][p];
	else return TComplex::Conjugate(QvectorSubGap2Middle[-n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::p(int n, int p)
{

	if(n>=0) return pvector[n][p];
	else return TComplex::Conjugate(pvector[n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::q(int n, int p)
{

	if(n>=0) return qvector[n][p];
	else return TComplex::Conjugate(qvector[n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::pGap0M(int n, int p)
{

	if(n>=0) return pvector0M[n][p];
	else return TComplex::Conjugate(pvector0M[n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::pGap0P(int n, int p)
{

	if(n>=0) return pvector0P[n][p];
	else return TComplex::Conjugate(pvector0P[n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::pGap4M(int n, int p)
{

	if(n>=0) return pvector4M[n][p];
	else return TComplex::Conjugate(pvector4M[n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::pGap4P(int n, int p)
{

	if(n>=0) return pvector4P[n][p];
	else return TComplex::Conjugate(pvector4P[n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::pGap8M(int n, int p)
{

	if(n>=0) return pvector8M[n][p];
	else return TComplex::Conjugate(pvector8M[n][p]);

}
//____________________________________________________________________
TComplex CorrelationCalculator::pGap8P(int n, int p)
{

	if(n>=0) return pvector8P[n][p];
	else return TComplex::Conjugate(pvector8P[n][p]);

}
//____________________________________________________________________
void CorrelationCalculator::ResetQ(const int nMaxHarm, const int nMaxPow)
{
	for(int i=0; i<nMaxHarm; i++)
	{
		for(int j=0; j<nMaxPow; j++)
		{
			Qvector[i][j] = TComplex(0.,0.);
		}
	}
}
//____________________________________________________________________
TComplex CorrelationCalculator::Two(int n1, int n2)
{

	TComplex formula = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap0(int n1, int n2)
{

	TComplex formula = QGap0M(n1,1)*QGap0P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap2(int n1, int n2)
{

	TComplex formula = QGap2M(n1,1)*QGap2P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap4(int n1, int n2)
{

	TComplex formula = QGap4M(n1,1)*QGap4P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap6(int n1, int n2)
{

	TComplex formula = QGap6M(n1,1)*QGap6P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap8(int n1, int n2)
{

	TComplex formula = QGap8M(n1,1)*QGap8P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap10(int n1, int n2)
{

	TComplex formula = QGap10M(n1,1)*QGap10P(n2,1);
	return formula;

}//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap12(int n1, int n2)
{

	TComplex formula = QGap12M(n1,1)*QGap12P(n2,1);
	return formula;

}

//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap14(int n1, int n2)
{

	TComplex formula = QGap14M(n1,1)*QGap14P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap16(int n1, int n2)
{

	TComplex formula = QGap16M(n1,1)*QGap16P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap18(int n1, int n2)
{

	TComplex formula = QGap18M(n1,1)*QGap18P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::Two_3SubLM(int n1, int n2)
{

	TComplex formula = QsubLeft(n1,1)*QsubMiddle(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::Two_3SubRM(int n1, int n2)
{

	TComplex formula = QsubMiddle(n1,1)*QsubRight(n2,1);
	return formula;

}
//
//____________________________________________________________________
TComplex CorrelationCalculator::Two_3SubLR(int n1, int n2)
{

	TComplex formula = QsubLeft(n1,1)*QsubRight(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::Two_3SubGap2LM(int n1, int n2)
{

	TComplex formula = QsubGap2Left(n1,1)*QsubGap2Middle(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::Two_3SubGap2RM(int n1, int n2)
{

	TComplex formula = QsubGap2Middle(n1,1)*QsubGap2Right(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap0M(int n1, int n2)
{

	TComplex formula = QGap0M(n1,1)*QGap0M(n2,1) - QGap0M(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap2M(int n1, int n2)
{

	TComplex formula = QGap2M(n1,1)*QGap2M(n2,1) - QGap2M(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap4M(int n1, int n2)
{

	TComplex formula = QGap4M(n1,1)*QGap4M(n2,1) - QGap4M(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap6M(int n1, int n2)
{

	TComplex formula = QGap6M(n1,1)*QGap6M(n2,1) - QGap6M(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap8M(int n1, int n2)
{

	TComplex formula = QGap8M(n1,1)*QGap8M(n2,1) - QGap8M(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap10M(int n1, int n2)
{

	TComplex formula = QGap10M(n1,1)*QGap10M(n2,1) - QGap10M(n1+n2,2);
	return formula;

}//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap12M(int n1, int n2)
{

	TComplex formula = QGap12M(n1,1)*QGap12M(n2,1) - QGap12M(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap0P(int n1, int n2)
{

	TComplex formula = QGap0P(n1,1)*QGap0P(n2,1) - QGap0P(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap2P(int n1, int n2)
{

	TComplex formula = QGap2P(n1,1)*QGap2P(n2,1) - QGap2P(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap4P(int n1, int n2)
{

	TComplex formula = QGap4P(n1,1)*QGap4P(n2,1) - QGap4P(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap6P(int n1, int n2)
{

	TComplex formula = QGap6P(n1,1)*QGap6P(n2,1) - QGap6P(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap8P(int n1, int n2)
{

	TComplex formula = QGap8P(n1,1)*QGap8P(n2,1) - QGap8P(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap10P(int n1, int n2)
{

	TComplex formula = QGap10P(n1,1)*QGap10P(n2,1) - QGap10P(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::TwoGap12P(int n1, int n2)
{

	TComplex formula = QGap12P(n1,1)*QGap12P(n2,1) - QGap12P(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::Three(int n1, int n2, int n3)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)
		- Q(n1,1)*Q(n2+n3,2)+2.*Q(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap0A(int n1, int n2, int n3)
{

	TComplex formula = QGap0M(n1,1)*QGap0P(n2,1)*QGap0P(n3,1)-QGap0M(n1,1)*QGap0P(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap0B(int n1, int n2, int n3)
{

	TComplex formula = QGap0P(n1,1)*QGap0M(n2,1)*QGap0M(n3,1)-QGap0P(n1,1)*QGap0M(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap2A(int n1, int n2, int n3)
{

	TComplex formula = QGap2M(n1,1)*QGap2P(n2,1)*QGap2P(n3,1)-QGap2M(n1,1)*QGap2P(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap2B(int n1, int n2, int n3)
{

	TComplex formula = QGap2P(n1,1)*QGap2M(n2,1)*QGap2M(n3,1)-QGap2P(n1,1)*QGap2M(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap4A(int n1, int n2, int n3)
{

	TComplex formula = QGap4M(n1,1)*QGap4P(n2,1)*QGap4P(n3,1)-QGap4M(n1,1)*QGap4P(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap4B(int n1, int n2, int n3)
{

	TComplex formula = QGap4P(n1,1)*QGap4M(n2,1)*QGap4M(n3,1)-QGap4P(n1,1)*QGap4M(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap6A(int n1, int n2, int n3)
{

	TComplex formula = QGap6M(n1,1)*QGap6P(n2,1)*QGap6P(n3,1)-QGap6M(n1,1)*QGap6P(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap6B(int n1, int n2, int n3)
{

	TComplex formula = QGap6P(n1,1)*QGap6M(n2,1)*QGap6M(n3,1)-QGap6P(n1,1)*QGap6M(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap8A(int n1, int n2, int n3)
{

	TComplex formula = QGap8M(n1,1)*QGap8P(n2,1)*QGap8P(n3,1)-QGap8M(n1,1)*QGap8P(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap8B(int n1, int n2, int n3)
{

	TComplex formula = QGap8P(n1,1)*QGap8M(n2,1)*QGap8M(n3,1)-QGap8P(n1,1)*QGap8M(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap10A(int n1, int n2, int n3)
{

	TComplex formula = QGap10M(n1,1)*QGap10P(n2,1)*QGap10P(n3,1)-QGap10M(n1,1)*QGap10P(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap10B(int n1, int n2, int n3)
{

	TComplex formula = QGap10P(n1,1)*QGap10M(n2,1)*QGap10M(n3,1)-QGap10P(n1,1)*QGap10M(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap12A(int n1, int n2, int n3)
{

	TComplex formula = QGap12M(n1,1)*QGap12P(n2,1)*QGap12P(n3,1)-QGap12M(n1,1)*QGap12P(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap12B(int n1, int n2, int n3)
{

	TComplex formula = QGap12P(n1,1)*QGap12M(n2,1)*QGap12M(n3,1)-QGap12P(n1,1)*QGap12M(n2+n3,2);
	return formula;

}


//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap0_subM(int n1, int n2, int n3)
{

	TComplex formula = QGap0M(n1,1)*QGap0M(n2,1)*QGap0M(n3,1)-QGap0M(n1+n2,2)*QGap0M(n3,1)-QGap0M(n2,1)*QGap0M(n1+n3,2)
		- QGap0M(n1,1)*QGap0M(n2+n3,2)+2.*QGap0M(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap0_subP(int n1, int n2, int n3)
{

	TComplex formula = QGap0P(n1,1)*QGap0P(n2,1)*QGap0P(n3,1)-QGap0P(n1+n2,2)*QGap0P(n3,1)-QGap0P(n2,1)*QGap0P(n1+n3,2)
		- QGap0P(n1,1)*QGap0P(n2+n3,2)+2.*QGap0P(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap2_subM(int n1, int n2, int n3)
{

	TComplex formula = QGap2M(n1,1)*QGap2M(n2,1)*QGap2M(n3,1)-QGap2M(n1+n2,2)*QGap2M(n3,1)-QGap2M(n2,1)*QGap2M(n1+n3,2)
		- QGap2M(n1,1)*QGap2M(n2+n3,2)+2.*QGap2M(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap2_subP(int n1, int n2, int n3)
{

	TComplex formula = QGap2P(n1,1)*QGap2P(n2,1)*QGap2P(n3,1)-QGap2P(n1+n2,2)*QGap2P(n3,1)-QGap2P(n2,1)*QGap2P(n1+n3,2)
		- QGap2P(n1,1)*QGap2P(n2+n3,2)+2.*QGap2P(n1+n2+n3,3);
	return formula;

}

//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap4_subM(int n1, int n2, int n3)
{

	TComplex formula = QGap4M(n1,1)*QGap4M(n2,1)*QGap4M(n3,1)-QGap4M(n1+n2,2)*QGap4M(n3,1)-QGap4M(n2,1)*QGap4M(n1+n3,2)
		- QGap4M(n1,1)*QGap4M(n2+n3,2)+2.*QGap4M(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap4_subP(int n1, int n2, int n3)
{

	TComplex formula = QGap4P(n1,1)*QGap4P(n2,1)*QGap4P(n3,1)-QGap4P(n1+n2,2)*QGap4P(n3,1)-QGap4P(n2,1)*QGap4P(n1+n3,2)
		- QGap4P(n1,1)*QGap4P(n2+n3,2)+2.*QGap4P(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap6_subM(int n1, int n2, int n3)
{

	TComplex formula = QGap6M(n1,1)*QGap6M(n2,1)*QGap6M(n3,1)-QGap6M(n1+n2,2)*QGap6M(n3,1)-QGap6M(n2,1)*QGap6M(n1+n3,2)
		- QGap6M(n1,1)*QGap6M(n2+n3,2)+2.*QGap6M(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap6_subP(int n1, int n2, int n3)
{

	TComplex formula = QGap6P(n1,1)*QGap6P(n2,1)*QGap6P(n3,1)-QGap6P(n1+n2,2)*QGap6P(n3,1)-QGap6P(n2,1)*QGap6P(n1+n3,2)
		- QGap6P(n1,1)*QGap6P(n2+n3,2)+2.*QGap6P(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap8_subM(int n1, int n2, int n3)
{

	TComplex formula = QGap8M(n1,1)*QGap8M(n2,1)*QGap8M(n3,1)-QGap8M(n1+n2,2)*QGap8M(n3,1)-QGap8M(n2,1)*QGap8M(n1+n3,2)
		- QGap8M(n1,1)*QGap8M(n2+n3,2)+2.*QGap8M(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap8_subP(int n1, int n2, int n3)
{

	TComplex formula = QGap8P(n1,1)*QGap8P(n2,1)*QGap8P(n3,1)-QGap8P(n1+n2,2)*QGap8P(n3,1)-QGap8P(n2,1)*QGap8P(n1+n3,2)
		- QGap8P(n1,1)*QGap8P(n2+n3,2)+2.*QGap8P(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap10_subM(int n1, int n2, int n3)
{

	TComplex formula = QGap10M(n1,1)*QGap10M(n2,1)*QGap10M(n3,1)-QGap10M(n1+n2,2)*QGap10M(n3,1)-QGap10M(n2,1)*QGap10M(n1+n3,2)
		- QGap10M(n1,1)*QGap10M(n2+n3,2)+2.*QGap10M(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap10_subP(int n1, int n2, int n3)
{

	TComplex formula = QGap10P(n1,1)*QGap10P(n2,1)*QGap10P(n3,1)-QGap10P(n1+n2,2)*QGap10P(n3,1)-QGap10P(n2,1)*QGap10P(n1+n3,2)
		- QGap10P(n1,1)*QGap10P(n2+n3,2)+2.*QGap10P(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap12_subM(int n1, int n2, int n3)
{

	TComplex formula = QGap12M(n1,1)*QGap12M(n2,1)*QGap12M(n3,1)-QGap12M(n1+n2,2)*QGap12M(n3,1)-QGap12M(n2,1)*QGap12M(n1+n3,2)
		- QGap12M(n1,1)*QGap12M(n2+n3,2)+2.*QGap12M(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::ThreeGap12_subP(int n1, int n2, int n3)
{

	TComplex formula = QGap12P(n1,1)*QGap12P(n2,1)*QGap12P(n3,1)-QGap12P(n1+n2,2)*QGap12P(n3,1)-QGap12P(n2,1)*QGap12P(n1+n3,2)
		- QGap12P(n1,1)*QGap12P(n2+n3,2)+2.*QGap12P(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::Three_3Sub(int n1, int n2, int n3)
{

	TComplex formula = QsubLeft(n1,1)*QsubMiddle(n2,1)*QsubRight(n3,1);
	return formula;
}

//____________________________________________________________________
TComplex CorrelationCalculator::Four(int n1, int n2, int n3, int n4)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
		- Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
		+ Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
		+ 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::FourGap0(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap0P(n1,1)*QGap0P(n2,1)*QGap0M(n3,1)*QGap0M(n4,1)-QGap0P(n1+n2,2)*QGap0M(n3,1)*QGap0M(n4,1)
		-QGap0P(n1,1)*QGap0P(n2,1)*QGap0M(n3+n4,2)+QGap0P(n1+n2,2)*QGap0M(n3+n4,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::FourGap0M(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap0M(n1,1)*QGap0M(n2,1)*QGap0M(n3,1)*QGap0M(n4,1)-QGap0M(n1+n2,2)*QGap0M(n3,1)*QGap0M(n4,1)-QGap0M(n2,1)*QGap0M(n1+n3,2)*QGap0M(n4,1)
		- QGap0M(n1,1)*QGap0M(n2+n3,2)*QGap0M(n4,1)+2.*QGap0M(n1+n2+n3,3)*QGap0M(n4,1)-QGap0M(n2,1)*QGap0M(n3,1)*QGap0M(n1+n4,2)
		+ QGap0M(n2+n3,2)*QGap0M(n1+n4,2)-QGap0M(n1,1)*QGap0M(n3,1)*QGap0M(n2+n4,2)+QGap0M(n1+n3,2)*QGap0M(n2+n4,2)
		+ 2.*QGap0M(n3,1)*QGap0M(n1+n2+n4,3)-QGap0M(n1,1)*QGap0M(n2,1)*QGap0M(n3+n4,2)+QGap0M(n1+n2,2)*QGap0M(n3+n4,2)
		+ 2.*QGap0M(n2,1)*QGap0M(n1+n3+n4,3)+2.*QGap0M(n1,1)*QGap0M(n2+n3+n4,3)-6.*QGap0M(n1+n2+n3+n4,4);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::FourGap0P(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap0P(n1,1)*QGap0P(n2,1)*QGap0P(n3,1)*QGap0P(n4,1)-QGap0P(n1+n2,2)*QGap0P(n3,1)*QGap0P(n4,1)-QGap0P(n2,1)*QGap0P(n1+n3,2)*QGap0P(n4,1)
		- QGap0P(n1,1)*QGap0P(n2+n3,2)*QGap0P(n4,1)+2.*QGap0P(n1+n2+n3,3)*QGap0P(n4,1)-QGap0P(n2,1)*QGap0P(n3,1)*QGap0P(n1+n4,2)
		+ QGap0P(n2+n3,2)*QGap0P(n1+n4,2)-QGap0P(n1,1)*QGap0P(n3,1)*QGap0P(n2+n4,2)+QGap0P(n1+n3,2)*QGap0P(n2+n4,2)
		+ 2.*QGap0P(n3,1)*QGap0P(n1+n2+n4,3)-QGap0P(n1,1)*QGap0P(n2,1)*QGap0P(n3+n4,2)+QGap0P(n1+n2,2)*QGap0P(n3+n4,2)
		+ 2.*QGap0P(n2,1)*QGap0P(n1+n3+n4,3)+2.*QGap0P(n1,1)*QGap0P(n2+n3+n4,3)-6.*QGap0P(n1+n2+n3+n4,4);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::FourGap2(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap2P(n1,1)*QGap2P(n2,1)*QGap2M(n3,1)*QGap2M(n4,1)-QGap2P(n1+n2,2)*QGap2M(n3,1)*QGap2M(n4,1)
		-QGap2P(n1,1)*QGap2P(n2,1)*QGap2M(n3+n4,2)+QGap2P(n1+n2,2)*QGap2M(n3+n4,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::FourGap4(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap4P(n1,1)*QGap4P(n2,1)*QGap4M(n3,1)*QGap4M(n4,1)-QGap4P(n1+n2,2)*QGap4M(n3,1)*QGap4M(n4,1)
		-QGap4P(n1,1)*QGap4P(n2,1)*QGap4M(n3+n4,2)+QGap4P(n1+n2,2)*QGap4M(n3+n4,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::FourGap6(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap6P(n1,1)*QGap6P(n2,1)*QGap6M(n3,1)*QGap6M(n4,1)-QGap6P(n1+n2,2)*QGap6M(n3,1)*QGap6M(n4,1)
		-QGap6P(n1,1)*QGap6P(n2,1)*QGap6M(n3+n4,2)+QGap6P(n1+n2,2)*QGap6M(n3+n4,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::FourGap8(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap8P(n1,1)*QGap8P(n2,1)*QGap8M(n3,1)*QGap8M(n4,1)-QGap8P(n1+n2,2)*QGap8M(n3,1)*QGap8M(n4,1)
		-QGap8P(n1,1)*QGap8P(n2,1)*QGap8M(n3+n4,2)+QGap8P(n1+n2,2)*QGap8M(n3+n4,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::FourGap10(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap10P(n1,1)*QGap10P(n2,1)*QGap10M(n3,1)*QGap10M(n4,1)-QGap10P(n1+n2,2)*QGap10M(n3,1)*QGap10M(n4,1)
		-QGap10P(n1,1)*QGap10P(n2,1)*QGap10M(n3+n4,2)+QGap10P(n1+n2,2)*QGap10M(n3+n4,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::FourGap12(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap12P(n1,1)*QGap12P(n2,1)*QGap12M(n3,1)*QGap12M(n4,1)-QGap12P(n1+n2,2)*QGap12M(n3,1)*QGap12M(n4,1)
		-QGap12P(n1,1)*QGap12P(n2,1)*QGap12M(n3+n4,2)+QGap12P(n1+n2,2)*QGap12M(n3+n4,2);
	return formula;

}
//____________________________________________________________________
TComplex CorrelationCalculator::Four_3SubMMLR(int n1, int n2, int n3, int n4)
{
	TComplex formula = QsubMiddle(n1,1)*QsubMiddle(n2,1)*QsubLeft(n3,1)*QsubRight(n4,1)-QsubMiddle(n1+n2,2)*QsubLeft(n3,1)*QsubRight(n4,1);
	return formula;
}

//____________________________________________________________________
TComplex CorrelationCalculator::Four_3SubLLMR(int n1, int n2, int n3, int n4)
{
	TComplex formula = QsubLeft(n1,1)*QsubLeft(n2,1)*QsubMiddle(n3,1)*QsubRight(n4,1)-QsubLeft(n1+n2,2)*QsubMiddle(n3,1)*QsubRight(n4,1);
	return formula;
}

//____________________________________________________________________
TComplex CorrelationCalculator::Four_3SubRRML(int n1, int n2, int n3, int n4)
{
	TComplex formula = QsubRight(n1,1)*QsubRight(n2,1)*QsubMiddle(n3,1)*QsubLeft(n4,1)-QsubRight(n1+n2,2)*QsubMiddle(n3,1)*QsubLeft(n4,1);
	return formula;
}

//____________________________________________________________________
TComplex CorrelationCalculator::Four_3SubGap2Evts(int n1, int n2, int n3, int n4)
{
	TComplex formula = QsubGap2Middle(n1,1)*QsubGap2Middle(n2,1)*QsubGap2Left(n3,1)*QsubGap2Right(n4,1)
		-QsubGap2Middle(n1+n2,2)*QsubGap2Left(n3,1)*QsubGap2Right(n4,1);
	return formula;
}
//___________________________________________________________________
TComplex CorrelationCalculator::Five(int n1, int n2, int n3, int n4, int n5)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)
		- Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)
		+ 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)
		+ Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)
		+ Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)
		- Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)
		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)
		- 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)
		+ Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)
		+ Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)
		- Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)
		+ Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)
		- 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)
		- 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)
		+ Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)
		+ Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)
		+ 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)
		+ 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)
		- 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)
		+ Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)
		+ Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)
		+ 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)
		+ 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)
		- 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)
		- 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)
		- 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)+24.*Q(n1+n2+n3+n4+n5,5);
	return formula;

}
//_____________________________________________________________________________
TComplex CorrelationCalculator::FiveGap0A(int n1, int n2, int n3, int n4, int n5) // ( + +,- - - )
{

	TComplex formula = TwoGap0M(n1, n2) * ThreeGap0_subP(n3, n4, n5);
	//TComplex formula = (QGap0M(n1,1)*QGap0M(n2,1) - QGap0M(n1+n2,2)) * (QGap0P(n3,1)*QGap0P(n4,1)*QGap0P(n5,1)-QGap0P(n3+n4,2)*QGap0P(n5,1)-QGap0P(n4,1)*QGap0P(n3+n5,2) - QGap0P(n3,1)*QGap0P(n4+n5,2)+2.*QGap0P(n3+n4+n5,3));
	return formula;

}

//_____________________________________________________________________________
TComplex CorrelationCalculator::FiveGap0A_2(int n1, int n2, int n3, int n4, int n5) // ( + +,- - - )
{

	//TComplex formula = TwoGap0M(n1, n2) * ThreeGap0_subP(n3, n4, n5);
	TComplex formula = (QGap0M(n1,1)*QGap0M(n2,1) - QGap0M(n1+n2,2)) * (QGap0P(n3,1)*QGap0P(n4,1)*QGap0P(n5,1)-QGap0P(n3+n4,2)*QGap0P(n5,1)-QGap0P(n4,1)*QGap0P(n3+n5,2) - QGap0P(n3,1)*QGap0P(n4+n5,2)+2.*QGap0P(n3+n4+n5,3));
	return formula;

}

//_____________________________________________________________________________
TComplex CorrelationCalculator::FiveGap0B(int n1, int n2, int n3, int n4, int n5) // (- - -, + +)
{

	TComplex formula = TwoGap0P(n1, n2) * ThreeGap0_subM(n3, n4, n5);
	return formula;

}
//_____________________________________________________________________________
TComplex CorrelationCalculator::FiveGap2A(int n1, int n2, int n3, int n4, int n5) // ( + +,- - - )
{

	TComplex formula = TwoGap2M(n1, n2) * ThreeGap2_subP(n3, n4, n5);
	return formula;

}

//_____________________________________________________________________________
TComplex CorrelationCalculator::FiveGap2B(int n1, int n2, int n3, int n4, int n5) // (- - -, + +)
{

	TComplex formula = TwoGap2P(n1, n2) * ThreeGap2_subM(n3, n4, n5);
	return formula;

}
//_____________________________________________________________________________
TComplex CorrelationCalculator::FiveGap4A(int n1, int n2, int n3, int n4, int n5) // ( + +,- - - )
{

	TComplex formula = TwoGap4M(n1, n2) * ThreeGap4_subP(n3, n4, n5);
	return formula;

}

//_____________________________________________________________________________
TComplex CorrelationCalculator::FiveGap4B(int n1, int n2, int n3, int n4, int n5) // (- - -, + +)
{

	TComplex formula = TwoGap4P(n1, n2) * ThreeGap4_subM(n3, n4, n5);
	return formula;

}
//_____________________________________________________________________________
TComplex CorrelationCalculator::FiveGap6A(int n1, int n2, int n3, int n4, int n5) // ( + +,- - - )
{

	TComplex formula = TwoGap6M(n1, n2) * ThreeGap6_subP(n3, n4, n5);
	return formula;

}

//_____________________________________________________________________________
TComplex CorrelationCalculator::FiveGap6B(int n1, int n2, int n3, int n4, int n5) // (- - -, + +)
{

	TComplex formula = TwoGap6P(n1, n2) * ThreeGap6_subM(n3, n4, n5);
	return formula;

}
//_____________________________________________________________________________
TComplex CorrelationCalculator::FiveGap8A(int n1, int n2, int n3, int n4, int n5) // ( + +,- - - )
{

	TComplex formula = TwoGap8M(n1, n2) * ThreeGap8_subP(n3, n4, n5);
	return formula;

}

//_____________________________________________________________________________
TComplex CorrelationCalculator::FiveGap8B(int n1, int n2, int n3, int n4, int n5) // (- - -, + +)
{

	TComplex formula = TwoGap8P(n1, n2) * ThreeGap8_subM(n3, n4, n5);
	return formula;

}
//_____________________________________________________________________________
TComplex CorrelationCalculator::FiveGap10A(int n1, int n2, int n3, int n4, int n5) // ( + +,- - - )
{

	TComplex formula = TwoGap10M(n1, n2) * ThreeGap10_subP(n3, n4, n5);
	return formula;

}

//_____________________________________________________________________________
TComplex CorrelationCalculator::FiveGap10B(int n1, int n2, int n3, int n4, int n5) // (- - -, + +)
{

	TComplex formula = TwoGap10P(n1, n2) * ThreeGap10_subM(n3, n4, n5);
	return formula;

}
//_____________________________________________________________________________
TComplex CorrelationCalculator::FiveGap12A(int n1, int n2, int n3, int n4, int n5) // ( + +,- - - )
{

	TComplex formula = TwoGap12M(n1, n2) * ThreeGap12_subP(n3, n4, n5);
	return formula;

}

//_____________________________________________________________________________
TComplex CorrelationCalculator::FiveGap12B(int n1, int n2, int n3, int n4, int n5) // (- - -, + +)
{

	TComplex formula = TwoGap12P(n1, n2) * ThreeGap12_subM(n3, n4, n5);
	return formula;

}
//___________________________________________________________________
TComplex CorrelationCalculator::Six(int n1, int n2, int n3, int n4, int n5, int n6)
{


	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)
		- Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)
		+ 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)
		+ Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)
		+ Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)
		- Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)
		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)
		- 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)
		+ Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)
		+ Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6,1)
		- Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)
		+ Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)
		- 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)
		- 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)
		+ Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)
		+ Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)
		+ 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6,1)
		+ 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)
		- 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)
		+ Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)
		+ Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)
		+ 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6,1)
		+ 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)
		- 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)
		- 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6,1)
		- 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)+24.*Q(n1+n2+n3+n4+n5,5)*Q(n6,1)
		- Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)
		+ Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)
		- 2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6,2)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)
		- Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6,2)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)
		- Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6,2)-2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6,2)
		+ Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)-Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6,2)
		- 2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6,2)-2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6,2)
		+ 6.*Q(n2+n3+n4+n5,4)*Q(n1+n6,2)-Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)
		+ Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)+Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)
		+ Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)-2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)
		+ Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)-Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6,2)
		+ Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)-Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)
		- 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6,2)+Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)
		- Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)-2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6,2)
		- 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)+6.*Q(n1+n3+n4+n5,4)*Q(n2+n6,2)
		+ 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)-2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6,3)
		- 2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6,3)-2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6,3)
		+ 4.*Q(n3+n4+n5,3)*Q(n1+n2+n6,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)
		+ Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)+Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)
		+ Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)-2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)
		+ Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)-Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6,2)
		+ Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)-Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)
		- 2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6,2)+Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)
		- Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)-2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6,2)
		- 2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)+6.*Q(n1+n2+n4+n5,4)*Q(n3+n6,2)
		+ 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6,3)
		- 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6,3)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6,3)
		+ 4.*Q(n2+n4+n5,3)*Q(n1+n3+n6,3)+2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)
		- 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)-2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6,3)
		- 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)+4.*Q(n1+n4+n5,3)*Q(n2+n3+n6,3)
		- 6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6,4)+6.*Q(n4+n5,2)*Q(n1+n2+n3+n6,4)
		- Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)+Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)
		+ Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)+Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)
		- 2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)+Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)
		- Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6,2)+Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)
		- Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)-2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6,2)
		+ Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)-Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)
		- 2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6,2)-2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)
		+ 6.*Q(n1+n2+n3+n5,4)*Q(n4+n6,2)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)
		- 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6,3)-2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6,3)
		- 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6,3)+4.*Q(n2+n3+n5,3)*Q(n1+n4+n6,3)
		+ 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)
		- 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6,3)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)
		+ 4.*Q(n1+n3+n5,3)*Q(n2+n4+n6,3)-6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6,4)
		+ 6.*Q(n3+n5,2)*Q(n1+n2+n4+n6,4)+2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)
		- 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)-2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6,3)
		- 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)+4.*Q(n1+n2+n5,3)*Q(n3+n4+n6,3)
		- 6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6,4)+6.*Q(n2+n5,2)*Q(n1+n3+n4+n6,4)
		- 6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)+6.*Q(n1+n5,2)*Q(n2+n3+n4+n6,4)
		+ 24.*Q(n5,1)*Q(n1+n2+n3+n4+n6,5)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)
		+ Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)
		+ Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)
		+ Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6,2)
		+ Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)
		- 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6,2)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)
		- Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6,2)
		- 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)+6.*Q(n1+n2+n3+n4,4)*Q(n5+n6,2)
		+ 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6,3)
		- 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6,3)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6,3)
		+ 4.*Q(n2+n3+n4,3)*Q(n1+n5+n6,3)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)
		- 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6,3)
		- 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)+4.*Q(n1+n3+n4,3)*Q(n2+n5+n6,3)
		- 6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6,4)+6.*Q(n3+n4,2)*Q(n1+n2+n5+n6,4)
		+ 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)
		- 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6,3)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)
		+ 4.*Q(n1+n2+n4,3)*Q(n3+n5+n6,3)-6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6,4)
		+ 6.*Q(n2+n4,2)*Q(n1+n3+n5+n6,4)-6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)
		+ 6.*Q(n1+n4,2)*Q(n2+n3+n5+n6,4)+24.*Q(n4,1)*Q(n1+n2+n3+n5+n6,5)
		+ 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)
		- 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)
		+ 4.*Q(n1+n2+n3,3)*Q(n4+n5+n6,3)-6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6,4)
		+ 6.*Q(n2+n3,2)*Q(n1+n4+n5+n6,4)-6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)
		+ 6.*Q(n1+n3,2)*Q(n2+n4+n5+n6,4)+24.*Q(n3,1)*Q(n1+n2+n4+n5+n6,5)
		- 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)+6.*Q(n1+n2,2)*Q(n3+n4+n5+n6,4)
		+ 24.*Q(n2,1)*Q(n1+n3+n4+n5+n6,5)+24.*Q(n1,1)*Q(n2+n3+n4+n5+n6,5)
		- 120.*Q(n1+n2+n3+n4+n5+n6,6);
	return formula;

}
//_____________________________________________________________________________
TComplex CorrelationCalculator::SixGap0(int n1, int n2, int n3, int n4, int n5, int n6)
{

	TComplex formula = ThreeGap0_subM(n1, n2, n3)*ThreeGap0_subP(n4, n5, n6);
	return formula;

}
//_____________________________________________________________________________
TComplex CorrelationCalculator::SixGap2(int n1, int n2, int n3, int n4, int n5, int n6)
{

	TComplex formula = ThreeGap2_subM(n1, n2, n3)*ThreeGap2_subP(n4, n5, n6);
	return formula;

}
//_____________________________________________________________________________
TComplex CorrelationCalculator::SixGap4(int n1, int n2, int n3, int n4, int n5, int n6)
{

	TComplex formula = ThreeGap4_subM(n1, n2, n3)*ThreeGap4_subP(n4, n5, n6);
	return formula;

}
//_____________________________________________________________________________
TComplex CorrelationCalculator::SixGap6(int n1, int n2, int n3, int n4, int n5, int n6)
{

	TComplex formula = ThreeGap6_subM(n1, n2, n3)*ThreeGap6_subP(n4, n5, n6);
	return formula;

}
//_____________________________________________________________________________
TComplex CorrelationCalculator::SixGap8(int n1, int n2, int n3, int n4, int n5, int n6)
{

	TComplex formula = ThreeGap8_subM(n1, n2, n3)*ThreeGap8_subP(n4, n5, n6);
	return formula;

}
//_____________________________________________________________________________
TComplex CorrelationCalculator::SixGap10(int n1, int n2, int n3, int n4, int n5, int n6)
{

	TComplex formula = ThreeGap10_subM(n1, n2, n3)*ThreeGap10_subP(n4, n5, n6);
	return formula;

}
//_____________________________________________________________________________
TComplex CorrelationCalculator::SixGap12(int n1, int n2, int n3, int n4, int n5, int n6)
{

	TComplex formula = ThreeGap12_subM(n1, n2, n3)*ThreeGap12_subP(n4, n5, n6);
	return formula;

}
//_________________________________________________________________________________
TComplex CorrelationCalculator::Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7)
{

	TComplex Correlation = {0, 0};
	int Narray[] = {n1, n2, n3, n4, n5, n6};

	for(int k=7; k-->0; )
	{// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

		int array[6] = {0,1,2,3,4,5};
		int iPerm = 0;
		int argument = 0;
		int count = 0;

		// k==6: there is just one combination, we can add it manually
		if(k==6){
			Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
				Six(n1, n2, n3, n4, n5, n6)*Q(n7, 7-k);
		}// k==6

		else if(k==5){
			do{
				iPerm += 1;
				if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
					count += 1;
					Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
						Five(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
								Narray[int(array[3])], Narray[int(array[4])])*
						Q(Narray[int(array[5])]+n7, 7-k);
				}
			}while(std::next_permutation(array, array+6));
		}// k==5

		else if(k==4){
			do{
				iPerm += 1;
				if(iPerm%2 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
							Four(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
									Narray[int(array[3])])*
							Q(Narray[int(array[4])]+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==4

		else if(k==3){
			do{
				iPerm += 1;
				if(iPerm%6 == 1){
					if(array[0] < array[1] && array[1] < array[2]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
							Three(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])])*
							Q(Narray[int(array[3])]+Narray[int(array[4])]+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==3

		else if(k==2){
			do{
				iPerm += 1;
				if(iPerm%24 == 1){
					if(array[0] < array[1]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
							Two(Narray[int(array[0])], Narray[int(array[1])])*
							Q(Narray[int(array[2])]+Narray[int(array[3])]+Narray[int(array[4])]
									+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==2

		else if(k == 1){
			Correlation = Correlation
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7, 7-k);
		}// k==1

		else if(k == 0){
			Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1+n2+n3+n4+n5+n6+n7, 7-k);
		}// k==0

		else{
			std::cout<<"invalid range of k"<<std::endl;
			return {0,0};
		}

	}// loop over k

	return Correlation;

}
//_____________________________________________________________________________
TComplex CorrelationCalculator::Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{

	TComplex Correlation = {0, 0};
	int Narray[] = {n1, n2, n3, n4, n5, n6, n7};

	for(int k=8; k-->0; )
	{// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

		int array[7] = {0,1,2,3,4,5,6};
		int iPerm = 0;
		int argument = 0;
		int count = 0;

		// k==7: there is just one combination, we can add it manually
		if(k==7){
			Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
				Seven(n1, n2, n3, n4, n5, n6, n7)*Q(n8, 8-k);
		}// k==7

		else if(k==6){
			do{
				iPerm += 1;
				if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4] && array[4] < array[5]){
					count += 1;
					Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
						Six(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
								Narray[int(array[3])], Narray[int(array[4])], Narray[int(array[5])])*
						Q(Narray[int(array[6])]+n8, 8-k);
				}
			}while(std::next_permutation(array, array+7));
		}// k==6

		else if(k==5){
			do{
				iPerm += 1;
				if(iPerm%2 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
							Five(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
									Narray[int(array[3])], Narray[int(array[4])])*
							Q(Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==5

		else if(k==4){
			do{
				iPerm += 1;
				if(iPerm%6 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
							Four(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])], Narray[int(array[3])])*
							Q(Narray[int(array[4])]+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==4

		else if(k==3){
			do{
				iPerm += 1;
				if(iPerm%24 == 1){
					if(array[0] < array[1] && array[1] < array[2]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
							Three(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])])*
							Q(Narray[int(array[3])]+Narray[int(array[4])]+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==3

		else if(k==2){
			do{
				iPerm += 1;
				if(iPerm%120 == 1){
					if(array[0] < array[1]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
							Two(Narray[int(array[0])], Narray[int(array[1])])*
							Q(Narray[int(array[2])]+Narray[int(array[3])]+Narray[int(array[4])]
									+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==2

		else if(k == 1){
			Correlation = Correlation
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n7, 1)*Q(n1+n2+n3+n4+n5+n6+n8, 8-k);
		}// k==1

		else if(k == 0){
			Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1+n2+n3+n4+n5+n6+n7+n8, 8-k);
		}// k==0

		else{
			std::cout<<"invalid range of k"<<std::endl;
			return {0,0};
		}

	}// loop over k

	return Correlation;

}

//_____________________________________________________________________________
TComplex CorrelationCalculator::EightGap0(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{

	TComplex formula = FourGap0M(n1, n2, n3, n4)*FourGap0P(n5, n6, n7, n8);
	return formula;

}

void CorrelationCalculator::FillQVector(TComplex _Qvector[20][20], double _Qcos[20][20], double _Qsin[20][20]) {
	for(int iharm=0; iharm<20; iharm++)
	{
		for(int ipow=0; ipow<20; ipow++)
		{
			_Qvector[iharm][ipow] = TComplex(_Qcos[iharm][ipow], _Qsin[iharm][ipow]);
		}
	}
}

