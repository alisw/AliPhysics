//This file is part of CRAB
//http://www.nscl.msu.edu/~pratt/freecodes/crab/home.html

/*--------------------------------------------------
complex.cxx 
Version 1.0.1

-- Some modification for speed
-- Some protections from overflow and underflow
-- Solving of equations of order 2,3,4:
     Solve2, Solve3, Solve4;
--------------------------------------------------*/

#ifndef __COMPLEX_CXX__
#define __COMPLEX_CXX__

#ifndef small_epsilon
#define small_epsilon 1e-12
#endif

#ifndef PI
#define PI 3.141592653589793238462643
#endif


#include <TMath.h>
#include <Riostream.h>
//#include <math.h>
#include <stdlib.h>

/////////////////////////// Helper Error Function

void ComplexError (char *msg)
{
	cout << endl <<  msg << endl;
	exit (1);
}

void ComplexWarning (char *msg)
{
	cout << endl <<  msg << endl;
}

class Complex
{
public:
	double Re,Im;

	Complex (double a=0, double b=0): Re(a),Im(b) {}

	double GetRe ()                  { return Re; }
	double GetIm ()                  { return Im; }
	void SetRe (double a)            { Re = a; }
	void SetIm (double a)            { Im = a; }
	void Set (double a, double b)    { Re = a; Im = b; }
	void SetPolar (double, double);
	void SetAbs (double);
	void SetArg (double);
	Complex operator !  () const;
	Complex operator !  ();
	Complex operator +  () const;
	Complex operator +  ();
	Complex operator -  () const;
	Complex operator -  ();
	Complex operator += (const Complex &);
	Complex operator += (double);
	Complex operator -= (const Complex &);
	Complex operator -= (double);
	Complex operator *= (const Complex &);
	Complex operator *= (double);
	Complex operator /= (const Complex &);
	Complex operator /= (double);

	friend ostream& operator << (ostream&, const Complex &);
	friend istream& operator >> (istream&, Complex &);
};

	const Complex ImUnit (0,1);

	void Zero (Complex &);
	Complex Polar (double, double);
        double  real  (const Complex &);
	double  imag (const Complex &);
	double  Arg  (const Complex &);
	double  phase(const Complex &);
	double  abs  (const Complex &);
    void Solve2 (Complex*, const Complex&, const Complex&);
    void Solve3 (Complex*, const Complex&, const Complex&, const Complex&);
    void Solve4 (Complex*, const Complex&, const Complex&, const Complex&, const Complex&);
    Complex Solve2 (const Complex&, const Complex&, int RootNumber);
    Complex Solve3 (const Complex&, const Complex&, const Complex&, int RootNumber);
    Complex Solve4 (const Complex&, const Complex&, const Complex&, const Complex&, int RootNumber);
	Complex conj (const Complex &);
	Complex sqr  (const Complex &);
	Complex sqrt (const Complex &, int=0);
	Complex pow  (Complex, int);
	Complex pow  (Complex, double);
	Complex pow  (Complex, Complex);
	Complex root (const Complex &, int, int=0);
	Complex exp  (const Complex &);
	Complex log  (const Complex &);
	Complex sin  (const Complex &);
	Complex cos  (const Complex &);
	Complex tan  (const Complex &);
	Complex cot  (const Complex &);
	Complex sec  (const Complex &);
	Complex csc  (const Complex &);
	Complex sinh (const Complex &);
	Complex cosh (const Complex &);
	Complex tanh (const Complex &);
	Complex coth (const Complex &);
	Complex sech (const Complex &);
	Complex csch (const Complex &);
	Complex asin (const Complex &, int=0);
	Complex acos (const Complex &, int=0);
	Complex atan (const Complex &);
	Complex acot (const Complex &);
	Complex asec (const Complex &, int=0);
	Complex acsc (const Complex &, int=0);
	Complex asinh(const Complex &, int=0);
	Complex acosh(const Complex &, int=0);
	Complex atanh(const Complex &);
	Complex acoth(const Complex &);
	Complex asech(const Complex &, int=0);
	Complex acsch(const Complex &, int=0);

	Complex operator ^  (const Complex &, int);
	Complex operator ^  (const Complex &, double);
	Complex operator ^  (const Complex &, const Complex &);

	double operator += (double&, const Complex &);
	double operator -= (double&, const Complex &);
	double operator *= (double&, const Complex &);
	double operator /= (double&, const Complex &);

	Complex operator +  (const Complex &, const Complex &);
	Complex operator +  (const Complex &, double);
	Complex operator +  (double, const Complex &);
	Complex operator -  (const Complex &, const Complex &);
	Complex operator -  (const Complex &, double);
	Complex operator -  (double, const Complex &);
	Complex operator *  (const Complex &, const Complex &);
	Complex operator *  (const Complex &, double);
	Complex operator *  (double, const Complex &);
	Complex operator /  (const Complex &, const Complex &);
	Complex operator /  (const Complex &, double);
	Complex operator /  (double, const Complex &);

	int operator == (const Complex &, double);
	int operator == (double, const Complex &);
	int operator == (const Complex &, const Complex &);
	int operator != (const Complex &, double);
	int operator != (double, const Complex &);
	int operator != (const Complex &, const Complex &);
	int operator >= (const Complex &, double);
	int operator >= (double, const Complex &);
	int operator >= (const Complex &, const Complex &);
	int operator <= (const Complex &, double);
	int operator <= (double, const Complex &);
	int operator <= (const Complex &, const Complex &);
	int operator > (const Complex &, double);
	int operator > (double, const Complex &);
	int operator > (const Complex &, const Complex &);
	int operator < (const Complex &, double);
	int operator < (double, const Complex &);
	int operator < (const Complex &, const Complex &);

////////////////////////////////////////////////

inline double sqr(double a)  {	return a*a; }

inline Complex Complex::operator + () const
{
	return *this;
}

inline Complex Complex::operator + ()
{
	return *this;
}

inline Complex Complex::operator - () const
{
	return Complex (-Re, -Im);
}

inline Complex Complex::operator - ()
{
	return Complex (-Re, -Im);
}

inline Complex Complex::operator ! () const  //  complex conjugate  //
{
	return Complex (Re, -Im);
}

inline Complex Complex::operator ! ()       //  complex conjugate  //
{
	return Complex (Re, -Im);
}

/////////////////////// +=,-=,*=,/=

double operator += (double &a, const Complex &b)
{
	if (TMath::Abs(b.Im) < small_epsilon) a+=b.Re;
	else
		ComplexError ("Error in double+=Complex: Complex is not double number");
	return a;
}

double operator -= (double &a, const Complex &b)
{
	if (TMath::Abs(b.Im) < small_epsilon) a-=b.Re;
	else
		ComplexError ("Error in double -=Complex: Complex is not double number");
	return a;
}

double operator *= (double &a, const Complex &b)
{
	if (TMath::Abs(b.Im) < small_epsilon) a*=b.Re;
	else
		ComplexError ("Error in double *=Complex: Complex is not double number");
	return a;
}

double operator /= (double &a, const Complex &b)
{
	if (TMath::Abs(b.Im) < small_epsilon) a/=b.Re;
	else
		ComplexError ("Error in double /=Complex: Complex is not double number");
	return a;
}

inline Complex Complex::operator += (const Complex &a)
{
	Re+=a.Re;
	Im+=a.Im;
	return *this;
}

inline Complex Complex::operator += (double a)
{
	Re+=a;
	return *this;
}

inline Complex Complex::operator -= (const Complex &a)
{
	Re-=a.Re;
	Im-=a.Im;
	return *this;
}

inline Complex Complex::operator -= (double a)
{
	Re-=a;
	return *this;
}

Complex Complex::operator *= (const Complex &a)
{
	double t1=Re*a.Re, t2=Im*a.Im;
	Im = (Re+Im) * (a.Re+a.Im) - t1 - t2;
	Re = t1 - t2;
	return *this;
}

inline Complex Complex::operator *= (double a)
{
	Re*=a;
	Im*=a;
	return *this;
}

Complex Complex::operator /= (const Complex &a)
{
    double t1, t2, temp;
    if (TMath::Abs(a.Re) >= TMath::Abs(a.Im))
    {
      t1   = a.Im / a.Re;
      t2   = a.Re + a.Im * t1;
      temp = (Re + Im * t1) / t2;
      Im   = (Im - Re * t1) / t2;
      Re   = temp;
    }
    else
    {
      t1   = a.Re / a.Im;
      t2   = a.Re * t1 + a.Im;
      temp = (Re * t1 + Im) / t2;
      Im   = (Im * t1 - Re) / t2;
      Re   = temp;
    }
	return *this;
}

inline Complex Complex::operator /= (double a)
{
	Re/=a;
	Im/=a;
	return *this;
}

/////////////////////////// +,-,*,/

inline Complex operator + (const Complex &a, const Complex &b)
{
	return Complex (a.Re+b.Re, a.Im+b.Im);
}

inline Complex operator + (const Complex &a, double b)
{
	return Complex (a.Re+b, a.Im);
}

inline Complex operator + (double b, const Complex &a)
{
	return Complex (a.Re+b, a.Im);
}

inline Complex operator - (const Complex &a, const Complex &b)
{
	return Complex (a.Re-b.Re, a.Im-b.Im);
}

inline Complex operator - (const Complex &a, double b)
{
	return Complex (a.Re-b, a.Im);
}

inline Complex operator - (double b, const Complex &a)
{
	return Complex (b-a.Re, -a.Im);
}

Complex operator * (const Complex &a, const Complex &b)
{
    double t1 = a.Re * b.Re;
    double t2 = a.Im * b.Im;
    return Complex (t1 - t2, (a.Re+a.Im) * (b.Re+b.Im) - t1 - t2);
}

inline Complex operator * (const Complex &a, double b)
{
	return Complex (a.Re*b, a.Im*b);
}

inline Complex operator * (double b, const Complex &a)
{
	return Complex (a.Re*b, a.Im*b);
}

inline Complex operator / (const Complex &a, double b)
{
	return Complex (a.Re/b, a.Im/b);
}

inline Complex operator / (const Complex &a, const Complex &b)
{
    double t1, t2;
    if (TMath::Abs(b.Re) >= TMath::Abs(b.Im))
    {
      t1   = b.Im / b.Re;
      t2   = b.Re + b.Im * t1;
      return Complex ((a.Re + a.Im * t1) / t2, (a.Im - a.Re * t1) / t2);
    }
    else
    {
      t1   = b.Re / b.Im;
      t2   = b.Re * t1 + b.Im;
      return Complex ((a.Re * t1 + a.Im) / t2, (a.Im * t1 - a.Re) / t2);
    }
}

inline Complex operator / (double b, const Complex &a)
{
	return (Complex (b,0)) / a;
}

//////////////////////// <<,>>

ostream& operator << (ostream &stream, const Complex &a)
{
	stream<<"   "<<a.Re<<"  "<<a.Im<<"  ";
	return stream;
}

istream& operator >> (istream &stream, Complex &a)
{
	stream>>a.Re>>a.Im;
	return stream;
}

/////////////////////// ==,!=

inline int operator == (const Complex &a, double b)
{
	return (a.Re==b) && (TMath::Abs(a.Im)<small_epsilon);
}

inline int operator == (double b, const Complex &a)
{
	return (a.Re==b) && (TMath::Abs(a.Im)<small_epsilon);
}

inline int operator == (const Complex &a, const Complex &b)
{
	return (a.Re==b.Re) && (a.Im==b.Im);
}

inline int operator != (const Complex &a, double b)
{
	return (a.Re!=b) || (TMath::Abs(a.Im)>small_epsilon);
}

inline int operator != (double b, const Complex &a)
{
	return (a.Re!=b) || (TMath::Abs(a.Im)>small_epsilon);
}

inline int operator != (const Complex &a, const Complex &b)
{
	return (a.Re!=b.Re) || (a.Im!=b.Im);
}

/////////////////////// >=,<=

inline int operator >= (const Complex &a, double b)
{
	return abs(a) >= b;
}

inline int operator >= (double b, const Complex &a)
{
	return b >= abs(a);
}

inline int operator >= (const Complex &a, const Complex &b)
{
	return abs(a) >= abs(b);
}

inline int operator <= (const Complex &a, double b)
{
	return abs(a) <= b;
}

inline int operator <= (double b, const Complex &a)
{
	return b <= abs(a);
}

inline int operator <= (const Complex &a, const Complex &b)
{
	return abs(a) <= abs(b);
}

/////////////////////// >,<

inline int operator > (const Complex &a, double b)
{
	return abs(a) > b;
}

inline int operator > (double b, const Complex &a)
{
	return b > abs(a);
}

inline int operator > (const Complex &a, const Complex &b)
{
	return abs(a) > abs(b);
}

inline int operator < (const Complex &a, double b)
{
	return abs(a) < b;
}

inline int operator < (double b, const Complex &a)
{
	return b < abs(a);
}

inline int operator < (const Complex &a, const Complex &b)
{
	return abs(a) < abs(b);
}

///////////////////////// Functions
double  real (const Complex &a)
{
	return a.Re;
}
double  imag (const Complex &a)
{
	return a.Im;
}
double abs (const Complex &a)
{
	if (a.Im == 0) return TMath::Abs(a.Re);
	if (a.Re == 0) return TMath::Abs(a.Im);
    double R = TMath::Abs(a.Re), I = TMath::Abs(a.Im);
	return (R >= I) ?
           (R * sqrt (1 + sqr(a.Im/a.Re))) :
           (I * sqrt (1 + sqr(a.Re/a.Im))) ;
}

double  Arg (const Complex &a)
{
	return atan2(a.Im,a.Re);
}

double  phase (const Complex &a)
{
	return atan2(a.Im,a.Re);
}

inline void Zero (Complex &a) {	a.Re=a.Im=0; }

Complex conj  (const Complex &a)              { return !a; }
Complex pow  (Complex a, int n)              { return a^n; }
Complex pow  (Complex a, double n)           { return a^n; }
Complex pow  (Complex a, Complex b)   { return a^b; }

inline Complex sqr (const Complex &a) { return a*a; }

Complex sqrt (const Complex &a, int flag)
{
	if ((a.Re>=0) && (TMath::Abs(a.Im)<small_epsilon))
       return flag ? -sqrt(a.Re) : sqrt(a.Re);
    double R = TMath::Abs(a.Re), I = TMath::Abs(a.Im);
	double w = (R >= I) ?
           sqrt (R/2 * (  1 + sqrt (1 + sqr(a.Im/a.Re)))):
           sqrt (I/2 * (R/I + sqrt (1 + sqr(a.Re/a.Im))));
    Complex c;
    if (a.Re >= 0)
    {
        c.Re = w;
        c.Im = a.Im / (2*w);
    }
    else
    {
        c.Re = I / (2*w);
        c.Im = (a.Im >= 0) ? w : -w;
    }
	return ((flag && (c.Re<0))  ||  (!flag && (c.Re>=0))) ? c : -c;
}

Complex operator ^ (const Complex &a, int n)
{
	Complex c(1,0);
	if (n==0) return c;
	if (n>0) for (int i=0;i<n;i++) c*=a;
			else for (int j=0;j>n;j--) c*=a;
	if (n>0) return c;
			else return 1/c;
}

Complex operator ^ (const Complex &a, double n)
{
	return exp(n*log(a));
}

Complex operator ^ (const Complex &a, const Complex &b)
{
	return exp(b*log(a));
}

Complex root (const Complex &z, int n, int k)
{
	double c=exp(log(abs(z))/n);
	double t=(Arg(z)+2*PI*k)/n;
	return Complex (c*cos(t), c*sin(t));
}

Complex exp (const Complex &a)
{
	double t=exp(a.Re);
	return Complex (t*cos(a.Im), t*sin(a.Im));
}

Complex log (const Complex &a)
{
	if (a==0)
		ComplexError("Error in function log(Complex): argument is 0");
	return Complex (log(abs(a)), Arg(a));
}

Complex sin (const Complex &a)
{
	return Complex (sin(a.Re)*cosh(a.Im), cos(a.Re)*sinh(a.Im));
}

Complex cos (const Complex &a)
{
	return Complex (cos(a.Re)*cosh(a.Im), -sin(a.Re)*sinh(a.Im));
}

Complex tan (const Complex &a)
{
	return sin(a)/cos(a);
}

Complex cot (const Complex &a)
{
	return cos(a)/sin(a);
}

Complex sec (const Complex &a)
{
	return 1/cos(a);
}

Complex csc (const Complex &a)
{
	return 1/sin(a);
}

Complex sinh (const Complex &a)
{
	return Complex (sinh(a.Re)*cos(a.Im), cosh(a.Re)*sin(a.Im));
}

Complex cosh (const Complex &a)
{
	return Complex (cosh(a.Re)*cos(a.Im), sinh(a.Re)*sin(a.Im));
}

Complex tanh (const Complex &a)
{
	return sinh(a)/cosh(a);
}

Complex coth (const Complex &a)
{
	return cosh(a)/sinh(a);
}

Complex sech (const Complex &a)
{
	return 1/cosh(a);
}

Complex csch (const Complex &a)
{
	return 1/sinh(a);
}
//////////////////////// Inverce trigonometric functions

Complex asin (const Complex &a, int flag)
{
	return -ImUnit * log(ImUnit*a + sqrt(1-sqr(a), flag));
}

Complex acos (const Complex &a, int flag)
{
	return -ImUnit * log(a + ImUnit*sqrt(1-sqr(a), flag));
}

Complex atan (const Complex &a)
{
	return ImUnit/2 * log((ImUnit+a)/(ImUnit-a));
}

Complex acot (const Complex &a)
{
	return ImUnit/2 * log((a-ImUnit)/(a+ImUnit));
}

Complex asec (const Complex &a, int flag)
{
	return acos(1/a, flag);
}

Complex acsc (const Complex &a, int flag)
{
	return asin(1/a, flag);
}

Complex asinh (const Complex &a, int flag)
{
	return log(a + sqrt(sqr(a)+1, flag));
}

Complex acosh (const Complex &a, int flag)
{
	return log(a + sqrt(sqr(a)-1, flag));
}

Complex atanh (const Complex &a)
{
	return log((1+a)/(1-a)) / 2;
}

Complex acoth (const Complex &a)
{
	 return log((a+1)/(a-1)) / 2;
}

Complex asech (const Complex &a, int flag)
{
	return acosh(1/a, flag);
}

Complex acsch (const Complex &a, int flag)
{
	return asinh(1/a, flag);
}

Complex Polar (double a, double b)
{
	return Complex (a*cos(b), a*sin(b));
}

void Complex::SetPolar (double a, double b)
{
	Re=a*cos(b);
	Im=a*sin(b);
}

void Complex::SetAbs (double a)
{
	double b=Arg(*this);
	Re=a*cos(b);
	Im=a*sin(b);
}

void Complex::SetArg (double b)
{
	double a=abs(*this);
	Re=a*cos(b);
	Im=a*sin(b);
}

void Solve2 (Complex* z, const Complex &b, const Complex &c)
   // finding z of equation:  z^2 + b*z + c = 0
{
  Complex t = sqrt(sqr(b)-4*c);
	Complex q = ((!b * t).Re >= 0) ? (-(b + t) / 2) : (-(b - t) / 2);
  z[0] = q;
  z[1] = c/q;
}

Complex Solve2 (const Complex &b, const Complex &c, int RootNumber)
{
	if ((RootNumber < 0) || (RootNumber > 1))
		ComplexError ("Error in Solve2: wrong root number");
	Complex z[2];
	Solve2 (z,b,c);
	return z[RootNumber];
}
    
void Solve3 (Complex* z, const Complex &a2, const Complex &a1, const Complex &a0)
   // finding z of equation:  z^3 + a2*z^2 + a1*z + a0 = 0
{
  Complex q, r, t, a, b, zero(0,0);
  q = (sqr(a2) - 3*a1) / 9;
  r = (2*(a2^3) - 9*a2*a1 + 27*a0) / 54;
  t = sqrt((r^2) - (q^3));
	a = ((!r * t) >=0.0) ? -((r+t)^(1.0/3)) : -((r-t)^(1.0/3));
  b = ((a == zero) ? zero : (q/a));
  z[0] = -(a+b)/2 - a2/3 + ImUnit*sqrt(3.0)*(a-b)/2;
  z[1] = -(a+b)/2 - a2/3 - ImUnit*sqrt(3.0)*(a-b)/2;
  z[2] = a + b - a2/3;
}

Complex Solve3 (const Complex &a2, const Complex &a1, const Complex &a0, int RootNumber)
{
	if ((RootNumber < 0) || (RootNumber > 2))
		ComplexError ("Error in Solve3: wrong root number");
	Complex z[3];
	Solve3 (z,a2,a1,a0);
	return z[RootNumber];
}

void Solve4 (Complex* z, const Complex &a3, const Complex &a2, const Complex &a1, const Complex &a0)
   // finding z of equation:  z^4 + a3*z^3 + a2*z^2 + a1*z + a0 = 0
{
  Complex u[3], t1, t2, t;
  Solve3 (u, -a2, (a1*a3 - 4*a0), -((a1^2) + a0*(a3^2) - 4*a0*a2));
	t = *u;
  t1 = sqrt((a3^2)/4 - a2 + t);
  t2 = sqrt((t^2)/4 - a0);
  Solve2 (      z, (a3/2 - t1), (t/2 + t2));
  Solve2 (&(z[2]), (a3/2 + t1), (t/2 - t2));
}

Complex Solve4 (const Complex &a3, const Complex &a2, const Complex &a1, const Complex &a0, int RootNumber)
{
	if ((RootNumber < 0) || (RootNumber > 3))
		ComplexError ("Error in Solve4: wrong root number");
	Complex z[4];
	Solve4 (z,a3,a2,a1,a0);
	return z[RootNumber];
}


#endif  // __COMPLEX_CXX__ //
