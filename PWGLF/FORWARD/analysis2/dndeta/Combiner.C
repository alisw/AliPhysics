#include <cmath>
#include <iterator>
#include <list>
#include <iosfwd>
// From http://www.slac.stanford.edu/~barlow/java/

struct Combiner
{
  /**
   * An experimental observation 
   */
  struct Result
  {
    /** Value @f$ x_i@f$ */
    double fX;
    /** Low error @f$ \sigma_i^-@f$ */
    double fEl;
    /** High error @f$ \sigma_i^+@f$ */
    double fEh;
    /** 
     * Create a single obersvation 
     * 
     * @param x  Value @f$ x_i@f$
     * @param el Low error @f$ \sigma_i^-@f$ on @f$ x_i@f$ 
     * @param eh High error @f$ \sigma_i^+@f$ on @f$ x_i@f$ 
     */
    Result(double x=0, double el=0, double eh=0)
      : fX(x), fEl(el), fEh(eh)
    {    
    }
    /** 
     * Virtual destructor 
     */
    virtual ~Result() {}
    /** 
     * Calculate 
     * 
     * @f[ 
     *   s_i = \sigma_i^+ \sigma_i^- / (\sigma_i^+ + \sigma_i^-)
     * @f]
     * 
     * 
     * @return @f$ s_i@f$
     */
    double S() const
    {
      return 2 * fEl * fEh / (fEl + fEh);
    }
    /** 
     * Calculate 
     * 
     * @f[ 
     *   s_i' = (\sigma_i^+ - \sigma_i^-) / (\sigma_i^+ + \sigma_i^-)
     * @f]
     * 
     * 
     * @return @f$ s_i'@f$
     */
    double Sprime() const
    {
      return (fEh - fEl) / (fEh + fEh);
    }
    /** 
     * Calculate 
     * 
     * @f[ 
     *   V_i = \sigma_i^+ \sigma_i^-
     * @f]
     * 
     * 
     * @return @f$ V_i@f$
     */
    double V() const
    {
      return fEl * fEh;
    }
    /** 
     * Calculate 
     * 
     * @f[ 
     *   V_i' = \sigma_i^+ - \sigma_i^-
     * @f]
     * 
     * 
     * @return @f$ V_i'@f$
     */
    double Vprime() const
    {
      return fEh - fEl;
    }
    /** 
     * Lower bound 
     * 
     * @return @f$ x_i - 3\sigma_i^+@f$ 
     */
    virtual double Low() const
    {
      return fX - 3 * fEl;
    }
    /** 
     * Upper bound 
     * 
     * @return @f$ x_i + 3\sigma_i^+@f$ 
     */
    virtual double High() const
    {
      return fX + 3 * fEh;
    }
  };
  
  /**
   * The final result 
   */
  struct Final : public Result
  {
    /** The final @f$\chi^2 @f$ */
    double fChi2;
    /** Lower bound to use */
    double fLower;
    /** Upper bound to use */
    double fUpper;
    /** 
     * The final result 
     * 
     * @param x     Best estimate of @f$ x@f$ 
     * @param el    Best estimate of @f$ \sigma^-@f$ 
     * @param eh    Best estimate of @f$ \sigma^+@f$ 
     * @param chi2  @f$\chi^2@f$ of best estimate 
     * @param low   The lower bound of the source data
     * @param high  The upper bound of the source data
     */
    Final(double x, double el, double eh,
	  double chi2, double low, double high)
      : Result(x, el, eh),
	fChi2(chi2),
	fLower(low),
	fUpper(high)
    {}
    /** 
     * Lower bound 
     * 
     * @return @f$ x_i - 3\sigma_i^+@f$ 
     */
    double Low() const
    {
      return fLower;
    }
    /** 
     * Upper bound 
     * 
     * @return @f$ x_i + 3\sigma_i^+@f$ 
     */
    double High() const
    {
      return fUpper;
    }
  };    

  /** 
   * A possible container of data 
   */
  struct List
  {
    /** The type of container */
    typedef std::list<Result> Container;
    /** Iterator */
    typedef Container::iterator iterator;
    /** Constant Iterator */
    typedef Container::const_iterator const_iterator;
    
    /** Our contianer */
    Container fData;
    /** 
     * Get iterator to the beginning 
     * 
     * @return Iterator pointing at the beginning 
     */
    iterator begin() { return fData.begin(); }
    /** 
     * Get iterator to the beginning 
     * 
     * @return Iterator pointing at the beginning 
     */
    const_iterator begin() const { return fData.begin(); }
    /** 
     * Get iterator to the ending 
     * 
     * @return Iterator pointing at the ending 
     */
    iterator end() { return fData.end(); }
    /** 
     * Get iterator to the ending 
     * 
     * @return Iterator pointing at the ending 
     */
    const_iterator end() const { return fData.end(); }

    /** 
     * Add an obervation 
     * 
     * @param r Observation
     */
    void Add(const Result& r)
    {
      fData.push_back(r);
    }
    /** 
     * Add an observation 
     * 
     * @param x  @f$ x_i@f$ 
     * @param el @f$ \sigma_i^-@f$  
     * @param eh @f$ \sigma_i^+@f$
     */
    void Add(double x, double el, double eh)
    {
      size_t n = fData.size();
      fData.resize(n+1);
      fData.back().fX = x;
      fData.back().fEl = el;
      fData.back().fEh = eh;
    }
  };
  typedef List::iterator iterator;
  typedef List::const_iterator const_iterator;
  
  /** 
   * Virtual destructor 
   */
  virtual ~Combiner() {}
  
  /** 
   * Calculate the weight 
   *
   * @param r Observation
   * 
   * @return @f$ W@f$ 
   */
  virtual double W(const Result& r) const = 0;
  /** 
   * Calculate the weight based on a guess of best @f$ x'@f$ 
   * 
   * @param guess Current guess @f$ x'@f$ 
   * @param r  	  Observation
   * 
   * @return @f$ W(x')@f$ 
   */
  virtual double StepW(double guess, const Result& r) const = 0;
  /** 
   * Calculate the bias. 
   * 
   * @return @f$\delta(x')@f$ 
   */
  virtual double StepOffset(double guess, const Result& r) const = 0;
  /** 
   * Calculate the contribution variance to the @f$\chi^2@f$ with
   * the guess @f$x'@f$. 
   *
   * @return @f$ v(x')@f$ 
   */
  virtual double TermVar(double guess, const Result& r) const = 0;
  /** 
   * Calculate the contribution variance to the @f$\chi^2@f$ with
   * the guess @f$ x'@f$. 
   *
   * @f[
   *   t_i(x') = (x' - x_i)^2 / v_i(x')
   * @f]
   *
   * where @f$ v_i(x')@f$ is the term variance
   *
   * @param guess @f$ x'@f$  
   * @param r     Obersvation 
   *
   * @return @f$ t(x')@f$ 
   */
  double Term(double guess, const Result& r) const
  {
    double var = TermVar(guess, r);
    if (var <= 0) return -1000;

    return std::pow(guess - r.fX, 2) / var;
  }
  /** 
   * Calculate the @f$ \chi^2(x')@f$ where @f$ x'@f$ is current guess
   * at the result.  
   * 
   * @param guess   Current guess @f$ x'@f$ 
   * @param chi2    Optional old @f$ \chi^2@f$ from best @f$ x@f$ value 
   * @param begin Iterator 
   * @param end   Iterator 
   * 
   * @return @f$ \chi^2(x')@f$
   */
  double F(double          guess,
	   double          chi2,
	   const_iterator& begin,
	   const_iterator& end) const
  {
    double s = -chi2;

    for (const_iterator i = begin; i != end; ++i) 
      s += Term(guess, *i);
    return s;
  }
  /** 
   * Try to find best error 
   * 
   * @param nIter Number of iterations 
   * @param begin Iterator 
   * @param end   Iterator 
   * @param sign  Direction (-1 is low, +1 is high)
   * @param best  Current best @f$ x@f$ value 
   * @param chi2  @f$ \chi^2@f$ of current best @f$ x@f$ value 
   * @param s     Summed weights in the direction 
   * 
   * @return The error in the chosen direction
   */
  double FindError(unsigned short   nIter,
		   const_iterator&  begin,
		   const_iterator&  end,
		   int              sign,
		   double           best,
		   double           chi2,
		   double           s)
  {
    // Step size 
    double delta = 0.1 * sign * s;

    // Iterations 
    for (unsigned short i = 0; i < nIter; i++) {
      // Calculate chi^2 with current guess 
      double got = F(best + sign * s, chi2, begin, end);

      if (std::fabs(got-1) < 1e-7)
	// We're close to 1 so get out
	break;

      // The next guess' chi^2 value e
      double guess = F(best + sign * s + delta, chi2, begin, end);

      // Where to search next 
      if ((got - 1) * (guess - 1) > 0) {
	if ((got - 1) / (guess - 1) < 1)
	  delta = -delta;
	else
	  s += sign * delta;
	continue;
      }

      // Update error value and decrease step size 
      s     += sign * delta * (1 - got) / (guess - got);
      delta /= 2;
    }
    return s;
  }
  /** 
   * Find best estimate of @f$ x@f$ 
   * 
   * @param nIter   Number of iterations 
   * @param begin Iterator 
   * @param end   Iterator 
   * @param lowest  Lower bound 
   * @param highest Upper bound 
   * 
   * @return  @f$ x@f$ 
   */
  double FindX(unsigned short   nIter,
	       const_iterator&  begin,
	       const_iterator&  end,
	       double           lowest,
	       double           highest)
  {
    // Starting values 
    double x    = (highest+lowest)/2;
    double oldX = -1e33;

    // Do the iterations 
    for (unsigned short i = 0; i < nIter; i++) {
      double sum    = 0;
      double sumw   = 0;
      double offset = 0;

      // Loop over observations 
      for (const_iterator j = begin; j != end; ++j) {
	double w =  StepW(x, *j);
	offset   += StepOffset(x,  *j);

	sum  += j->fX * w;
	sumw += w;
      }
      x = (sum - offset) / sumw;

      if (std::fabs(x - oldX) < (highest-lowest) * 1e-6) break;
      oldX = x;
    }
    return x;
  }
  /** 
   * Do the calculation 
   * 
   * @param begin Iterator 
   * @param end   Iterator 
   * @param nIter How many iterations to do. 
   * 
   * @return The best estimate of @f$ x@f$ and associated errors 
   */
  Final Calculate(const_iterator& begin,
		  const_iterator& end,
		  unsigned short nIter=50)
  {
    double lowest  = +1e33;
    double highest = -1e33;
    double sumLow  = 0;
    double sumHigh = 0;

    // Find boundaries and sum weights
    for (const_iterator i = begin; i != end; ++i) {
      lowest  = std::min(i->Low(),  lowest);
      highest = std::max(i->High(), highest);
      sumLow  = 1./std::pow(i->fEl, 2);
      sumHigh = 1./std::pow(i->fEh, 2);
    }
    // Summed weights 
    double sLow  = 1. / std::sqrt(sumLow);
    double sHigh = 1. / std::sqrt(sumHigh);
    
    // Now do the calculations
    double bestX    = FindX(nIter, begin, end, lowest, highest);
    double bestChi2 = F(bestX, 0, begin, end);
    double bestLow  = FindError(nIter,begin,end,-1,bestX,bestChi2,sLow);
    double bestHigh = FindError(nIter,begin,end,+1,bestX,bestChi2,sHigh);

    return Final(bestX, bestLow, bestHigh, bestChi2, lowest, highest);
  }
};

/** 
 * Put-to operator of an observation 
 * 
 * @param o Output stream 
 * @param r Observation 
 * 
 * @return The output stream 
 */
std::ostream& operator<<(std::ostream& o, const Combiner::Result& r)
{
  return o << r.fX << "\t-" << r.fEl << "\t+" << r.fEh;
}
/** 
 * Put-to operator of an observation 
 * 
 * @param o Output stream 
 * @param r Observation 
 * 
 * @return The output stream 
 */
std::ostream& operator<<(std::ostream& o, const Combiner::Final& r)
{
  return o << static_cast<const Combiner::Result&>(r)  << "\t" << r.fChi2;
}

/**
 * A combiner that uses a linear @f$\sigma@f$ approximation 
 */
struct LinearSigmaCombiner : public Combiner
{
  /** 
   * Calculate the weight 
   * 
   * @f[ 
   *   w = 1/2 (s + x s')^3 / s
   * @f]
   *
   * @param r Observation
   * 
   * @return @f$ W@f$ 
   */
  double W(const Result& r) const
  {
    double s = r.S();
    return .5 * std::pow(s + r.fX * r.Sprime(), 3) / s;
  }
  /** 
   * Calculate the weight based on a guess of best @f$ x'@f$ 
   * 
   * @f[ 
   *   w(x') = s / [s + s' (x' - x)]^3
   * @f]
   *
   * @param guess Current guess @f$ x'@f$ 
   * @param r	  Observation
   * 
   * @return @f$ W(x')@f$ 
   */
  double StepW(double guess, const Result& r) const
  {
    double s = r.S();
    return s / std::pow(s+r.Sprime()*(guess - r.fX),3);
  }
  /** 
   * Calculate the bias. 
   * 
   * @return 0
   */
  double StepOffset(double, const Result&) const
  {
    return 0;
  }
  /** 
   * Calculate the contribution variance to the @f$\chi^2@f$ with
   * the guess @f$x'@f$. 
   *
   * @f[
   *   v(x') = [s + s' (x' - x)]^2
   * @f]
   * 
   * @param guess Current guess @f$ x'@f$ 
   * @param r	  Observation
   * 
   * @return @f$ v(x')@f$ 
   */
  double TermVar(double guess, const Result& r) const
  {
    return std::pow(r.S() + r.Sprime() * (guess - r.fX),2);
  }
  /** 
   * Return the likely-hood function value at @f$ x'@f$:
   *
   * @f[
   *   L(x') = \left[(x'-x) / (s + s'(x'-x))\right]^2
   * @f] 
   * 
   * where 
   * @f[
   *   s = 2\sigma^+\sigma^-/(\sigma^++\sigma^-)
   * @f]
   * @f[
   *   s' = (\sigma^+-\sigma^-)/(\sigma^++\sigma^-)
   * @f]
   * 
   * @param guess @f$ x'@f$ 
   * @param x     @f$ x@f$  
   * @param el    @f$ \sigma^-@f$  
   * @param eh    @f$ \sigma^+@f$  
   * 
   * @return 
   */
  static double L(double guess, double x, double el, double eh)
  {
    double d  = (guess-x);
    double s  = 2 * el * eh / (el + eh);
    double sp = (eh - el) / (eh + eh);
    return std::pow(d / (s+sp*d), 2);
  }
  /** 
   * Wrap likely-hood function for ROOT 
   * 
   * @param xp Pointer to independent variables 
   * @param pp Pointer to parameters 
   * 
   * @return Likely-hood function evaluate at @a xp
   */
  static double WrapL(double* xp, double* pp)
  {
    return L(xp[0], pp[0], pp[1], pp[2]);
  }
};

/**
 * A combiner that uses a linear variance approximation 
 */
struct LinearVarianceCombiner : public Combiner
{
  /** 
   * Calculate the weight 
   * 
   * @f[
   *   w = (V + x V')^2 / (2 V + x V')
   * @f] 
   *
   * @param r Observation
   * 
   * @return @f$ W@f$ 
   */
  double W(const Result& r) const
  {
    double v  = r.V();
    double vp = r.Vprime();
    return std::pow(v + r.fX * vp, 2) / (2 * v + r.fX * vp);
  }
  /** 
   * Calculate the weight based on a guess of best @f$ x'@f$ 
   * 
   * @f[
   *   W(x') = v / [V + V' (x' - x)]^2
   * @f] 
   *
   * @param guess Current guess @f$ x'@f$ 
   * @param r  	  Observation
   * 
   * @return @f$ W(x')@f$ 
   */
  double StepW(double guess, const Result& r) const
  {
    double v = r.V();
    return v / std::pow(v+r.Vprime()*(guess - r.fX), 2);
  }
  /** 
   * Calculate the bias. 
   * 
   * @f[
   *   \delta(x') = 1/2 V' [(x'-x) / (V + V'(x' - x))]^2
   * @f] 
   * 
   * @param guess Current guess @f$ x'@f$ 
   * @param r	  Observation
   * 
   * @return @f$\delta(x')@f$ 
   */
  double StepOffset(double guess, const Result& r) const
  {
    double vp = r.Vprime();
    return 0.5 * vp * std::pow((guess-r.fX)/(r.V()+vp*(guess-r.fX)),2);
  }
  /** 
   * Calculate the contribution variance to the @f$\chi^2@f$ with
   * the guess @f$x'@f$. 
   *
   * @f[ 
   *   V(x') = V + V' (x' - x)
   * @f] 
   * 
   * @param guess Current guess @f$ x'@f$ 
   * @param r	  Observation
   * 
   * @return @f$ v(x')@f$ 
   */
  double TermVar(double guess, const Result& r) const
  {
    return r.V() + r.Vprime() * (guess - r.fX);
  }
  /** 
   * Return the likely-hood function value at @f$ x'@f$:
   *
   * @f[
   *   L(x') = (x'-x)^2 / (V + V'(x'-x))
   * @f] 
   * 
   * where 
   * @f[
   *   v = \sigma^\sigma^-
   * @f]
   * @f[
   *   v' = \sigma^+-\sigma^-
   * @f]
   * 
   * @param guess @f$ x'@f$ 
   * @param x     @f$ x@f$  
   * @param el    @f$ \sigma^-@f$  
   * @param eh    @f$ \sigma^+@f$  
   * 
   * @return 
   */
  static double L(double guess, double x, double el, double eh)
  {
    double d  = (guess-x);
    double v  = eh * el;
    double vp = eh - el;
    return std::pow(d,2) / (v+vp*d);
  }
  /** 
   * Wrap likely-hood function for ROOT 
   * 
   * @param xp Pointer to independent variables 
   * @param pp Pointer to parameters 
   * 
   * @return Likely-hood function evaluate at @a xp
   */
  static double WrapL(double* xp, double* pp)
  {
    return L(xp[0], pp[0], pp[1], pp[2]);
  }

};

#include <TF1.h>
#include <TList.h>
#include <TLine.h>
#include <TH1.h>

struct DrawResult
{
  static TF1* MakeF(const Combiner::Result& r,
		    Int_t                   j,
		    Bool_t                  isSigma)
  {
    TF1* f = new TF1(Form("f%02d", j),
		     isSigma ? LinearSigmaCombiner::WrapL 
		     : LinearVarianceCombiner::WrapL,
		     r.Low(), r.High(), 3);
    f->SetParNames("x", "#sigma^{-}", "#sigma^{+}");
    f->SetParameters(r.fX, r.fEl, r.fEh);
    f->SetLineStyle((j%3)+2);
    f->SetLineColor(kBlack);
    f->SetLineWidth(2);
    return f;
  }
  static TLine* MakeL(TF1* f)
  {
    Double_t m = f->GetParameter(0);
    TLine* l = new TLine(m-f->GetParameter(1), 1,
			 m+f->GetParameter(2), 1);
    // l->SetName(Form("l%s", f->GetName()));
    l->SetLineColor(f->GetLineColor());
    l->SetLineStyle(f->GetLineStyle());
    l->SetLineWidth(f->GetLineWidth());
    return l;
  }
    
  static void Draw(Combiner::const_iterator b,
		   Combiner::const_iterator e,
		   Combiner::Result         r,
		   Bool_t                   isSigma)
  {
    TList fs; fs.SetOwner(false);
    Int_t j = 0;
    for (Combiner::const_iterator i = b; i != e; i++) {
      TF1* f = MakeF(*i, j, isSigma);
      f->SetRange(r.Low(), r.High());
      fs.Add(f, j == 0 ? "" : "same");
      fs.Add(MakeL(f));
      j++;
    }
    TF1* fr = MakeF(r, j, isSigma);
    fr->SetLineColor(kRed+2);
    fr->SetLineStyle(1);
    fs.Add(fr);
    fs.Add(MakeL(fr));
    TIter next(&fs);
    TObject* o = 0;
    j = 0;
    while((o = next())) { o->Draw(j == 0 ? "" : "same"); j++; }
    // fs.Draw();
    static_cast<TF1*>(fs.First())->GetHistogram()
      ->GetYaxis()->SetRangeUser(-.1, 5.1);
  }
};

#include <iostream>
#include <algorithm>

void RunTest(Combiner::const_iterator b,
	     Combiner::const_iterator e,
	     Combiner::Final          s,
	     Combiner::Final          v)
{
  std::ostream_iterator<Combiner::Result> out(std::cout, "\n");
  std::copy(b, e, out);

  LinearSigmaCombiner sc;
  Combiner::Final sf = sc.Calculate(b, e);
  Combiner::Final sd(sf.fX-s.fX, sf.fEl-s.fEl,sf.fEh-s.fEh,
		     sf.fChi2-s.fChi2,0,0);
  std::cout << "Linear sigma:    " << sf  << "\n"
	    << " Expected:       " << s   << "\n"
	    << " Difference:     " << sd  << std::endl;
  DrawResult::Draw(b, e, sf, true);
  // return;
    
  LinearVarianceCombiner sv;
  Combiner::Final vf = sv.Calculate(b, e);
  Combiner::Final vd(vf.fX-v.fX, vf.fEl-v.fEl,vf.fEh-v.fEh,
		     vf.fChi2-v.fChi2,0,0);
  std::cout << "Linear variance: " << vf  << "\n"
	    << " Expected:       " << v   << "\n"
	    << " Difference:     " << vd  << std::endl;

  DrawResult::Draw(b, e, vf, false);
}

	     
void Test1()
{
  Combiner::List l;
  l.Add(4, 1.682, 2.346);
  l.Add(5, 1.912, 2.581);
  Combiner::Final s(4.49901, 1.33301, 1.66701, 0.115, 0, 11);
  Combiner::Final v(4.50001, 1.33701, 1.66901, 0.113, 0, 11);

  RunTest(l.begin(), l.end(), s, v);
}

void Test2()
{
  Combiner::List l;
  l.Add(6.32064,	0.567382,	0.379042);
  l.Add(6.15549,	0.159504,	0.170811);

  Combiner::Final s(6.2, 0.2, 0.2, 1, 0, 11);

  RunTest(l.begin(), l.end(), s, s);
}
void Test3()
{
  Combiner::List l;
  l.Add(6.20449,	0.451232,	0.495192);
  l.Add(6.16188,	0.165785,	0.176712);

  Combiner::Final s(6.18, 0.3, 0.3, 1, 0, 11);

  RunTest(l.begin(), l.end(), s, s);
  
}

#if 0
int
main()
{
  Test1();

  return 0;
}

#endif
