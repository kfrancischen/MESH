/* errcodes.h
 * Part of quadpack++.
 */

#ifndef GSL_ERRCODES_H
#define GSL_ERRCODES_H

enum { 
	GSL_SUCCESS  = 0, 
	GSL_FAILURE  = -1,
	GSL_CONTINUE = -2,  /* iteration has not converged */
	GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
	GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
	GSL_EFAULT   = 3,   /* invalid pointer */
	GSL_EINVAL   = 4,   /* invalid argument supplied by user */
	GSL_EFAILED  = 5,   /* generic failure */
	GSL_EFACTOR  = 6,   /* factorization failed */
	GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
	GSL_ENOMEM   = 8,   /* malloc failed */
	GSL_EBADFUNC = 9,   /* problem with user-supplied function */
	GSL_ERUNAWAY = 10,  /* iterative process is out of control */
	GSL_EMAXITER = 11,  /* exceeded max number of iterations */
	GSL_EZERODIV = 12,  /* tried to divide by zero */
	GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
	GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
	GSL_EUNDRFLW = 15,  /* underflow */
	GSL_EOVRFLW  = 16,  /* overflow  */
	GSL_ELOSS    = 17,  /* loss of accuracy */
	GSL_EROUND   = 18,  /* failed because of roundoff error */
	GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
	GSL_ENOTSQR  = 20,  /* matrix not square */
	GSL_ESING    = 21,  /* apparent singularity detected */
	GSL_EDIVERGE = 22,  /* integral or series is divergent */
	GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
	GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
	GSL_ECACHE   = 25,  /* cache limit exceeded */
	GSL_ETABLE   = 26,  /* table limit exceeded */
	GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
	GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
	GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
	GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
	GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
	GSL_EOF      = 32   /* end of file */
} ;



/** \file errcodes.h
 \brief Error handling
 
 Adaptive quadrature routine(s) in quadpack++ throw an exception when a 
 numerical difficulty is encountered. For continuity with the original QUADPACK
 routnes and their GSL implemenation, an integer error code is also returned.
 Usage might be along the lines
 \code
 int errcode;
 try {
    errcode = Work.qag(f, a, b, epsabs, epsrel, limit, &result, &abserr);
 }
 catch (const char* reason) {
    std::cerr << reason << std::endl;
    ...
 }
 \endcode
 The table of error codes from the GSL is used here without modification.
 */

/* function.h
 * Part of the quadpack++.
 */


//! C++ version of <tt>gsl_function</tt> passed to quadrature routines.
template<typename Real>
class FtnBase {
public:
	virtual Real operator() (Real x) =0;
};

/** \brief C++ version of <tt>gsl_function</tt> with constructors. 

 User-defined functions for numerical routines may require additional parameters
 depending on the application. However function-pointers, when passed as 
 arguments to C-routines, depend on their "signature" type. To accomodate a
 universal signature, the GSL uses the struct <tt>gsl_function</tt> with the
 member <tt>void* params</tt> for all possible parameter/structure types. 
 
 The approach in quadpack++ is analgous except that the parameter type is 
 templated. The universality is accomodated in that only the base type,
 \ref FtnBase, is passed as an argument to routines. 
 
 <b>Note:</b> the parameter \ref params_ is must be a pointer instead of a
 refererence since this object is declared outside of the scope.
 */
template<typename Real, class param_t>
class Function : public FtnBase<Real> {
public:
	//! Signature defining function \ref FtnBase::operator()
	typedef Real defn_t(Real, param_t*);
	
	//! Reference to definition
	defn_t& function_;
	//! Pointer to current parameters.
	param_t* params_;

	//! Overload the \ref FtnBase::operator() with \ref function_.
	virtual Real operator() (Real x) 
	{
		return function_(x, params_); 
	}
	
	//! Constructor without parameters.
	Function(defn_t& function) : function_(function), params_(0) {}
	//! General constructor.
	Function(defn_t& function, param_t* params) : function_(function), params_(params) {}
	~Function() {}
};

/* machar.hpp
 *
 * Class that provides key floating-point parameters based on the
 * original MACHAR routine.
 *
 * Copyright (C) 2010 Jerry Gagelman
 * Copyright (C) 2006 John Burkardt
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \brief Key prarameters for the \a Real floating-point type.

 The class constructor adapts source code from
 <http://people.sc.fsu.edu/~jburkardt/c_src/machar/machar.html>.
 The original algorithm is described in
 - William Cody, <em>Algorithm 665: MACHAR, a subroutine to dynamically 
 determine machine parameters,</em> ACM Transactions on Mathematical Software,
 Volume 14, (1988) pp. 303-311.
 */
template <class Real>
class Machar
{
protected:
	Real eps_;	///< Machine epsilon
	Real xmin_;	///< Underflow threshold
	Real xmax_;	///< Overflow threshold
	
	Real inline abs(Real x) { return (x < Real(0)) ? -(x) : x; }
	Real inline max(Real a, Real b) { return a > b ? a : b; }
	Real inline min(Real a, Real b) { return a < b ? a : b; }
	
public:
	Machar();
	~Machar() {}
	
	Real inline get_eps() { return eps_; }
	Real inline get_xmin() { return xmin_; }
	Real inline get_xmax() { return xmax_; }
};

/* Naming conventions from the original FORTRAN program have been maintained.
 * A brief description is as follows: 
 *
 C  IBETA   - the radix for the floating-point representation.
 C
 C  IT      - the number of base IBETA digits in the floating-point
 C            significand.
 C
 C  IRND    - 0 if floating-point addition chops,
 C            1 if floating-point addition rounds, but not in the
 C              IEEE style,
 C            2 if floating-point addition rounds in the IEEE style,
 C            3 if floating-point addition chops, and there is
 C              partial underflow,
 C            4 if floating-point addition rounds, but not in the
 C              IEEE style, and there is partial underflow,
 C            5 if floating-point addition rounds in the IEEE style,
 C              and there is partial underflow.
 C
 C  NGRD    - the number of guard digits for multiplication with
 C            truncating arithmetic.  It is
 C            0 if floating-point arithmetic rounds, or if it
 C              truncates and only  IT  base  IBETA digits
 C              participate in the post-normalization shift of the
 C              floating-point significand in multiplication;
 C            1 if floating-point arithmetic truncates and more
 C              than  IT  base  IBETA  digits participate in the
 C              post-normalization shift of the floating-point
 C              significand in multiplication.
 C
 C  MACHEP  - the largest negative integer such that
 C            1.0 + FLOAT(IBETA)**MACHEP != 1.0, except that
 C            MACHEP is bounded below by  -(IT+3).
 C
 C  NEGEP   - the largest negative integer such that
 C            1.0 - FLOAT(IBETA)**NEGEP != 1.0, except that
 C            NEGEP is bounded below by  -(IT+3).
 C
 C  IEXP    - the number of bits (decimal places if IBETA = 10)
 C            reserved for the representation of the exponent
 C            (including the bias or sign) of a floating-point
 C            number.
 C
 C  MINEXP  - the largest in magnitude negative integer such that
 C            FLOAT(IBETA)**MINEXP is positive and normalized.
 C
 C  MAXEXP  - the smallest positive power of  BETA  that overflows.
 C
 C  EPS     - the smallest positive floating-point number such
 C            that  1.0+EPS != 1.0. In particular, if either
 C            IBETA = 2  or  IRND = 0, EPS = FLOAT(IBETA)**MACHEP.
 C            Otherwise,  EPS = (FLOAT(IBETA)**MACHEP)/2.
 C
 C  EPSNEG  - A small positive floating-point number such that
 C            1.0-EPSNEG .NE. 1.0. In particular, if IBETA = 2
 C            or  IRND = 0, EPSNEG = FLOAT(IBETA)**NEGEP.
 C            Otherwise,  EPSNEG = (IBETA**NEGEP)/2.  Because
 C            NEGEP is bounded below by -(IT+3), EPSNEG may not
 C            be the smallest number that can alter 1.0 by
 C            subtraction.
 C
 C  XMIN    - the smallest non-vanishing normalized floating-point
 C            power of the radix, i.e.,  XMIN = FLOAT(IBETA)**MINEXP.
 C
 C  XMAX    - the largest finite floating-point number.  In
 C            particular  XMAX = (1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP
 C            Note - on some machines  XMAX  will be only the
 C            second, or perhaps third, largest number, being
 C            too small by 1 or 2 units in the last digit of
 C            the significand. 
 */
template <class Real>
Machar<Real>::Machar()
{
	long int 
	ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp;
	Real epsneg;
	
	Real one = (Real)1.0;
	Real two = one + one;
	Real zero = (Real)0.0;
	
	/* Determine IBETA and BETA. */
	Real a = one;
	Real temp, temp1;
	do {
		a = a + a;
		temp = a + one;
		temp1 = temp - a;
	} while ( temp1 - one == zero );
	
	Real b = one;
	do {
		b = b + b;
		temp = a + b;
	} while ( temp - a == zero );
	
	Real beta = temp - a;
	ibeta = (long int)beta;
	Real betah = beta / two;
	Real betain = one / beta;
	
	/* Determine IRND, IT. */
	it = 0;
	b = one;
	do {
		it = it + 1;
		b = b * beta;
		temp = b + one;
		temp1 = temp - b;
	} while (temp1 - one == zero);
	
	irnd = 0;
	temp = a + betah;
	temp1 = temp - a;
	if ( temp1 != zero )
		irnd = 1;
	
	Real tempa = a + beta;
	temp = tempa + betah;
	
	if (irnd == 0 && (temp - tempa != zero ))
		irnd = 2;
	
	/* Determine NEGEP, EPSNEG. */
	negep = it + 3;
	a = one;
	for (long int i = 1; i <= negep; ++i )
		a = a * betain;
	
	tempa = a;
	temp = one - a;
	while (temp - one == zero ) {
		a = a * beta;
		negep = negep - 1;
		temp = one - a;
	}
	
	negep = -negep;
	epsneg = a;
	
	if ( ibeta != 2 && irnd != 0 ) {
		a = ( a * (one + a) ) / two;
		temp = one - a;
		if ( temp - one != zero )
			epsneg = a;
	}
	
	/* Determine MACHEP, EPS. */
	machep = -it - 3;
	a = tempa;
	temp = one + a;
	while ( temp - one == zero ) {
		a = a * beta;
		machep = machep + 1;
		temp = one + a;
	}
	
	eps_ = a;
	
	/* Determine NGRD. */
	ngrd = 0;
	temp = one + eps_;
	if ( irnd == 0 && ( temp * one - one ) != zero )
		ngrd = 1;
	
	/* Determine IEXP, MINEXP and XMIN.
	 * Loop to determine largest I such that (1/BETA) ** (2**(I))
	 * does not underflow.  Exit from loop is signaled by an underflow.
	 */
	long int i = 0, k = 1, nxres = 0, iz = 0, mx = 0;
	Real y, z = betain, t = one + eps_;
		
	for ( ; ; ) {
		y = z;
		z = y * y;
		/* Check for underflow */
		a = z * one;
		temp = z * t;
		if ( ( a + a == zero ) || z > y ) {
			break;
		}
		temp1 = temp * betain;
		if ( temp1 * beta == z ) {
			iexp = i + 1;
			mx = k + k;
			break;
		}
		i = i + 1;
		k = k + k;
	}

	if ( ibeta == 10 ) { 
		/* For decimal machines only */
		iexp = 2;
		iz = ibeta;
		while ( k >= iz ) {
			iz = iz * ibeta;
			iexp = iexp + 1;
		}
		mx = iz + iz - 1;
	}
		
	/* Loop to determine MINEXP, XMIN.
	 * Exit from loop is signaled by an underflow.
	 */
	for ( ; ; )
	{
		xmin_ = y;
		y = y * betain;
		a = y * one;
		temp = y * t;
		if ( ( a + a == zero ) || y >=  xmin_ )
			break;
		k = k + 1;
		temp1 = temp * betain;
		if ( temp1 * beta  == y ) {
			nxres = 3;
			xmin_ = y;
			break;
		}
	}
	
	minexp = -k;
	
	/* Determine MAXEXP, XMAX. */
	if ( mx <= ( k + k - 3 ) && ibeta != 10 ) {
		mx = mx + mx;
		iexp = iexp + 1;
	}
	
	maxexp = mx + minexp;
	
	/* Adjust IRND to reflect partial underflow. */
	irnd = irnd + nxres;
	
	/* Adjust for IEEE-style machines. */
	if ( ( irnd >= 2 ) || ( irnd == 5 ) )
		maxexp = maxexp - 2;
	
	/* Adjust for non-IEEE machines with partial underflow. */
	if ( ( irnd == 3 ) || ( irnd == 4 ) ) {
		maxexp = maxexp - it;
	}
	
	/* Adjust for machines with implicit leading bit in binary
	 * significand and machines with radix point at extreme
	 * right of significand.
	 */
	if ( ( ibeta == 2 ) && ( maxexp + minexp == 0 ) )
		maxexp = maxexp - 1;
	
	if ( maxexp + minexp > 20 ) 
		maxexp = maxexp - 1;
	
	if ( a != y )
		maxexp = maxexp - 2;
	
	xmax_ = one - epsneg;
	if ( xmax_ * one != xmax_ )
		xmax_ = one - beta * epsneg;
	
	xmax_ /= ( beta * beta * beta * xmin_ );
	i = maxexp + minexp + 3;
	
	if ( i > 0 ) {
		for (long int j = 1; j <= i; j++ ) {
			if ( ibeta == 2 )
				xmax_ = xmax_ + xmax_;
			else
				xmax_ = xmax_ * beta;
		}
	}
}
/* gauss-kronrod.hpp
 *
 * Gauss-Kronrod quadrature using templated floating-point type;
 * based on "qk*.c" codes from <http://www.gnu.org/software/gsl/>.
 *
 * Copyright (C) 2010, 2011 Jerry Gagelman
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \brief Gauss-Kronrod quadrature rule(s) plus error-estimate data for 
 adaptive quadrature routines.
 
 The Gauss-Kronrod abscissae consist of 2m+1 points \f$x_1 < \cdots < x_{2m+1}\f$
 in the interval (-1, 1) used for a low-order and a high-order quadrature rule:
 \f[ Q_m^G f = \sum_{k=1}^m a_k f(x_{2k}), \qquad
 Q_m^{GK} f = \sum_{k=1}^{2m+1} b_k f(x_k). 
 \f]
 
 The weights and abscissae are available through member functions, however 
 they are stored according to compact QUADPACK convention. Due to symmetry,
 the positive abscissae \f$x_{2m+1}\f$, \f$x_{2m}\f$,..., \f$x_m\f$ are returned 
 as values <em>xgk(0)</em>, <em>xgk(1)</em>,..., <em>xgk(m+1)</em> respectively 
 of the member function \ref xgk(). Note the reverse order. The corresponding
 weights \f$b_{2m+1}\f$, \f$b_{2m}\f$,..., \f$b_m\f$ are returned by the 
 respective values of \ref wgk(). The weights \f$a_m\f$, \f$a_{m-1}\f$,...,
 corresponding to the even-indexed \f$x_{2m}\f$, \f$x_{2m-2}\f$, ...., are
 returned by the values of \ref wg() in their reverse order.
 
 <h3>Computational details</h3>
 
 The even-indexed abscissae \f$x_2\f$, ..., \f$x_{2m}\f$ are the zeros of the m-th
 Legendre polynomial \f$P_m\f$. The odd indexed points are zeros of a polynomial
 that is represented as a Chebyshev sum,
 \f[ E_{m+1} = T_{m+1} + c_{m-1}T_{m-1} + c_{m-3} T_{m-3} + \cdots,
 \f]
 whose coefficients are defined by explicit formulae in
 - Giovanni Monegato, <em>Some remarks on the construction of extended Gaussian 
 quadrature rules,</em> Math. Comp., Vol. 32 (1978) pp. 247-252. 
 <a href="http://www.jstor.org/stable/2006272">[jstor]</a>.
 
 The zeros of both of these polynomials are computed by Newton's method. Upper 
 bounds for their round-off errors, as functions of machine epsilon, are 
 incorporated in the stopping criteria for for the root finders.

 The weights \f$a_1\f$, ..., \f$a_m\f$ are Gauss-Legendre weights. 
 The \f$b_k\f$ are given by the formulae
 \f[ b_{2k} = a_k + \frac{2 p_m}{(2m+1)t_{m+1} P_m'(x_{2k}) E_{m+1}(x_{2k})},
 \qquad k = 1,\ldots,m, \f] 
 and
 \f[ b_{2k+1} = \frac{2 p_m}{(2m+1) t_{m+1} P_m(x_{2k+1}) E_{m+1}'(x_{2k+1})},
 \qquad k = 0, \ldots, m, \f]
 where \f$p_m\f$ and \f$t_{m+1}\f$ are the leading coefficients of the 
 polynomials \f$P_m\f$ and \f$T_{m+1}\f$ respectively. These are from
 - Giovanni Monegato, <em>A note on extended Gaussian quadrature rules,</em> 
 Math. Comp., Vol. 30 (1976) pp. 812-817. 
 <a href="http://www.jstor.org/stable/2005400">[jstor]</a>. 
 */
template <class Real>
class GaussKronrod : public Machar<Real>
{
private:
	size_t m_;  // Gauss-Legendre degree
	size_t n_;  // size of Gauss-Kronrod arrays
	Real*  xgk_;  // Gauss-Kronrod abscissae
	Real*  wg_;   // Gauss-Legendre weights
	Real*  wgk_;  // Gauss-Kronrod weights
	
	Real *coefs;  // Chebyshev coefficients 
	Real *zeros;  // zeros of Legendre polynomial
	Real *fv1, *fv2;  // scratch space for error estimator
	
	Real rescale_error(Real err, const Real result_abs, 
							 const Real result_asc);
	
	void legendre_zeros();
	void chebyshev_coefs();
	void gauss_kronrod_abscissae();
	void gauss_kronrod_weights();
	
	Real legendre_err(int deg, Real x, Real& err);
	Real legendre_deriv(int deg, Real x);	
	Real chebyshev_series(Real x, Real& err);
	Real chebyshev_series_deriv(Real x);
	
public:
	//! Initializes class for (2m+1)-point Gauss-Kronrod quadrature.
	GaussKronrod(size_t m = 10);
	~GaussKronrod();
	
	//! Approximates \f$\int_a^b f\,dx\f$ using the Gauss-Kronrod rule.
	void qk(FtnBase<Real>& f, Real a, Real b,
			  Real& result, Real& abserr, Real& resabs, Real& resasc);
	
	//! Size of arrays of Gauss-Kronrod abscissae and weights.
	size_t size() { return n_; };
	
	//! Array of Gauss-Kronrod abscissae in (0, 1); QUADPACK convention.
	Real xgk(int k) 
	{ 
		return (0 <= k && k < n_) ? xgk_[k] : Real(0); 
	}
	
	//! Array of corresponding Gauss-Kronrod weights; QUADPACK convention. 
	Real wgk(int k)
	{ 
		return (0 <= k && k < n_) ? wgk_[k] : Real(0); 
	}

	//! Gauss-Legendre weights for odd indexed abscissae; QUADPACK convention.
	Real wg(int k) 
	{ 
		return (0 <= k && k < n_/2) ? wg_[k] : Real(0); 
	}
};

template <class Real>
GaussKronrod<Real>::GaussKronrod(size_t m) : Machar<Real>()
{
	m_ = m;
	n_ = m_ + 1;
	xgk_ = new Real[n_];
	wg_  = new Real[n_ / 2];
	wgk_ = new Real[n_];	
	coefs = new Real[n_ + 1];
	zeros = new Real[m_ + 2];
	fv1 = new Real[n_];
	fv2 = new Real[n_];
	
	legendre_zeros();
	chebyshev_coefs();
	gauss_kronrod_abscissae();
	gauss_kronrod_weights();
}

template <class Real>
GaussKronrod<Real>::~GaussKronrod()
{
	delete[] xgk_;
	delete[] wgk_;
	delete[] wg_;
	delete[] coefs;
	delete[] zeros;
	delete[] fv1;
	delete[] fv2;
}

/** 
 Computes the zeros of the Legendre polynomial \f$P_m\f$. Upon exit, these
 are stored consecutively as elements of the array \c zeros[] <b>indexed by 
 1,..., m.</b>
 */
template <class Real>
void GaussKronrod<Real>::legendre_zeros()
{
	Real* temp = new Real[m_+1];
	zeros[0] = Real(-1);
	zeros[1] = Real(1);
	Real delta, epsilon;
	
	for (int k = 1; k <= m_; ++k) {
		// loop to locate zeros of P_k interlacing z_0,...,z_k
		for (int j = 0; j < k; ++j) {
			// Newton's method for P_k :
			// initialize solver at midpoint of (z_j, z_{j+1})
			delta = 1;
			Real x_j = (zeros[j] + zeros[j+1]) / 2;
			Real P_k = legendre_err(k, x_j, epsilon);
			while (this->abs(P_k) > epsilon &&
					 this->abs(delta) > this->eps_) 
			{
				delta = P_k / legendre_deriv(k, x_j);
				x_j -= delta;
				P_k = legendre_err(k, x_j, epsilon);
			}
			temp[j] = x_j;
		}
		
		// copy roots tmp_0,...,tmp_{k-1} to z_1,...z_k:
		zeros[k+1] = zeros[k];
		for (int j = 0; j < k; ++j)
			zeros[j+1] = temp[j];
		
	}
	delete[] temp;
}

/** 
 Computes coefficients of polynomial \f$E_{m+1}\f$ in the array \c coefs[].
 */
template <class Real>
void GaussKronrod<Real>::chebyshev_coefs()
{
	size_t ell = (m_ + 1)/2;
	Real* alpha = new Real[ell+1];
	Real* f = new Real[ell+1];
	
	/* Care must be exercised in initalizing the constants in the definitions.
	 * Compilers interpret expressions like "(2*k + 1.0)/(k + 2.0)" as floating
	 * point precision, before casting to Real.
	 */
	f[1] = Real(m_+1)/Real(2*m_ + 3);
	alpha[0] = Real(1); // coefficient of T_{m+1}
	alpha[1] = -f[1];
	
	for (int k = 2; k <= ell; ++k) {
		f[k] = f[k-1] * (2*k - 1) * (m_ + k) / (k * (2*m_ + 2*k + 1));
		alpha[k] = -f[k];
		for (int i = 1; i < k; ++i)
			alpha[k] -= f[i] * alpha[k-i];
	}
	
	for (int k = 0; k <= ell; ++k) {
		coefs[m_ + 1 - 2*k] = alpha[k];
		if (m_  >= 2*k)
			coefs[m_ - 2*k] = Real(0);
	}
	
	delete[] alpha;
	delete[] f;
}

/** 
 Computes Gauss-Legendre weights \c wg_[] and Gauss-Kronrod weights \c wgk_[]. 
 */
template <class Real>
void GaussKronrod<Real>::gauss_kronrod_weights()
{
	Real err;
	/* Gauss-Legendre weights:
	 */
	for (int k = 0; k < n_ / 2; ++k) 
	{
		Real x = xgk_[2*k + 1];
		wg_[k] = (Real(-2) / 
					 ((m_ + 1) * legendre_deriv(m_, x) * legendre_err(m_+1, x, err)));
	}
	
	/* The ratio of leading coefficients of P_n and T_{n+1} is computed
	 * from the recursive formulae for the respective polynomials.
	 */
	Real F_m = Real(2) / Real(2*m_ + 1);
	for (int k = 1; k <= m_; ++k)
		F_m *= (Real(2*k) / Real(2*k - 1));
	
	/* Gauss-Kronrod weights:  
	 */
	for (size_t k = 0; k < n_; ++k) 
	{
		Real x = xgk_[k];
		if (k % 2 == 0) 
		{
			wgk_[k] = F_m / (legendre_err(m_, x, err) * chebyshev_series_deriv(x));
		}
		else 
		{
			wgk_[k] = (wg_[k/2] + 
						  F_m / (legendre_deriv(m_, x) * chebyshev_series(x, err)));
		}
	}	
}

/**
 Computes the zeros of the polynomial \f$E_{m+1}\f$, using the fact that these
 interlace the zeros of the Legendre polynomial \f$P_m\f$, which are stored in 
 the array \c zeros[]. Appropriate elements of \c zeros[] are then copied
 into \c xgk_[].
 */
template <class Real>
void GaussKronrod<Real>::gauss_kronrod_abscissae()
{
	Real delta, epsilon;
	
	for (int k = 0; k < n_ / 2; ++k) 
	{
		delta = 1;
		// Newton's method for E_{n+1} :
		Real x_k = (zeros[m_-k] + zeros[m_+1-k])/Real(2);
		Real E = chebyshev_series(x_k, epsilon);
		while (this->abs(E) > epsilon &&
				 this->abs(delta) > this->eps_) 
		{
			delta = E / chebyshev_series_deriv(x_k);
			x_k -= delta;
			E = chebyshev_series(x_k, epsilon);
		}
		xgk_[2*k] = x_k;
		// copy adjacent Legendre-zero into the array:
		if (2*k+1 < n_)
			xgk_[2*k+1] = zeros[m_-k];
	}
}

/**
 Recursive definition of the Legendre polynomials
 \f[ (k+1) P_{k+1}(x) = (2k+1) x P_k(x) - k P_{k-1}(x), 
 \f]
 and estimate of the rounding error,
 \f[ E_{k+1} = \frac{(2k+1)|x|E_k + kE_{k-1}}{2(k+1)},
 \f]
 are from the routine <tt>gsl_sf_legendre_Pl_e</tt> distributed with GSL.
 */
template <class Real>
Real GaussKronrod<Real>::legendre_err(int n, Real x, Real& err)
{
	if (n == 0) {
		err = Real(0);
		return Real(1);
	}
	else if (n == 1) {
		err = Real(0);
		return x;
	}
	
	Real P0 = Real(1), P1 = x, P2;
	Real E0 = this->eps_; 
	Real E1 = this->abs(x) * this->eps_; 
	for (int k = 1; k < n; ++k) 
	{
		P2 = ((2*k + 1) * x * P1 - k * P0) / (k + 1);
		err = ((2*k + 1) * this->abs(x) * E1 + k * E0) / (2*(k + 1));
		P0 = P1; P1 = P2;
		E0 = E1; E1 = err;
	}
	return P2;	
}

/**
 Three-term recursion identity for the Legendre derivatives:
 \f[ P_{k+1}'(x) = (2k+1) P_k(x) + P_{k-1}'(x).
 \f]
 */
template <class Real>
Real GaussKronrod<Real>::legendre_deriv(int n, Real x)
{
	if (n == 0)   
		return Real(0);
	else if (n == 1)
		return Real(1);
	
	Real P0 = Real(1), P1 = x, P2;
	Real dP0 = Real(0), dP1 = Real(1), dP2;
	for (int k = 1; k < n; ++k) 
	{
		P2 = ((2*k + 1) * x * P1 - k * P0) / (k + Real(1));
		dP2 = (2*k + 1) * P1 + dP0;
		P0 = P1; P1 = P2;
		dP0 = dP1; dP1 = dP2;
	}
	return dP2;
}

/**
 Evaluation of the polynomial \f$E_{m+1}\f$ is using the Clenshaw method and
 (truncation) error estimate is taken from the routine 
 <tt>gsl_cheb_eval_err</tt> distributed with GSL.
*/
template <class Real>
Real GaussKronrod<Real>::chebyshev_series(Real x, Real& err)
{
	Real d1(0), d2(0);
	Real absc = this->abs(coefs[0]); // final term for truncation error
	Real y2 = 2 * x; // linear term for Clenshaw recursion
	
	for (int k = n_; k >= 1; --k) {
		Real temp = d1;
		d1 = y2 * d1 - d2 + coefs[k];
      d2 = temp;
		absc += this->abs(coefs[k]);
	}
	
	err = absc * this->eps_;
	return x * d1 - d2 + coefs[0]/2;
}

/** 
 Derivatives of Chebyshev polynomials satisfy the identity \f$T_n' = nU_{n-1}\f$,
 where the \f$U_k\f$ are Chebyshev polynomials of the second kind. The derivative
 of the polynomial \f$E_{m+1}\f$ is implemented using the Clenshaw algorithm for
 the latter polynomials.
 */
template <class Real>
Real
GaussKronrod<Real>::chebyshev_series_deriv(Real x)
{
	Real d1(0), d2(0);
	Real y2 = 2 * x; // linear term for Clenshaw recursion

	for (int k = n_; k >= 2; --k) {
		Real temp = d1;
		d1 = y2 * d1 - d2 + k * coefs[k];
      d2 = temp;
	}

	return y2 * d1 - d2 + coefs[1];
}

/**
 QUADPACK's nonlinear formula for the absolute error.
 */
template <class Real>
Real GaussKronrod<Real>::rescale_error (Real err, const Real result_abs, 
													 const Real result_asc)
{
	err = this->abs(err);
	
	if (result_asc != Real(0) && err != Real(0))
	{
		// cast 1.5 as Real number
		Real exponent = Real(3)/Real(2);
		Real scale = pow((200 * err / result_asc), exponent);
		
		if (scale < Real(1))
		{
			err = result_asc * scale ;
		}
		else 
		{
			err = result_asc ;
		}
	}
	
	if (result_abs > this->xmin_ / (50 * this->eps_))
	{
      Real min_err = 50 * this->eps_ * result_abs ;
		
      if (min_err > err) 
		{
			err = min_err ;
		}
	}
	
	return err ;
}

template <class Real>
void GaussKronrod<Real>::qk(FtnBase<Real>& f, Real a, Real b,
									 Real& result, Real& abserr, 
									 Real& resabs, Real& resasc)
{
	const Real center = (a + b) / 2;
	const Real half_length = (b - a) / 2;
	const Real abs_half_length = this->abs(half_length);
	// const Real f_center = f.function(center, f.params);
	const Real f_center = f(center);
	
	Real result_gauss = Real(0);
	Real result_kronrod = f_center * wgk_[n_ - 1];
	Real result_abs = this->abs(result_kronrod);
	Real result_asc = Real(0);
	Real mean = Real(0), err = Real(0);
	
	int j;
	
	 if (n_ % 2 == 0)
	 {
		 result_gauss = f_center * wg_[n_/2 - 1];
	 }
	
	for (j = 0; j < (n_ - 1) / 2; j++)
	{
      int jtw = j * 2 + 1;        /* j=1,2,3 jtw=2,4,6 */
      Real abscissa = half_length * xgk_[jtw];
//      Real fval1 = f.function( center - abscissa , f.params);
//      Real fval2 = f.function( center + abscissa , f.params);
		Real fval1 = f(center - abscissa);
		Real fval2 = f(center + abscissa);
      Real fsum = fval1 + fval2;
      fv1[jtw] = fval1;
      fv2[jtw] = fval2;
      result_gauss += wg_[j] * fsum;
      result_kronrod += wgk_[jtw] * fsum;
      result_abs += wgk_[jtw] * (this->abs(fval1) + this->abs(fval2));
	}
	
	for (j = 0; j < n_ / 2; j++)
	{
      int jtwm1 = j * 2;
      Real abscissa = half_length * xgk_[jtwm1];
//      Real fval1 = f.function( center - abscissa , f.params);
//      Real fval2 = f.function( center + abscissa , f.params);
 		Real fval1 = f(center - abscissa);
		Real fval2 = f(center + abscissa);
		fv1[jtwm1] = fval1;
      fv2[jtwm1] = fval2;
      result_kronrod += wgk_[jtwm1] * (fval1 + fval2);
      result_abs += wgk_[jtwm1] * (this->abs(fval1) + this->abs(fval2));
	};
	
	mean = result_kronrod / 2;
	
	result_asc = wgk_[n_ - 1] * this->abs(f_center - mean);
	
	for (j = 0; j < n_ - 1; j++)
	{
      result_asc += wgk_[j] * (this->abs(fv1[j] - mean) + 
										 this->abs(fv2[j] - mean));
	}
	
	/* scale by the width of the integration region */
	
	err = (result_kronrod - result_gauss) * half_length;
	
	result_kronrod *= half_length;
	result_abs *= abs_half_length;
	result_asc *= abs_half_length;
	
	result = result_kronrod;
	resabs = result_abs;
	resasc = result_asc;
	abserr = rescale_error (err, result_abs, result_asc);
}

/* workspace.hpp
 *
 * Adaptive quadrature for a templated floating-point type based on the
 * routine QAG implemented in <http://www.gnu.org/software/gsl/>.
 *
 * Copyright (C) 2010 Jerry Gagelman
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \brief The <tt>gsl_integration_workspace</tt> structure with member
 function \ref qag() for adaptive quadrature. 
 */
template <class Real>
class Workspace : public GaussKronrod<Real>
{	
public:
	Workspace();
	//! Initialize workspace for \a limit refinement intervals.
	Workspace(size_t limit);
	//! Initialize for \a limit refinement intervals and (2m+1)-point quadrature.
	Workspace(size_t limit, size_t m);
	~Workspace();
	
	/** \brief Implementation of <tt>gsl_integration_qag</tt> using class data.
	 
	 The arguments correspond to the <a href=
	 "http://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html">
	 original GSL</a> function except that <em>limit</em> and <em>key</em> are
	 use the workspace size and Gauss-Kronrod rule initializing the class. 
	 The stopping criterion is
	 <center>
	 |\a result \f$-\int f\,dx\f$| \f$\leq\f$
	 max{\a epsabs, \a epsrel\f$|\int f\,dx|\f$},
	 </center>
	 and \a abserr is the internal estimate of the right hand side. 
	 */
	int qag(FtnBase<Real>& f, Real a, Real b,
			  Real epsabs, Real epsrel, Real& result, Real& abserr);
	
private:
	void	allocate(size_t limit);
	
	/* data from gsl_integration_workspace */
	size_t limit;
	size_t size;
	size_t nrmax;
	size_t i_work;
	size_t maximum_level;
	Real   *alist;
	Real   *blist;
	Real   *rlist;
	Real   *elist;
	size_t *order;
	size_t *level;
	
	/* auxillary functions for adaptive quadrature */
	void  append_interval(Real a1, Real b1, Real area1, Real error1);
	void  initialise(Real a, Real b);
	void  set_initial_result(Real result, Real error);
	void  qpsrt(); 
	void  sort_results();
	void  retrieve(Real& a, Real& b, Real& r, Real& e);
	int   subinterval_too_small (Real a1, Real a2, Real b2);
	Real  sum_results();
	void  update(Real a1, Real b1, Real area1, Real error1,
					 Real a2, Real b2, Real area2, Real error2);
	
};

template <class Real>
void
Workspace<Real>::allocate(size_t n)
{
	alist = new Real[n];
	blist = new Real[n];
	rlist = new Real[n];
	elist = new Real[n];
	order = new size_t[n];
	level = new size_t[n];
	
	size = 0;
	limit = n;
	maximum_level = 0;	
}

template <class Real>
Workspace<Real>::Workspace() : GaussKronrod<Real>()
{
	allocate(1);
}

template <class Real>
Workspace<Real>::Workspace(size_t n) : GaussKronrod<Real>()
{
	if (n == 0) n = 1;	// ensure that workspace has positive size
	allocate(n);
}

template <class Real>
Workspace<Real>::Workspace(size_t n, size_t m) : GaussKronrod<Real>(m)
{
	if (n == 0) n = 1;	// ensure that workspace has positive size
	allocate(n);
}

template <class Real>
Workspace<Real>::~Workspace()
{
	delete[] alist;
	delete[] blist;
	delete[] rlist;
	delete[] elist;
	delete[] order;
	delete[] level;
}

template <class Real>
void
Workspace<Real>::append_interval (Real a1, Real b1, Real area1, Real error1)
{
	alist[size] = a1;
	blist[size] = b1;
	rlist[size] = area1;
	elist[size] = error1;
	order[size] = size;
	level[size] = 0;
	
	size++;
}	

template <class Real>
void 
Workspace<Real>::initialise (Real a, Real b)
{
	size = 0;
	nrmax = 0;
	i_work = 0;
	alist[0] = a;
	blist[0] = b;
	rlist[0] = Real(0);
	elist[0] = Real(0);
	order[0] = 0;
	level[0] = 0;
	
	maximum_level = 0;
}

template <class Real>
void 
Workspace<Real>::sort_results ()
{
	size_t i;
	
	for (i = 0; i < size; i++)
	{
      size_t i1 = order[i];
      Real e1 = elist[i1];
      size_t i_max = i1;
      size_t j;
		
      for (j = i + 1; j < size; j++)
		{
			size_t i2 = order[j];
			Real e2 = elist[i2];
			
			if (e2 >= e1)
			{
				i_max = i2;
				e1 = e2;
			}
		}
		
      if (i_max != i1)
		{
			order[i] = order[i_max];
			order[i_max] = i1;
		}
	}
	
	i_work = order[0] ;
}

template <class Real>
void
Workspace<Real>::qpsrt ()
{
	const size_t last = size - 1;
	
	Real errmax ;
	Real errmin ;
	int i, k, top;
	
	size_t i_nrmax = nrmax;
	size_t i_maxerr = order[i_nrmax] ;
	
	/* Check whether the list contains more than two error estimates */
	
	if (last < 2) 
	{
      order[0] = 0 ;
      order[1] = 1 ;
      i_work = i_maxerr ;
      return ;
	}
	
	errmax = elist[i_maxerr] ;
	
	/* This part of the routine is only executed if, due to a difficult
	 integrand, subdivision increased the error estimate. In the normal
	 case the insert procedure should start after the nrmax-th largest
	 error estimate. */
	
	while (i_nrmax > 0 && errmax > elist[order[i_nrmax - 1]]) 
	{
      order[i_nrmax] = order[i_nrmax - 1] ;
      i_nrmax-- ;
	} 
	
	/* Compute the number of elements in the list to be maintained in
	 descending order. This number depends on the number of
	 subdivisions still allowed. */
	
	if(last < (limit/2 + 2)) 
	{
      top = last ;
	}
	else
	{
      top = limit - last + 1;
	}
	
	/* Insert errmax by traversing the list top-down, starting
	 comparison from the element elist(order(i_nrmax+1)). */
	
	i = i_nrmax + 1 ;
	
	/* The order of the tests in the following line is important to
	 prevent a segmentation fault */
	
	while (i < top && errmax < elist[order[i]])
	{
      order[i-1] = order[i] ;
      i++ ;
	}
	
	order[i-1] = i_maxerr ;
	
	/* Insert errmin by traversing the list bottom-up */
	
	errmin = elist[last] ;
	
	k = top - 1 ;
	
	while (k > i - 2 && errmin >= elist[order[k]])
	{
      order[k+1] = order[k] ;
      k-- ;
	}
	
	order[k+1] = last ;
	
	/* Set i_max and e_max */
	
	i_maxerr = order[i_nrmax] ;
	
	i_work = i_maxerr ;
	nrmax = i_nrmax ;
}

template <class Real>
void
Workspace<Real>::set_initial_result (Real result, Real error)
{
	size = 1;
	rlist[0] = result;
	elist[0] = error;
}

template <class Real>
void
Workspace<Real>::retrieve (Real& a, Real& b, Real& r, Real& e)
{
	a = alist[i_work] ;
	b = blist[i_work] ;
	r = rlist[i_work] ;
	e = elist[i_work] ;
}

template <class Real>
Real
Workspace<Real>::sum_results ()
{
	size_t k;
	Real result_sum = Real(0);
	
	for (k = 0; k < size; k++)
	{
      result_sum += rlist[k];
	}
	return result_sum;
}

template <class Real>
int
Workspace<Real>::subinterval_too_small (Real a1, Real a2, Real b2)
{
	const Real e = this->eps_;
	const Real u = this->xmin_;
	
	Real tmp = (1 + 100 * e) * (this->abs(a2) + 1000 * u);
	int status = (this->abs(a1) <= tmp && this->abs(b2) <= tmp);
	return status;
}

template <class Real>
void 
Workspace<Real>::update (Real a1, Real b1, Real area1, Real error1,
								 Real a2, Real b2, Real area2, Real error2)
{
	const size_t i_max = i_work ;
	const size_t i_new = size ;
	
	const size_t new_level = level[i_max] + 1;
	
	/* append the newly-created intervals to the list */
	
	if (error2 > error1)
	{
      alist[i_max] = a2;        /* blist[maxerr] is already == b2 */
      rlist[i_max] = area2;
      elist[i_max] = error2;
      level[i_max] = new_level;
      
      alist[i_new] = a1;
      blist[i_new] = b1;
      rlist[i_new] = area1;
      elist[i_new] = error1;
      level[i_new] = new_level;
	}
	else
	{
      blist[i_max] = b1;        /* alist[maxerr] is already == a1 */
      rlist[i_max] = area1;
      elist[i_max] = error1;
      level[i_max] = new_level;
      
      alist[i_new] = a2;
      blist[i_new] = b2;
      rlist[i_new] = area2;
      elist[i_new] = error2;
      level[i_new] = new_level;
	}
	
	size++;
	
	if (new_level > maximum_level)
	{
      maximum_level = new_level;
	}
	
	qpsrt () ;
}

/* -------------------------------------------------------------------------- */
template <class Real>
int Workspace<Real>::qag(FtnBase<Real>& f, Real a, Real b,
							Real epsabs, Real epsrel, Real& result, Real& abserr)
{
	Real area, errsum;
	Real result0, abserr0, resabs0, resasc0;
	Real tolerance;
	size_t iteration = 0;
	int roundoff_type1 = 0, roundoff_type2 = 0, error_type = 0;
	
	Real round_off;     
	
	/* Initialize results */
	
	initialise (a, b);
	
	result = Real(0);
	abserr = Real(0);
	
	if (epsabs <= Real(0) && epsrel < 50 * this->eps_)
	{
      static const char* 
		message = "tolerance cannot be acheived with given epsabs and epsrel";
		throw message;
		return GSL_EBADTOL;
	}
	
	/* perform the first integration */
	
	this->qk(f, a, b, result0, abserr0, resabs0, resasc0);
	
	set_initial_result (result0, abserr0);
	
	/* Test on accuracy */
	
	tolerance = this->max (epsabs, epsrel * this->abs (result0));
	
	/* need IEEE rounding here to match original quadpack behavior
	 *
	 * round_off = GSL_COERCE_DBL (50 * GSL_DBL_EPSILON * resabs0);
	 */
	round_off = 50 * this->eps_ * resabs0;
	
	if (abserr0 <= round_off && abserr0 > tolerance)
	{
      result = result0;
      abserr = abserr0;
		
      static const char* 
		message = "cannot reach tolerance because of roundoff error on first attempt";
		throw message;
		return GSL_EROUND;
	}
	else if ((abserr0 <= tolerance && abserr0 != resasc0) 
				|| abserr0 == Real(0))
	{
      result = result0;
      abserr = abserr0;
		
      return GSL_SUCCESS;
	}
	else if (limit == 1)
	{
      result = result0;
      abserr = abserr0;
		
      static const char* 
		message = "a maximum of one iteration was insufficient";
		throw message;
		return GSL_EMAXITER;
	}
	
	area = result0;
	errsum = abserr0;
	
	iteration = 1;
	
	do
	{
      Real a1, b1, a2, b2;
      Real a_i, b_i, r_i, e_i;
      Real area1 = Real(0), area2 = Real(0), area12 = Real(0);
      Real error1 = Real(0), error2 = Real(0), error12 = Real(0);
      Real resasc1, resasc2;
      Real resabs1, resabs2;
		
      /* Bisect the subinterval with the largest error estimate */
		
      retrieve(a_i, b_i, r_i, e_i);
		
      a1 = a_i; 
      b1 = (a_i + b_i) / Real(2);
      a2 = b1;
      b2 = b_i;
		
      this->qk(f, a1, b1, area1, error1, resabs1, resasc1);
      this->qk(f, a2, b2, area2, error2, resabs2, resasc2);
		
      area12 = area1 + area2;
      error12 = error1 + error2;
		
      errsum += (error12 - e_i);
      area += area12 - r_i;
		
		/*
		 * ISSUE: Should checking for roundoff errors using QUADPACK's 
		 * tolerances for double precision be disabled for multiple
		 * precision arithmetic.
		 *
		 */
		 
      if (resasc1 != error1 && resasc2 != error2)
		{
			Real delta = r_i - area12;
			
			if (this->abs(delta) <= 1.0e-5 * this->abs(area12) 
				 && error12 >= 0.99 * e_i)
			{
				roundoff_type1++;
			}
			if (iteration >= 10 && error12 > e_i)
			{
				roundoff_type2++;
			}
		}
		
		tolerance = this->max (epsabs, epsrel * this->abs(area));
		
      if (errsum > tolerance)
		{
			if (roundoff_type1 >= 6 || roundoff_type2 >= 20)
			{
				error_type = 2;   /* round off error */
			}
			
			/* set error flag in the case of bad integrand behaviour at
			 a point of the integration range */
			
			if (subinterval_too_small (a1, a2, b2))
			{
				error_type = 3;
			}
		}
		
      update (a1, b1, area1, error1, a2, b2, area2, error2);
		
      retrieve(a_i, b_i, r_i, e_i);
		
      iteration++;
		
	}
	while (iteration < limit && !error_type && errsum > tolerance);
	
	result = sum_results();
	abserr = errsum;
	
	if (errsum <= tolerance)
	{
      return GSL_SUCCESS;
	}
	else if (error_type == 2)
	{
      static const char* 
		message = "roundoff error prevents tolerance from being achieved";
		throw message;
		return GSL_EROUND;
	}
	else if (error_type == 3)
	{
      static const char* 
		message = "bad integrand behavior found in the integration interval";
		throw message;
		return GSL_ESING;
	}
	else if (iteration == limit)
	{
      static const char* 
		message = "maximum number of subdivisions reached";
		throw message;
		return GSL_EMAXITER;
	}
	else
	{
      static const char* 
		message = "could not integrate function";
		throw message;
		return GSL_EFAILED;
	}
}

#endif // GSL_ERRCODES_H

