/* errcodes.h
 * Part of quadpack++.
 */

#ifndef GSL_ERRCODES_H
#define GSL_ERRCODES_H

namespace QuadPack
{

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


	inline const char* const returnCodeToString(int error)
	{
	    switch(error)
	    {
	        case GSL_SUCCESS:
	            return "GSL_SUCCESS";
	        case GSL_FAILURE:
	            return "GSL_FAILURE";
	        case GSL_CONTINUE:
	            return "GSL_CONTINUE: iteration has not converged";
	        case GSL_EDOM:
	            return "GSL_EDOM: input domain error, e.g sqrt(-1)";
	        case GSL_ERANGE:
	            return "GSL_ERANGE: output range error, e.g. exp(1e100)";
	        case GSL_EFAULT:
	            return "GSL_EFAULT: invalid pointer";
	        case GSL_EINVAL:
	            return "GSL_EINVAL: invalid argument supplied by user";
	        case GSL_EFAILED:
	            return "GSL_EFAILED: generic failure";
	        case GSL_EFACTOR:
	            return "GSL_EFACTOR: factorization failed";
	        case GSL_ESANITY:
	            return "GSL_ESANITY: sanity check failed - shouldn't happen";
	        case GSL_ENOMEM:
	            return "GSL_ENOMEM: malloc failed";
	        case GSL_EBADFUNC:
	            return "GSL_EBADFUNC: problem with user-supplied function";
	        case GSL_ERUNAWAY:
	            return "GSL_ERUNAWAY: iterative process is out of control";
	        case GSL_EMAXITER:
	            return "GSL_EMAXITER: exceeded max number of iterations";
	        case GSL_EZERODIV:
	            return "GSL_EZERODIV: tried to divide by zero";
	        case GSL_EBADTOL:
	            return "GSL_EBADTOL: user specified an invalid tolerance";
	        case GSL_ETOL:
	            return "GSL_ETOL: failed to reach the specified tolerance";
	        case GSL_EUNDRFLW:
	            return "GSL_EUNDRFLW: underflow";
	        case GSL_EOVRFLW:
	            return "GSL_EOVRFLW: overflow ";
	        case GSL_ELOSS:
	            return "GSL_ELOSS: loss of accuracy";
	        case GSL_EROUND:
	            return "GSL_EROUND: failed because of roundoff error";
	        case GSL_EBADLEN:
	            return "GSL_EBADLEN: matrix, vector lengths are not conformant";
	        case GSL_ENOTSQR:
	            return "GSL_ENOTSQR: matrix not square";
	        case GSL_ESING:
	            return "GSL_ESING: apparent singularity detected";
	        case GSL_EDIVERGE:
	            return "GSL_EDIVERGE: integral or series is divergent";
	        case GSL_EUNSUP:
	            return "GSL_EUNSUP: requested feature is not supported by the hardware";
	        case GSL_EUNIMPL:
	            return "GSL_EUNIMPL: requested feature not (yet) implemented";
	        case GSL_ECACHE:
	            return "GSL_ECACHE: cache limit exceeded";
	        case GSL_ETABLE:
	            return "GSL_ETABLE: table limit exceeded";
	        case GSL_ENOPROG:
	            return "GSL_ENOPROG: iteration is not making progress towards solution";
	        case GSL_ENOPROGJ:
	            return "GSL_ENOPROGJ: jacobian evaluations are not improving the solution";
	        case GSL_ETOLF:
	            return "GSL_ETOLF: cannot reach the specified tolerance in F";
	        case GSL_ETOLX:
	            return "GSL_ETOLX: cannot reach the specified tolerance in X";
	        case GSL_ETOLG:
	            return "GSL_ETOLG: cannot reach the specified tolerance in gradient";
	        case GSL_EOF:
	            return "GSL_EOF: end of file";
	    }
	    return "Unkown";
	}

}

#endif // GSL_ERRCODES_H

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

