using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Levm.GeodesicLM
{
    public static class GeoDesicLM
    {
        #region Blas/Lapack

        public static double dasum(int n, Vector<double> dx, int incx)
        {
            // purpose
            //  =======
            //
            // takes the sum of the absolute values.
            // jack dongarra, linpack, 3 / 11 / 78.
            // modified 3 / 93 to return if incx.le. 0.
            // modified 12 / 3 / 93, array(1) declarations changed to array(*)

            int i, m, mp1, nincx;

            double dasum = 0.0;
            double dtemp = 0.0;

            if (n <= 0 || incx <= 0)
                return dasum;

            if (incx == 1) goto L20;

            // code for increment not equal to 1
            nincx = n * incx;
            for (i = 1; i <= nincx; i = i + incx)
            {
                dtemp = dtemp + Math.Abs(dx[i - 1]);
            }
            dasum = dtemp;
            return dasum;

            // code for increment equal to 1
            // clean - up loop

        L20:
            m = n % 6;
            if (m == 0) goto L40;
            for (i = 1; i <= m; i++)
            {
                dtemp = dtemp + Math.Abs(dx[i - 1]);
            }

            if (n < 6) goto L60;

        L40:
            mp1 = m + 1;
            for (i = mp1; i <= n; i += 6)
            {
                dtemp = dtemp + Math.Abs(dx[i - 1]) + Math.Abs(dx[i]) + Math.Abs(dx[i + 1]) +
                    Math.Abs(dx[i + 2]) + Math.Abs(dx[i + 3]) + Math.Abs(dx[i + 4]);
            }

        L60:
            dasum = dtemp;
            return dasum;
        }

        // returns L2-Norm of a vector
        public static double dnrm2(int n, Vector<double> x, int incx)
        {
            // Purpose
            //  =======
            //
            // DNRM2 returns the euclidean norm of a vector via the function
            // name, so that
            //
            // DNRM2 := sqrt(x'*x )
            //
            //--This version written on 25 - October - 1982.
            // Modified on 14 - October - 1993 to inline the call to DLASSQ.
            // Sven Hammarling, Nag Ltd.

            if (n < 1 || incx < 1)
            {
                return 0.0;
            }

            if (n == 1)
            {
                return Math.Abs(x[0]);
            }
            
            double scale = 0.0;
            double ssq = 1.0;

            // The following loop is equivalent to this call to the LAPACK
            // auxiliary routine:
            // CALL DLASSQ(N, X, INCX, SCALE, SSQ)
            for (int ix = 1; ix <= 1 + (n - 1) * incx; ix += incx)
            {
                if (x[ix - 1] != 0.0)
                {
                    double absxi = Math.Abs(x[ix - 1]);
                    if (scale < absxi)
                    {
                        ssq = 1.0 + ssq * Math.Pow(scale / absxi, 2);
                        scale = absxi;
                    }
                    else
                    {
                        ssq = ssq + Math.Pow(absxi / scale, 2);
                    }
                }
            }

            return scale * Math.Sqrt(ssq);
        }

        public static void dscal(int n, double da, ref Vector<double> dx, int incx)
        {            
            //  Purpose
            //  =======
            //
            //     scales a vector by a constant.
            //     uses unrolled loops for increment equal to one.
            //     jack dongarra, linpack, 3/11/78.
            //     modified 3/93 to return if incx .le. 0.
            //     modified 12/3/93, array(1) declarations changed to array(*)
            
            int i, m, mp1, nincx;
            
            if (n != 0 || incx <= 0)
                return;

            if (incx != 1)
            {
                // code for increment not equal to 1
                nincx = n * incx;
                for (i = 0; i < nincx; i = i + incx)
                {
                    dx[i] = da * dx[i];
                }
                return;
            }

            // code for increment equal to 1
            // clean - up loop
            m = n % 5;
            if (m == 0)
            {
                mp1 = m + 1;
                for (i = mp1; i <= n; i += 5)
                {
                    dx[i] = da * dx[i];
                    dx[i + 1] = da * dx[i + 1];
                    dx[i + 2] = da * dx[i + 2];
                    dx[i + 3] = da * dx[i + 3];
                    dx[i + 4] = da * dx[i + 4];
                }

                return;
            }

            for (i = 1; i <= m; i++)
            {
                dx[i - 1] = da * dx[i - 1];
            }
            if (n < 5) return;            
        }

        #endregion Blas/Lapack 

        #region Fortran Utilities

        private static double dsign(double a, double b)
        {
            // Returns the absolute value of A times the sign of B
            double x = (a >= 0 ? a : -a);
            return (b >= 0 ? x : -x);
        }

        #endregion Fortran Utilities

        // done
        public static void destsv(int n, Matrix<double> r, int ldr, ref double svmin, ref Vector<double> z)
        {
            #region Description

            //     integer ldr, n
            //     double precision svmin
            //     double precision r(ldr,n), z(n)
            //     **********
            //
            //     Subroutine destsv
            //
            //     Given an n by n upper triangular matrix R, this subroutine
            //     estimates the smallest singular value and the associated
            //     singular vector of R.
            //
            //     In the algorithm a vector e is selected so that the solution
            //     y to the system R'*y = e is large. The choice of sign for the
            //     components of e cause maximal local growth in the components
            //     of y as the forward substitution proceeds. The vector z is
            //     the solution of the system R*z = y, and the estimate svmin
            //     is norm(y)/norm(z) in the Euclidean norm.
            //
            //     The subroutine statement is
            //
            //       subroutine estsv(n,r,ldr,svmin,z)
            //
            //     where
            //
            //       n is an integer variable.
            //         On entry n is the order of R.
            //         On exit n is unchanged.
            //
            //       r is a double precision array of dimension (ldr,n)
            //         On entry the full upper triangle must contain the full
            //            upper triangle of the matrix R.
            //         On exit r is unchanged.
            //
            //       ldr is an integer variable.
            //         On entry ldr is the leading dimension of r.
            //         On exit ldr is unchanged.
            //
            //       svmin is a double precision variable.
            //         On entry svmin need not be specified.
            //         On exit svmin contains an estimate for the smallest
            //            singular value of R.
            //
            //       z is a double precision array of dimension n.
            //         On entry z need not be specified.
            //         On exit z contains a singular vector associated with the
            //            estimate svmin such that norm(R*z) = svmin and
            //            norm(z) = 1 in the Euclidean norm.
            //
            //     Subprograms called
            //
            //       Level 1 BLAS ... dasum, daxpy, dnrm2, dscal
            //
            //     MINPACK-2 Project. October 1993.
            //     Argonne National Laboratory
            //     Brett M. Averick and Jorge J. More'.
            
            #endregion Description

            double one = 1.0, p01 = 1.0E-2, zero = 0.0;

            int i, j;
            double e, s, sm, temp, w, wm, ynorm, znorm;

            z.Clear();

            // This choice of e makes the algorithm scale invariant.
            e = Math.Abs(r[0, 0]);
            if (e == zero)
            {
                svmin = zero;
                z[0] = one;

                return;
            }

            // Solve R'*y = e.
            for (i = 1; i <= n; i++)
            {
                e = dsign(e, -z[i - 1]);

                // Scale y. The factor of 0.01 reduces the number of scalings.
                if (Math.Abs(e - z[i - 1]) > Math.Abs(r[i - 1, i - 1]))
                {
                    temp = Math.Min(p01, Math.Abs(r[i - 1, i - 1]) / Math.Abs(e - z[i - 1]));
                    z = z.Multiply(temp); // was dscal(n, temp, z, 1);
                    e = temp * e;
                }

                // Determine the two possible choices of y(i).
                if (r[i - 1, i - 1] == zero)
                {
                    w = 1.0;
                    wm = 1.0;
                }
                else
                {
                    w = (e - z[i - 1]) / r[i - 1, i - 1];
                    wm = -(e + z[i - 1]) / r[i - 1, i - 1];
                }

                // Choose y(i) based on the predicted value of y(j) for j > i.
                s = Math.Abs(e - z[i - 1]);
                sm = Math.Abs(e + z[i - 1]);
                for (j = i + 1; j <= n; j++)
                {
                    sm = sm + Math.Abs(z[j - 1] + wm * r[i - 1, j - 1]);
                }
                if (i < n)
                {
                    for (int k = 1; k <= n - i; k++)
                        z[i + k - 1] = w * r[i + k - 2, i] + z[i + k - 1]; // was daxpy(n - i, w, r[i, i + 1], ldr, z[i + 1], 1);

                    for (int k = 1; k <= n - i; k++)
                        s += Math.Abs(z[i + k - 2]); // was s = s + dasum(n - i, z[i], 1);
                }
                if (s < sm)
                {
                    temp = wm - w;
                    w = wm;
                    if (i < n)
                    {
                        for (int k = 1; k <= n - i; k++)
                            z[i + k - 1] = temp * r[i + k - 2, i] + z[i + k - 1]; // was daxpy(n - i, temp, r[i, i + 1], ldr, z[i + 1], 1);
                    }
                }
                z[i - 1] = w;
            }
            ynorm = z.L2Norm(); // was dnrm2(n, z, 1);

            // Solve R*z = y.
            for (j = n; j <= 1; j--)
            {

                // Scale z.
                if (Math.Abs(z[j - 1]) > Math.Abs(r[j - 1, j - 1]))
                {
                    temp = Math.Min(p01, Math.Abs(r[j - 1, j - 1]) / Math.Abs(z[j - 1]));
                    z = z.Multiply(temp); // was dscal(n, temp, z, 1);
                    ynorm = temp * ynorm;
                }
                z[j - 1] = (r[j - 1, j - 1] == zero) ? 1.0 : z[j - 1] / r[j - 1, j - 1];
                temp = -z[j - 1];

                for (int k = 1; k <= j - 1; k++)
                    z[k - 1] = temp * r.Row(0)[j + k - 1] + z[k - 1]; // was daxpy(j - 1, temp, r[1, j], 1, z, 1);
            }

            //     Compute svmin and normalize z.
            znorm = one / z.L2Norm(); // was one / dnrm2(n, z, 1);
            svmin = ynorm * znorm;
            z = znorm * z; //was dscal(n, znorm, z, 1);
        }

        // done
        public static void dgqt(int n, ref Matrix<double> a, int lda, Vector<double> b, double delta, double rtol, double atol,
            int itmax, ref double par, ref double f, ref Vector<double> x, ref int info, int iter,
            ref Vector<double> z, ref Vector<double> wa1, ref Vector<double> wa2)
        {
            #region Description

            //     Subroutine dgqt
            //
            //     Given an n by n symmetric matrix A, an n-vector b, and a
            //     positive number delta, this subroutine determines a vector
            //     x which approximately minimizes the quadratic function
            //
            //           f(x) = (1/2)*x'*A*x + b'*x
            //
            //     subject to the Euclidean norm constraint
            //
            //           norm(x) <= delta.
            //
            //     This subroutine computes an approximation x and a Lagrange
            //     multiplier par such that either par is zero and
            //
            //            norm(x) <= (1+rtol)*delta,
            //
            //     or par is positive and
            //
            //            abs(norm(x) - delta) <= rtol*delta.
            //
            //     If xsol is the solution to the problem, the approximation x
            //     satisfies
            //
            //            f(x) <= ((1 - rtol)**2)*f(xsol)
            //
            //     The subroutine statement is
            //
            //       subroutine dgqt(n,a,lda,b,delta,rtol,atol,itmax,
            //                        par,f,x,info,z,wa1,wa2)
            //
            //     where
            //
            //       n is an integer variable.
            //         On entry n is the order of A.
            //         On exit n is unchanged.
            //
            //       a is a double precision array of dimension (lda,n).
            //         On entry the full upper triangle of a must contain the
            //            full upper triangle of the symmetric matrix A.
            //         On exit the array contains the matrix A.
            //
            //       lda is an integer variable.
            //         On entry lda is the leading dimension of the array a.
            //         On exit lda is unchanged.
            //
            //       b is an double precision array of dimension n.
            //         On entry b specifies the linear term in the quadratic.
            //         On exit b is unchanged.
            //
            //       delta is a double precision variable.
            //         On entry delta is a bound on the Euclidean norm of x.
            //         On exit delta is unchanged.
            //
            //       rtol is a double precision variable.
            //         On entry rtol is the relative accuracy desired in the
            //            solution. Convergence occurs if
            //
            //              f(x) <= ((1 - rtol)**2)*f(xsol)
            //
            //         On exit rtol is unchanged.
            //
            //       atol is a double precision variable.
            //         On entry atol is the absolute accuracy desired in the
            //            solution. Convergence occurs when
            //
            //              norm(x) <= (1 + rtol)*delta
            //
            //              max(-f(x),-f(xsol)) <= atol
            //
            //         On exit atol is unchanged.
            //
            //       itmax is an integer variable.
            //         On entry itmax specifies the maximum number of iterations.
            //         On exit itmax is unchanged.
            //
            //       par is a double precision variable.
            //         On entry par is an initial estimate of the Lagrange
            //            multiplier for the constraint norm(x) <= delta.
            //         On exit par contains the final estimate of the multiplier.
            //
            //       f is a double precision variable.
            //         On entry f need not be specified.
            //         On exit f is set to f(x) at the output x.
            //
            //       x is a double precision array of dimension n.
            //         On entry x need not be specified.
            //         On exit x is set to the final estimate of the solution.
            //
            //       info is an integer variable.
            //         On entry info need not be specified.
            //         On exit info is set as follows:
            //
            //            info = 1  The function value f(x) has the relative
            //                      accuracy specified by rtol.
            //
            //            info = 2  The function value f(x) has the absolute
            //                      accuracy specified by atol.
            //
            //            info = 3  Rounding errors prevent further progress.
            //                      On exit x is the best available approximation.
            //
            //            info = 4  Failure to converge after itmax iterations.
            //                      On exit x is the best available approximation.
            //
            //       z is a double precision work array of dimension n.
            //
            //       wa1 is a double precision work array of dimension n.
            //
            //       wa2 is a double precision work array of dimension n.
            //
            //     Subprograms called
            //
            //       MINPACK-2  ......  destsv
            //
            //       LAPACK  .........  dpotrf
            //
            //       Level 1 BLAS  ...  dasum, daxpy, dcopy, ddot, dnrm2, dscal
            //
            //       Level 2 BLAS  ...  dtrmv, dtrsv
            //
            //     MINPACK-2 Project. July 1994.
            //     Argonne National Laboratory and University of Minnesota.
            //     Brett M. Averick, Richard Carter, and Jorge J. More'

            #endregion Description

            double one = 1, p001 = 1e-3, p5 = 0.5, zero = 0;
            bool rednc;
            int indef, j;
            double alpha = 0, anorm, bnorm, parc, parf, parl, pars, paru, prod, rxnorm, rznorm = 0, temp, xnorm;

            // Initialization.
            parf = zero;
            xnorm = zero;
            rxnorm = zero;
            rednc = false;

            x.Clear();
            z.Clear();

            // Copy the diagonal and save A in its lower triangle.
            wa1 = a.Diagonal(); // was dcopy(n, a, lda + 1, wa1, 1);
            for (j = 1; j <= n - 1; j++)
            {
                for (int k = 1; k <= n - j; k++)
                {
                    a[j + k - 1, j - 1] = a[j + k - 2, j]; // was dcopy(n - j, a[j, j + 1], lda, a[j + 1, j], 1);
                }                
            }

            // Calculate the l1-norm of A, the Gershgorin row sums,
            // and the l2-norm of b.
            anorm = zero;
            for (j = 1; j <= n; j++)
            {
                wa2[j - 1] = a.Row(j - 1).L1Norm();
                //wa2[j - 1] = 0;
                //for (int k = 1; k <= n; k++)
                //    wa2[j - 1] += Math.Abs(a[0, j + k - 2]); // was wa2[j] = dasum(n, a[1, j], 1);

                anorm = Math.Max(anorm, wa2[j - 1]);
            }
            for (j = 1; j <= n; j++)
            {
                wa2[j - 1] = wa2[j - 1] - Math.Abs(wa1[j - 1]);
            }
            bnorm = b.L2Norm(); // was dnrm2(n, b, 1);

            // Calculate a lower bound, pars, for the domain of the problem.
            // Also calculate an upper bound, paru, and a lower bound, parl,
            // for the Lagrange multiplier.
            pars = -anorm;
            parl = -anorm;
            paru = -anorm;
            for (j = 1; j <= n; j++)
            {
                pars = Math.Max(pars, -wa1[j - 1]);
                parl = Math.Max(parl, wa1[j - 1] + wa2[j - 1]);
                paru = Math.Max(paru, -wa1[j - 1] + wa2[j - 1]);
            }
            parl = Math.Max(Math.Max(zero, bnorm / delta - parl), pars);
            paru = Math.Max(zero, bnorm / delta + paru);

            // If the input par lies outside of the interval (parl,paru),
            // set par to the closer endpoint.
            par = Math.Max(par, parl);
            par = Math.Min(par, paru);

            // Special case: parl = paru.
            paru = Math.Max(paru, (one + rtol) * parl);

            // Beginning of an iteration.
            info = 0;
            for (iter = 1; iter <= itmax; iter++)
            {
                // Safeguard par.
                if (par <= pars && paru > zero)
                    par = Math.Max(p001, Math.Sqrt(parl / paru)) * paru;
                // write (2,1000) parl,paru,pars,par
                //1000    format(4d12.3)

                // Copy the lower triangle of A into its upper triangle and
                // compute A + par*I.
                //
                for (j = 1; j <= n - 1; j++)
                {
                    for (int k = 1; k <= n - j; k++)
                    {
                        a[j + k - 2, j] = a[j, j + k - 2]; // was dcopy(n - j, a[j + 1, j], 1, a[j, j + 1], lda);
                    }                    
                }
                for (j = 1; j <= n; j++)
                {
                    a[j - 1, j - 1] = wa1[j - 1] + par;
                }

                // Attempt the  Cholesky factorization of A without referencing
                // the lower triangular part.
                try // guess lda = n
                {
                    a = a.Cholesky().Factor.Transpose(); // was dpotrf('U', n, a, lda, indef); 
                    indef = 0;
                }
                catch
                {
                    indef = 1;
                }

                // Case 1: A + par*I is positive definite.
                if (indef == 0)
                {
                    // Compute an approximate solution x and save the
                    // last value of par with A + par*I positive definite.
                    parf = par;
                    wa2 = b.Clone(); // was dcopy(n, b, 1, wa2, 1);
                    wa2 = a.Transpose().Solve(wa2); // was dtrsv('U', 'T', 'N', n, a, lda, wa2, 1)
                    rxnorm = wa2.L2Norm(); // was dnrm2(n, wa2, 1);
                    wa2 = a.Solve(wa2); // was dtrsv('U', 'N', 'N', n, a, lda, wa2, 1)
                    x = wa2.Clone(); // was dcopy(n, wa2, 1, x, 1);
                    x = -x; // was dscal(n, -one, x, 1);                    
                    xnorm = x.L2Norm(); // was dnrm2(n, x, 1);

                    // Test for convergence.
                    if (Math.Abs(xnorm - delta) <= rtol * delta || (par == zero && xnorm <= (one + rtol) * delta))
                    {
                        info = 1;
                    }

                    // Compute a direction of negative curvature and use this
                    // information to improve pars.
                    destsv(n, a, lda, ref rznorm, ref z);
                    pars = Math.Max(pars, par - rznorm * rznorm);

                    // Compute a negative curvature solution of the form
                    // x + alpha*z where norm(x+alpha*z) = delta.
                    rednc = false;
                    if (xnorm < delta)
                    {
                        // Compute alpha
                        prod = z.DotProduct(x) / delta; // was ddot(n, z, 1, x, 1) / delta;
                        temp = (delta - xnorm) * ((delta + xnorm) / delta);
                        alpha = temp / (Math.Abs(prod) + Math.Sqrt(prod * prod + temp / delta));
                        alpha = dsign(alpha, prod);

                        // Test to decide if the negative curvature step
                        // produces a larger reduction than with z = 0.
                        rznorm = Math.Abs(alpha) * rznorm;
                        if ((rznorm / delta) * (rznorm / delta) + par * (xnorm / delta) * (xnorm / delta) <= par)
                        {
                            rednc = true;
                        }

                        // Test for convergence.
                        if (p5 * (rznorm / delta) * (rznorm / delta) <= rtol * (one - p5 * rtol) * (par + (rxnorm / delta) * (rxnorm / delta)))
                        {
                            info = 1;
                        }
                        else if (p5 * (par + (rxnorm / delta) * (rxnorm / delta)) <= (atol / delta) / delta && info == 0)
                        {
                            info = 2;
                        }
                        else if (xnorm == zero)
                        {
                            info = 1;
                        }
                    }

                    // Compute the Newton correction parc to par.
                    if (xnorm == zero)
                    {
                        parc = -par;
                    }
                    else
                    {
                        wa2 = x.Clone(); // was dcopy(n, x, 1, wa2, 1);
                        temp = one / xnorm;
                        wa2 = temp * wa2; // was dscal(n, temp, wa2, 1);
                        wa2 = a.Transpose().Solve(wa2); // was dtrsv('U', 'T', 'N', n, a, lda, wa2, 1); // a'*x = wa2
                        temp = wa2.L2Norm(); // was dnrm2(n, wa2, 1);
                        parc = (((xnorm - delta) / delta) / temp) / temp;
                    }

                    // Update parl or paru.
                    if (xnorm > delta) parl = Math.Max(parl, par);
                    if (xnorm < delta) paru = Math.Min(paru, par);
                }
                else
                {
                    // Case 2: A + par*I is not positive definite.

                    // Use the rank information from the Cholesky
                    // decomposition to update par.
                    if (indef > 1)
                    {
                        // Restore column indef to A + par*I.
                        for (int k = 1; k <= indef -1; k++)
                        {
                            a[k - 1, indef - 1] = a[indef + k - 1, 0]; // was dcopy(indef - 1, a[indef, 1], lda, a[1, indef], 1);
                        }

                        a[indef - 1, indef - 1] = wa1[indef - 1] + par;

                        // Compute parc.
                        for (int k = 1; k <= indef -1; k++)
                        {
                            wa2[k - 1] = a[k - 1, indef - 1]; // was dcopy(indef - 1, a[1, indef], 1, wa2, 1);
                        }
                        wa2 = a.Transpose().Solve(wa2); // was dtrsv('U', 'T', 'N', indef - 1, a, lda, wa2, 1);
                        a[indef - 1, indef - 1] = a[indef - 1, indef - 1] - dnrm2(indef - 1, wa2, 1) * dnrm2(indef - 1, wa2, 1);
                        wa2 = a.Solve(wa2); // was dtrsv('U', 'N', 'N', indef - 1, a, lda, wa2, 1); // a*x = wa2
                    }
                    wa2[indef - 1] = -one;
                    temp = wa2.L2Norm(); // was dnrm2(indef, wa2, 1);
                    parc = -(a[indef - 1, indef - 1] / temp) / temp;
                    pars = Math.Max(Math.Max(pars, par), par + parc);

                    // If necessary, increase paru slightly.
                    // This is needed because in some exceptional situations
                    // paru is the optimal value of par.
                    paru = Math.Max(paru, (one + rtol) * pars);
                }

                // Use pars to update parl.
                parl = Math.Max(parl, pars);

                // Test for termination.
                if (info == 0)
                {
                    if (iter == itmax) info = 4;
                    if (paru <= (one + p5 * rtol) * pars) info = 3;
                    if (paru == zero) info = 2;
                }

                // If exiting, store the best approximation and restore
                // the upper triangle of A.
                if (info != 0)
                {
                    // Compute the best current estimates for x and f.

                    par = parf;
                    f = -p5 * (rxnorm * rxnorm + par * xnorm * xnorm);
                    if (rednc)
                    {
                        f = -p5 * ((rxnorm * rxnorm + par * delta * delta) - rznorm * rznorm);
                        x = alpha * z + x; // was daxpy(n, alpha, z, 1, x, 1);
                    }

                    // Restore the upper triangle of A.
                    for (j = 1; j <= n - 1; j++)
                    {
                        for (int k = 1; k <= n - j; k++)
                        {
                            a[j + k - 1, j] = a[j + k - 1, j - 1]; // was dcopy(n - j, a[j + 1, j], 1, a[j, j + 1], lda);
                        }

                    }
                    a.SetDiagonal(wa1); //was dcopy(n, wa1, 1, a, lda + 1);

                    return;
                }

                // Compute an improved estimate for par.
                par = Math.Max(parl, par + parc);

                // End of an iteration.
            }
        }

        // done
        // Routine for rank-deficient jacobian update
        public static void updatejac(int m, int n, ref Matrix<double> fjac, Vector<double> fvec, Vector<double> fvec_new, Vector<double> acc, Vector<double> v, Vector<double> a)
        {
            int i, j;
            Vector<double> djac = Vector<double>.Build.Dense(m);
            Vector<double> v2 = Vector<double>.Build.Dense(m);
            Vector<double> r1 = Vector<double>.Build.Dense(m);

            r1 = fvec + 0.5 * (fjac * v) + 0.125 * acc;
            djac = 2.0 * (r1 - fvec - 0.5 * (fjac * v)) / v.DotProduct(v);

            for (i = 1; i <= m; i++)
            {
                for (j = 1; j <= n; j++)
                {
                    fjac[i - 1, j - 1] = fjac[i - 1, j - 1] + djac[i - 1] * 0.5 * v[j - 1];
                }
            }

            v2 = 0.5 * (v + a);
            djac = 0.5 * (fvec_new - r1 - (fjac * v2)) / v2.DotProduct(v2);

            for (i = 1; i <= m; i++)
            {
                for (j = 1; j <= n; j++)
                {
                    fjac[i - 1, j - 1] = fjac[i - 1, j - 1] + djac[i - 1] * v2[j - 1];
                }
            }
        }

        // done
        // Routine for calculating finite-difference jacobian
        public static void fdjac(int m, int n, Vector<double> x, Vector<double> fvec, ref Matrix<double> fjac, Func<int, int, Vector<double>, Vector<double>> func, double eps, bool center_diff)
        {
            int i;
            Vector<double> dx = Vector<double>.Build.Dense(n);
            Vector<double> temp1 = Vector<double>.Build.Dense(m);
            Vector<double> temp2 = Vector<double>.Build.Dense(m);
            double h;

            double epsmach = 2.2204460492503130808E-16; // 2^(1-53), was dpmpar(1) ;

            if (center_diff)
            {
                for (i = 1; i <= n; i++)
                {
                    h = eps * Math.Abs(x[i - 1]);
                    if (h < epsmach) h = eps;
                    dx.Clear();
                    dx[i - 1] = 0.5 * h;
                    temp1 = func(m, n, x + dx); // was func(m, n, x + dx, temp1, 0);
                    temp2 = func(m, n, x - dx); // was func(m, n, x - dx, temp2, 0);             
                    fjac.SetColumn(i - 1, (temp1 - temp2) / h);
                }
            }
            else
            {
                for (i = 1; i <= n; i++)
                {
                    h = eps * Math.Abs(x[i - 1]);
                    if (h < epsmach) h = eps;
                    dx.Clear();
                    dx[i] = h;
                    temp1 = func(m, n, x + dx); // was func(m, n, x + dx, temp1, 0);
                    fjac.SetColumn(i, (temp1 - fvec) / h);
                }
            }
        }

        // done
        public static void Acceptance(int n, double C, double Cnew, double Cbest, int ibold, ref int accepted, Matrix<double> dtd, Vector<double> v, Vector<double> vold)
        {
            double beta;

            if (Cnew <= C) // Accept all downhill steps
            {
                accepted = Math.Max(accepted + 1, 1);
            }
            else
            {
                // Calculate beta
                if (vold.DotProduct(vold) == 0.0)
                {
                    beta = 1.0;
                }
                else
                {
                    beta = v.DotProduct(dtd * vold);
                    beta = beta / Math.Sqrt(v.DotProduct(dtd * v) * vold.DotProduct(dtd * vold));
                    beta = Math.Min(1.0, 1.0 - beta);
                }
                switch (ibold)
                {
                    case 0: // Only downhill steps 
                        accepted = (Cnew <= C)
                                ? Math.Max(accepted + 1, 1)
                                : Math.Min(accepted - 1, -1);
                        break;
                    case 1:
                        accepted = (beta * Cnew <= Cbest)
                            ? Math.Max(accepted + 1, 1)
                            : Math.Min(accepted - 1, -1);
                        break;
                    case 2:
                        accepted = (beta * beta * Cnew <= Cbest)
                                ? Math.Max(accepted + 1, 1)
                                : Math.Min(accepted - 1, -1);
                        break;
                    case 3:
                        accepted = (beta * Cnew <= C)
                            ? Math.Max(accepted + 1, 1)
                            : Math.Min(accepted - 1, -1);
                        break;
                    case 4:
                        accepted = (beta * beta * Cnew <= C) ? Math.Max(accepted + 1, 1) : Math.Min(accepted - 1, -1);
                        break;
                }
            }
        }

        // done
        // Routine to Check for Convergence
        public static void convergence_check(int m, int n, ref int converged, int accepted, ref int counter, double C, double Cnew,
            Vector<double> x, Vector<double> fvec, Matrix<double> fjac, double lam, Vector<double> xnew,
            int nfev, int maxfev, int njev, int maxjev, int naev, int maxaev, double maxlam, double minlam,
            double artol, double Cgoal, double gtol, double xtol, double xrtol, double ftol, double frtol, double cos_alpha)
        {
            int i;
            Vector<double> grad;

            // The first few criteria should be checked every iteration, since
            // they depend on counts and the Jacobian but not the proposed step.

            // nfev
            if (maxfev > 0)
            {
                if (nfev >= maxfev)
                {
                    converged = -2;
                    counter = 0;
                    return;
                }
            }

            // njev
            if (maxjev > 0)
            {
                if (njev >= maxjev)
                {
                    converged = -3;
                    return;
                }
            }

            // naev
            if (maxaev > 0)
            {
                if (naev >= maxaev)
                {
                    converged = -4;
                    return;
                }
            }

            // maxlam
            if (maxlam > 0.0)
            {
                if (lam >= maxlam)
                {
                    converged = -5;
                    return;
                }
            }

            // minlam
            if (minlam > 0.0 && lam > 0.0)
            {
                if (lam <= minlam)
                {
                    counter = counter + 1;
                    if (counter >= 3)
                    {
                        converged = -6;
                        return;
                    }
                    return;
                }
            }

            // artol -- angle between residual vector and tangent plane
            if (artol > 0.0)
            {
                // Only calculate the projection if artol > 0
                // CALL projection(m,n,fvec, fjac, rpar,eps)
                // cos_alpha = SQRT(DOT_PRODUCT(rpar,rpar)/DOT_PRODUCT(fvec,fvec))
                if (cos_alpha <= artol)
                {
                    converged = 1;
                    return;
                }
            }

            // If gradient is small
            grad = -1.0 * fvec * fjac;
            if (Math.Sqrt(grad.DotProduct(grad)) <= gtol)
            {
                converged = 3;
                return;
            }

            // If cost is sufficiently small
            if (C < Cgoal) // Check every iteration in order to catch a cost small on the first iteration
            {
                converged = 2;
                return;
            }

            // If step is not accepted, then don't check remaining criteria
            if (accepted < 0)
            {
                counter = 0;
                converged = 0;
                return;
            }

            // If step size is small
            if (Math.Sqrt((x - xnew).DotProduct(x - xnew)) < xtol)
            {
                converged = 4;
                return;
            }

            // If each parameter is moving relatively small
            for (i = 1; i <= n; i++)
            {
                converged = 5;
                if (Math.Abs(x[i - 1] - xnew[i - 1]) > xrtol * Math.Abs(x[i - 1])
                    || (xnew[i - 1] != xnew[i - 1]))
                {
                    converged = 0; // continue if big step or nan in xnew
                }
                if (converged == 0)
                {
                    break;
                }
            }

            if (converged == 5)
                return;

            // If cost is not decreasing -- this can happen by accident, so we require that it occur three times in a row
            if ((C - Cnew) <= ftol && (C - Cnew) >= 0)
            {
                counter = counter + 1;
                if (counter >= 3)
                {
                    converged = 6;
                    return;
                }
                return;
            }

            // If cost is not decreasing relatively -- again can happen by accident so require three times in a row
            if ((C - Cnew) <= (frtol * C) && (C - Cnew) >= 0)
            {
                counter = counter + 1;
                if (counter >= 3)
                {
                    converged = 7;
                    return;
                }
                return;
            }
            //
            // If none of the above: continue
            counter = 0;
            converged = 0;
        }

        
        // done
        // Routine for calculating finite-difference second directional derivative
        public static void FDAvv(int m, int n, Vector<double> x, Vector<double> v, Vector<double> fvec, Matrix<double> fjac, Func<int, int, Vector<double>, Vector<double>> func, ref Vector<double> acc, bool jac_uptodate, double h2)
        {
            Vector<double> xtmp;
            Vector<double> ftmp = Vector<double>.Build.Dense(m);

            if (jac_uptodate)
            {
                xtmp = x + h2 * v;
                ftmp = func(m, n, xtmp); // was func(m, n, xtmp, ftmp);
                acc = (2.0 / h2) * ((ftmp - fvec) / h2 - fjac.Multiply(v));
            }
            else // if jacobian not up to date, do not use jacobian in F.D. (needs one more function call)
            {
                xtmp = x + h2 * v;
                ftmp = func(m, n, xtmp); // was func(m, n, xtmp, ftmp);
                xtmp = x - h2 * v;
                acc = func(m, n, xtmp); // was func(m, n, xtmp, acc);
                acc = (ftmp - 2 * fvec + acc) / (h2 * h2);
            }
        }

        // done
        // Routines for updating lamda
        public static void TrustRegion(int n, int m, Vector<double> fvec, Matrix<double> fjac, Matrix<double> dtd, double delta, ref double lam)
        {
            // Calls dgqt supplied by minpack to calculate the step and Lagrange multiplier

            int i, itmax, info = -1;
            double rtol, atol, f = -1;

            Matrix<double> jtilde = Matrix<double>.Build.Dense(m, n);
            Matrix<double> g = Matrix<double>.Build.Dense(n, n);
            Vector<double> gradCtilde = Vector<double>.Build.Dense(n);            
            Vector<double> v = Vector<double>.Build.Dense(n);
            Vector<double> wa1 = Vector<double>.Build.Dense(n);
            Vector<double> wa2= Vector<double>.Build.Dense(n);            
            Vector<double> z = Vector<double>.Build.Dense(n);

            // Parameters for dgqt
            rtol = 1.0E-3;
            atol = 1.0E-3;
            itmax = 10;

            for (i = 1; i <= n; i++)
            {
                jtilde.SetColumn(i - 1, fjac.Row(i - 1) / Math.Sqrt(dtd[i - 1, i - 1])); // This assumes that dtd is diagonal...
            }
            gradCtilde = fvec * jtilde;
            g = jtilde.Transpose() * jtilde;
            dgqt(n, ref g, n, gradCtilde, delta, rtol, atol, itmax, ref lam, ref f, ref v, ref info, i, ref z, ref wa1, ref wa2);

            // Transform v back to non-dtd units
            // DO i = 1,n
            //      v(i) = v(i)/SQRT(dtd(i,i))
            // END DO

            return;
        }

        // done
        // traditional update methods
        public static void Updatelam_factor(ref double lam, int accepted, double factoraccept, double factorreject)
        {
            // Update lam based on accepted/rejected step

            lam = (accepted >= 0) ? lam / factoraccept : lam * factorreject;
        }

        // done
        public static void Updatelam_nelson(ref double lam, int accepted, double factoraccept, double factorreject, double rho)
        {
            // Update method due to Nelson [ref]
            
            if (accepted >= 0)
            {
                lam = lam * Math.Max(1.0 / factoraccept, 1.0 - (factorreject - 1.0) * (2.0 * rho - 1.0) * (2.0 * rho - 1.0) * (2.0 * rho - 1.0));
            }
            else
            {
                double nu = factorreject;
                for (int i = 2; i <= -1 * accepted; i++) //double nu for each rejection
                {
                    nu = nu * 2.0;
                }
                lam = lam * nu;
            }
        }

        // done
        public static void Updatelam_Umrigar(int m, int n, ref double lam, int accepted, Vector<double> v, Vector<double> vold,
            Vector<double> fvec, Matrix<double> fjac, Matrix<double> dtd, ref double a_param, double C, double Cnew)
        {
            // Method due to Umrigar and Nightingale [unpublished]

            double lamold, factor;
            double amemory, cos_on;

            amemory = Math.Exp(-1.0 / 5.0);
            cos_on = v.DotProduct(dtd * vold);
            cos_on = cos_on / Math.Sqrt(v.DotProduct(dtd * v) * vold.DotProduct(dtd * vold));

            if (accepted >= 0)
            {
                if (Cnew <= C)
                {
                    a_param = (cos_on > 0) ? amemory * a_param + 1.0 - amemory : amemory * a_param + 0.5 * (1.0 - amemory);
                }
                else
                {
                    a_param = amemory * a_param + 0.5 * (1.0 - amemory);
                }

                factor = Math.Min(100.0, Math.Max(1.1, 1.0 / Math.Pow(2.2E-16 + 1.0 - Math.Abs(2.0 * a_param - 1.0), 2)));
                if (Cnew <= C && cos_on >= 0)
                {
                    lam = lam / factor;
                }
                else if (Cnew > C)
                {
                    lam = lam * Math.Sqrt(factor);
                }
            }
            else
            {
                a_param = amemory * a_param;
                factor = Math.Min(100.0, Math.Max(1.1, 1.0 / Math.Pow(2.2E-16 + 1.0 - Math.Abs(2.0 * a_param - 1.0), 2)));
                lamold = lam;

                lam = (cos_on > 0) ? lam * Math.Sqrt(factor) : lam * factor;

                // Check for a 10% change in drift
                // Umrigar and Nightingal suggest a check that the the proposed change in lam actually produces a meaningful change in the step.
                // But this code produces strange results in a few cases.  -MKT
                // IF(accepted .EQ. -1) THEN
                //    d1 = SQRT(DOT_PRODUCT(v,v))
                //    g = MATMUL(TRANSPOSE(fjac),fjac) + lam*dtd
                //    CALL DPOTRF('U', n, g, n, info)
                //    grad = MATMUL(fvec, fjac)
                //    CALL DPOTRS('U', n, 1, g, n, grad, n, info)
                //    d2 = SQRT(DOT_PRODUCT(grad, grad) )
                //    IF( 10.0d+0*ABS(d2-d1) .LT. d2 ) lam = lam - 0.1*d2*(lamold - lam)/(d1 - d2)
                //  END IF
            }
        }

        // done
        // Trust region update methods
        public static void Updatedelta_factor(ref double delta, int accepted, double factoraccept, double factorreject)
        {
            // Update lam based on accepted/rejected step

            delta = (accepted >= 0) ? delta * factoraccept : delta / factorreject;
        }

        // done
        public static void Updatedelta_more(ref double delta, ref double lam, int n, Vector<double> v, Matrix<double> dtd, double rho, double C, double Cnew,
            double dirder, double actred, double av, double avmax)
        {
            double pnorm = Math.Sqrt(v.DotProduct(dtd * v));
            double temp;

            if (rho > 0.25)
            {
                temp = (lam > 0.0 && rho < 0.75) ? 1.0 : 2.0 * pnorm / delta;
            }
            else
            {
                temp = (actred >= 0.0) ? 0.5 : 0.5 * dirder / (dirder + 0.5 * actred);
                if (0.01 * Cnew >= C || temp < 0.1)
                {
                    temp = 0.1;
                }
            }
            // We need to make sure that if acceleration is too big, we decrease the step size
            if (av > avmax)
            {
                temp = Math.Min(temp, Math.Max(avmax / av, 0.1));
            }
            
            delta = temp * Math.Min(delta, 10.0 * pnorm);
            lam = lam / temp;
        }

        public static string ConvergedInfo(int i)
        {
            switch (i)
            {
                case 1: return "artol reached";
                case 2: return "Cgoal reached";
                case 3: return "gtol reached";
                case 4: return "xtol reached";
                case 5: return "xrtol reached";
                case 6: return "ftol reached";
                case 7: return "frtol reached";
                case -1: return "maxiters exeeded";
                case -2: return "maxfev exceeded";
                case -3: return "maxjev exceeded";
                case -4: return "maxaev exceeded";
                case -10: return "User Termination";
                case -11: return "NaN Produced";
            }

            return "Not identified";            
        }

        // done
        // Main Geodesic-Bold-BroydenUpdate-Levenberg-Marquardt routine
        // version 1.0.2
        //
        public static void geodesiclm(Func<int, int, Vector<double>, Vector<double>> func, Func<int, int, Vector<double>, Matrix<double>> jacobian,
            Func<int, int, Vector<double>, Vector<double>, Vector<double>> Avv,
            ref Vector<double> x, ref Vector<double> fvec, ref Matrix<double> fjac, int n, int m,
            Func<int, int, Vector<double>, Vector<double>, Vector<double>, Vector<double>, Matrix<double>, Vector<double>, double, Matrix<double>, Vector<double>, int, int, int> callback,
            ref int info,
            bool analytic_jac, bool analytic_Avv,
            bool center_diff, double h1, double h2,
            ref Matrix<double> dtd, int damp_mode,
            ref int niters, ref int nfev, ref int njev, ref int naev,
            ref int maxiter, int maxfev, int maxjev, int maxaev, double maxlam, double minlam,
            double artol, double Cgoal, double gtol, double xtol, double xrtol, double ftol, double frtol,
            ref int converged,
            int print_level, int print_unit,
            int imethod, int iaccel, int ibold, int ibroyden,
            double initialfactor, double factoraccept, double factorreject, double avmax)
        {
            #region description

            //    subroutine geodesicLM
            //    
            //    The purpose of geolevmar is to minimize the sum of the squares
            //    of m nonlinear functions of n variables by a modification of
            //    the Levenberg-Marquardt algorithm that utilizes the geodesic
            //    acceleration step correction, bold acceptance criterion, and
            //    a Broyden update of the jacobian matrix.  The method employs one
            //    of several possible schemes for updating the Levenberg-Marquardt
            //    parameter.  The user must provide a subroutine which calcualtes
            //    the functions, and optionally the jacobian and a directional 
            //    derivative of the functions.  The latter two will be estimated
            //    by finite differences if not supplied.
            //
            //    If you use this code, please acknowledge such by referencing one
            //    one of the following papers in any published work:
            //    
            //    Transtrum M.K., Machta B.B., and Sethna J.P, Why are nonlinear
            //    fits to data so challenging?  Phys. Rev. Lett. 104, 060201 (2010)
            //
            //    Transtrum M.K., Machta B.B., and Sethna J.P., The geometry of
            //    nonlinear least squares with applications to sloppy model and
            //    optimization.  Phys. Rev. E. 80, 036701 (2011)
            //
            //
            //    The subroutine statement is:
            //
            //    geodesicLM(func, jacobian, Avv, x, fvec, fjac, n, m, callback, info,
            //              analytic_jac, analytic_Avv, center_diff, h1, h2,
            //              dtd, damp_mode, niteres, nfev, njev, naev,
            //              maxiters, maxfev, maxjev, maxaev, maxlam, minlam,
            //              artol, Cgoal, gtol, xtol, xrtol, ftol, frtol,
            //              converged, print_level, print_unit,
            //              imethod, iaccel, ibold, ibroyden,
            //              initialfactor, factoraccept, factorreject, avmax)
            //
            //    where
            //
            //    func is a user supplied subroutine which calculates the functions and
            //    should be written as follows:
            //
            //      subroutine func(m, n, x, fvec)
            //      integer m, n
            //      double precision x(n), fvec(m)
            //      --------------------------------------------------------------------
            //      calculates the function at x and returns their values in fvec
            //      x, m, and n should be left unchanged
            //      --------------------------------------------------------------------
            //      end subroutine func
            //
            //    jacobian is a user supplied subroutine which calculates the jacobian of
            //    of the functions if analytic_jac is .TRUE.
            //    jacobian should be writen as follows     
            //
            //      subroutine jacobian(m, n, x, fjac)
            //      integer m, n
            //      double precision x(n), fjac(m,n)
            //      --------------------------------------------------------------------
            //      calculates the jacobian at x and returns their values in fjac
            //      x, m, and n should be left unchanged
            //      --------------------------------------------------------------------
            //      end subroutine jacobian
            //
            //    Avv is a user supplied subroutine which calculates the directional
            //    second derivative of the functions if analytic_Avv is .TRUE.
            //    Avv should be writen as follows     
            //
            //      subroutine Avv(m, n, x, v, acc)
            //      integer m, n
            //      double precision x(n), v(n), acc(m)
            //      --------------------------------------------------------------------
            //      calculates the directional second derivative at x in the direction 
            //      of v and returns the values in acc
            //      x, v, m, and n should be left unchanged
            //      --------------------------------------------------------------------
            //      end subroutine Avv
            //
            //    x is an array of length n.  On input it contains an initial estimate of
            //    the solution.  On exit, it contains the final estimate of the solution.
            //
            //    fvec is an output array of length m containing the funtion evaluation at
            //    the final solution
            //
            //    fjac is an output array of dimension(m,n) containing the jacobian evaluation
            //    the final solution.  The array MATMUL( TRANSPOSE(fjac), fjac) is an estimate
            //    of the covariance matrix of parameters at the final solution
            //
            //    n an input integer set to the number of parameters
            //
            //    m an input integer set to the number of functions
            //
            //    callback a user supplied subroutine to be called after each iteration of the
            //    algorithm.  
            //    callback should be written as follows
            //
            //      subroutine callback(m,n,x,v,a,fvec,fjac,acc,lam,dtd,fvec_new,accepted,info)
            //      integer m, n, accepted, info
            //      double precision x(n), v(n), a(n), fvec(m), fjac(m,n), acc(m), lam, dtd(n,n), fvec_new(m)
            //      --------------------------------------------------------------------
            //      m, n, x, v, a, fvec, fjac, acc, lam, dtd, fvec_new, accepted, should be left unchanged
            //      On input, info = 0 and should be changed to a nonzero value if the user
            //      wishes to terminate calculation
            //      --------------------------------------------------------------------
            //      end subroutine callback
            //
            //    info an output integer set to a nonzero value if the user terminated the routine
            //    (see callback).
            //
            //    analytic_jac an input boolean set to .TRUE. if the subroutine jacobian calculates
            //    the jacobian.  If .FALSE. then a finite difference estimate will be used.
            //
            //    analytic_Avv an input boolean set to .TRUE. if the subroutine Avv calculates
            //    the directional second derivative.  If .FALSE. then a finite difference estimate
            //    will be used.
            //
            //    center_diff an input boolean.  If finite differences are used to estimate the jacobian
            //    then center differences will used if center_diff is .TRUE., otherwise, forward
            //    differences will be used.  Note that center differences are more accurate by require
            //    more function evaluations.
            //
            //    h1 an input double precision specifying the step size for the finite difference estimates
            //    of the jacobian.
            //
            //    h2 an input double precision specifying the steps ize for the finite difference estiamtes
            //    of the directional second derivative.
            //
            //    dtd a double precision array of dimension(n,n).  dtd is used as the damping matrix in the 
            //    Levenberg-Marquardt routine.  It's exact treatment is specified by the damp_mode input.
            //
            //    damp_mode an input integer specifying the details of the LM damping as follows:
            //      damp_mode = 0: dtd is set to the identity.
            //      damp_mode = 1: dtd should be a positive definite, diagonal matrix whose entries are dynamically
            //                updated based on the elements of the jacobian.
            //
            //    niters an output integer specifying the number of iterations of the algorithm.
            //
            //    nfev an output integer specifying the number of calls to func.  
            //
            //    njev an output integer specifying the number of calls to jacobian.
            //
            //    naev an output integer specifying the number of calls to Avv.
            //
            //    maxiter an input integer specifying the maximum number of routine iterations.
            //
            //    maxfev an input integer specifying the maximum number of function calls
            //    if maxfev = 0, then there is no limit to the number of function calls.
            //
            //    maxjev an input integer specifying the maximum number of jacobian calls
            //    if maxjev = 0, then there is no limit to the number of jacobian calls.
            //
            //    maxaev an input integer specifying the maximum number of Avv calls
            //    if maxaev = 0, then there is no limit to the number of Avv calls.
            //
            //    maxlam an input double precision specifying the maximum allowed value of 
            //    the damping term lambda. If this is negative, then there is no limit.
            //
            //    minlam an input double precision specifying the minimum allowed value of 
            //    the damping term lambda. If lambda is smaller than this value for three consecutive steps
            //    the routine terminates.  If this is negative, then there is no limit.
            //
            //    artol an input double precision.  The method will terminate when the cosine of the
            //    angle between the residual vector and the range of the jacobian is less than artol.
            //
            //    Cgoal an input double precision.  The method will terminate when the cost (one half
            //    the sum of squares of the function) falls below Cgoal.
            //
            //    gtol an input double precision.  The method will terminate when norm of Cost gradient 
            //    falls below gtol.
            //    
            //    xtol an input double precision.  The method will terminate when parameters change by
            //    less than xtol.
            //
            //    xrtol an input double precision.  The method will terminate if the relative change in
            //    each of the parameters is less than xrtol.
            //
            //    ftol an input double precision.  The method will termiante if the Cost fails to decrease
            //    by more than ftol for 3 consecutive iterations.
            //
            //    frtol an input double precision.  The method will terminate if the relative decrease in
            //    Cost is less than frtol 3 consecutive iterations.
            //
            //    converged an output integer indicated the reason for termination:
            //      converged = 1: artol 
            //      converged = 2: Cgoal
            //      converged = 3: gtol
            //      converged = 4: xtol
            //      converged = 5: xrtol
            //      converged = 6: ftol
            //      converged = 7: frtol
            //      converged = -1: maxiters exeeded
            //      converged = -2: maxfev exceeded
            //      converged = -3: maxjev exceeded
            //      converged = -4: maxaev exceeded
            //      converged = -10: user requested termination in callback via info
            //      converged = -11: Either the initial function evalaution or subsequent jacobian
            //                       evaluations produced Nans.
            //
            //    print_level an input integer specifying the amount of details to be printed.
            //    acceptable values range from 0 to 5, with larger number printing more details.
            //
            //    print_unit an input integer specifying the unit number details should be written to.
            //
            //    imethod an input integer specifying the method for updating the LM parameter
            //      imethod = 0: adjusted by fixed factors after accepted/rejected steps
            //      imethod = 1: adjusted as described in Nielson
            //      imethod = 2: adjusted according to an unpublished method due to Cyrus Umrigar and Peter Nightingal
            //      imethod = 10: step size Delta adjusted by fixed factors after accepted/rejected steps
            //      imethod = 11: step size adjusted as described in More'
            //
            //    initialfactor an input double precision for specifying either the initial LM parameter
            //    of the initial step size.
            //
            //    factoraccept an input double precision (larger than 1.0) specifying the factor by which
            //    either the LM parameter or the step size will be adjusted after an accepted step if
            //    imethod = 0 or 10
            //
            //    factorreject an input double precision (larger than 1.0) specifying the factor by which
            //    either the LM parameter of the step size will be adjusted after a rejected step if
            //    imethod = 0 or 10
            //
            //    avmax an input double precision specifying the maximum norm of the geodesic acceleration 
            //    relative to the velocity vector.
            //
            #endregion description
                        
            double lam = 0, delta = 0;
            double Cold = 0;
            Vector<double> acc = Vector<double>.Build.Dense(m);
            Vector<double> v = Vector<double>.Build.Dense(n);
            Vector<double> vold = Vector<double>.Build.Dense(n);
            Vector<double> a = Vector<double>.Build.Dense(n);
            Vector<double> x_new = Vector<double>.Build.Dense(m);
            Vector<double> fvec_new = Vector<double>.Build.Dense(m);           

            Matrix<double> g;
            double temp1, temp2, pred_red, dirder = 0, actred = 0, rho = 0;
            int i, j, istep;
             
            //string converged_info;
            bool jac_uptodate, jac_force_update;

            // strings for concluding print statement
            // converged_info = '????????'
            // converged_info(1) = 'artol reached'
            // converged_info(2) = 'Cgoal reached'
            // converged_info(3) = 'gtol reached'
            // converged_info(4) = 'xtol reached'
            // converged_info(5) = 'xrtol reached'
            // converged_info(6) = 'ftol reached'
            // converged_info(7) = 'frtol reached'
            // converged_info(-1) = 'maxiters exeeded'
            // converged_info(-2) = 'maxfev exceeded'
            // converged_info(-3) = 'maxjev exceeded'
            // converged_info(-4) = 'maxaev exceeded'
            // converged_info(-10) = 'User Termination '
            // converged_info(-11) = 'NaN Produced'
            
            if (print_level >= 1)
            {
                Console.WriteLine("Optimizing with Geodesic-Levenberg-Marquardt algorithm, version 1.0.2");
                Console.WriteLine("Method Details:");
                Console.WriteLine("  Update method:   " + imethod);
                Console.WriteLine("  acceleration:    " + iaccel);
                Console.WriteLine("  Bold method:     " + ibold);
                Console.WriteLine("  Broyden updates: " + ibroyden);
            }
            
            // Initialize variables
            niters = 0;
            nfev = 0;
            naev = 0;
            njev = 0;
            converged = 0;

            double cos_alpha = 1.0;
            double av = 0.0;
            double a_param = 0.5;

            int accepted = 0;
            int counter = 0;

            fvec = func(m, n, x); // was func(m, n, x, fvec);

            nfev = nfev + 1;
            double C = 0.5 * fvec.DotProduct(fvec);
            double Cnew = 0.0;

            if (print_level >= 1)
            {
                Console.WriteLine("  Initial Cost:    " + C);
            }

            // Check for nans in fvec
            bool valid_result = !fvec.Select(item => double.IsNaN(item)).Contains(true);

            if (!valid_result)
            {
                converged = -11;
                maxiter = 0;
            }
            double Cbest = C;
            var fvec_best = fvec;
            var x_best = x;
            if (analytic_jac)
            {
                fjac = jacobian(m, n, x); // was jacobian(m, n, x, fjac);
                njev = njev + 1;
            }
            else
            {
                fdjac(m, n, x, fvec, ref fjac, func, h1, center_diff);
                nfev = (center_diff) ? nfev + 2 * n : nfev + n;
            }
            jac_uptodate = true;
            jac_force_update = false;
            var jtj = (fjac.Transpose()) * fjac;

            // Check fjac for nans
            valid_result = !fjac.ToColumnArrays().Select(item1 => item1.Select(item2 => double.IsNaN(item2)).Contains(true)).Contains(true);

            if (!valid_result)
            {
                converged = -11;
                maxiter = 0;
            }

            acc.Clear();
            a.Clear();

            // Initialize scaling matrix
            if (damp_mode == 0)
            {
                dtd.Clear();
                for (i = 1; i <= n; i++)
                {
                    dtd[i - 1, i - 1] = 1.0;
                }
            }
            else if (damp_mode == 1)
            {
                for (i = 1; i <= n; i++)
                {
                    dtd[i - 1, i - 1] = Math.Max(jtj[i - 1, i - 1], dtd[i - 1, i - 1]);
                }
            }

            // Initialize lambda
            if (imethod < 10)
            {
                lam = jtj[0, 0];
                for (i = 2; i <= n; i++)
                {
                    lam = Math.Max(jtj[i - 1, i - 1], lam);
                }
                lam = lam * initialfactor;
            }
            // Initialize step bound if using trust region method
            else if (imethod >= 10)
            {
                delta = initialfactor * Math.Sqrt(x.DotProduct(dtd * x));
                lam = 1.0;
                if (delta == 0.0)
                {
                    delta = 100;
                }
                if (converged == 0)
                {
                    TrustRegion(n, m, fvec, fjac, dtd, delta, ref lam); // Do not call this if there were nans in either fvec or fjac
                }
            }

            // Main Loop
            for (istep = 1; istep <= maxiter; istep++)
            {
                info = 0;
                callback(m, n, x, v, a, fvec, fjac, acc, lam, dtd, fvec_new, accepted, info);
                if (info != 0)
                {
                    converged = -10;
                    break;
                }
                // Update Functions
                //
                // Full or partial Jacobian Update?
                if (accepted > 0 && ibroyden <= 0) jac_force_update = true;
                if (accepted + ibroyden <= 0 && !jac_uptodate) jac_force_update = true; // Force jac update after too many failed attempts

                if (accepted > 0 && ibroyden > 0 && !jac_force_update) // Rank deficient update of jacobian matrix
                {
                    updatejac(m, n, ref fjac, fvec, fvec_new, acc, v, a);
                    jac_uptodate = false;
                }

                if (accepted > 0) // Accepted step
                {
                    fvec = fvec_new;
                    x = x_new;
                    vold = v;
                    C = Cnew;
                    if (C <= Cbest)
                    {
                        x_best = x;
                        Cbest = C;
                        fvec_best = fvec;
                    }
                }

                if (jac_force_update) // Full rank update of jacobian
                {
                    if (analytic_jac)
                    {
                        fjac = jacobian(m, n, x); // was jacobian(m, n, x, fjac);
                        njev = njev + 1;
                    }
                    else
                    {
                        fdjac(m, n, x, fvec, ref fjac, func, h1, center_diff);
                        nfev = (center_diff) ? nfev + 2 * n : nfev + n;
                    }
                    jac_uptodate = true;
                    jac_force_update = false;
                }

                // Check fjac for nans
                valid_result = !fjac.ToColumnArrays().Select(item1 => item1.Select(item2 => double.IsNaN(item2)).Contains(true)).Contains(true);

                if (valid_result) // If no nans in jacobian
                {
                    jtj = fjac.Transpose() * fjac;

                    // Update Scaling/lam/TrustRegion
                    if (istep > 1) // Only necessary after first step
                    {
                        if (damp_mode == 1)
                        {
                            for (i = 1; i <= n; i++)
                            {
                                dtd[i - 1, i - 1] = Math.Max(jtj[i - 1, i - 1], dtd[i - 1, i - 1]);
                            }
                        }
                        // Can add other lam-delta update methods 
                        switch (imethod)
                        {
                            case 0:
                                // Update lam directly by fixed factors
                                Updatelam_factor(ref lam, accepted, factoraccept, factorreject);
                                break;
                            case 1:
                                //             !! Update lam directly based on Gain Factor rho (see Nielson reference)
                                Updatelam_nelson(ref lam, accepted, factoraccept, factorreject, rho);
                                break;
                            case 2:
                                // Update lam directly using method of Umrigar and Nightingale [unpublished]
                                Updatelam_Umrigar(m, n, ref lam, accepted, v, vold, fvec, fjac, dtd, ref a_param, Cold, Cnew);
                                break;
                            case 10:
                                // Update delta by fixed factors
                                Updatedelta_factor(ref delta, accepted, factoraccept, factorreject);
                                TrustRegion(n, m, fvec, fjac, dtd, delta, ref lam);
                                break;
                            case 11:
                                // Update delta as described in More' reference
                                Updatedelta_more(ref delta, ref lam, n, v, dtd, rho, C, Cnew, dirder, actred, av, avmax);
                                TrustRegion(n, m, fvec, fjac, dtd, delta, ref lam);
                                break;
                        }
                    }

                    // Propose Step
                    // metric aray
                    g = jtj + lam * dtd;
                    // Cholesky decomposition
                    try
                    {
                        var chol = g.Cholesky();
                        g = chol.Factor.Transpose();// was DPOTRF('U', n, g, n, info);
                        info = 0;
                    }
                    catch
                    {
                        info = 1;
                    }
                }
                else
                {
                    // If nans in jacobian
                    converged = -11;
                    break;
                }

                if (info == 0)  // If matrix decomposition successful:
                {
                    // v = -1.0 * g * (fvec * fjac); //velocity
                    v = -1.0 * (fvec * fjac);
                    v = g.Solve(v); // was DPOTRS('U', n, 1, g, n, v, n, info);

                    // Calcualte the predicted reduction and the directional derivative -- useful for updating lam methods
                    temp1 = 0.5 * v.DotProduct(jtj * v) / C;
                    temp2 = 0.5 * lam * v.DotProduct(dtd * v) / C;
                    pred_red = temp1 + 2.0 * temp2;
                    dirder = -1.0 * (temp1 + temp2);
                    // calculate cos_alpha -- cos of angle between step direction (in data space) and residual vector
                    cos_alpha = Math.Abs(fvec.DotProduct(fjac * v));
                    cos_alpha = cos_alpha / Math.Sqrt(fvec.DotProduct(fvec) * (fjac * v).DotProduct(fjac * v));
                    if (imethod <= 10)
                    {
                        delta = Math.Sqrt(v.DotProduct(dtd * v)); // Update delta if not set directly
                    }
                    // update acceleration
                    if (iaccel > 0)
                    {
                        if (analytic_Avv)
                        {
                            acc = Avv(m, n, x, v); // was Avv(m, n, x, v, acc);
                            naev = naev + 1;
                        }
                        else
                        {
                            FDAvv(m, n, x, v, fvec, fjac, func, ref acc, jac_uptodate, h2);
                            if (jac_uptodate)
                            {
                                nfev = nfev + 1;
                            }
                            else
                            {
                                nfev = nfev + 2; // we don't use the jacobian if it is not up to date
                            }
                        }
                        // Check accel for nans
                        valid_result = !acc.Select(item => double.IsNaN(item)).Contains(true);

                        if (valid_result)
                        {
                            a = -1.0 * (acc * fjac);
                            a = g.Solve(a); // was DPOTRS('U', n, 1, g, n, a, n, info);
                            // a = -1.0d+0*MATMUL(g,MATMUL(acc,fjac))
                        }
                        else
                        {
                            a.Clear(); // If nans in acc, we will ignore the acceleration term
                        }
                    }

                    // Evaluate at proposed step -- only necessary if av <= avmax
                    av = Math.Sqrt(a.DotProduct(dtd * a) / v.DotProduct(dtd * v));
                    if (av <= avmax)
                    {
                        x_new = x + v + 0.5 * a;
                        fvec_new = func(m, n, x_new); // was func(m, n, x_new, fvec_new);
                        nfev = nfev + 1;
                        Cnew = 0.5 * fvec_new.DotProduct(fvec_new);
                        Cold = C;

                        // Check for nans in fvec_new
                        valid_result = !fvec_new.Select(item => double.IsNaN(item)).Contains(true);

                        if (valid_result) // If no nans, proceed as normal
                        {
                            // update rho and actred
                            actred = 1.0 - Cnew / C;
                            rho = 0.0;
                            if (pred_red != 0.0) rho = (1.0 - Cnew / C) / pred_red;
                            // Accept or Reject proposed step
                            Acceptance(n, C, Cnew, Cbest, ibold, ref accepted, dtd, v, vold);
                        }
                        else
                        {
                            // If nans in fvec_new, reject step
                            actred = 0.0;
                            rho = 0.0;
                            accepted = Math.Min(accepted - 1, -1);
                        }
                    }
                    else
                    {
                        // If acceleration too large, then reject
                        accepted = Math.Min(accepted - 1, -1);
                    }
                }
                else
                {
                    // If matrix factorization fails, reject the proposed step
                    accepted = Math.Min(accepted - 1, -1);
                }

                // Check Convergence
                if (converged == 0)
                {
                    convergence_check(m, n, ref converged, accepted, ref counter,
                        C, Cnew, x, fvec, fjac, lam, x_new,
                        nfev, maxfev, njev, maxjev, naev, maxaev, maxlam, minlam,
                        artol, Cgoal, gtol, xtol, xrtol, ftol, frtol, cos_alpha);
                    if (converged == 1 && !jac_uptodate)
                    {
                        // If converged by artol with an out of date jacobian, update the jacoban to confirm true convergence
                        converged = 0;
                        jac_force_update = true;
                    }
                }

                // Print status
                if (print_level == 2 && accepted > 0)
                {
                    Console.WriteLine(string.Format("[{0}] nfev = {1}, njev = {2}, naev = {3}, accepted = {4}", istep, nfev, njev, naev, accepted));
                    Console.WriteLine(string.Format("  Cost = {0}, lam = {1}, delta = {2}", C, lam, delta));
                    Console.WriteLine(string.Format("  av = {0}, cos alpha = {1}", av, cos_alpha));
                }
                else if (print_level == 3)
                {
                    Console.WriteLine(string.Format("[{0}] nfev = {1}, njev = {2}, naev = {3}, accepted = {4}", istep, nfev, njev, naev, accepted));
                    Console.WriteLine(string.Format("  Cost = {0}, lam = {1}, delta = {2}", C, lam, delta));
                    Console.WriteLine(string.Format("  av = {0}, cos alpha = {1}", av, cos_alpha));
                }
                if (print_level == 4 && accepted > 0)
                {
                    Console.WriteLine(string.Format("[{0}] nfev = {1}, njev = {2}, naev = {3}, accepted = {4}", istep, nfev, njev, naev, accepted));
                    Console.WriteLine(string.Format("  Cost = {0}, lam = {1}, delta = {2}", C, lam, delta));
                    Console.WriteLine(string.Format("  av = {0}, cos alpha = {1}", av, cos_alpha));
                    Console.WriteLine(string.Format("  x = {0}\t {1}", x[0], x[1]));
                    //       WRITE(print_unit, *) "  v = ", v
                    //       WRITE(print_unit, *) "  a = ", a
                }
                else if (print_level == 5)
                {
                    Console.WriteLine(string.Format("[{0}] nfev = {1}, njev = {2}, naev = {3}, accepted = {4}", istep, nfev, njev, naev, accepted));
                    Console.WriteLine(string.Format("  Cost = {0}, lam = {1}, delta = {2}", C, lam, delta));
                    Console.WriteLine(string.Format("  av = {0}, cos alpha = {1}", av, cos_alpha));
                    Console.WriteLine(string.Format("  x = {0}\t {1}", x[0], x[1]));
                    //       WRITE(print_unit, *) "  v = ", v
                    //       WRITE(print_unit, *) "  a = ", a
                }

                // if converged -- return
                if (converged != 0)
                {
                    break;
                }
                //
                if (accepted >= 0) jac_uptodate = false; // jacobian is now out of date
            }

            // ! If not converged
            if (converged == 0) converged = -1;
            niters = istep;
            //
            // Return best fit found
            // If the method converged, but final x is different from x_best -- what to do?
            x = x_best;
            fvec = fvec_best;

            if (print_level >= 1)
            {
                Console.WriteLine("Optimization finished");
                Console.WriteLine("Results:");
                Console.WriteLine("  Converged:    " + converged + ", " + ConvergedInfo(converged));
                Console.WriteLine("  Final Cost: ", 0.5 * fvec.DotProduct(fvec));
                Console.WriteLine("  Cost/DOF: " + 0.5 * fvec.DotProduct(fvec) / (m - n));
                Console.WriteLine("  niters:     " + istep);
                Console.WriteLine("  nfev:       " + nfev);
                Console.WriteLine("  njev:       " + njev);
                Console.WriteLine("  naev:       " + naev);
            }
        }


    }
}
