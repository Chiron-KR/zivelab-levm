using Levm.GeodesicLM;
using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Levm.GeodesicLM.Test
{
    class Program
    {
        #region Define Rosenbrock

        private static Vector<double> Rosenbrock(int m, int n, Vector<double> x)
        {
            Vector<double> fvec = Vector<double>.Build.Dense(2);

            fvec[0] = ((1.0 - x[0]) * (1.0 - x[0]) + 105.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]));
            fvec[1] = ((1.0 - x[0]) * (1.0 - x[0]) + 105.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]));

            return fvec;
        }

        private static Matrix<double> RosenbrockJacobian(int m, int n, Vector<double> x)
        {

            Matrix<double> fjac = Matrix<double>.Build.Dense(2, 2);

            fjac[0, 0] = (-2.0 + 2.0 * x[0] - 4.0 * 105.0 * (x[1] - x[0] * x[0]) * x[0]); ;
            fjac[0, 1] = (2.0 * 105.0 * (x[1] - x[0] * x[0]));
            fjac[1, 0] = (-2.0 + 2.0 * x[0] - 4.0 * 105.0 * (x[1] - x[0] * x[0]) * x[0]); ;
            fjac[1, 1] = (2.0 * 105.0 * (x[1] - x[0] * x[0]));

            return fjac;
        }

        private static Vector<double> RosenbrockAvv(int m, int n, Vector<double> x, Vector<double> v)
        {
            Vector<double> acc = Vector<double>.Build.Dense(2);

            acc[0] = 0;
            acc[1] = -200 * v[0] * v[0];

            return acc;
        }

        #endregion

        #region Define Norris

        // https://www.itl.nist.gov/div898/strd/lls/data/Norris.shtml

        private static double[] norris_p = { -0.262323073774029, 1.00211681802045 }; // best parameters
        private static double[] norris_y = {
            0.1,    338.8,  118.1,  888,    9.2,    228.1,  668.5,  998.5,  449.1,  778.9,
            559.2,  0.3,    0.1,    778.1,  668.8,  339.3,  448.9,  10.8,   557.7,  228.3,
            998,    888.8,  119.6,  0.3,    0.6,    557.6,  339.3,  888,    998.5,  778.9,
            10.2,   117.6,  228.9,  668.4,  449.2,  0.2
        };
        private static double[] norris_x = {
            0.2,    337.4,  118.2,  884.6,  10.1,   226.5,  666.3,  996.3,  448.6,  777,
            558.2,  0.4,    0.6,    775.5,  666.9,  338,    447.5,  11.6,   556,    228.1,
            995.8,  887.6,  120.2,  0.3,    0.3,    556.8,  339.1,  887.2,  999,    779,
            11.1,   118.3,  229.2,  669.1,  448.9,  0.5
        };
        private static Vector<double> norris(int m, int n, Vector<double> p)
        {
            int i;
            double xi, yi;

            Vector<double> y = Vector<double>.Build.Dense(m);

            for (i = 0; i < m; i++)
            {
                xi = norris_x[i];
                yi = norris_y[i];

                y[i] = (p[0] + p[1] * xi) - yi;
            }

            return y;
        }
        private static Matrix<double> norrisPrime(int m, int n, Vector<double> p)
        {
            int i, j;
            double x;

            Matrix<double> jac = Matrix<double>.Build.Dense(m, n);
            for (i = 0; i < m; i++)
            {
                x = norris_x[i];

                jac[i, 0] = 1.0;
                jac[i, 1] = x;
            }

            return jac;
        }

        # endregion

        private static int callback(int m, int n, Vector<double> x, Vector<double> v, Vector<double> a,
            Vector<double> fvec, Matrix<double> fjac, Vector<double> acc, double lam,
            Matrix<double> dtd, Vector<double> fvec_new, int accepted, int info)
        {
            //implicit none
            //integer m, n, accepted, info
            //double precision x(n), v(n), a(n), fvec(m), fjac(m, n), acc(m), lam, dtd(n, n)
            //double precision fvec_new(m)

            //StringBuilder sb = new StringBuilder();

            //for (int i = 0; i < m; ++i)
            //    sb.AppendFormat("{0, 12:0.0000000}", x[i]);

            //Console.WriteLine(sb);

            return 0;
        }

        static void Main(string[] args)
        {
            int m, n, info = 0, damp_mode, niters = 0, nfev = 0, njev = 0, naev = 0;
            int maxiter, maxfev, maxjev, maxaev, converged = 0, print_level, print_unit;
            int imethod, iaccel, ibold, ibroyden;

            double h1, h2;
            double maxlam, minlam;
            double artol, Cgoal, gtol, xtol, xrtol, ftol, frtol;
            double initialfactor, factoraccept, factorreject, avmax;

            bool analytic_jac, analytic_Avv, center_diff;

            Vector<double> x = Vector<double>.Build.Dense(2);
            Vector<double> fvec = Vector<double>.Build.Dense(2);
            Matrix<double> fjac = Matrix<double>.Build.Dense(2, 2);

            // dimensions of problem
            m = 2;
            n = 2;

            // starting point
            x = Vector<double>.Build.DenseOfArray(new double[] { -1.2, 1 });

            analytic_jac = true;
            analytic_Avv = true;

            //Irrelevant since we use analytic derivatives
            center_diff = false;
            h1 = 1.0E-8;
            h2 = 1.0E-1;

            // damping matrix
            Matrix<double> dtd = Matrix<double>.Build.DenseOfArray(new double[,] { { 1, 0 }, { 0, 1 } });

            // dynamically update dtd
            damp_mode = 1;

            // stopping criterion
            maxiter = 10000;
            maxfev = 0;
            maxjev = 0;
            maxaev = 0;
            maxlam = 1E7;
            minlam = -1.0;
            // this is the preferred stopping criterion for most problems in which m > n

            // I recommend using 0.001.
            // since m = n in this case we do not use it(indicated by a negative value)
            artol = -1.0;
            // instead, we only use the Cgoal as a stopping criterion
            Cgoal = 1.0E-8;
            gtol = -1.0;
            xtol = -1.0;
            xrtol = -1.0;
            ftol = -1.0;
            frtol = -1.0;

            // setup printing of algorithm status
            print_level = 5;
            print_unit = 6;

            // flags to determine algorithm
            imethod = 0;// 11;
            iaccel = 0;
            ibold = 0;
            ibroyden = 0;
            initialfactor = 1; // initial lamda
            factoraccept = 5; // 10 for standard LM
            factorreject = 2; // 10 for standard LM
            avmax = 0.75;

            //GeoDesicLM.geodesiclm(Rosenbrock, RosenbrockJacobian, RosenbrockAvv,
            //    ref x, ref fvec, ref fjac, n, m,
            //    callback,
            //    ref info,
            //    analytic_jac, analytic_Avv,
            //    center_diff, h1, h2, ref dtd, damp_mode,
            //    ref niters, ref nfev, ref njev, ref naev,
            //    ref maxiter, maxfev, maxjev, maxaev, maxlam, minlam,
            //    artol, Cgoal, gtol, xtol, xrtol, ftol, frtol,
            //    ref converged, print_level, print_unit,
            //    imethod, iaccel, ibold, ibroyden,
            //    initialfactor, factoraccept, factorreject, avmax);

            //Console.WriteLine("To exit applicaiton, type anykey.");
            //Console.ReadKey(true);

            // this took 266 iterations and 104 jacobian evaluations for me.

            // now turn on geodesic acceleration and repeat

            //iaccel = 1;

            //x = Vector<double>.Build.DenseOfArray(new double[] { -1, 1 });

            //GeoDesicLM.geodesiclm(Rosenbrock, RosenbrockJacobian, RosenbrockAvv,
            //    x, ref fvec, ref fjac, n, m,
            //    callback,
            //    ref info,
            //    analytic_jac, analytic_Avv,
            //    center_diff, h1, h2, ref dtd, damp_mode,
            //    ref niters, ref nfev, ref njev, ref naev,
            //    ref maxiter, maxfev, maxjev, maxaev, maxlam, minlam,
            //    artol, Cgoal, gtol, xtol, xrtol, ftol, frtol,
            //    ref converged, print_level, print_unit,
            //    imethod, iaccel, ibold, ibroyden,
            //    initialfactor, factoraccept, factorreject, avmax);

            //this took 2 iterations and 2 jacobian evaluations for me.

            m = 36; // number of data
            n = 2; // number of parameters

            x = Vector<double>.Build.DenseOfArray(new double[] { 1, 1 });
            analytic_jac = true; analytic_Avv = false;

            fvec = Vector<double>.Build.Dense(m);
            fjac = Matrix<double>.Build.Dense(m, n);
            dtd = Matrix<double>.Build.DenseDiagonal(n, 1.0);

            // dynamically update dtd
            damp_mode = 1;

            // stopping criterion
            maxiter = 100000;
            maxfev = 0;
            maxjev = 0;
            maxaev = 0;
            maxlam = 1E8;
            minlam = 1E-8;

            // tolerances
            artol = 1.0E-8;
            Cgoal = 1.0E-8;
            gtol = 1.0E-8;
            xtol = 1.0E-8;
            xrtol = 1.0E-8;
            ftol = 1.0E-8;
            frtol = 1.0E-8;

            // flags to determine algorithm
            //      imethod = 0: adjusted by fixed factors after accepted/rejected steps
            //      imethod = 1: adjusted as described in Nielson
            //      imethod = 2: adjusted according to an unpublished method due to Cyrus Umrigar and Peter Nightingal
            //      imethod = 10: step size Delta adjusted by fixed factors after accepted/rejected steps
            //      imethod = 11: step size adjusted as described in More'
            imethod = 2; // 2 >> 0; 10 : zero lamda, 1, 10, 11: maxlamda
            iaccel = 0; // 0 or 1
            ibold = 0; // can be 0, 1, 2, 3, 4 
            ibroyden = 0; // 0 or 1
            initialfactor = 0.001; // initial lamda
            factoraccept = 2.5; // 2.5 for standard LM
            factorreject = 10; // 10 for standard LM
            avmax = 0.75;

            GeoDesicLM.geodesiclm(norris, norrisPrime, null,
                ref x, ref fvec, ref fjac, n, m,
                callback,
                ref info,
                analytic_jac, analytic_Avv,
                center_diff, h1, h2, ref dtd, damp_mode,
                ref niters, ref nfev, ref njev, ref naev,
                ref maxiter, maxfev, maxjev, maxaev, maxlam, minlam,
                artol, Cgoal, gtol, xtol, xrtol, ftol, frtol,
                ref converged, print_level, print_unit,
                imethod, iaccel, ibold, ibroyden,
                initialfactor, factoraccept, factorreject, avmax);

            Console.WriteLine("To exit applicaiton, type anykey.");
            Console.ReadKey(true);
        }

        
    }
}
