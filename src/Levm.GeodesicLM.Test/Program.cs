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
            initialfactor = 0.001; // initial lamda
            factoraccept = 5; // 10 for standard LM
            factorreject = 2; // 10 for standard LM
            avmax = 0.75;

            GeoDesicLM.geodesiclm(Rosenbrock, RosenbrockJacobian, RosenbrockAvv,
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
        }

        
    }
}
