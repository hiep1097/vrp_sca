package vrp.solver.algorithm.gwo;

import vrp.solver.algorithm.f_xj;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

public class GWO_pro
{
    double r1;
    double r2;
    int N;
    int D;
    int maxiter;
    double alfa[];
    double beta[];
    double delta[];
    double Lower[];
    double Upper[];
    f_xj ff;
    double XX[][];
    double X1[][];
    double X2[][];
    double X3[][];
    double a[];
    double A1[];
    double C1[];
    double A2[];
    double C2[];
    double A3[];
    double C3[];

    double[][] Result;
    double[][] arrRandomBestVal;
    double a_gwo;
    double initArray[];

    public GWO_pro(f_xj iff, double iLower[], double iUpper[], int imaxiter, int iN, double a_gwo, double initArray[])
    {
        maxiter = imaxiter;
        ff = iff;
        Lower = iLower;
        Upper = iUpper;
        N = iN;
        D = Upper.length;
        a = new double[D];
        XX = new double[N][D];
        alfa = new double[D];
        beta = new double[D];
        delta = new double[D];
        A1 = new double[D];
        C1 = new double[D];
        A2 = new double[D];
        C2 = new double[D];
        A3 = new double[D];
        C3 = new double[D];
        X1 = new double[N][D];
        X2 = new double[N][D];
        X3 = new double[N][D];
        this.a_gwo = a_gwo;
        this.initArray = initArray;
        arrRandomBestVal = new double[maxiter][D];
    }

    double[][] sort_and_index(double[][] XXX) {
        double[] yval = new double[N];

        for(int i = 0; i < N; i++) {
            yval[i] = ff.func(XXX[i]);
        }

        ArrayList<Double> nfit = new ArrayList<Double>();

        for(int i = 0; i < N; i++) {
            nfit.add(yval[i]);
        }

        ArrayList<Double> nstore = new ArrayList<Double>(nfit);
        Collections.sort(nfit);

        double[] ret = new double[nfit.size()];
        Iterator<Double> iterator = nfit.iterator();

        int ii = 0;

        while(iterator.hasNext()) {
            ret[ii] = iterator.next().doubleValue();
            ii++;
        }

        int[] indexes = new int[nfit.size()];

        for(int n = 0; n < nfit.size(); n++) {
            indexes[n] = nstore.indexOf(nfit.get(n));
        }

        double[][] B = new double[N][D];

        for(int i = 0; i < N; i++) {
            for(int j = 0; j < D; j++) {
                B[i][j] = XXX[indexes[i]][j];
            }
        }

        return B ;
    }

    void init(double arr[]) {
        if (arr == null){
            for(int i = 0; i < N; i++) {
                for(int j = 0; j < D; j++) {
                    XX[i][j] = Lower[j] + (Upper[j] - Lower[j]) * Math.random();
                }
            }
        } else {
            for(int i = 0; i < N; i++) {
                for(int j = 0; j < D; j++) {
                    XX[i][j] = arr[j];
                }
            }
        }

        XX=sort_and_index(XX);

        for(int i = 0; i < D; i++) {
            alfa[i] = XX[0][i];
        }

        for(int i = 0; i < D; i++) {
            beta[i] = XX[1][i];
        }

        for(int i = 0; i < D; i++) {
            delta[i] = XX[2][i];
        }
    }

    double[][] simplebounds(double s[][])
    {
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < D; j++) {
                if(s[i][j] < Lower[j]) {
                    s[i][j] = Lower[j] * ((Upper[j] - Lower[j]) * Math.random());
                }

                if(s[i][j] > Upper[j]) {
                    s[i][j] = Lower[j] * ((Upper[j] - Lower[j]) * Math.random());
                }
            }
        }
        return s;
    }

    double[][] solution(){
        init(initArray);
        int iter = 1;
        while(iter <= maxiter)
        {
            for(int j = 0; j < D; j++) {
                a[j] = a_gwo -((double)iter * (a_gwo / (double)maxiter));
            }

            for(int i = 0; i < N; i++)
            {
                for(int j = 0; j < D; j++)
                {
                    r1 = Math.random();
                    r2 = Math.random();
                    for(int ii = 0; ii < D; ii++) {
                        A1[ii] = 2.0 * a[ii] * r1 - a[ii];
                    }

                    for(int ii = 0; ii < D; ii++) {
                        C1[ii] = 2.0 * r2;
                    }

                    X1[i][j] = alfa[j] - A1[j] * (Math.abs(C1[j] * alfa[j] - XX[i][j]));
                    X1 = simplebounds(X1);
                    r1 = Math.random();
                    r2 = Math.random();

                    for(int ii = 0; ii < D; ii++){
                        A2[ii] = 2.0 * a[ii] * r1 - a[ii];
                    }

                    for(int ii = 0; ii < D; ii++) {
                        C2[ii] = 2.0*r2;
                    }

                    X2[i][j] = beta[j] - A2[j] * (Math.abs(C2[j] * beta[j] - XX[i][j]));
                    X2 = simplebounds(X2);
                    r1 = Math.random();
                    r2 = Math.random();

                    for(int ii = 0; ii < D; ii++) {
                        A3[ii] = 2.0 * a[ii] * r1 - a[ii];
                    }

                    for(int ii = 0; ii < D; ii++) {
                        C3[ii] = 2.0 * r2;
                    }

                    X3[i][j] = delta[j] - A3[j] * (Math.abs(C3[j] * delta[j] - XX[i][j]));
                    X3 = simplebounds(X3);
                    XX[i][j] = (X1[i][j] + X2[i][j] + X3[i][j]) / 3.0;
                }
            }
            XX = simplebounds(XX);
            XX = sort_and_index(XX);

            for(int i = 0; i < D; i++) {
                XX[N-1][i] = XX[0][i];
            }

            for(int i = 0; i < D; i++) {
                alfa[i] = XX[0][i];
            }

            for(int i = 0; i < D; i++) {
                beta[i] = XX[1][i];
            }

            for(int i = 0; i < D; i++) {
                delta[i] = XX[2][i];
            }

            for(int i = 0; i < D; i++) {
                arrRandomBestVal[iter-1][i] = XX[0][i];
            }
            System.out.println("Iteration: "+iter);
            if (ff.func(arrRandomBestVal[iter-1]) == 100000){
                if (iter == 1){
                    while (ff.func(arrRandomBestVal[iter-1]) == 100000){
                        for(int i = 0; i < D; i++) {
                            arrRandomBestVal[iter-1][i] = Lower[i] + (Upper[i] - Lower[i]) * Math.random();;
                        }
                    }
                } else {
                    for(int i = 0; i < D; i++) {
                        arrRandomBestVal[iter-1][i] = arrRandomBestVal[iter-2][i];
                    }
                }
            }
            System.out.println("Best score: "+ff.func(arrRandomBestVal[iter-1]));
            iter++;
        }

        double[][] out = new double[2][D];

        for(int i = 0; i < D; i++){
            out[1][i] = alfa[i];
        }

        out[0][0] = ff.func(alfa);
        return out;
    }

    public void execute(){
        Result = solution();
    }

    public double[][] getArrayRandomResult(){
        return arrRandomBestVal;
    }
}
