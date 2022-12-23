package vrp.solver.algorithm.sca;

import vrp.solver.algorithm.f_xj;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

public class HSCA {
    double infinity = 10E+50;
    int N, N_GWO, N_SCA;
    int dim;
    int maxiter;
    double Lower[];
    double Upper[];
    f_xj ff;
    double X[][], X_GWO[][], X_SCA[][];
    double Destination_position[];
    double Destination_fitness;
    double a_sca, r1, r2, r3, r4;
    double Objective_values[];
    double[][] Result;
    double[][] arrRandomBestVal;
    double alfa[], beta[], delta[];
    double X1[][], X2[][], X3[][];
    double a[], A1[], C1[], A2[], C2[], A3[], C3[];
    double percent;


    public HSCA(f_xj iff, double iLower[], double iUpper[], int imaxiter, int iN) {
        maxiter = imaxiter;
        ff = iff;
        Lower = iLower;
        Upper = iUpper;
        N = iN;
        dim = Upper.length;
        X = new double[N][dim];
        Destination_position = new double[dim];
        Destination_fitness = infinity;
        Objective_values = new double[N];
        arrRandomBestVal = new double[maxiter][dim];
        percent = 1;
        N_GWO = (int) (percent * N);
        N_SCA = N - N_GWO;
        System.out.println("N_GWO = "+N_GWO);
        System.out.println("N_SCA = "+N_SCA);
        a = new double[dim];
        alfa = new double[dim];
        beta = new double[dim];
        delta = new double[dim];
        A1 = new double[dim];
        C1 = new double[dim];
        A2 = new double[dim];
        C2 = new double[dim];
        A3 = new double[dim];
        C3 = new double[dim];
        X1 = new double[N][dim];
        X2 = new double[N][dim];
        X3 = new double[N][dim];
    }

    void init() {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < dim; j++) {
                X[i][j] = Lower[j] + (Upper[j] - Lower[j]) * Math.random();
            }
        }

        X = sort_and_index(X);
        for(int i = 0; i < dim; i++) {
            alfa[i] = X[0][i];
        }

        for(int i = 0; i < dim; i++) {
            beta[i] = X[1][i];
        }

        for(int i = 0; i < dim; i++) {
            delta[i] = X[2][i];
        }
        Destination_fitness = ff.func(X[0]);
    }

    double[][] solution() {
        init();
        int iter = 1; //start from the second iteration since the first iteration was dedicated to calculating the fitness
        while (iter <= maxiter) {
            //SCA
//            a_sca = 2.0;
//            double b = 0.5;
//            //r1 = a - iter*(a/(double)maxiter); // r1 decreases linearly from a to 0
//            r1 = a_sca * Math.sin((1 - iter / (double) maxiter) * (Math.PI / 2)) + b;
//
//            // Update the position of solutions with respect to destination
//            for (int i = N_GWO; i < N; i++) {
//                for (int j = 0; j < dim; j++) {
//                    // Update r2, r3, and r4 for Eq. (3.3)
//                    r2 = (2.0 * Math.PI) * Math.random();
//                    r3 = 2.0 * Math.random();
//                    r4 = Math.random();
//                    // Eq. (3.3)
//                    if (r4 < 0.5) {
//                        // Eq. (3.1)
//                        X[i][j] = X[i][j] + (r1 * Math.sin(r2) * Math.abs(r3 * Destination_position[j] - X[i][j]));
//                    } else {
//                        // Eq. (3.2)
//                        X[i][j] = X[i][j] + (r1 * Math.cos(r2) * Math.abs(r3 * Destination_position[j] - X[i][j]));
//                    }
//                }
//            }
//
//            for (int i = N_GWO; i < N; i++) {
//                for (int j = 0; j < dim; j++) {
//                    if (X[i][j] > Upper[j]) {
//                        X[i][j] = Upper[j];
//                        //X[i][j] = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
//                    }
//                    if (X[i][j] < Lower[j]) {
//                        X[i][j] = Lower[j];
//                        //X[i][j] = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
//                    }
//                }
//
//                //OBL
////                double[] x_OBL = new double[dim];
////                double[] x_MCS = new double[dim];
////                double[] x_MCS_OBL = new double[dim];
////
////                for (int j = 0; j < dim; j++) {
////                    x_OBL[j] = (1.0 - (double) iter * (1.0 - 0.1) / (double) maxiter) * (Lower[j] + Upper[j] - X[i][j]);
////                }
////
////                //MCS of default
////                for (int j = 0; j < dim; j++) {
////                    x_MCS[j] = X[i][j];
////                    double random = Math.random();
////                    if (random < 0.15) {
////                        x_MCS[j] = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
////                    }
////                }
////
////                //MCS of OBL
////                for (int j = 0; j < dim; j++) {
////                    x_MCS_OBL[j] = x_OBL[j];
////                    double random = Math.random();
////                    if (random < 0.15) {
////                        x_MCS_OBL[j] = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
////                    }
////                }
//
//                double minValue = infinity;
//                double minPosititon[] = new double[dim];
//
//                double defaultFitness = ff.func(X[i]);
////                double oblFitness = ff.func(x_OBL);
////                double mcsDefaultFitness = ff.func(x_MCS);
////                double mcsOblFitness = ff.func(x_MCS_OBL);
//
//                if (minValue > defaultFitness) {
//                    minValue = defaultFitness;
//                    for (int j = 0; j < dim; j++) {
//                        minPosititon[j] = X[i][j];
//                    }
//                }
//
////                if (minValue > oblFitness) {
////                    minValue = oblFitness;
////                    for (int j = 0; j < dim; j++) {
////                        minPosititon[j] = x_OBL[j];
////                    }
////                }
////
////                if (minValue > mcsDefaultFitness) {
////                    minValue = mcsDefaultFitness;
////                    for (int j = 0; j < dim; j++) {
////                        minPosititon[j] = x_MCS[j];
////                    }
////                }
////
////                if (minValue > mcsOblFitness) {
////                    minValue = mcsOblFitness;
////                    for (int j = 0; j < dim; j++) {
////                        minPosititon[j] = x_MCS_OBL[j];
////                    }
////                }
//
//                Objective_values[i] = minValue;
//                // Update the destination if there is a better solution
//                if (Objective_values[i] < Destination_fitness) {
//                    for (int j = 0; j < dim; j++) {
//                        Destination_position[j] = minPosititon[j];
//                    }
//                    Destination_fitness = Objective_values[i];
//                }
//            }
//
//            //GWO
            for (int j = 0; j < dim; j++) {
                a[j] = 2.0 - ((double) iter * (2.0 / (double) maxiter));
            }

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < dim; j++) {
                    r1 = Math.random();
                    r2 = Math.random();
                    for (int ii = 0; ii < dim; ii++) {
                        A1[ii] = 2.0 * a[ii] * r1 - a[ii];
                    }

                    for (int ii = 0; ii < dim; ii++) {
                        C1[ii] = 2.0 * r2;
                    }

                    X1[i][j] = alfa[j] - A1[j] * (Math.abs(C1[j] * alfa[j] - X[i][j]));
                    X1 = simplebounds(X1);
                    r1 = Math.random();
                    r2 = Math.random();

                    for (int ii = 0; ii < dim; ii++) {
                        A2[ii] = 2.0 * a[ii] * r1 - a[ii];
                    }

                    for (int ii = 0; ii < dim; ii++) {
                        C2[ii] = 2.0 * r2;
                    }

                    X2[i][j] = beta[j] - A2[j] * (Math.abs(C2[j] * beta[j] - X[i][j]));
                    X2 = simplebounds(X2);
                    r1 = Math.random();
                    r2 = Math.random();

                    for (int ii = 0; ii < dim; ii++) {
                        A3[ii] = 2.0 * a[ii] * r1 - a[ii];
                    }

                    for (int ii = 0; ii < dim; ii++) {
                        C3[ii] = 2.0 * r2;
                    }

                    X3[i][j] = delta[j] - A3[j] * (Math.abs(C3[j] * delta[j] - X[i][j]));
                    X3 = simplebounds(X3);
                    X[i][j] = (X1[i][j] + X2[i][j] + X3[i][j]) / 3.0;
                }
            }

            X = simplebounds(X);
            X = sort_and_index(X);

            for (int i = 0; i < dim; i++) {
                X[N - 1][i] = X[0][i];
            }

            for(int i = 0; i < dim; i++) {
                alfa[i] = X[0][i];
            }

            for(int i = 0; i < dim; i++) {
                beta[i] = X[1][i];
            }

            for(int i = 0; i < dim; i++) {
                delta[i] = X[2][i];
            }

            // Update the destination if there is a better solution
            if (ff.func(alfa) < Destination_fitness) {
                for (int j = 0; j < dim; j++) {
                    Destination_position[j] = alfa[j];
                }
                Destination_fitness = ff.func(alfa);
            }

            for (int j = 0; j < dim; j++) {
                arrRandomBestVal[iter-1][j] = alfa[j];
            }
            System.out.println("Iteration: "+iter);
            System.out.println("Best score: "+ff.func(alfa));
            iter++;
        }

        double[][] out = new double[2][dim];

        for (int i = 0; i < dim; i++) {
            out[1][i] = alfa[i];
        }

        out[0][0] = ff.func(alfa);
        return out;
    }

    public void execute() {
        Result = solution();
    }

    void toStringnew() throws IOException {
        double[][] in = solution();
        System.out.println("Optimized value = " + in[0][0]);
        for (int i = 0; i < dim; i++) {
            System.out.print("x[" + i + "] = " + in[1][i] + "\t");
        }
        System.out.println();
    }

    public double getRes() throws IOException {
        double[][] in = solution();
        return in[0][0];
    }

    public void toStringNew(String sMessage) throws IOException {
        System.out.println(sMessage + Result[0][0]);

        for (int i = 0; i < dim; i++) {
            System.out.println("x[" + i + "] = " + Result[1][i]);
        }

        System.out.println("----------------------------------------");
    }

    public double[] getBestArray() {
        return Result[1];
    }

    public double[][] getArrayRandomResult() {
        return arrRandomBestVal;
    }

    double[][] simplebounds(double s[][]) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < dim; j++) {
                if (s[i][j] < Lower[j]) {
                    s[i][j] = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
                }

                if (s[i][j] > Upper[j]) {
                    s[i][j] = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
                }
            }
        }
        return s;
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

        double[][] B = new double[N][dim];

        for(int i = 0; i < N; i++) {
            for(int j = 0; j < dim; j++) {
                B[i][j] = XXX[indexes[i]][j];
            }
        }

        return B ;
    }
}
