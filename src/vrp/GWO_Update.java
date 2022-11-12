package vrp;

import java.io.IOException;

public class GWO_Update {
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
    double X1;
    double X2;
    double X3;
    double BESTVAL[];
    double iterdep[];
    double a;
    double A1;
    double C1;
    double A2;
    double C2;
    double A3;
    double C3;
    double Alpha_score, Beta_score, Delta_score;
    double inf = 10E+50;

    double[][] Result;
    double[][] arrRandomBestVal;

    public GWO_Update(f_xj iff, double iLower[], double iUpper[], int imaxiter, int iN) {
        maxiter = imaxiter;
        ff = iff;
        Lower = iLower;
        Upper = iUpper;
        N = iN;
        D = Upper.length;
        XX = new double[N][D];
        alfa = new double[D];
        beta = new double[D];
        delta = new double[D];
        Alpha_score = Beta_score = Delta_score = inf;
        BESTVAL = new double[maxiter];
        iterdep = new double[maxiter];

        arrRandomBestVal = new double[maxiter][D];
    }

    void init() {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < D; j++) {
                XX[i][j] = Lower[j] + (Upper[j] - Lower[j]) * Math.random();
            }
        }
    }

    double[][] solution() {
        init();
        int iter = 1;
        while (iter < maxiter) {
            for (int i = 0; i < N; i++){
                double fitnessValue = ff.func(XX[i]);
                if (fitnessValue < Alpha_score){
                    Alpha_score = fitnessValue;
                    for (int j = 0; j < D; j++) {
                        alfa[j] = XX[i][j];
                    }
                }
                if (fitnessValue > Alpha_score && fitnessValue < Beta_score){
                    Beta_score = fitnessValue;
                    for (int j = 0; j < D; j++) {
                        beta[j] = XX[i][j];
                    }
                }
                if (fitnessValue > Alpha_score && fitnessValue > Beta_score && fitnessValue < Delta_score){
                    Delta_score = fitnessValue;
                    for (int j = 0; j < D; j++) {
                        delta[j] = XX[i][j];
                    }
                }
            }

            a = 2.0 - ((double) iter * (2.0 / (double) maxiter));

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < D; j++) {
                    r1 = Math.random();
                    r2 = Math.random();
                    A1 = 2.0 * a * r1 - a;
                    C1 = 2.0 * r2;
                    X1 = alfa[j] - A1 * (Math.abs(C1 * alfa[j] - XX[i][j]));
                    if (X1 < Lower[j] || X1 > Upper[j]) X1 = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());

                    r1 = Math.random();
                    r2 = Math.random();
                    A2 = 2.0 * a * r1 - a;
                    C2 = 2.0 * r2;
                    X2 = beta[j] - A2 * (Math.abs(C2 * beta[j] - XX[i][j]));
                    if (X2 < Lower[j] || X2 > Upper[j]) X2 = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());

                    r1 = Math.random();
                    r2 = Math.random();
                    A3 = 2.0 * a * r1 - a;
                    C3 = 2.0 * r2;
                    X3 = delta[j] - A3 * (Math.abs(C3 * delta[j] - XX[i][j]));
                    if (X3 < Lower[j] || X3 > Upper[j]) X3 = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
                    XX[i][j] = (X1 + X2 + X3) / 3.0;
                }
            }
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < D; j++) {
                    if (XX[i][j] < Lower[j]) {
                        XX[i][j] = Lower[j];
                    }
                    if (XX[i][j] > Upper[j]) {
                        XX[i][j] = Upper[j];
                    }
                }
            }

            BESTVAL[iter] = Alpha_score;
            arrRandomBestVal[iter] = alfa;
            System.out.println("Iter: "+iter);
            System.out.println(BESTVAL[iter]);
            iter++;
        }

        double[][] out = new double[2][D];

        for (int i = 0; i < D; i++) {
            out[1][i] = alfa[i];
        }

        out[0][0] = ff.func(alfa);
        return out;
    }

    public void execute() {
        Result = solution();
    }

    public double[][] getArrayRandomResult() {
        return arrRandomBestVal;
    }
}
