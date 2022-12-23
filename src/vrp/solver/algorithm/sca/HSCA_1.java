package vrp.solver.algorithm.sca;

import vrp.solver.algorithm.f_xj;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

public class HSCA_1
{
    double r1_gwo;
    double r2_gwo;
    int N, N_GWO, N_SCA;
    int dim;
    int maxiter;
    double alfa[];
    double beta[];
    double delta[];
    double Lower[];
    double Upper[];
    f_xj ff;
    double X[][];
    double X_GWO[][];
    double X_SCA[][];
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

    double Destination_position[];
    double Destination_fitness;
    double a_sca, r1_sca, r2_sca, r3_sca, r4_sca;
    double Objective_values[];


    double[][] Result;
    double[][] arrRandomBestVal;
    double [] worstArr;
    double infinity = 10E+50;
    double percent;

    public HSCA_1(f_xj iff,double iLower[],double iUpper[],int imaxiter,int iN)
    {
        maxiter = imaxiter;
        ff = iff;
        Lower = iLower;
        Upper = iUpper;
        N = iN;
        dim = Upper.length;
        a = new double[dim];
        X = new double[N][dim];
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
        arrRandomBestVal = new double[maxiter][dim];
        Destination_position = new double[dim];
        Destination_fitness = infinity;
        Objective_values = new double[N];
        percent = 1;
        N_GWO = (int) (percent * N);
        N_SCA = N - N_GWO;
        X_GWO = new double[N_GWO][dim];
        X_SCA = new double[N_SCA][dim];
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

    void init() {
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < dim; j++) {
                X[i][j] = Lower[j] + (Upper[j] - Lower[j]) * Math.random();
            }
        }

        X=sort_and_index(X);

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

    double[][] simplebounds(double s[][]) {
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < dim; j++) {
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
        init();
        int iter = 1;
        while(iter <= maxiter) {
            //tach quan the thanh GWO va SCA
            for (int i=0; i<N_GWO; i++){
                for (int j=0; j<dim; j++){
                    X_GWO[i][j] = X[i][j];
                }
            }
            for (int i=N_GWO; i<N; i++){
                for (int j=0; j<dim; j++){
                    X_SCA[i-N_GWO][j] = X[i][j];
                }
            }

            //GWO
            for(int j = 0; j < dim; j++) {
                a[j] = 2.0 -((double)iter * (2.0 / (double)maxiter));
            }

            for(int i = 0; i < N_GWO; i++) {
                for(int j = 0; j < dim; j++) {
                    r1_gwo = Math.random();
                    r2_gwo = Math.random();
                    for(int ii = 0; ii < dim; ii++) {
                        A1[ii] = 2.0 * a[ii] * r1_gwo - a[ii];
                    }

                    for(int ii = 0; ii < dim; ii++) {
                        C1[ii] = 2.0 * r2_gwo;
                    }

                    X1[i][j] = alfa[j] - A1[j] * (Math.abs(C1[j] * alfa[j] - X[i][j]));
                    X1 = simplebounds(X1);
                    r1_gwo = Math.random();
                    r2_gwo = Math.random();

                    for(int ii = 0; ii < dim; ii++){
                        A2[ii] = 2.0 * a[ii] * r1_gwo - a[ii];
                    }

                    for(int ii = 0; ii < dim; ii++) {
                        C2[ii] = 2.0*r2_gwo;
                    }

                    X2[i][j] = beta[j] - A2[j] * (Math.abs(C2[j] * beta[j] - X[i][j]));
                    X2 = simplebounds(X2);
                    r1_gwo = Math.random();
                    r2_gwo = Math.random();

                    for(int ii = 0; ii < dim; ii++) {
                        A3[ii] = 2.0 * a[ii] * r1_gwo - a[ii];
                    }

                    for(int ii = 0; ii < dim; ii++) {
                        C3[ii] = 2.0 * r2_gwo;
                    }

                    X3[i][j] = delta[j] - A3[j] * (Math.abs(C3[j] * delta[j] - X[i][j]));
                    X3 = simplebounds(X3);
                    X_GWO[i][j] = (X1[i][j] + X2[i][j] + X3[i][j]) / 3.0;
                }
            }

            //SCA
            a_sca = 2.0;
            double b = 0.5;
            //r1 = a - iter*(a/(double)maxiter); // r1 decreases linearly from a to 0
            r1_sca = a_sca*Math.sin((1-iter/(double) maxiter)*(Math.PI/2))+b;

            // Update the position of solutions with respect to destination
            for(int i = 0; i < N_SCA; i++) {
                for(int j = 0; j < dim; j++) {
                    // Update r2, r3, and r4 for Eq. (3.3)
                    r2_sca=(2.0*Math.PI)*Math.random();
                    r3_sca=2.0*Math.random();
                    r4_sca=Math.random();
                    // Eq. (3.3)
                    if (r4_sca < 0.5){
                        // Eq. (3.1)
                        X_SCA[i][j] = X_SCA[i][j]+(r1_sca*Math.sin(r2_sca)*Math.abs(r3_sca*Destination_position[j]-X_SCA[i][j]));
                    } else {
                        // Eq. (3.2)
                        X_SCA[i][j] = X_SCA[i][j]+(r1_sca*Math.cos(r2_sca)*Math.abs(r3_sca*Destination_position[j]-X_SCA[i][j]));
                    }
                }
            }

            for (int i=0; i < N_SCA; i++) {
                for (int j=0; j<dim; j++){
                    if (X_SCA[i][j] > Upper[j]){
                        X_SCA[i][j] = Upper[j];
                        //X[i][j] = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
                    }
                    if (X_SCA[i][j] < Lower[j]){
                        X_SCA[i][j] = Lower[j];
                        //X[i][j] = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
                    }
                }

                //OBL
//                double [] x_OBL = new double[dim];
//                double [] x_MCS = new double[dim];
//                double [] x_MCS_OBL = new double[dim];
//
//                for (int j=0; j<dim; j++){
//                    x_OBL[j] = (1.0 - (double) iter *(1.0 - 0.1)/ (double) maxiter) * (Lower[j]+ Upper[j] - X[i][j]);
//                }
//
//                //MCS of default
//                for (int j=0; j<dim; j++){
//                    x_MCS[j] = X[i][j];
//                    double random = Math.random();
//                    if (random < 0.15){
//                        x_MCS[j] = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
//                    }
//                }
//
//                //MCS of OBL
//                for (int j=0; j<dim; j++){
//                    x_MCS_OBL[j] = x_OBL[j];
//                    double random = Math.random();
//                    if (random < 0.15){
//                        x_MCS_OBL[j] = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
//                    }
//                }
//
                double minValue = infinity;
                double minPosititon[] = new double[dim];
                double defaultFitness = ff.func(X_SCA[i]);
//                double oblFitness = ff.func(x_OBL);
//                double mcsDefaultFitness = ff.func(x_MCS);
//                double mcsOblFitness = ff.func(x_MCS_OBL);
//
                if (minValue > defaultFitness){
                    minValue = defaultFitness;
                    for (int j=0; j<dim; j++){
                        minPosititon[j] = X[i][j];
                    }
                }
//
//                if (minValue > oblFitness){
//                    minValue = oblFitness;
//                    for (int j=0; j<dim; j++){
//                        minPosititon[j] = x_OBL[j];
//                    }
//                }
//
//                if (minValue > mcsDefaultFitness){
//                    minValue = mcsDefaultFitness;
//                    for (int j=0; j<dim; j++){
//                        minPosititon[j] = x_MCS[j];
//                    }
//                }
//
//                if (minValue > mcsOblFitness){
//                    minValue = mcsOblFitness;
//                    for (int j=0; j<dim; j++){
//                        minPosititon[j] = x_MCS_OBL[j];
//                    }
//                }

                Objective_values[i] = minValue;
                // Update the destination if there is a better solution
                if (Objective_values[i]<Destination_fitness){
                    for (int j=0; j<dim; j++){
                        Destination_position[j] = minPosititon[j];
                    }
                    Destination_fitness=Objective_values[i];
                }
            }

            X_GWO = simplebounds(X_GWO);
            //Hop nhat GWO va SCA
            for (int i=0; i<N_GWO; i++){
                for (int j=0; j<dim; j++){
                    X[i][j]=X_GWO[i][j];
                }
            }
            for (int i=N_GWO; i<N; i++){
                for (int j=0; j<dim; j++){
                    X[i][j]=X_SCA[i-N_GWO][j];
                }
            }

            X = sort_and_index(X);

            for(int i = 0; i < dim; i++) {
                X[N-1][i] = X[0][i];
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

            for(int i = 0; i < dim; i++) {
                arrRandomBestVal[iter-1][i] = X[0][i];
            }
            System.out.println("Iteration: "+iter);
            if (ff.func(arrRandomBestVal[iter-1]) == 100000){
                if (iter == 1){
                    while (ff.func(arrRandomBestVal[iter-1]) == 100000){
                        for(int i = 0; i < dim; i++) {
                            arrRandomBestVal[iter-1][i] = Lower[i] + (Upper[i] - Lower[i]) * Math.random();;
                        }
                    }
                } else {
                    for(int i = 0; i < dim; i++) {
                        arrRandomBestVal[iter-1][i] = arrRandomBestVal[iter-2][i];
                    }
                }
            }
            System.out.println("Best score: "+ff.func(arrRandomBestVal[iter-1]));
            iter++;
        }

        double[][] out = new double[2][dim];

        for(int i = 0; i < dim; i++){
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
