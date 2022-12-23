package vrp.solver.algorithm.sca;

import vrp.solver.algorithm.f_xj;

import java.io.IOException;

public class SCA {
    double infinity = 10E+50;
    int N;
    int dim;
    int maxiter;
    double Lower[];
    double Upper[];
    f_xj ff;
    double X[][];
    double Destination_position[];
    double Destination_fitness;
    double a, r1, r2, r3, r4;
    double Objective_values[];
    double[][] Result;
    double[][] arrRandomBestVal;

    public SCA(f_xj iff, double iLower[], double iUpper[], int imaxiter, int iN) {
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
    }

    //double abc[] = new double[]{0.04,0.26,0.02,0.19,0.28,0.13,0.09,0.21,0.16,0.07,0.11,0.3,0.1,0.01,0.27,0.14,0.03,0.06,0.22,0.15,0.05,0.08,0.2,0.24,0.25,0.17,0.12,0.29,0.18,0.23};
    double abc[] = new double[]{0.05,0.2,0.08,0.29,0.19,0.17,0.13,0.24,0.01,0.1,0.15,0.22,0.12,0.07,0.21,0.11,0.04,0.09,0.03,0.02,0.06,0.14,0.25,0.27,0.28,0.3,0.16,0.18,0.23,0.26};

    void init() {
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < dim; j++) {
//                X[i][j]=abc[j];
//                if (i != 0){
//                    X[i][j] = Lower[j] + (Upper[j] - Lower[j]) * Math.random();
//                } else{
//
//                }
                X[i][j] = Lower[j] + (Upper[j] - Lower[j]) * Math.random();
            }
        }

        //Calculate the fitness of the first set and find the best one
        for (int i=0; i<N; i++){
            Objective_values[i] = ff.func(X[i]);
            if (i == 0){
                for (int j=0; j<dim; j++){
                    Destination_position[j] = X[i][j];
                }
                Destination_fitness=Objective_values[i];

            } else if (Objective_values[i] < Destination_fitness) {
                for (int j=0; j<dim; j++){
                    Destination_position[j] = X[i][j];
                }
                Destination_fitness=Objective_values[i];
            }
        }
        for (int j=0; j<dim; j++){
            arrRandomBestVal[0][j] = Destination_position[j];
        }
        System.out.println("Iter: 1");
        System.out.println(Destination_fitness);
    }

    double[][] simplebounds(double s[][]) {
        for(int i = 0; i < N; i++) {
        for(int j = 0; j < dim; j++) {
            if(s[i][j] < Lower[j]) {
                s[i][j] = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
            }

            if(s[i][j] > Upper[j]) {
                s[i][j] = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
            }
        }
    }
        return s;
}

    double[][] solution()  {
        init();
        int iter = 2; //start from the second iteration since the first iteration was dedicated to calculating the fitness
        while(iter <= maxiter) {
            a = 2.0;
            r1 = a - iter*(a/(double)maxiter); // r1 decreases linearly from a to 0

            // Update the position of solutions with respect to destination
            for(int i = 0; i < N; i++) {
                for(int j = 0; j < dim; j++) {
                    // Update r2, r3, and r4 for Eq. (3.3)
                    r2=(2.0*Math.PI)*Math.random();
                    r3=2.0*Math.random();
                    r4=Math.random();
                    // Eq. (3.3)
                    if (r4 < 0.5){
                        // Eq. (3.1)
                        X[i][j] = X[i][j]+(r1*Math.sin(r2)*Math.abs(r3*Destination_position[j]-X[i][j]));
                    } else {
                        // Eq. (3.2)
                        X[i][j] = X[i][j]+(r1*Math.cos(r2)*Math.abs(r3*Destination_position[j]-X[i][j]));
                    }
                }
            }

            for (int i=0; i < N; i++) {
                for (int j=0; j<dim; j++){
                    if (X[i][j] > Upper[j]){
                        X[i][j] = Upper[j];
                        //X[i][j] = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
                    }
                    if (X[i][j] < Lower[j]){
                        X[i][j] = Lower[j];
                        //X[i][j] = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
                    }
                }
                Objective_values[i] = ff.func(X[i]);
                // Update the destination if there is a better solution
                if (Objective_values[i]<Destination_fitness){
                    for (int j=0; j<dim; j++){
                        Destination_position[j] = X[i][j];
                    }
                    Destination_fitness=Objective_values[i];
                }
            }
            for (int j=0; j<dim; j++){
                arrRandomBestVal[iter-1][j] = Destination_position[j];
            }
            System.out.println("Iter: "+iter);
            System.out.println(Destination_fitness);
            iter++;
        }

        double[][] out = new double[2][dim];

        for(int i = 0; i < dim; i++){
            out[1][i] = Destination_position[i];
        }

        out[0][0] = Destination_fitness;
        return out;
    }

    public void execute() {
        Result = solution();
    }

    void toStringnew() throws IOException {
        double[][] in=solution();
        System.out.println("Optimized value = "+in[0][0]);
        for(int i=0;i<dim;i++) {
            System.out.print("x["+i+"] = "+in[1][i]+"\t");
        }
        System.out.println();
    }

    public double getRes() throws IOException {
        double[][] in=solution();
        return in[0][0];
    }

    public void toStringNew(String sMessage) throws IOException {
        System.out.println(sMessage + Result[0][0]);

        for(int i = 0; i < dim;i++) {
            System.out.println("x["+i+"] = "+ Result[1][i]);
        }

        System.out.println("----------------------------------------");
    }

    public double[] getBestArray()
    {
        return Result[1];
    }

    public double[][] getArrayRandomResult(){
        return arrRandomBestVal;
    }
}
