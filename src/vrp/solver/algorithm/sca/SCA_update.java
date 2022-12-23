package vrp.solver.algorithm.sca;

import vrp.solver.algorithm.Algorithm;
import vrp.solver.algorithm.f_xj;

import java.io.IOException;

public class SCA_update implements Algorithm  {
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

    public SCA_update(f_xj iff, double iLower[], double iUpper[], int imaxiter, int iN) {
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

    void init() {
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < dim; j++) {
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
            double b = 0.5;
            //r1 = a - iter*(a/(double)maxiter); // r1 decreases linearly from a to 0
            r1 = a*Math.sin((1-iter/(double) maxiter)*(Math.PI/2))+b;

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

                //OBL
                double [] x_OBL = new double[dim];
                double [] x_MCS = new double[dim];
                double [] x_MCS_OBL = new double[dim];

                for (int j=0; j<dim; j++){
                    x_OBL[j] = (1.0 - (double) iter *(1.0 - 0.1)/ (double) maxiter) * (Lower[j]+ Upper[j] - X[i][j]);
                }

                //MCS of default
                for (int j=0; j<dim; j++){
                    x_MCS[j] = X[i][j];
                    double random = Math.random();
                    if (random < 0.15){
                        x_MCS[j] = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
                    }
                }

                //MCS of OBL
                for (int j=0; j<dim; j++){
                    x_MCS_OBL[j] = x_OBL[j];
                    double random = Math.random();
                    if (random < 0.15){
                        x_MCS_OBL[j] = Lower[j] + ((Upper[j] - Lower[j]) * Math.random());
                    }
                }

                double minValue = infinity;
                double minPosititon[] = new double[dim];

                double defaultFitness = ff.func(X[i]);
                double oblFitness = ff.func(x_OBL);
                double mcsDefaultFitness = ff.func(x_MCS);
                double mcsOblFitness = ff.func(x_MCS_OBL);

                if (minValue > defaultFitness){
                    minValue = defaultFitness;
                    for (int j=0; j<dim; j++){
                        minPosititon[j] = X[i][j];
                    }
                }

                if (minValue > oblFitness){
                    minValue = oblFitness;
                    for (int j=0; j<dim; j++){
                        minPosititon[j] = x_OBL[j];
                    }
                }

                if (minValue > mcsDefaultFitness){
                    minValue = mcsDefaultFitness;
                    for (int j=0; j<dim; j++){
                        minPosititon[j] = x_MCS[j];
                    }
                }

                if (minValue > mcsOblFitness){
                    minValue = mcsOblFitness;
                    for (int j=0; j<dim; j++){
                        minPosititon[j] = x_MCS_OBL[j];
                    }
                }

                Objective_values[i] = minValue;
                // Update the destination if there is a better solution
                if (Objective_values[i]<Destination_fitness){
                    for (int j=0; j<dim; j++){
                        Destination_position[j] = minPosititon[j];
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
