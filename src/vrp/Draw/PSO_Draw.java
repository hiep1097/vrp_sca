package vrp.Draw;

import vrp.solver.algorithm.f_xj;
import vrp.solver.algorithm.pso.Particles;
import vrp.solver.algorithm.pso.Swarm;

import java.io.IOException;

public class PSO_Draw {
    double r1, r2, c1, c2, w, wMax, wMin;
    int N;
    int dim;
    int maxIter;
    double lb[];
    double ub[];
    double vMax[];
    double vMin[];
    f_xj ff;
    double X[][];
    double[][] Result;
    double[][] arrRandomBestVal;
    Swarm swarm;
    double [] worstArr;

    //for draw
    double X_1[];   //gia tri x1 cua search agent dau tien sau moi lan lap
    double X_2[];   //gia tri x2 cua search agent dau tien sau moi lan lap

    public PSO_Draw(f_xj iff, double[] Lower, double[] Upper, int imaxIter, int iN) {
        maxIter = imaxIter;
        ff = iff;
        lb = Lower;
        ub = Upper;
        vMax = new double[lb.length];
        vMin = new double[ub.length];
        N = iN;
        dim = ub.length;
        X = new double[N][dim];
        arrRandomBestVal = new double[maxIter][dim];
        worstArr = new double[dim];

        X_1 = new double[maxIter];
        X_2 = new double[maxIter];
    }

    void init() {
        wMax = 0.9;
        wMin = 0.2;
        c1 = 2;
        c2 = 2;
        for (int i=0; i<lb.length; i++){
            vMax[i] = (ub[i]-lb[i])*0.2;
            vMin[i] = -vMax[i];
        }

        swarm = new Swarm(N, dim);
        for(int i = 0; i < N; i++) {
            swarm.particles[i] = new Particles(dim);
            for(int j = 0; j < dim; j++) {
                swarm.particles[i].X[j] = lb[j] + (ub[j] - lb[j]) * nextRand();
            }
        }

    }

    double[][] solution() {
        init();
        int iter = 1;
        while(iter <= maxIter) {
            // Calcualte the objective value
            for(int i = 0; i < N; i++) {
                double [] currentX = new double[dim];
                for (int j=0; j<dim; j++){
                    currentX[j] = swarm.particles[i].X[j];
                }
                swarm.particles[i].O = ff.func(currentX);

                if (swarm.particles[i].O > ff.func(worstArr)){
                    for (int j = 0; j<dim; j++){
                        worstArr[j] = swarm.particles[i].X[j];
                    }
                }

                //Update the PBEST
                if (swarm.particles[i].O < swarm.particles[i].PBEST.O){
                    for (int j=0; j<dim; j++){
                        swarm.particles[i].PBEST.X[j] = swarm.particles[i].X[j];
                    }
                    swarm.particles[i].PBEST.O = swarm.particles[i].O;
                }

                //Update the GBEST
                if (swarm.particles[i].O < swarm.GBEST.O){
                    for (int j=0; j<dim; j++){
                        swarm.GBEST.X[j] = swarm.particles[i].X[j];
                    }
                    swarm.GBEST.O = swarm.particles[i].O;
                }
            }

            //Update the X and V vectors
            w = wMax - iter * ((wMax - wMin) / maxIter);

            for(int i = 0; i < N; i++) {
                for (int j=0; j<dim; j++){
                    r1 = nextRand();
                    r2 = nextRand();
                    swarm.particles[i].V[j] = w*swarm.particles[i].V[j] + c1*r1*(swarm.particles[i].PBEST.X[j] - swarm.particles[i].X[j])
                            + c2*r2*(swarm.GBEST.X[j] - swarm.particles[i].X[j]);

                    //Check velocities
                    if (swarm.particles[i].V[j] > vMax[j]){
                        swarm.particles[i].V[j] = vMax[j];
                    }
                    if (swarm.particles[i].V[j] < vMin[j]){
                        swarm.particles[i].V[j] = vMin[j];
                    }

                    swarm.particles[i].X[j] = swarm.particles[i].X[j] + swarm.particles[i].V[j];

                    //Check positions
                    if (swarm.particles[i].X[j] > ub[j]){
                        swarm.particles[i].X[j] = lb[j] + ((ub[j] - lb[j]) * Math.random());
                    }
                    if (swarm.particles[i].X[j] < lb[j]){
                        swarm.particles[i].X[j] = lb[j] + ((ub[j] - lb[j]) * Math.random());
                    }
                }
            }

            for (int i=0; i<dim; i++){
                arrRandomBestVal[iter-1][i] = swarm.GBEST.X[i];
            }
            System.out.println("Iteration: "+iter);
            System.out.println("Best score: "+swarm.GBEST.O);
            if (iter==1){
                X[0][0] = 0.8;
                X[0][1] = 0.8;
            }
            X_1[iter-1] = swarm.particles[0].X[0];
            X_2[iter-1] = swarm.particles[0].X[1];
            iter++;
        }

        double[][] out = new double[2][dim];

        for (int i = 0; i < dim; i++){
            out[1][i] = swarm.GBEST.X[i];
        }

        out[0][0] = swarm.GBEST.O;
        return out;
    }

    public void execute()  {
        Result = solution();
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

    public double[] getWorstArray() {
        return worstArr;
    }

    public double getRes() throws IOException {
        ExcelUtils.fillX1X2ToExcelForDraw(X_1, X_2, maxIter, 9);
        return swarm.GBEST.O;
    }

    double nextRand(){
//        return 0.7;
//        position++;
//        return randomm[position-1];
        return Math.random();
    }
}
