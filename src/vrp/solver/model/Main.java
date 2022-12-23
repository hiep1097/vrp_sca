package vrp.solver.model;

import vrp.solver.algorithm.gwo.GWO;
import vrp.solver.algorithm.sca.HSCA;
import vrp.solver.algorithm.sca.HSCA_1;

public class Main {

    public static void main(String[] args) {
        int maxiter = 200;
        int numOfAgents = 50;

        fVRP fVRP = new fVRP();

        long startTime = System.currentTimeMillis();
        GWO gwo = new GWO(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
        gwo.execute();
        long endTime = System.currentTimeMillis();
        long totalTime = endTime - startTime;

        for (int i=0; i<gwo.getArrayRandomResult().length; i++){
            System.out.println("Lan chay thu "+(i+1)+ ":");
            test vrp = new test();
            vrp.Execute(gwo.getArrayRandomResult()[i]);
            vrp.printRoute();
        }

        System.out.println("Thoi gian chay: "+ (totalTime / 1000.0) + " sec");
    }
}
