package vrp.solver.model;

import vrp.solver.algorithm.gwo.GWO;
import vrp.solver.algorithm.gwo.GWO_Update;

public class Main {

    public static void main(String[] args) {
        int maxiter = 100;
        int numOfAgents = 50;

        fVRP fVRP = new fVRP();

        GWO_Update gwo = new GWO_Update(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
        gwo.execute();

        for (int i=0; i<gwo.getArrayRandomResult().length; i++){
            System.out.println("Lan chay thu "+(i+1)+ ":");
            test vrp = new test();
            vrp.Execute(gwo.getArrayRandomResult()[i]);
            vrp.printRoute();
        }
    }
}
