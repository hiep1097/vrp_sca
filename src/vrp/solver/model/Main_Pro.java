package vrp.solver.model;

import vrp.solver.algorithm.gwo.GWO_pro;

public class Main_Pro {

    public static void main(String[] args) {
        fVRP fVRP = new fVRP();

        int maxiter = 200;
        int numOfAgents = 60;
        double a_gwo = 2.0;
        GWO_pro gwo = new GWO_pro(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents, a_gwo, null);
        gwo.execute();

        int maxiter1 = 150;
        int numOfAgents1 = 60;
        double a_gwo1 = 1.5;
        double initArray1[] = gwo.getArrayRandomResult()[199];
        GWO_pro gwo1 = new GWO_pro(fVRP, fVRP.Lower, fVRP.Upper, maxiter1, numOfAgents1, a_gwo1, initArray1);
        gwo1.execute();

        int maxiter2 = 100;
        int numOfAgents2 = 60;
        double a_gwo2 = 1.0;
        double initArray2[] = gwo1.getArrayRandomResult()[149];
        GWO_pro gwo2 = new GWO_pro(fVRP, fVRP.Lower, fVRP.Upper, maxiter2, numOfAgents2, a_gwo2, initArray2);
        gwo2.execute();

        int maxiter3 = 50;
        int numOfAgents3 = 60;
        double a_gwo3 = 0.5;
        double initArray3[] = gwo2.getArrayRandomResult()[99];
        GWO_pro gwo3 = new GWO_pro(fVRP, fVRP.Lower, fVRP.Upper, maxiter3, numOfAgents3, a_gwo3, initArray3);
        gwo3.execute();
        for (int i=0; i<gwo.getArrayRandomResult().length; i++){
            System.out.println("Lan chay thu "+(i+1)+ ":");
            test vrp = new test();
            vrp.Execute(gwo.getArrayRandomResult()[i]);
            vrp.printRoute();
        }

        for (int i=0; i<gwo.getArrayRandomResult().length; i++){
            System.out.println("1 Lan chay thu "+(i+1)+ ":");
            test vrp = new test();
            vrp.Execute(gwo.getArrayRandomResult()[i]);
            vrp.printRoute();
        }

        for (int i=0; i<gwo1.getArrayRandomResult().length; i++){
            System.out.println("2 Lan chay thu "+(i+1)+ ":");
            test vrp = new test();
            vrp.Execute(gwo1.getArrayRandomResult()[i]);
            vrp.printRoute();
        }

        for (int i=0; i<gwo2.getArrayRandomResult().length; i++){
            System.out.println("3 Lan chay thu "+(i+1)+ ":");
            test vrp = new test();
            vrp.Execute(gwo2.getArrayRandomResult()[i]);
            vrp.printRoute();
        }

        for (int i=0; i<gwo3.getArrayRandomResult().length; i++){
            System.out.println("4 Lan chay thu "+(i+1)+ ":");
            test vrp = new test();
            vrp.Execute(gwo3.getArrayRandomResult()[i]);
            vrp.printRoute();
        }
    }
}
