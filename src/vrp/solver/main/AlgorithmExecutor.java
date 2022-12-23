package vrp.solver.main;

import vrp.solver.algorithm.Algorithm;
import vrp.solver.algorithm.alo.ALO;
import vrp.solver.algorithm.da.DA;
import vrp.solver.algorithm.gwo.GWO;
import vrp.solver.algorithm.pso.PSO;
import vrp.solver.algorithm.sca.HSCA;
import vrp.solver.algorithm.sca.SCA;
import vrp.solver.algorithm.sca.SCA_update;

import java.io.File;
import java.util.Scanner;

public class AlgorithmExecutor {
    public static final String ALO = "ALO";
    public static final String PSO = "PSO";
    public static final String DA = "DA";
    public static final String SCA = "SCA";
    public static final String HSCA = "HSCA";


    public static void executeAlgorithm(String algorithmName, String customerFileName, String chiadeuxe, int maxiter, int numOfAgents) throws Exception {

        File f = new File("Data/"+customerFileName+".txt");
        Scanner scan = new Scanner(f);

        fVRP fVRP = new fVRP(customerFileName, chiadeuxe);

        Algorithm algorithm = null;

        switch (algorithmName){
            case ALO:
                algorithm = new ALO(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
                break;
            case PSO:
                algorithm = new PSO(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
                break;
            case DA:
                algorithm = new DA(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
                break;
            case SCA:
                if (numOfAgents < 20){
                    algorithm = new SCA(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
                } else {
                    algorithm = new PSO(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
                }
                break;
            case HSCA:
                if (numOfAgents < 20){
                    algorithm = new SCA_update(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
                } else {
                    algorithm = new HSCA(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
                }
                break;
            default:
                throw new Exception("Algorithm name wrong!");
        }

        if (algorithmName.equals(ALO))
        algorithm = new ALO(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);

        long startTime = System.currentTimeMillis();

        algorithm.execute();
        long endTime = System.currentTimeMillis();
        long totalTime = endTime - startTime;

        for (int i=0; i<algorithm.getArrayRandomResult().length; i++){
            System.out.println("Lan chay thu "+(i+1)+ ":");
            VRP vrp = new VRP(customerFileName, chiadeuxe);
            vrp.Execute(algorithm.getArrayRandomResult()[i]);
            vrp.printRoute();
        }

        System.out.println("Thoi gian chay: "+ (totalTime / 1000.0) + " sec");



    }
}
