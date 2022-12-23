package vrp.solver.main;

public class Main {

    public static void main(String[] args) throws Exception {
        /*
        #################################################
            - Valid algorithmName:
            HSCA
            SCA
            PSO
            DA
            ALO
        #################################################
            - Valid customerFileName:
            5customer
            8customer
            25customer
            25customerM
            30customer
            30customerM
        #################################################
        */
        String algorithmName = "HSCA";
        String customerFileName = "25customer";
        String chiadeuxe = "NO";  // YES or NO
        int maxiter = 200;
        int numOfAgents = 50;

        AlgorithmExecutor.executeAlgorithm(algorithmName, customerFileName, chiadeuxe, maxiter, numOfAgents);
    }
}
