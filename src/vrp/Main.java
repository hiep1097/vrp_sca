package vrp;

public class Main {

    public static void main(String[] args) {
        int maxiter = 50;
        int numOfAgents = 20;

        fVRP fVRP = new fVRP();

        SCA gwo = new SCA(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
        gwo.execute();

        for (int i=0; i<gwo.getArrayRandomResult().length; i++){
            System.out.println("Lan chay thu "+(i+1)+ ":");
            test vrp = new test();
            vrp.Execute(gwo.getArrayRandomResult()[i]);
            vrp.printRoute();
        }
    }


}
