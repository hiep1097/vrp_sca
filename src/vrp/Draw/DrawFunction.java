package vrp.Draw;
import vrp.solver.main.fVRP;

public class DrawFunction {
    static String customerFileName = "30customerM";
    public static void main(String[] args) throws Exception {
        int maxiter = 200;
        int numOfAgents = 50;
        fVRP fVRP = new fVRP(customerFileName);
//        PSO_Draw pso_draw = new PSO_Draw(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
//        pso_draw.solution();
//        pso_draw.getRes();

        DA_Draw da_draw = new DA_Draw(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
        da_draw.solution();
        da_draw.getRes();
    }
}
