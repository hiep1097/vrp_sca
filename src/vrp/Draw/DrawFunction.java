package vrp.Draw;
import vrp.solver.model.fVRP;

public class DrawFunction {
    public static void main(String[] args) throws Exception {
        int maxiter = 50;
        int numOfAgents = 20;
        fVRP fVRP = new fVRP();
//        SCA_Draw result_SCA = new SCA_Draw(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
//        result_SCA.getRes();

        SCA_Update_Draw result_SCA_Update = new SCA_Update_Draw(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
        result_SCA_Update.getRes();
    }
}
