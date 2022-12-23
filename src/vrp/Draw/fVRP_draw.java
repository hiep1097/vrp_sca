package vrp.Draw;

import vrp.solver.algorithm.f_xj;
import vrp.solver.model.test;

public class fVRP_draw extends f_xj {
    public double[] Lower;
    public double[] Upper;

    public fVRP_draw(){
        test_draw ff = new test_draw();

        Lower = new double[ff.N];
        Upper = new double[ff.N];

        for(int i = 0; i < ff.N; i++) {
            Lower[i] = 0;
            Upper[i] = 1;
        }

        ff = null;
        System.gc();
    }

    public double func(double x[]){
        test_draw ff = new test_draw();

        double res = ff.Execute(x);

        ff = null;
        System.gc();

        return Math.abs(res);
    }
}
