package vrp.solver.model;

import vrp.solver.algorithm.f_xj;

public class fVRP extends f_xj {
    public double[] Lower;
    public double[] Upper;

    public fVRP(){
        test ff = new test();

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
        test ff = new test();

        double res = ff.Execute(x);

        ff = null;
        System.gc();

        return Math.abs(res);
    }
}
