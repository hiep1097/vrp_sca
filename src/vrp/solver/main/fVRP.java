package vrp.solver.main;

import vrp.solver.algorithm.f_xj;

import java.io.FileNotFoundException;

public class fVRP extends f_xj {
    public double[] Lower;
    public double[] Upper;
    String customerFileName;
    String chiadeuxe = "";

    public fVRP(String customerFileName) {
        this.customerFileName = customerFileName;
        VRP ff = new VRP(customerFileName, chiadeuxe);

        Lower = new double[ff.N];
        Upper = new double[ff.N];

        for(int i = 0; i < ff.N; i++) {
            Lower[i] = 0;
            Upper[i] = 1;
        }

        ff = null;
        System.gc();
    }

    public fVRP(String customerFileName, String chiadeuxe) {
        this.customerFileName = customerFileName;
        this.chiadeuxe = chiadeuxe;
        VRP ff = new VRP(customerFileName, chiadeuxe);

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
        VRP ff = new VRP(customerFileName, chiadeuxe);

        double res = ff.Execute(x);

        ff = null;
        System.gc();

        return Math.abs(res);
    }
}
