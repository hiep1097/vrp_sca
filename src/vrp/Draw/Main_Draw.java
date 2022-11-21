package vrp.Draw;

import vrp.solver.algorithm.alo.ALO;
import vrp.solver.algorithm.da.DA;
import vrp.solver.algorithm.gwo.GWO;
import vrp.solver.algorithm.pso.PSO;
import vrp.solver.algorithm.sca.SCA;
import vrp.solver.algorithm.sca.SCA_update;
import vrp.solver.model.fVRP;

public class Main_Draw {
    static double GWO_avg[] = new double[20];
    static double DA_avg[] = new double[20];
    static double PSO_avg[] = new double[20];
    static double ALO_avg[] = new double[20];
    static double SCA_avg[] = new double[20];
    static double SCA_Update_avg[] = new double[20];
    static int times = 20; //so lan chay
    static int maxiter = 50;
    static int numOfAgents = 20;

    public static void main(String[] args) throws Exception {
        SCA(times);
//        SCA_update(times);
//        DA(times);
//        GWO(times);
//        PSO(times);
//        ALO(times);
    }

    public static void SCA(int times) throws Exception {
        double avg[] = new double[maxiter];
        double sum[] = new double[maxiter];
        for (int i = 0; i < times; i++) {
            System.out.println("Times: "+i);
            fVRP fVRP = new fVRP();
            SCA result = new SCA(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
            result.execute();
            double position[][] = result.getArrayRandomResult();
            for (int j = 0; j < maxiter; j++) {
                sum[j] = sum[j] + fVRP.func(position[j]);
            }
        }

        for (int j = 0; j < maxiter; j++) {
            avg[j] = sum[j] / (double) times;
            System.out.println("SCA Avg Fitness at iteration " + (j + 1) + ": " + avg[j]);
        }

        ExcelUtils.fillAvgToExcel(maxiter, avg, 4);
    }

    public static void SCA_update(int times) throws Exception {
        double avg[] = new double[maxiter];
        double sum[] = new double[maxiter];
        for (int i = 0; i < times; i++) {
            System.out.println("Times: "+i);
            fVRP fVRP = new fVRP();
            SCA_update result = new SCA_update(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
            result.execute();
            double position[][] = result.getArrayRandomResult();
            for (int j = 0; j < maxiter; j++) {
                sum[j] = sum[j] + fVRP.func(position[j]);
            }
        }

        for (int j = 0; j < maxiter; j++) {
            avg[j] = sum[j] / (double) times;
            System.out.println("SCA_update Avg Fitness at iteration " + (j + 1) + ": " + avg[j]);
        }

        ExcelUtils.fillAvgToExcel(maxiter, avg, 10);
    }

    public static void DA(int times) throws Exception {
        double avg[] = new double[maxiter];
        double sum[] = new double[maxiter];
        for (int i = 0; i < times; i++) {
            System.out.println("Times: "+i);
            fVRP fVRP = new fVRP();
            DA result = new DA(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
            result.execute();
            double position[][] = result.getArrayRandomResult();
            for (int j = 0; j < maxiter; j++) {
                sum[j] = sum[j] + fVRP.func(position[j]);
            }
        }

        for (int j = 0; j < maxiter; j++) {
            avg[j] = sum[j] / (double) times;
            System.out.println("DA Avg Fitness at iteration " + (j + 1) + ": " + avg[j]);
        }

        ExcelUtils.fillAvgToExcel(maxiter, avg, 16);
    }

    public static void GWO(int times) throws Exception {
        double avg[] = new double[maxiter];
        double sum[] = new double[maxiter];
        for (int i = 0; i < times; i++) {
            System.out.println("Times: "+i);
            fVRP fVRP = new fVRP();
            GWO result = new GWO(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
            result.execute();
            double position[][] = result.getArrayRandomResult();
            for (int j = 0; j < maxiter; j++) {
                sum[j] = sum[j] + fVRP.func(position[j]);
            }
        }

        for (int j = 0; j < maxiter; j++) {
            avg[j] = sum[j] / (double) times;
            System.out.println("GWO Avg Fitness at iteration " + (j + 1) + ": " + avg[j]);
        }

        ExcelUtils.fillAvgToExcel(maxiter, avg, 22);
    }

    public static void PSO(int times) throws Exception {
        double avg[] = new double[maxiter];
        double sum[] = new double[maxiter];
        for (int i = 0; i < times; i++) {
            System.out.println("Times: "+i);
            fVRP fVRP = new fVRP();
            PSO result = new PSO(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
            result.execute();
            double position[][] = result.getArrayRandomResult();
            for (int j = 0; j < maxiter; j++) {
                sum[j] = sum[j] + fVRP.func(position[j]);
            }
        }

        for (int j = 0; j < maxiter; j++) {
            avg[j] = sum[j] / (double) times;
            System.out.println("PSO Avg Fitness at iteration " + (j + 1) + ": " + avg[j]);
        }

        ExcelUtils.fillAvgToExcel(maxiter, avg, 28);
    }

    public static void ALO(int times) throws Exception {
        double avg[] = new double[maxiter];
        double sum[] = new double[maxiter];
        for (int i = 0; i < times; i++) {
            System.out.println("Times: "+i);
            fVRP fVRP = new fVRP();
            ALO result = new ALO(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
            result.execute();
            double position[][] = result.getArrayRandomResult();
            for (int j = 0; j < maxiter; j++) {
                sum[j] = sum[j] + fVRP.func(position[j]);
            }
        }

        for (int j = 0; j < maxiter; j++) {
            avg[j] = sum[j] / (double) times;
            System.out.println("ALO Avg Fitness at iteration " + (j + 1) + ": " + avg[j]);
        }

        ExcelUtils.fillAvgToExcel(maxiter, avg, 34);
    }
}
