package vrp.Draw;

import vrp.solver.algorithm.alo.ALO;
import vrp.solver.algorithm.da.DA;
import vrp.solver.algorithm.gwo.GWO;
import vrp.solver.algorithm.pso.PSO;
import vrp.solver.algorithm.sca.HSCA;
import vrp.solver.algorithm.sca.SCA;
import vrp.solver.algorithm.sca.SCA_update;
import vrp.solver.model.fVRP;
import vrp.solver.model.test;

public class Main_Draw {
    static int times = 20; //so lan chay
    static int maxiter = 1000;
    static int numOfAgents = 60;

    public static void main(String[] args) throws Exception {
        GWO(times);
//        PSO(times);
//        DA(times);

//        ALO(times);
//        SCA(times);
//        SCA_update(times);
    }

    public static void SCA(int times) throws Exception {
        double avg[] = new double[maxiter];
        double sum[] = new double[maxiter];
        double fitness_optimize[] = new double[times];
        for (int i = 0; i < times; i++) {
            System.out.println("Times: "+i);
            fVRP fVRP = new fVRP();
            PSO result = new PSO(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
            result.execute();
            double position[][] = result.getArrayRandomResult();
            for (int j = 0; j < maxiter; j++) {
                sum[j] = sum[j] + fVRP.func(position[j]);
                System.out.println("Lan chay thu "+(j+1)+ ":");
                test vrp = new test();
                vrp.Execute(position[j]);
                vrp.printRoute();
            }
            fitness_optimize[i] = fVRP.func(result.getArrayRandomResult()[maxiter-1]);
        }

        for (int j = 0; j < maxiter; j++) {
            avg[j] = sum[j] / (double) times;
            System.out.println("SCA Avg Fitness at iteration " + (j + 1) + ": " + avg[j]);
        }

        ExcelUtils.fillAvgToExcel(maxiter, avg, 4);
        System.out.println("Optimize fitness in each time:");
        for (int i=0; i<times; i++){
            System.out.println(fitness_optimize[i]);
        }
    }

    public static void SCA_update(int times) throws Exception {
        double avg[] = new double[maxiter];
        double sum[] = new double[maxiter];
        double fitness_optimize[] = new double[times];
        for (int i = 0; i < times; i++) {
            System.out.println("Times: "+i);
            fVRP fVRP = new fVRP();
            SCA_update result = new SCA_update(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
            result.execute();
            double position[][] = result.getArrayRandomResult();
            for (int j = 0; j < maxiter; j++) {
                sum[j] = sum[j] + fVRP.func(position[j]);
                System.out.println("Lan chay thu "+(j+1)+ ":");
                test vrp = new test();
                vrp.Execute(position[j]);
                vrp.printRoute();
            }
            fitness_optimize[i] = fVRP.func(result.getArrayRandomResult()[maxiter-1]);
        }

        for (int j = 0; j < maxiter; j++) {
            avg[j] = sum[j] / (double) times;
            System.out.println("SCA_update Avg Fitness at iteration " + (j + 1) + ": " + avg[j]);
        }

        ExcelUtils.fillAvgToExcel(maxiter, avg, 10);
        System.out.println("Optimize fitness in each time:");
        for (int i=0; i<times; i++){
            System.out.println(fitness_optimize[i]);
        }
    }

    public static void DA(int times) throws Exception {
        double avg[] = new double[maxiter];
        double sum[] = new double[maxiter];
        double fitness_optimize[] = new double[times];
        for (int i = 0; i < times; i++) {
            System.out.println("Times: "+i);
            fVRP fVRP = new fVRP();
            DA result = new DA(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
            result.execute();
            double position[][] = result.getArrayRandomResult();
            for (int j = 0; j < maxiter; j++) {
                sum[j] = sum[j] + fVRP.func(position[j]);
                System.out.println("Lan chay thu "+(j+1)+ ":");
                test vrp = new test();
                vrp.Execute(position[j]);
                vrp.printRoute();
            }
            fitness_optimize[i] = fVRP.func(result.getArrayRandomResult()[maxiter-1]);
        }

        for (int j = 0; j < maxiter; j++) {
            avg[j] = sum[j] / (double) times;
            System.out.println("DA Avg Fitness at iteration " + (j + 1) + ": " + avg[j]);
        }

        ExcelUtils.fillAvgToExcel(maxiter, avg, 16);
        System.out.println("Optimize fitness in each time:");
        for (int i=0; i<times; i++){
            System.out.println(fitness_optimize[i]);
        }
    }

    public static void GWO(int times) throws Exception {
        double avg[] = new double[maxiter];
        double sum[] = new double[maxiter];
        double fitness_optimize[] = new double[times];
        for (int i = 0; i < times; i++) {
            System.out.println("Times: "+i);
            fVRP fVRP = new fVRP();
            GWO result = new GWO(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
            result.execute();
            double position[][] = result.getArrayRandomResult();
            for (int j = 0; j < maxiter; j++) {
                sum[j] = sum[j] + fVRP.func(position[j]);
            }
            fitness_optimize[i] = fVRP.func(result.getArrayRandomResult()[maxiter-1]);
            test vrp = new test();
            vrp.Execute(position[maxiter-1]);
            vrp.printRoute();
        }

        for (int j = 0; j < maxiter; j++) {
            avg[j] = sum[j] / (double) times;
            System.out.println("GWO Avg Fitness at iteration " + (j + 1) + ": " + avg[j]);
        }

        ExcelUtils.fillAvgToExcel(maxiter, avg, 22);
        System.out.println("Optimize fitness in each time:");
        for (int i=0; i<times; i++){
            System.out.println(fitness_optimize[i]);
        }
    }

    public static void PSO(int times) throws Exception {
        double avg[] = new double[maxiter];
        double sum[] = new double[maxiter];
        double fitness_optimize[] = new double[times];
        for (int i = 0; i < times; i++) {
            System.out.println("Times: "+i);
            fVRP fVRP = new fVRP();
            PSO result = new PSO(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
            result.execute();
            double position[][] = result.getArrayRandomResult();
            for (int j = 0; j < maxiter; j++) {
                sum[j] = sum[j] + fVRP.func(position[j]);
                System.out.println("Lan chay thu "+(j+1)+ ":");
                test vrp = new test();
                vrp.Execute(position[j]);
                vrp.printRoute();
            }
            fitness_optimize[i] = fVRP.func(result.getArrayRandomResult()[maxiter-1]);
        }

        for (int j = 0; j < maxiter; j++) {
            avg[j] = sum[j] / (double) times;
            System.out.println("PSO Avg Fitness at iteration " + (j + 1) + ": " + avg[j]);
        }

        ExcelUtils.fillAvgToExcel(maxiter, avg, 28);
        System.out.println("Optimize fitness in each time:");
        for (int i=0; i<times; i++){
            System.out.println(fitness_optimize[i]);
        }
    }

    public static void ALO(int times) throws Exception {
        double avg[] = new double[maxiter];
        double sum[] = new double[maxiter];
        double fitness_optimize[] = new double[times];
        for (int i = 0; i < times; i++) {
            System.out.println("Times: "+i);
            fVRP fVRP = new fVRP();
            ALO result = new ALO(fVRP, fVRP.Lower, fVRP.Upper, maxiter, numOfAgents);
            result.execute();
            double position[][] = result.getArrayRandomResult();
            for (int j = 0; j < maxiter; j++) {
                sum[j] = sum[j] + fVRP.func(position[j]);
                System.out.println("Lan chay thu "+(j+1)+ ":");
                test vrp = new test();
                vrp.Execute(position[j]);
                vrp.printRoute();
            }
            fitness_optimize[i] = fVRP.func(result.getArrayRandomResult()[maxiter-1]);
        }

        for (int j = 0; j < maxiter; j++) {
            avg[j] = sum[j] / (double) times;
            System.out.println("ALO Avg Fitness at iteration " + (j + 1) + ": " + avg[j]);
        }

        ExcelUtils.fillAvgToExcel(maxiter, avg, 34);
        System.out.println("Optimize fitness in each time:");
        for (int i=0; i<times; i++){
            System.out.println(fitness_optimize[i]);
        }
    }
}
