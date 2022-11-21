package vrp.solver.model;

import java.io.File;
import java.util.*;

public class VRP{
    public double infinity = 100000;
    public List<Customer> customers = new ArrayList<>();
    public List<Vehicle> vehicles = new ArrayList<>();
    List<CustomerSchedule> customerSchedules;
    public int N;  //number of customer
    public int V;  //number of vehicle
    public double c[][];
    public int x[][];

    public double RES = 0;
    public double[] x_rand;

    public VRP() {
        readData("Data/30customer.txt");
        x = new int[N+1][N+1];
        x_rand = new double[N];
        //print_matrix(c, N+1);
//        System.out.println("So khach hang: "+N);
//        for (Customer customer: customers){
//            System.out.println("Name: "+customer.getName() + " - Demand: "+customer.getDemand());
//        }
//        System.out.println();
//        System.out.println("So xe: "+V);
//        for (Vehicle vehicle: vehicles){
//            System.out.println("Name: "+vehicle.getName() + " - Capacity: "+vehicle.getCapacity());
//        }
    }

    public double Execute(double x[]) {
        for (int i=0; i<x.length; i++){
            x_rand[i] = x[i];
        }
        ExecuteAlgorithm();
        return RES;
    }

    public void printRoute(){
        System.out.println("Result: "+ RES);
        System.out.println("Rout is: ");
        for (int k=0; k<V; k++){
            System.out.print("0->");
            for (int i = 0; i < N; i++){
                if (customerSchedules.get(i).vehicleName.equals(vehicles.get(k).name)){
                    System.out.print(customerSchedules.get(i).customerName+"->");
                }
            }
            System.out.print("0");
            System.out.println();
        }
    }

    public void ExecuteAlgorithm(){
        //Thực hiện khởi tạo quần thể sắp lịch ban đầu từ kết quả của thuật toán GWO thông qua x_rand[]
        customerSchedules = new ArrayList<CustomerSchedule>();
        for(int i = 0; i < customers.size(); i++){
            CustomerSchedule customerSchedule = new CustomerSchedule();
            customerSchedule.customerName = customers.get(i).name;
            customerSchedule.customerValue = (int) (x_rand[i]*10000);
            customerSchedules.add(customerSchedule);
        }

        // Dựa và giá trị value của dãy Random/Thuật toán GWO, sắp sếp lại thứ tự đổ tại các công trường
        customerSchedules.sort(Comparator.comparingDouble(o -> o.customerValue));

        RES = infinity;
        double sum = 0;

        int k = 0;
        int i = 0;
        while (k < V){
            Vehicle vehicle = vehicles.get(k);
            boolean ok = false;
            int count = 0;
            while (i < N &&
                    vehicle.capacity-getCustomerByName(customerSchedules.get(i).customerName).getDemand() >= 0){
                ok = true;
                Customer customer = getCustomerByName(customerSchedules.get(i).customerName);
                vehicle.capacity = vehicle.capacity - customer.demand;
                customerSchedules.get(i).vehicleName = vehicle.name;
                if (count == 0){
                    sum = sum + c[0][customer.position];
                } else {
                    Customer customerBefore = getCustomerByName(customerSchedules.get(i-1).customerName);
                    sum = sum + c[customer.position][customerBefore.position];
                }
                //System.out.println("gia tri kkkkkkk: "+k);
                //System.out.println("gia tri iiiiiii: "+i);
                count++;
                i++;
            }
            if (ok) {
                Customer customer = getCustomerByName(customerSchedules.get(i-1).customerName);
                sum = sum + c[0][customer.position];
            } else {
                RES = infinity;
                //System.out.println("Gia tri 1 ressssssss: "+RES);
                return;
            }
            k++;
        }
        if (k == V && i == N){
            RES = sum;
        } else {
            RES = infinity;
        }
        //System.out.println("Gia tri ressssssss: "+RES);
    }

    public void readData(String fileName) {
        try {
            File f = new File(fileName);
            Scanner scan = new Scanner(f);
            N = scan.nextInt();
            c = new double[N+1][N+1];    //distance between customers

            for (int i=1; i<=N; i++){
                int demand = scan.nextInt();
                customers.add(new Customer("Customer "+i, demand, i));
            }

            for (int i=0; i<=N; i++){
                for (int j=0; j<=N; j++){
                    c[i][j] = scan.nextDouble();
                    if (c[i][j] == -1.0){
                        c[i][j] = infinity;
                    }
                }
            }

            V = scan.nextInt();
            for (int i=1; i<=V; i++){
                int capacity = scan.nextInt();
                vehicles.add(new Vehicle("Vehicle "+(char)+(97+i-1), capacity));
            }
        } catch (Exception e){
            System.out.println("Error when read data: "+e.getMessage());
        }


    }

    public void print_matrix(int a[][], int length){
        for (int i=0; i<length; i++){
            for (int j=0; j<length; j++) System.out.print(a[i][j] + "\t\t");
            System.out.println();
        }
    }

    public Customer getCustomerByName(String name){
        return customers.stream()
                .filter(customer -> customer.name.equals(name))
                .findFirst()
                .orElse(null);
    }
}

