package vrp;

import java.io.File;
import java.util.*;

public class testtt {
    public double infinity = 10E+50;
    public List<Customer> customers = new ArrayList<>();
    public List<Vehicle> vehicles = new ArrayList<>();
    List<CustomerSchedule> customerSchedules;
    public int N;  //number of customer
    public int V;  //number of vehicle
    public double c[][];
    public int x[][];

    public double RES = 0;
    public double[] x_rand;

    public testtt() {
        readData("Data/30customer.txt");
        x = new int[N + 1][N + 1];
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
//        System.out.println("x_randd:");
        for (int i = 0; i < x.length; i++) {
            x_rand[i] = x[i];
//            System.out.print(x_rand[i] + " ");
        }
//        System.out.println();
        ExecuteAlgorithm();
        return RES;
    }

    public void printRoute() {
        System.out.println("Result: " + RES);
        if (RES != infinity) {
            System.out.println("Rout is: ");
            for (int k = 0; k < V; k++) {
                System.out.print("0->");
                for (int i = 0; i < N; i++) {
                    if (customerSchedules.get(i).vehicleName.equals(vehicles.get(k).name)) {
                        System.out.print(customerSchedules.get(i).customerName + "->");
                    }
                }
                System.out.print("0");
                System.out.println();
            }
        } else {
            System.out.println("RES is INFINITY!");
        }

    }

    public void ExecuteAlgorithm() {
        customerSchedules = new ArrayList<CustomerSchedule>();
        for (int i = 0; i < customers.size(); i++) {
            CustomerSchedule customerSchedule = new CustomerSchedule();
            customerSchedule.customerName = customers.get(i).name;
            customerSchedule.customerValue = (x_rand[i]);
            customerSchedules.add(customerSchedule);
        }
        customerSchedules.sort(Comparator.comparingDouble(o -> o.customerValue));

        RES = infinity;
        double sum = 0;

        for (int i = 0; i < V; i++) {
            Vehicle vehicle = vehicles.get(i);
            int count = 0;
            int positionBefore = 0;
            for (int j = 0; j < customerSchedules.size(); j++) {
                if ((i == 0 && customerSchedules.get(j).customerValue == 0.0) ||
                        (customerSchedules.get(j).customerValue > (double) i && customerSchedules.get(j).customerValue <= (double) (i + 1))) {
                    if (vehicle.capacity - getCustomerByName(customerSchedules.get(j).customerName).getDemand() >= 0) {
                        Customer customer = getCustomerByName(customerSchedules.get(j).customerName);
                        vehicle.capacity = vehicle.capacity - customer.getDemand();
                        customerSchedules.get(j).vehicleName = vehicle.name;
                        if (count == 0) {
                            sum = sum + c[0][customer.position];
                        } else {
                            Customer customerBefore = getCustomerByName(customerSchedules.get(positionBefore).customerName);
                            sum = sum + c[customer.position][customerBefore.position];
                        }
                        positionBefore = j;
                        count++;
//                        System.out.println("Summmmmm: "+sum);
//                        System.out.println("iiiiiiiiii: "+i +" Countttt: "+count);
                    } else {
//                        System.out.println("aasfafafaf");
                        return;
                    }
                }
            }
            Customer customer = getCustomerByName(customerSchedules.get(positionBefore).customerName);
            sum = sum + c[0][customer.position];
        }
        System.out.println("Doneeeee");
        RES = sum;
    }

    public void readData(String fileName) {
        try {
            File f = new File(fileName);
            Scanner scan = new Scanner(f);
            N = scan.nextInt();
            c = new double[N + 1][N + 1];    //distance between customers

            for (int i = 1; i <= N; i++) {
                int demand = scan.nextInt();
                customers.add(new Customer("Customer " + i, demand, i));
            }

            for (int i = 0; i <= N; i++) {
                for (int j = 0; j <= N; j++) {
                    c[i][j] = scan.nextDouble();
                    if (c[i][j] == -1.0) {
                        c[i][j] = infinity;
                    }
                }
            }

            V = scan.nextInt();
            for (int i = 1; i <= V; i++) {
                int capacity = scan.nextInt();
                vehicles.add(new Vehicle("Vehicle " + (char) +(97 + i - 1), capacity));
            }
        } catch (Exception e) {
            System.out.println("Error when read data: " + e.getMessage());
        }


    }

    public void print_matrix(int a[][], int length) {
        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) System.out.print(a[i][j] + "\t\t");
            System.out.println();
        }
    }

    public Customer getCustomerByName(String name) {
        return customers.stream()
                .filter(customer -> customer.name.equals(name))
                .findFirst()
                .orElse(null);
    }
}

