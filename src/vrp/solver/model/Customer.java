package vrp.solver.model;

public class Customer {
    public String name;
    public int demand;
    public int position;

    public Customer() {
    }

    public Customer(String name, int demand, int position) {
        this.name = name;
        this.demand = demand;
        this.position = position;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public int getDemand() {
        return demand;
    }

    public void setDemand(int demand) {
        this.demand = demand;
    }

    public int getPosition() {
        return position;
    }

    public void setPosition(int position) {
        this.position = position;
    }

}
