package vrp.solver.algorithm.pso;

public class Swarm {
    public Particles[] particles;
    public BestResult GBEST;
    public Swarm(int numOfParticles, int dim){
        particles = new Particles[numOfParticles];
        GBEST = new BestResult(dim);
    }
}
