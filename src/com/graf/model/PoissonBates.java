package com.graf.model;



import net.finmath.functions.NormalDistribution;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.randomnumbers.MersenneTwister;
import net.finmath.stochastic.RandomVariableInterface;


public class PoissonBates {
	

	//Poisson RV like in exercise 10
	public static int generatePoisson(double lambda, MersenneTwister mersenne){
		int k=0;
		double p=1;
		double u;
		
		do{
			k++;
			u = mersenne.nextDouble();
			p *= u;
		}
		while(p > Math.exp(-lambda));
		
		return k-1;		
		
		
	}
	
	public static RandomVariableInterface getPoisson (double timeStep, int numberOfPaths, double k, double delta, RandomVariableInterface lambda, int seed){
		
		MersenneTwister mersenne = new MersenneTwister(seed);
		double[] jumpSize = new double[numberOfPaths];
	
		for(int path = 0; path < numberOfPaths; path++ ){
			for(int numberOfJumps = 0; numberOfJumps < generatePoisson(lambda.get(path) * timeStep, mersenne); numberOfJumps++){
//				jumpSize[path] += Math.log(Math.exp(Math.log(1-k) - 0.5*delta*delta + delta*NormalDistribution.inverseCumulativeDistribution(mersenne.nextDouble())));
				jumpSize[path] += Math.log(1-k) - 0.5*delta*delta + delta*NormalDistribution.inverseCumulativeDistribution(mersenne.nextDouble());
			}
		}
			
		return new RandomVariable(0, jumpSize);
			
	}
		

}
