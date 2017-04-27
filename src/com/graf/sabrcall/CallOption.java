package com.graf.sabrcall;

import net.finmath.functions.NormalDistribution;

public class CallOption {

	
	private final double alpha; //initial volatility
	private final double beta;
	private final double initialForward;
	private final double rho;	
	
	private final double v;
	private final double strike;
	private final double maturity;
	private final double riskFreeRate;
	

	//j is set to 2 to extract the characteristic function

	//Two Factor Bates Model
	public CallOption(
			double alpha,
			double beta,
			double initialForward,
			double rho,
			double v,
			double strike,
			double maturity,
			double riskFreeRate) {
		super();
		this.alpha			= alpha;
		this.beta			= beta;
		this.initialForward	= initialForward;
		this.rho			= rho;
		this.v				= v;
		this.strike			= strike;
		this.maturity		= maturity;
		this.riskFreeRate	= riskFreeRate;
	}
	
	public double getPrice(){
		
		double z = (v/alpha) * Math.pow(initialForward*strike, 0.5*(1-beta))*Math.log(initialForward/strike);
		
		double x = Math.log( ( Math.sqrt(1-2*rho*z + z*z) + z - rho)/(1-rho) );
		
		double help1 = (1-beta)*Math.log(initialForward/strike); //to shorten the formula
		
		double taylor = 1+help1*help1/24 + help1*help1*help1*help1/1920;
		
		double impliedVolatility = (alpha/(Math.pow(initialForward*strike, 0.5*(1-beta))*taylor)
							)*
							(z/x)*
							(1+ maturity*((1-beta)*(1-beta)*alpha*alpha/(24*Math.pow(initialForward*strike, 1-beta))
										  + rho*beta*v*alpha/(4*Math.pow(initialForward*strike, 0.5*(1-beta)))
										  + (2-3*rho*rho)*v*v/24
										 )
							)
							;
		
		double d1 = (Math.log(initialForward/strike) + 0.5*impliedVolatility*impliedVolatility*maturity)/(impliedVolatility*Math.sqrt(maturity));
		
		double d2 = (Math.log(initialForward/strike) - 0.5*impliedVolatility*impliedVolatility*maturity)/(impliedVolatility*Math.sqrt(maturity));
		
		double discountFactor = Math.exp(-riskFreeRate * maturity);
		
		return discountFactor*(initialForward*NormalDistribution.cumulativeDistribution(d1) - strike * NormalDistribution.cumulativeDistribution(d2));
	}
}
