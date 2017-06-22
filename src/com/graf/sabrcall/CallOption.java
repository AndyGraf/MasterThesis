package com.graf.sabrcall;

import net.finmath.functions.NormalDistribution;

public class CallOption {

	
	private final double alpha; 
	private final double beta;
	private final double initialForward;
	private final double rho;	
	
	private final double v;
	private final double strike;
	private final double maturity;
	private final double riskFreeRate;
	
	private final double volatility;
	
	//for the model
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
		this.volatility		= -1;
	}
	
	//for given implied volatility
	public CallOption(
			double initialForward,
			double volatility,
			double maturity,
			double strike,
			double discountCurve) {
		super();
		this.alpha			= 0;
		this.beta			= 0;
		this.initialForward	= initialForward;
		this.rho			= 0;
		this.v				= 0;
		this.strike			= strike;
		this.maturity		= maturity;
		this.riskFreeRate	= discountCurve;
		this.volatility		= volatility;
	}
	
	
	public double getPrice(){
		
		double impliedVolatility;
		
		if(volatility < 0){
		
			double z = (v/alpha) * Math.pow(initialForward*strike, 0.5*(1-beta))*Math.log(initialForward/strike);
		
			double x = Math.log( ( Math.sqrt(1-2*rho*z + z*z) + z - rho)/(1-rho) );
		
			double help1 = (1-beta)*Math.log(initialForward/strike); //to shorten the formula
		
			double taylor = 1+help1*help1/24 + help1*help1*help1*help1/1920;
		
			impliedVolatility = (alpha/(Math.pow(initialForward*strike, 0.5*(1-beta))*taylor)
							)*
							(z/x)*
							(1+ maturity*((1-beta)*(1-beta)*alpha*alpha/(24*Math.pow(initialForward*strike, 1-beta))
										  + rho*beta*v*alpha/(4*Math.pow(initialForward*strike, 0.5*(1-beta)))
										  + (2-3*rho*rho)*v*v/24
										 )
							)
							;
		}
		else{impliedVolatility = volatility;}
		
		double d1 = (Math.log(initialForward/strike) + 0.5*impliedVolatility*impliedVolatility*maturity)/(impliedVolatility*Math.sqrt(maturity));
		
		double d2 = (Math.log(initialForward/strike) - 0.5*impliedVolatility*impliedVolatility*maturity)/(impliedVolatility*Math.sqrt(maturity));
		
		double discountFactor;
		
		if(volatility < 0){
			discountFactor = Math.exp(-riskFreeRate * maturity);
		}else{discountFactor = riskFreeRate;}
		
		return discountFactor*(initialForward*NormalDistribution.cumulativeDistribution(d1) - strike * NormalDistribution.cumulativeDistribution(d2));
	}
}
