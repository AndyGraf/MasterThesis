package com.graf.fouriertransform;

import net.finmath.functions.NormalDistribution;

public class SABRVolatility {

	

	
	static public double getVolatility(
			double alpha,	//initial volatility
			double beta,
			double rho,
			double v,	//vol of vol
			double initialForward,
			double strike,
			double maturity) {
		
		double z = (v/alpha) * Math.pow(initialForward*strike, 0.5*(1-beta))*Math.log(initialForward/strike);
		
		double x = Math.log( ( Math.sqrt(1-2*rho*z + z*z) + z - rho)/(1-rho) );
		
		double help1 = (1-beta)*Math.log(initialForward/strike); //to shorten the formula of volatility
		
		double volatility = (alpha/(Math.pow(initialForward*strike, 0.5*(1-beta))*
									(1+help1*help1/24 + help1*help1*help1*help1/1920)
									)
							)*(z/x)*
							(1+ maturity*((1-beta)*(1-beta)*alpha*alpha/(24*Math.pow(initialForward*strike, 1-beta))
										  + rho*beta*v*alpha/(4*Math.pow(initialForward*strike, 0.5*(1-beta)))
										  + (2-3*rho*rho)*v*v/24
										 )
							)
							;
		
		
		return volatility;
	}
}
