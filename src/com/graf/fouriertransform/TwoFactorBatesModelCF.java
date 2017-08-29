package com.graf.fouriertransform;

import org.apache.commons.math3.complex.Complex;

import net.finmath.fouriermethod.CharacteristicFunctionInterface;
import net.finmath.fouriermethod.models.ProcessCharacteristicFunctionInterface;


public class TwoFactorBatesModelCF implements ProcessCharacteristicFunctionInterface {

	
	private final double[] alpha;
	private final double[] beta;
	private final double[] sigma;
	private final double[] rho;	
	
	private final double[] lambda; //3 constants
	private final double k;
	private final double delta;
	private final double[] volatility;
	private final double initialValue;
	private final double riskFreeRate;
	private final double discountFactor;
	
	private final int numberOfFactors;
	
	//Two Factor Bates Model
	public TwoFactorBatesModelCF(
			double[] alpha,
			double[] beta,
			double[] sigma,
			double[] rho,
			double[] lambda,
			double k,
			double delta,
			double[] volatility,
			double initialPrice,
			double riskFreeRate,
			double discountFactor) {
		super();
		this.alpha			= alpha;
		this.beta			= beta;
		this.sigma			= sigma;
		this.rho			= rho;
		this.lambda			= lambda;
		this.k				= k;
		this.delta			= delta;
		this.volatility		= volatility;
		this.initialValue	= initialPrice;
		this.riskFreeRate	= riskFreeRate;
		this.discountFactor	= discountFactor;
		
		this.numberOfFactors = 2;
	}
	
	//Two-factor Bates with non-vector input
	public TwoFactorBatesModelCF(
			double alphaOne,
			double alphaTwo,
			double betaOne,
			double betaTwo,
			double sigmaOne,
			double sigmaTwo,
			double rhoOne,
			double rhoTwo,
			double lambdaZero,
			double lambdaOne,
			double lambdaTwo,
			double k,
			double delta,
			double volatilityOne,
			double volatilityTwo,
			double initialPrice,
			double riskFreeRate,
			double discountFactor) {
		super();
		this.alpha			= new double[]{alphaOne, alphaTwo};
		this.beta			= new double[]{betaOne, betaTwo};
		this.sigma			= new double[]{sigmaOne, sigmaTwo};
		this.rho			= new double[]{rhoOne, rhoTwo};
		this.lambda			= new double[]{lambdaZero, lambdaOne, lambdaTwo};
		this.k				= k;
		this.delta			= delta;
		this.volatility		= new double[]{volatilityOne, volatilityTwo};
		this.initialValue	= initialPrice;
		this.riskFreeRate	= riskFreeRate;
		this.discountFactor	= discountFactor;
		
		this.numberOfFactors = 2;
	}
	
	//One Factor Bates Model
	public TwoFactorBatesModelCF(
			double alpha,
			double beta,
			double sigma,
			double rho,
			double lambda,
			double lambdaOne,
			double k,
			double delta,
			double volatility,
			double initialPrice,
			double riskFreeRate,
			double discountFactor) {
		super();
		this.alpha			= new double[]{alpha};
		this.beta			= new double[]{beta};
		this.sigma			= new double[]{sigma};
		this.rho			= new double[]{rho};
		this.lambda			= new double[]{lambda, lambdaOne};
		this.k				= k;
		this.delta			= delta;
		this.volatility		= new double[]{volatility};
		this.initialValue	= initialPrice;
		this.riskFreeRate	= riskFreeRate;
		this.discountFactor	= discountFactor;
		
		this.numberOfFactors = 1;
	}

	/* (non-Javadoc)
	 * @see net.finmath.fouriermethod.models.ProcessCharacteristicFunctionInterface#apply(double)
	 */
	@Override
	public CharacteristicFunctionInterface apply(final double time) {
        return new CharacteristicFunctionInterface() {
            @Override
            public Complex apply(Complex argument) {
            	
                Complex iargument = argument.multiply(Complex.I);
                
                Complex[] gamma = new Complex[numberOfFactors];
                
                Complex[] A = new Complex[numberOfFactors];
                
                Complex[] B = new Complex[numberOfFactors];
                
                Complex C = iargument
                				.multiply(iargument)
                				.add(iargument.multiply(-1))
                				.multiply(0.5*delta*delta)			
                				.exp()
                				.multiply(new Complex(1+k).pow(iargument))
                				.add(-1)							
                				.add(iargument.multiply(-k));
                
                
               
                for(int i = 0; i < numberOfFactors; i++){
                	
                	if(sigma[i] != 0){
                	
                		gamma[i] = iargument
                				.multiply(rho[i]*sigma[i])
                				.subtract(beta[i])
                				.pow(2)
                				.subtract(
                					iargument.multiply(iargument)
                					.add(iargument.multiply(-1))
                					.multiply(0.5)
                					.add(C.multiply(lambda[i+1]))
                					.multiply(2*sigma[i]*sigma[i])
                					)
                				.sqrt()
                				;
                	
                	
                		A[i] =
                			iargument
                			.multiply(rho[i] * sigma[i])
                			.subtract(beta[i])
                			.subtract(gamma[i])
                			.multiply((-alpha[i]*time)/(sigma[i]*sigma[i]))
                			.subtract(iargument
                					.multiply(rho[i]*sigma[i])
                					.subtract(beta[i])
                					.subtract(gamma[i])
                					.multiply(new Complex(1).divide(gamma[i].multiply(time).exp())
                							.subtract(1)
                							.divide(gamma[i])
                							)
                					.multiply(0.5)
                					.add(new Complex(1).divide(gamma[i].multiply(time).exp()))
                					.log()
                					.add(gamma[i].multiply(time))
                					.multiply((2*alpha[i])/(sigma[i]*sigma[i]))
                					)
                			;
                	
                		B[i] = iargument
                			.multiply(iargument)
                			.add(iargument.multiply(-1))
                			.multiply(0.5)
                			.add(C.multiply(lambda[i+1]))
                			.multiply(-2)
                			.divide(iargument
                					.multiply(rho[i] * sigma[i])
                					.subtract(beta[i])
                					.add(gamma[i]
                							.multiply(new Complex(1).divide(gamma[i].multiply(time).exp())
                									.add(1)
                									.divide(new Complex(1).divide(gamma[i].multiply(time).exp())
                											.subtract(1)
                											)
                									)
                							)
                					)
                			;
                	}else{
                		B[i] = iargument
                    		.multiply(iargument)
                    		.add(iargument.multiply(-1))
                    		.multiply(0.5)
                    		.add(C.multiply(lambda[i+1]))
                    		.divide(-beta[i])
                    		.multiply(Math.exp(-beta[i]*time)-1)
                    		;
                		
                		
                		A[i] = iargument
                        	.multiply(iargument)
                        	.add(iargument.multiply(-1))
                        	.multiply(0.5)
                        	.add(C.multiply(lambda[i+1]))
                			.multiply(alpha[i]/-beta[i])
                			.multiply((-1/beta[i]) * Math.exp(-beta[i]*time) - time + (1/beta[i]))
                			;
                		
                		
                	}
                		
                	
                }
                

//                if(numberOfFactors == 2){
//                	return	A[0]
//                				.add(A[1])
//                				.add(B[0].multiply(volatility[0]))
//                				.add(B[1].multiply(volatility[1]))
//                				.add(C.multiply(time*lambda[0]))
//                				.add(iargument.multiply(Math.log(initialPrice)+time*riskFreeRate))
//                				.exp()
//                				.multiply(discountFactor);
//                }
//                else{
//                    return	A[0]
//                    			.add(B[0].multiply(volatility[0]))
//                    			.add(C.multiply(time*lambda[0]))
//                    			.add(iargument.multiply(Math.log(initialPrice)+time*riskFreeRate))
//                    			.exp()
//                    			.multiply(discountFactor);
//                }
				Complex characteristicFunction =
						A[0]
								.add(B[0].multiply(volatility[0]))
								.add(C.multiply(time*lambda[0]))
								.add(iargument.multiply(Math.log(initialValue)+time*riskFreeRate))
								.add(-riskFreeRate*time);

				if(numberOfFactors == 2) {
					characteristicFunction = characteristicFunction
							.add(A[1])
							.add(B[1].multiply(volatility[1]));
				}

				characteristicFunction = characteristicFunction.exp().multiply(discountFactor);;

				return characteristicFunction;
				
            };
        };
	}
}
