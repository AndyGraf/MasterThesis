package com.graf.fouriertransform;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Scanner;

import net.finmath.exception.CalculationException;
import net.finmath.fouriermethod.products.AbstractProductFourierTransform;
import net.finmath.fouriermethod.products.EuropeanOption;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.optimizer.LevenbergMarquardt;
import net.finmath.optimizer.SolverException;

//  Calibration of the One- and Two-factor Bates and Heston models as well as the SABR model 
//  towards the swaption implied volatility surface.
//
//  The function net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility used in this class has been adjusted.
//  optionValue has been adjusted to 2*optionValue: 
//  double volatilityUpperBound = (2*optionValue + Math.abs(forward-optionStrike)) / Math.sqrt(optionMaturity) / payoffUnit;
//  found in line (845) of version 1.9
//
//  The Levenberg-Marquardt algorithm uses the damping of Levenberg
//  (JTJ + lambda Id) = JT error
//  alphaElement += 1 * lambda;
//  instead of the originally implemented version of Marquardt
//  (JTJ + lambda diag(JTJ)) = JT error
//	alphaElement *= 1 + lambda;
//  found in line (688) of version 1.6 net.finmath.optimizer.LevenbergMarquardt

public class LevenbergCalibrationBachelierVolatility {
	
	//maximum number of iterations for the LM-algorithm
	static int maxIteration = 150;
	
	static DecimalFormat df = new DecimalFormat("#.####");
	static DecimalFormat error = new DecimalFormat("#.#########");
    
	//defining the options
	static int[] maturities = {1, 2, 5, 10, 15, 20};
	static double[] strikes = {-1,-0.5, -0.25, 0, 0.25, 0.5, 1, 1.5, 2};
	static int tenor;
	static double shift;
	
	
	static double [] targetValues = new double[maturities.length*strikes.length];
	static double [] targetPrices = new double[maturities.length*strikes.length];
	static double [] weights = new double[maturities.length*strikes.length];
	static double [] values = new double[maturities.length*strikes.length];

	static AbstractProductFourierTransform[][] europeanMatrix = new AbstractProductFourierTransform[maturities.length][strikes.length];
	
    //jump parameters for the heston models
    static double lambdaZero	= 0;
    static double lambdaOne		= 0;
    static double lambdaTwo		= 0;
    static double k				= 0;
    static double delta			= 0;
    
    //we use a SABR model without displacement
    static double sabrDisplacement = 0;
	

	public static void main(String[] args) throws IOException {

		
		SABRdata sheetdata = new SABRdata();
		
		System.out.println( "Model to calibrate (input number): \n" + 
							"1) Two-factor Bates \n" + 
							"2) Two-factor Heston \n" + 
							"3) One-factor Bates \n" + 
							"4) One-factor Heston \n" + 
							"5) SABR \n" +
							"6) maturity-wise SABR"
							);
		Scanner in = new Scanner(System.in);
		int input = in.nextInt();
		
		//setting the market data
		tenor = 10;
		if(tenor == 2){shift = 2.65;}else if(tenor == 5){shift = 1.6;}else if(tenor == 10){shift = 1.5;}else{shift = 0;};
		double[][] volatilities = sheetdata.getSmileData(tenor);
		for(int w = 0; w < weights.length; w++){weights[w] = 1;};


		//calculate the target values
		System.out.println("Target swaption prices in % of the notional principal amount with a tenor of " + tenor + " years:\n");
		System.out.println("Maturty" + "\t\t\t\t\t\t\t\t Strikes\n" + "\t-1%      \t-0.5%    \t-0.25%    \t0%      \t0.25%      \t0.5%      \t1%      \t1.5%      \t2%\n");
		for(int i = 0; i < maturities.length; i++){
			System.out.print(maturities[i] + "  \t");
			for(int j = 0; j < strikes.length; j++){
				
				
				targetPrices[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.blackScholesGeneralizedOptionValue(sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,volatilities[i][j],maturities[i], strikes[j]+shift,sheetdata.getSwapAnnuity(maturities[i], tenor));
				
				//calculate target bachelier implied volatilities
				targetValues[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(
						sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
						maturities[i],
						strikes[j]+shift,
						sheetdata.getSwapAnnuity(maturities[i], tenor),
						targetPrices[i*strikes.length + j]);
				
				europeanMatrix[i][j] = new EuropeanOption(maturities[i], strikes[j]+shift);
				System.out.print(df.format(targetValues[i*strikes.length + j]) + "     \t");

			}
			System.out.println();
		}
	  	
		
		switch(input){
			case 1:
				double[] twoBatesParameters = twoFactorBates(sheetdata);
				printTwoFactorBates(twoBatesParameters, sheetdata);
				break;
			case 2:
				double[] twoHestonParameters = twoFactorHeston(sheetdata);
				printTwoFactorHeston(twoHestonParameters, sheetdata);
				break;
			case 3:
				double[] oneBatesParameters = oneFactorBates(sheetdata);
				printOneFactorBates(oneBatesParameters, sheetdata);
				break;
			case 4:
				double[] oneHestonParameters = oneFactorHeston(sheetdata);
				printOneFactorHeston(oneHestonParameters, sheetdata);
				break;
			case 5:
				double[] sabrParameters = sabr(sheetdata);
				System.out.println("Root mean squared error: " + error.format(Math.sqrt((printSABR(sabrParameters, sheetdata)/(maturities.length*strikes.length)))));
				break;
			case 6:
				double err = 0;
				System.out.println("\nSwaption prices calculated with a maturity-wise calibration and the error in brackets:\n");
				System.out.println("Maturty" + "\t\t\t\t\t\t\t\t Strikes\n" + "\t-1%      \t-0.5%    \t-0.25%    \t0%      \t0.25%      \t0.5%      \t1%      \t1.5%      \t2%\n");
				double[] mTargetValues = new double[strikes.length];
				for(int i = 0; i<maturities.length; i++){
					int[] maturity = {maturities[i]};
					for(int j = 0; j < strikes.length; j++){
						mTargetValues[j] = targetValues[i*strikes.length + j];
					}
					double[] mSabrParameters = sabr(sheetdata, maturity, mTargetValues);
					

//					 optional print for the maturity-wise parameters
//					System.out.println(
//							"\n\talpha = \t\t" + Math.exp(mSabrParameters[0]) + 
//							";\n\tbeta  = \t\t" + Math.exp(-mSabrParameters[1]*mSabrParameters[1]) + 
//							";\n\trho = \t\t\t" + Math.tanh(0.5*mSabrParameters[2]) + 
//							";\n\tv = \t\t\t" + Math.exp(mSabrParameters[3]) + ";"
//					);
					
					err += printSABR(mSabrParameters, sheetdata, maturity, mTargetValues);
				}
				System.out.println("Root mean squared error: " + error.format(Math.sqrt((err/(maturities.length*strikes.length)))));
				break;
		}

		

	  	

		

		

		

	}
	
	
	
	//calibration of the 5 models with Levenberg-Marquardt, returning the best fit parameters
	
	
	
	public static double[] twoFactorBates(SABRdata sheetdata){
		
		//initial parameters  
		
        double alphaOne     = 0.01;
        double alphaTwo		= 0.04;
        double betaOne      = 0.91;
        double betaTwo      = 1.76;
        double sigmaOne     = 0.582;
        double sigmaTwo     = 0.346;
        double tanhRhoOne   = -0.848;
        double tanhRhoTwo  	= -0.402;
        double rhoOne		= Math.log( ( 1.0 + tanhRhoOne) / ( 1 - tanhRhoOne));
        double rhoTwo		= Math.log( ( 1.0 + tanhRhoTwo) / ( 1 - tanhRhoTwo));

        
        double lambdaZero   = 0.1143;
        double lambdaOne    = 1.56;
        double lambdaTwo    = 0.28;
        double k       		= -0.057;
        double delta   		= 0.102;
        
        double volatilityOne   = 0.10963;
        double volatilityTwo   = 0.001352;   

        
        double[] initialParameters = new double[]{			
            	Math.log(alphaOne),
            	Math.log(alphaTwo),
    			Math.log(betaOne),
    			Math.log(betaTwo),
    			Math.log(sigmaOne),
    			Math.log(sigmaTwo),
    			rhoOne,
    			rhoTwo,
    			Math.log(lambdaZero),
    			Math.log(lambdaOne),
    			Math.log(lambdaTwo),
    			k,
    			delta,
    			Math.log(volatilityOne),
    			Math.log(volatilityTwo)
            };
        
		
        
	 	LevenbergMarquardt optimizer = new LevenbergMarquardt() {

	 		@Override
		 	public void setValues(double[] parameters, double[] values) {
		 		for(int i = 0; i < maturities.length; i++){
					for(int j = 0; j < strikes.length; j++){
						try {
							values[i*strikes.length + j] = 	net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(
									sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
									maturities[i],
									strikes[j]+shift,
									sheetdata.getSwapAnnuity(maturities[i], tenor),
									europeanMatrix[i][j].getValue(
																new TwoFactorBatesModelCF(
																Math.exp(parameters[0]),
																Math.exp(parameters[1]),
																Math.exp(parameters[2]),
																Math.exp(parameters[3]),																	
																Math.exp(parameters[4]),
																Math.exp(parameters[5]),
																Math.tanh(0.5*parameters[6]),
																Math.tanh(0.5*parameters[7]),
																Math.exp(parameters[8]),
																Math.exp(parameters[9]),
																Math.exp(parameters[10]),
																parameters[11],
																Math.abs(parameters[12]),
																Math.exp(parameters[13]),
																Math.exp(parameters[14]),
																sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																0,
																sheetdata.getSwapAnnuity(maturities[i], tenor))
															)
									)
							;
						} catch (IOException e) {
							e.printStackTrace();
						} catch (CalculationException e) {
							e.printStackTrace();
						}
					}
		 		}
			}
		};
	  	optimizer.setInitialParameters(initialParameters);
	  	optimizer.setWeights(weights);
	  	optimizer.setMaxIteration(maxIteration);
	  	
	  	optimizer.setTargetValues(targetValues);
	  
	  	try {
			optimizer.run();
		} catch (SolverException e) {
			e.printStackTrace();
		}
	  
	  	System.out.println("\n\n\tThe solver required " + optimizer.getIterations() + " iterations.");
		
		return optimizer.getBestFitParameters();
		
	}
	
	public static double[] twoFactorHeston(SABRdata sheetdata){
		
		//initial parameters 
		
        double alphaOne     = 0.028;
        double alphaTwo		= 0.130;
        double betaOne      = 0.0;
        double betaTwo      = 5.58;
        double sigmaOne     = 1.039;
        double sigmaTwo     = 0.667;
        double tanhRhoOne   = -0.775;
        double tanhRhoTwo  	= -0.382;
        double rhoOne		= Math.log( ( 1.0 + tanhRhoOne) / ( 1 - tanhRhoOne));
        double rhoTwo		= Math.log( ( 1.0 + tanhRhoTwo) / ( 1 - tanhRhoTwo));
        
        double volatilityOne   = 0.00963;
        double volatilityTwo   = 0.01352;      
        
        double[] initialParameters = new double[]{			
            	Math.log(alphaOne),
            	Math.log(alphaTwo),
    			Math.log(betaOne),
    			Math.log(betaTwo),
    			Math.log(sigmaOne),
    			Math.log(sigmaTwo),
    			rhoOne,
    			rhoTwo,
    			Math.log(volatilityOne),
    			Math.log(volatilityTwo)
            };
        
        LevenbergMarquardt optimizer = new LevenbergMarquardt() {

		 	public void setValues(double[] parameters, double[] values) {
		 		for(int i = 0; i < maturities.length; i++){
					for(int j = 0; j < strikes.length; j++){
						try {
							values[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(
									sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
									maturities[i],
									strikes[j]+shift,
									sheetdata.getSwapAnnuity(maturities[i], tenor),
									europeanMatrix[i][j].getValue(
																new TwoFactorBatesModelCF(
																Math.exp(parameters[0]),
																Math.exp(parameters[1]),
																Math.exp(parameters[2]),
																Math.exp(parameters[3]),																	
																Math.exp(parameters[4]),
																Math.exp(parameters[5]),
																Math.tanh(0.5*parameters[6]),
																Math.tanh(0.5*parameters[7]),
																lambdaZero,
																lambdaOne,
																lambdaTwo,
																k,
																delta,
																Math.exp(parameters[8]),
																Math.exp(parameters[9]),
																sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																0,
																sheetdata.getSwapAnnuity(maturities[i], tenor))
															)
									)
							;
						} catch (IOException e) {
							e.printStackTrace();
						} catch (CalculationException e) {
							e.printStackTrace();
						}
					}
		 		}
			}
		};
	  	optimizer.setInitialParameters(initialParameters);
	  	optimizer.setWeights(weights);
	  	optimizer.setMaxIteration(maxIteration);
	  	
	  	optimizer.setTargetValues(targetValues);
	  
	  	try {
			optimizer.run();
		} catch (SolverException e) {
			e.printStackTrace();
		}
	  	
	  	System.out.println("\n\n\tThe solver required " + optimizer.getIterations() + " iterations.");
	  
	  	return optimizer.getBestFitParameters();
        
	}
	
	public static double[] oneFactorBates(SABRdata sheetdata){
		
		//initial parameters  
        
        double alphaOne     = 0.049;
        double betaOne      = 2.45;
        double sigmaOne     = 0.378;
        double tanhRhoOne   = -0.545;
        double rhoOne		= Math.log( ( 1.0 + tanhRhoOne) / ( 1 - tanhRhoOne));
        
        double lambdaZero   = 0;
        double lambdaOne	= 27.19;
        double k       		= -0.095;
        double delta   		= 0.109;
        
        double volatilityOne   = 0.00963; 
        
        double[] initialParameters = new double[]{			
            	Math.log(alphaOne),
    			Math.log(betaOne),
    			Math.log(sigmaOne),
    			rhoOne,
    			Math.log(lambdaZero),
    			Math.log(lambdaOne),
    			k,
    			delta,
    			Math.log(volatilityOne)
            };
        
        LevenbergMarquardt optimizer = new LevenbergMarquardt() {

		 	public void setValues(double[] parameters, double[] values) {
		 		for(int i = 0; i < maturities.length; i++){
					for(int j = 0; j < strikes.length; j++){
						try {
							values[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(
									sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
									maturities[i],
									strikes[j]+shift,
									sheetdata.getSwapAnnuity(maturities[i], tenor),
									europeanMatrix[i][j].getValue(
																new TwoFactorBatesModelCF(
																Math.exp(parameters[0]),
																Math.exp(parameters[1]),
																Math.exp(parameters[2]),
																Math.tanh(0.5*parameters[3]),
																Math.exp(parameters[4]),
																Math.exp(parameters[5]),
																parameters[6],
																Math.abs(parameters[7]),
																Math.exp(parameters[8]),
																sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																0,
																sheetdata.getSwapAnnuity(maturities[i], tenor))
															)
									)
							;
						} catch (IOException e) {
							e.printStackTrace();
						} catch (CalculationException e) {
							e.printStackTrace();
						}
					}
		 		}
			}
		};
	  	optimizer.setInitialParameters(initialParameters);
	  	optimizer.setWeights(weights);
	  	optimizer.setMaxIteration(maxIteration);
	  	
	  	optimizer.setTargetValues(targetValues);
	  
	  	try {
			optimizer.run();
		} catch (SolverException e) {
			e.printStackTrace();
		}
	  
	  	System.out.println("\n\n\tThe solver required " + optimizer.getIterations() + " iterations.");
	  	
	  	return optimizer.getBestFitParameters();
        
	}
	
	
	public static double[] oneFactorHeston(SABRdata sheetdata){
		
		//initial parameters  
        
        double alpha		= 0.1;
        double beta			= 1.49;
        double sigma		= 0.742;
        double tanhRho		= -0.571;
        double rho			= Math.log( ( 1.0 + tanhRho) / ( 1 - tanhRho));

        double volatilityOne   = 0.00963;          
        
        
        double[] initialParameters = new double[]{			
            	Math.log(alpha),
    			Math.log(beta),
    			Math.log(sigma),
    			rho,
    			Math.log(volatilityOne)
            };
        
	 	LevenbergMarquardt optimizer = new LevenbergMarquardt() {

		 	public void setValues(double[] parameters, double[] values) {
		 		for(int i = 0; i < maturities.length; i++){
					for(int j = 0; j < strikes.length; j++){
						try {
							values[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(
									sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
									maturities[i],
									strikes[j]+shift,
									sheetdata.getSwapAnnuity(maturities[i], tenor),
									europeanMatrix[i][j].getValue(
																new TwoFactorBatesModelCF(
																Math.exp(parameters[0]),
																Math.exp(parameters[1]),
																Math.exp(parameters[2]),
																Math.tanh(0.5*parameters[3]),
																lambdaZero,
																lambdaOne,
																k,
																delta,
																Math.exp(parameters[4]),
																sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																0,
																sheetdata.getSwapAnnuity(maturities[i], tenor))
															)
									)
							;
						} catch (IOException e) {
							e.printStackTrace();
						} catch (CalculationException e) {
							e.printStackTrace();
						}
					}
		 		}
			}
		};
	  	optimizer.setInitialParameters(initialParameters);
	  	optimizer.setWeights(weights);
	  	optimizer.setMaxIteration(maxIteration);
	  	
	  	optimizer.setTargetValues(targetValues);
	  
	  	try {
			optimizer.run();
		} catch (SolverException e) {
			e.printStackTrace();
		}
	  	
	  	System.out.println("\n\n\tThe solver required " + optimizer.getIterations() + " iterations.");
	  
	  	return optimizer.getBestFitParameters();
        
	}
	
	public static double[] sabr(SABRdata sheetdata) throws IOException{
		
		return sabr(sheetdata, maturities, targetValues);
	}
	
	public static double[] sabr(SABRdata sheetdata, int [] maturities, double[] targetVolatilities) throws IOException{
		
		//initial parameters 
		
        double alpha	    = 0.2;
        double beta 		= 1;
        double rho      	= -0.3;
        double v   			= 0.25;

		

        
		double[] initialParameters = new double[]{			
        		Math.log(alpha),
        		beta,
        		rho,
        		Math.log(v),

            };
		
        
        LevenbergMarquardt optimizer = new LevenbergMarquardt() {

		 	public void setValues(double[] parameters, double[] values) {
		 		for(int i = 0; i < maturities.length; i++){
					for(int j = 0; j < strikes.length; j++){
							try {
								values[i*strikes.length + j] = AnalyticFormulas.sabrBerestyckiNormalVolatilityApproximation(
																	Math.exp(parameters[0]),
																	Math.exp(-parameters[1]*parameters[1]),
																	Math.tanh(0.5*parameters[2]),
																	Math.exp(parameters[3]),
																	sabrDisplacement,
																	sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																	strikes[j]+shift,
																	maturities[i]);
							} catch (IOException e) {
								e.printStackTrace();
							}
					}
		 		}
			}
		};
	  	optimizer.setInitialParameters(initialParameters);
	  	optimizer.setWeights(weights);
	  	optimizer.setMaxIteration(maxIteration);
	  	
	  	optimizer.setTargetValues(targetVolatilities);
	  
	  	try {
			optimizer.run();
		} catch (SolverException e) {
			e.printStackTrace();
		}
	  
	  	
	  	if(maturities.length >1)System.out.println("\n\n\tThe solver required " + optimizer.getIterations() + " iterations.");
	  	
	  	return optimizer.getBestFitParameters();
	}
	
	
	
	
	
	
	//print methods to show the volatility surface with the calibrated parameters
	
	
	public static void printTwoFactorBates(double[] bestParameters, SABRdata sheetdata) throws IOException{
		System.out.println( "\n\tThe best fit paramters for the Two-factor Bates model are: \n\n" +
				"\talphaOne = \t\t" + Math.exp(bestParameters[0]) + 
				"\n\talphaTwo = \t\t" + Math.exp(bestParameters[1]) + 
				"\n\tbetaOne = \t\t" + Math.exp(bestParameters[2]) + 
				"\n\tbetaTwo = \t\t" + Math.exp(bestParameters[3]) + 
				"\n\tsigmaOne = \t\t" + Math.exp(bestParameters[4]) + 
				"\n\tsigmaTwo = \t\t" + Math.exp(bestParameters[5]) + 
				"\n\trhoOne = \t\t" + Math.tanh(0.5*bestParameters[6]) + 
				"\n\trhoTwo = \t\t" + Math.tanh(0.5*bestParameters[7]) + 
				"\n\tlambdaZero = \t\t" + Math.exp(bestParameters[8]) + 
				"\n\tlambdaOne = \t\t" + Math.exp(bestParameters[9])  +
				"\n\tlambdaTwo = \t\t" + Math.exp(bestParameters[10]) + 
				"\n\tk = \t\t\t" + bestParameters[11] + 
				"\n\tdelta = \t\t" + Math.abs(bestParameters[12]) +			
				"\n\tvolatilityOne = \t" + Math.exp(bestParameters[13]) + 
				"\n\tvolatilityTwo = \t" + Math.exp(bestParameters[14])
						);


		System.out.println("\nSwaption prices calculated with the calibrated parameters and the error in brackets:\n");
		System.out.println("Maturty" + "\t\t\t\t\t\t\t\t Strikes\n" + "\t-1%      \t-0.5%    \t-0.25%    \t0%      \t0.25%      \t0.5%      \t1%      \t1.5%      \t2%\n");
		
		double squaredError = 0;
		
		for(int i = 0; i < maturities.length; i++){
			System.out.print(maturities[i] + "  \t");
			for(int j = 0; j < strikes.length; j++){
				try {
					values[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(
							sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
							maturities[i],
							strikes[j]+shift,
							sheetdata.getSwapAnnuity(maturities[i], tenor),
							europeanMatrix[i][j].getValue(
															new TwoFactorBatesModelCF(
																	Math.exp(bestParameters[0]),
																	Math.exp(bestParameters[1]),
																	Math.exp(bestParameters[2]),
																	Math.exp(bestParameters[3]),
																	Math.exp(bestParameters[4]),
																	Math.exp(bestParameters[5]),
																	Math.tanh(0.5*bestParameters[6]),
																	Math.tanh(0.5*bestParameters[7]),
																	Math.exp(bestParameters[8]),
																	Math.exp(bestParameters[9]),
																	Math.exp(bestParameters[10]),
																	bestParameters[11],
																	Math.abs(bestParameters[12]),
																	Math.exp(bestParameters[13]),
																	Math.exp(bestParameters[14]),
																	sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																	0,
																	sheetdata.getSwapAnnuity(maturities[i], tenor))
													));
				} catch (CalculationException e) {
					e.printStackTrace();
				}
				squaredError +=(values[i*strikes.length + j]-targetValues[i*strikes.length + j])*(values[i*strikes.length + j]-targetValues[i*strikes.length + j]);
				System.out.print(df.format(values[i*strikes.length + j]) + "       \t");
			}
			System.out.print("\n  \t");
			for(int j = 0; j < strikes.length; j++){
				System.out.print("(" + df.format(Math.abs((values[i*strikes.length + j]-targetValues[i*strikes.length + j]))) + ")     \t");
			}
			System.out.print("\n\n");
		}
		System.out.println("Root mean squared error: " + error.format(Math.sqrt((squaredError/(maturities.length*strikes.length)))));

	}
	

	public static void printTwoFactorHeston(double[] bestParameters, SABRdata sheetdata) throws IOException{
		System.out.println( "\n\tThe best fit paramters for the Two-factor Heston model are: \n\n" +
				"\talphaOne = \t\t" + Math.exp(bestParameters[0]) + 
				"\n\talphaTwo = \t\t" + Math.exp(bestParameters[1]) + 
				"\n\tbetaOne = \t\t" + Math.exp(bestParameters[2]) + 
				"\n\tbetaTwo = \t\t" + Math.exp(bestParameters[3]) + 
				"\n\tsigmaOne = \t\t" + Math.exp(bestParameters[4]) + 
				"\n\tsigmaTwo = \t\t" + Math.exp(bestParameters[5]) + 
				"\n\trhoOne = \t\t" + Math.tanh(0.5*bestParameters[6]) + 
				"\n\trhoTwo = \t\t" + Math.tanh(0.5*bestParameters[7]) + 
				"\n\tvolatilityOne = \t" + Math.exp(bestParameters[8]) + 
				"\n\tvolatilityTwo = \t" + Math.exp(bestParameters[9])
						);


		System.out.println("\nSwaption prices calculated with the calibrated parameters and the error in brackets:\n");
		System.out.println("Maturty" + "\t\t\t\t\t\t\t\t Strikes\n" + "\t-1%      \t-0.5%    \t-0.25%    \t0%      \t0.25%      \t0.5%      \t1%      \t1.5%      \t2%\n");
		
		double squaredError = 0;
		
		for(int i = 0; i < maturities.length; i++){
			System.out.print(maturities[i] + "  \t");
			for(int j = 0; j < strikes.length; j++){
				try {
					values[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(
							sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
							maturities[i],
							strikes[j]+shift,
							sheetdata.getSwapAnnuity(maturities[i], tenor),
							europeanMatrix[i][j].getValue(
															new TwoFactorBatesModelCF(
																	Math.exp(bestParameters[0]),
																	Math.exp(bestParameters[1]),
																	Math.exp(bestParameters[2]),
																	Math.exp(bestParameters[3]),
																	Math.exp(bestParameters[4]),
																	Math.exp(bestParameters[5]),
																	Math.tanh(0.5*bestParameters[6]),
																	Math.tanh(0.5*bestParameters[7]),
																	lambdaZero,
																	lambdaOne,
																	lambdaTwo,
																	k,
																	delta,
																	Math.exp(bestParameters[8]),
																	Math.exp(bestParameters[9]),
																	sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																	0,
																	sheetdata.getSwapAnnuity(maturities[i], tenor))
													));
				} catch (CalculationException e) {
					e.printStackTrace();
				}
				squaredError +=(values[i*strikes.length + j]-targetValues[i*strikes.length + j])*(values[i*strikes.length + j]-targetValues[i*strikes.length + j]);
				System.out.print(df.format(values[i*strikes.length + j]) + "       \t");
			}
			System.out.print("\n  \t");
			for(int j = 0; j < strikes.length; j++){
				System.out.print("(" + df.format(Math.abs((values[i*strikes.length + j]-targetValues[i*strikes.length + j]))) + ")     \t");
			}
			System.out.print("\n\n");
		}
		System.out.println("Root mean squared error: " + error.format(Math.sqrt((squaredError/(maturities.length*strikes.length)))));

	}
	
	public static void printOneFactorBates(double[] bestParameters, SABRdata sheetdata) throws IOException{
		System.out.println( "\n\tThe best fit paramters for the One-factor Bates model are: \n\n" +
				"\talpha = \t\t" + Math.exp(bestParameters[0]) + 
				"\n\tbeta = \t\t\t" + Math.exp(bestParameters[1]) + 
				"\n\tsigma = \t\t" + Math.exp(bestParameters[2]) + 
				"\n\trho = \t\t\t" + Math.tanh(0.5*bestParameters[3]) + 
				"\n\tlambdaZero = \t\t" + Math.exp(bestParameters[4]) + 
				"\n\tlambdaOne = \t\t" + Math.exp(bestParameters[5]) + 
				"\n\tk = \t\t\t" + bestParameters[6] + 
				"\n\tdelta = \t\t" + Math.abs(bestParameters[7]) +
				"\n\tvolatility = \t\t" + Math.exp(bestParameters[8]) 
						);


		System.out.println("\nSwaption prices calculated with the calibrated parameters and the error in brackets:\n");
		System.out.println("Maturty" + "\t\t\t\t\t\t\t\t Strikes\n" + "\t-1%      \t-0.5%    \t-0.25%    \t0%      \t0.25%      \t0.5%      \t1%      \t1.5%      \t2%\n");
		
		double squaredError = 0;
		
		for(int i = 0; i < maturities.length; i++){
			System.out.print(maturities[i] + "  \t");
			for(int j = 0; j < strikes.length; j++){
				try {
					values[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(
							sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
							maturities[i],
							strikes[j]+shift,
							sheetdata.getSwapAnnuity(maturities[i], tenor),
							europeanMatrix[i][j].getValue(
														new TwoFactorBatesModelCF(
																Math.exp(bestParameters[0]),
																Math.exp(bestParameters[1]),
																Math.exp(bestParameters[2]),
																Math.tanh(0.5*bestParameters[3]),
																Math.exp(bestParameters[4]),
																Math.exp(bestParameters[5]),
																bestParameters[6],
																Math.abs(bestParameters[7]),
																Math.exp(bestParameters[8]),//[i],
																sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																0,
																sheetdata.getSwapAnnuity(maturities[i], tenor))
													));
				} catch (CalculationException e) {
					e.printStackTrace();
				}
				squaredError +=(values[i*strikes.length + j]-targetValues[i*strikes.length + j])*(values[i*strikes.length + j]-targetValues[i*strikes.length + j]);
				System.out.print(df.format(values[i*strikes.length + j]) + "       \t");
			}
			System.out.print("\n  \t");
			for(int j = 0; j < strikes.length; j++){
				System.out.print("(" + df.format(Math.abs((values[i*strikes.length + j]-targetValues[i*strikes.length + j]))) + ")     \t");
			}
			System.out.print("\n\n");
		}
		System.out.println("Root mean squared error: " + error.format(Math.sqrt((squaredError/(maturities.length*strikes.length)))));

	}
	
	public static void printOneFactorHeston(double[] bestParameters, SABRdata sheetdata) throws IOException{
		System.out.println( "\n\tThe best fit paramters for the One-factor Heston model are: \n\n" +
				"\talpha = \t\t" + Math.exp(bestParameters[0]) + 
				"\n\tbeta = \t\t\t" + Math.exp(bestParameters[1]) + 
				"\n\tsigma = \t\t" + Math.exp(bestParameters[2]) + 
				"\n\trho = \t\t\t" + Math.tanh(0.5*bestParameters[3]) + 
				"\n\tvolatility = \t\t" + Math.exp(bestParameters[4]) 
						);
						


		System.out.println("\nSwaption prices calculated with the calibrated parameters and the error in brackets:\n");
		System.out.println("Maturty" + "\t\t\t\t\t\t\t\t Strikes\n" + "\t-1%      \t-0.5%    \t-0.25%    \t0%      \t0.25%      \t0.5%      \t1%      \t1.5%      \t2%\n");

		double squaredError = 0;
		
		for(int i = 0; i < maturities.length; i++){
			System.out.print(maturities[i] + "  \t");
			for(int j = 0; j < strikes.length; j++){
				try {
					values[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(
							sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
							maturities[i],
							strikes[j]+shift,
							sheetdata.getSwapAnnuity(maturities[i], tenor),
							europeanMatrix[i][j].getValue(
														new TwoFactorBatesModelCF(
																Math.exp(bestParameters[0]),
																Math.exp(bestParameters[1]),
																Math.exp(bestParameters[2]),
																Math.tanh(0.5*bestParameters[3]),
																lambdaZero,
																lambdaOne,
																k,
																delta,
																Math.exp(bestParameters[4]),
																sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																0,
																sheetdata.getSwapAnnuity(maturities[i], tenor))
													));
				} catch (CalculationException e) {
					e.printStackTrace();
				}
				squaredError +=(values[i*strikes.length + j]-targetValues[i*strikes.length + j])*(values[i*strikes.length + j]-targetValues[i*strikes.length + j]);
				System.out.print(df.format(values[i*strikes.length + j]) + "       \t");
			}
			System.out.print("\n  \t");
			for(int j = 0; j < strikes.length; j++){
				System.out.print("(" + df.format(Math.abs((values[i*strikes.length + j]-targetValues[i*strikes.length + j]))) + ")     \t");
			}
			System.out.print("\n\n");
		}
		System.out.println("Root mean squared error: " + error.format(Math.sqrt((squaredError/(maturities.length*strikes.length)))));

	}
	
	public static double printSABR(double[] bestParameters, SABRdata sheetdata) throws IOException{
		System.out.println( "\n\tThe best fit paramters for the SABR model are: \n\n" +
				"\talpha = \t\t" + Math.exp(bestParameters[0]) + 
				"\n\tbeta  = \t\t" + Math.exp(-bestParameters[1]*bestParameters[1]) + 
				"\n\trho = \t\t\t" + Math.tanh(0.5*bestParameters[2]) + 
				"\n\tv = \t\t\t" + Math.exp(bestParameters[3])
						);
		
		System.out.println("\nSwaption prices calculated with the calibrated parameters and the error in brackets:\n");
		System.out.println("Maturty" + "\t\t\t\t\t\t\t\t Strikes\n" + "\t-1%      \t-0.5%    \t-0.25%    \t0%      \t0.25%      \t0.5%      \t1%      \t1.5%      \t2%\n");
	
		
		return printSABR(bestParameters, sheetdata, maturities, targetValues);
	}
	
	public static double printSABR(double[] bestParameters, SABRdata sheetdata, int[] maturities, double[] targetValues) throws IOException{
		
		double squaredError = 0;
		for(int i = 0; i < maturities.length; i++){
			System.out.print(maturities[i] + "  \t");
			for(int j = 0; j < strikes.length; j++){
				values[i*strikes.length + j] = AnalyticFormulas.sabrBerestyckiNormalVolatilityApproximation(
													Math.exp(bestParameters[0]),
													Math.exp(-bestParameters[1]*bestParameters[1]),
													Math.tanh(0.5*bestParameters[2]),
													Math.exp(bestParameters[3]),
													sabrDisplacement,
													sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
													strikes[j]+shift,
													maturities[i]);
 				
				squaredError +=(values[i*strikes.length + j]-targetValues[i*strikes.length + j])*(values[i*strikes.length + j]-targetValues[i*strikes.length + j]);
				System.out.print(df.format(values[i*strikes.length + j]) + "       \t");
			}
			System.out.print("\n  \t");
			for(int j = 0; j < strikes.length; j++){
				System.out.print("(" + df.format(Math.abs((values[i*strikes.length + j]-targetValues[i*strikes.length + j]))) + ")     \t");
			}
			System.out.print("\n\n");
		}
		return squaredError;

	}
	

	
	
	
}
