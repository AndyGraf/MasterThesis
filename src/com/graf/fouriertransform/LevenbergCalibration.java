package com.graf.fouriertransform;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Scanner;

import net.finmath.exception.CalculationException;
import net.finmath.fouriermethod.products.AbstractProductFourierTransform;
import net.finmath.fouriermethod.products.EuropeanOption;
import net.finmath.optimizer.LevenbergMarquardt;
import net.finmath.optimizer.SolverException;

public class LevenbergCalibration {
	
	//maximum number of iterations for the LM-algorithm
	static int maxIteration = 100;
	
	static DecimalFormat df = new DecimalFormat("#.###");
    
	//defining the options
	static int[] maturities = {1, 2, 5, 10, 15, 20};
	static double[] strikes = {-0.01,-0.005, -0.0025, 0, 0.0025, 0.005, 0.01, 0.015, 0.02};
	static int tenor;
	static double shift;
	
	
	static double [] targetValues = new double[maturities.length*strikes.length];
	static double [] targetVolatilities = new double[maturities.length*strikes.length];
	static double [] weights = new double[maturities.length*strikes.length];
	static double [] values = new double[maturities.length*strikes.length];

	static AbstractProductFourierTransform[][] europeanMatrix = new AbstractProductFourierTransform[maturities.length][strikes.length];
	
    
    static double lambdaZero	= 0;
    static double lambdaOne		= 0;
    static double lambdaTwo		= 0;
    static double k				= 0;
    static double delta			= 0;
    
    //set beta equal to 1 since the market volatilities are log-normal
    static double sabrBeta	    	= 1;
	

	public static void main(String[] args) throws IOException {

		
		SABRdata sheetdata = new SABRdata();
		
		System.out.println( "Model to calibrate: \n" + 
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
		if(tenor == 2){shift = 0.0265;}else if(tenor == 5){shift = 0.016;}else if(tenor == 10){shift = 0.015;}else{shift = 0;};
		double[][] volatilities = sheetdata.getSmileData(tenor);
		for(int w = 0; w < weights.length; w++){weights[w] = 1;};


		//calculate the target values
		System.out.println("Target swaption prices in % of the notional principal amount with a tenor of " + tenor + ":\n");
		System.out.println("Maturty" + "\t\t\t\t\t\t\t\t Strikes\n" + "\t-1%      \t-0.5%    \t-0.25%    \t0%      \t0.25%      \t0.5%      \t1%      \t1.5%      \t2%\n");
		for(int i = 0; i < maturities.length; i++){
			System.out.print(maturities[i] + "  \t");
			for(int j = 0; j < strikes.length; j++){
				targetVolatilities[i*strikes.length + j] = volatilities[i][j];
				targetValues[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.blackScholesGeneralizedOptionValue(sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,targetVolatilities[i*strikes.length + j],maturities[i], strikes[j]+shift,sheetdata.getSwapAnnuity(maturities[i], tenor));
				europeanMatrix[i][j] = new EuropeanOption(maturities[i], strikes[j]+shift);
				System.out.print(df.format(targetValues[i*strikes.length + j]*100) + "     \t");
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
				System.out.println("Root mean squared error: " + df.format(Math.sqrt((printSABR(sabrParameters, sheetdata)/(maturities.length*strikes.length)))));
				break;
			case 6:
				double error = 0;
				System.out.println("\nSwaption prices calculated with a maturity-wise calibration and the error in brackets:\n");
				System.out.println("Maturty" + "\t\t\t\t\t\t\t\t Strikes\n" + "\t-1%      \t-0.5%    \t-0.25%    \t0%      \t0.25%      \t0.5%      \t1%      \t1.5%      \t2%\n");
				double[] mTargetValues = new double[strikes.length];
				double[] mTargetVolatilities = new double[strikes.length];
				for(int i = 0; i<maturities.length; i++){
					int[] maturity = {maturities[i]};
					for(int j = 0; j < strikes.length; j++){
						mTargetValues[j] = targetValues[i*strikes.length + j];
						mTargetVolatilities[j] = volatilities[i][j];
					}
					double[] mSabrParameters = sabr(sheetdata, maturity, mTargetVolatilities);
					error += printSABR(mSabrParameters, sheetdata, maturity, mTargetValues);
				}
				System.out.println("Root mean squared error: " + df.format(Math.sqrt((error/(maturities.length*strikes.length)))));
				break;
		}

		

	  	

		

		

		

	}
	
	
	
	
	public static double[] twoFactorBates(SABRdata sheetdata){
		
		//initial parameters  
		
        double alphaOne     = 0.01;
        double alphaTwo		= 0.04;
        double betaOne      = 0.91;
        double betaTwo      = 1.76;
        double sigmaOne     = 0.582;
        double sigmaTwo     = 0.346;
        double rhoOne       = -0.848;
        double rhoTwo       = -0.402;
        
        double lambdaZero   = 0.1143;
        double lambdaOne    = 81.56;
        double lambdaTwo    = 0.28;
        double k       		= -0.057;
        double delta   		= 0.102;
        
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
    			Math.log(lambdaZero),
    			Math.log(lambdaOne),
    			Math.log(lambdaTwo),
    			Math.log(-k),
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
							values[i*strikes.length + j] = europeanMatrix[i][j].getValue(
																new TwoFactorBatesModelCF(
																Math.exp(parameters[0]),
																Math.exp(parameters[1]),
																Math.exp(parameters[2]),
																Math.exp(parameters[3]),																	
																Math.exp(parameters[4]),
																Math.exp(parameters[5]),
																parameters[6],
																parameters[7],
																Math.exp(parameters[8]),
																Math.exp(parameters[9]),
																Math.exp(parameters[10]),
																-Math.exp(parameters[11]),
																parameters[12],
																Math.exp(parameters[13]),
																Math.exp(parameters[14]),
																sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																0,
																sheetdata.getSwapAnnuity(maturities[i], tenor))
															)
							;
						} catch (IOException e) {
							e.printStackTrace();
						} catch (CalculationException e) {
							// TODO Auto-generated catch block
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
        double rhoOne       = -0.775;
        double rhoTwo       = -0.382;
        
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
							values[i*strikes.length + j] = europeanMatrix[i][j].getValue(
																new TwoFactorBatesModelCF(
																Math.exp(parameters[0]),
																Math.exp(parameters[1]),
																Math.exp(parameters[2]),
																Math.exp(parameters[3]),																	
																Math.exp(parameters[4]),
																Math.exp(parameters[5]),
																parameters[6],
																parameters[7],
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
							;
						} catch (IOException e) {
							e.printStackTrace();
						} catch (CalculationException e) {
							// TODO Auto-generated catch block
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
	  
	  	return optimizer.getBestFitParameters();
        
	}
	
	public static double[] oneFactorBates(SABRdata sheetdata){
		
		//initial parameters  
        
        double alphaOne     = 0.049;
        double betaOne      = 2.45;
        double sigmaOne     = 0.378;
        double rhoOne       = -0.545;
        
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
							values[i*strikes.length + j] = europeanMatrix[i][j].getValue(
																new TwoFactorBatesModelCF(
																Math.exp(parameters[0]),
																Math.exp(parameters[1]),
																Math.exp(parameters[2]),
																parameters[3],
																Math.exp(parameters[4]),
																Math.exp(parameters[5]),
																parameters[6],
																parameters[7],
																Math.exp(parameters[8]),
																sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																0,
																sheetdata.getSwapAnnuity(maturities[i], tenor))
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
	  
	  	return optimizer.getBestFitParameters();
        
	}
	
	
	public static double[] oneFactorHeston(SABRdata sheetdata){
		
		//initial parameters  
        
        double alpha		= 0.1;
        double beta			= 1.49;
        double sigma		= 0.742;
        double rho			= -0.571;

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
							values[i*strikes.length + j] = europeanMatrix[i][j].getValue(
																new TwoFactorBatesModelCF(
																Math.exp(parameters[0]),
																Math.exp(parameters[1]),
																Math.exp(parameters[2]),
																parameters[3],
																lambdaZero,
																lambdaOne,
																k,
																delta,
																Math.exp(parameters[4]),
																sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																0,
																sheetdata.getSwapAnnuity(maturities[i], tenor))
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
	  
	  	return optimizer.getBestFitParameters();
        
	}
	
	public static double[] sabr(SABRdata sheetdata) throws IOException{
		
		return sabr(sheetdata, maturities, targetVolatilities);
	}
	
	public static double[] sabr(SABRdata sheetdata, int [] maturities, double[] targetVolatilities) throws IOException{
		
		//initial parameters 
		
        double alpha	    = 0.926;
        double rho      	= -0.52;
        double sigma   		= 0.073;
		

        
		double[] initialParameters = new double[]{			
        		alpha,
        		rho,
        		Math.log(sigma)
            };
		
        
        LevenbergMarquardt optimizer = new LevenbergMarquardt() {

		 	public void setValues(double[] parameters, double[] values) {
		 		for(int i = 0; i < maturities.length; i++){
					for(int j = 0; j < strikes.length; j++){
							try {
								values[i*strikes.length + j] = SABRVolatility.getVolatility(
																	parameters[0],
																	sabrBeta,
																	parameters[1],
																	Math.exp(parameters[2]),
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
	  
	  	return optimizer.getBestFitParameters();
	}
	
	
	
	public static void printTwoFactorBates(double[] bestParameters, SABRdata sheetdata) throws IOException{
		System.out.println( "\n\tThe best fit paramters for the Two-factor Bates model are: \n\n" +
				"\talphaOne = \t\t" + Math.exp(bestParameters[0]) + 
				"\n\talphaTwo = \t\t" + Math.exp(bestParameters[1]) + 
				"\n\tbetaOne = \t\t" + Math.exp(bestParameters[2]) + 
				"\n\tbetaTwo = \t\t" + Math.exp(bestParameters[3]) + 
				"\n\tsigmaOne = \t\t" + Math.exp(bestParameters[4]) + 
				"\n\tsigmaTwo = \t\t" + Math.exp(bestParameters[5]) + 
				"\n\trhoOne = \t\t" + bestParameters[6] + 
				"\n\trhoTwo = \t\t" + bestParameters[7] + 
				"\n\tlambdaZero = \t\t" + Math.exp(bestParameters[8]) + 
				"\n\tlambdaOne = \t\t" + Math.exp(bestParameters[9])  +
				"\n\tlambdaTwo = \t\t" + Math.exp(bestParameters[10]) + 
				"\n\tk = \t\t\t" + -Math.exp(bestParameters[11]) + 
				"\n\tdelta = \t\t" + bestParameters[12] +
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
					values[i*strikes.length + j] = europeanMatrix[i][j].getValue(
															new TwoFactorBatesModelCF(
																	Math.exp(bestParameters[0]),
																	Math.exp(bestParameters[1]),
																	Math.exp(bestParameters[2]),
																	Math.exp(bestParameters[3]),
																	Math.exp(bestParameters[4]),
																	Math.exp(bestParameters[5]),
																	bestParameters[6],
																	bestParameters[7],
																	Math.exp(bestParameters[8]),
																	Math.exp(bestParameters[9]),
																	Math.exp(bestParameters[10]),
																	-Math.exp(bestParameters[11]),
																	bestParameters[12],
																	Math.exp(bestParameters[13]),
																	Math.exp(bestParameters[14]),
																	sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																	0,
																	sheetdata.getSwapAnnuity(maturities[i], tenor))
													);
				} catch (CalculationException e) {
					e.printStackTrace();
				}
				squaredError +=100*(values[i*strikes.length + j]-targetValues[i*strikes.length + j])*100*(values[i*strikes.length + j]-targetValues[i*strikes.length + j]);
				System.out.print(df.format(values[i*strikes.length + j]*100) + "       \t");
			}
			System.out.print("\n  \t");
			for(int j = 0; j < strikes.length; j++){
				System.out.print("(" + df.format(100*Math.abs((values[i*strikes.length + j]-targetValues[i*strikes.length + j]))) + ")     \t");
			}
			System.out.print("\n\n");
		}
		System.out.println("Root mean squared error: " + df.format(Math.sqrt((squaredError/(maturities.length*strikes.length)))));

	}
	

	public static void printTwoFactorHeston(double[] bestParameters, SABRdata sheetdata) throws IOException{
		System.out.println( "\n\tThe best fit paramters for the Two-factor Heston model are: \n\n" +
				"\talphaOne = \t\t" + Math.exp(bestParameters[0]) + 
				"\n\talphaTwo = \t\t" + Math.exp(bestParameters[1]) + 
				"\n\tbetaOne = \t\t" + Math.exp(bestParameters[2]) + 
				"\n\tbetaTwo = \t\t" + Math.exp(bestParameters[3]) + 
				"\n\tsigmaOne = \t\t" + Math.exp(bestParameters[4]) + 
				"\n\tsigmaTwo = \t\t" + Math.exp(bestParameters[5]) + 
				"\n\trhoOne = \t\t" + bestParameters[6] + 
				"\n\trhoTwo = \t\t" + bestParameters[7] + 
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
					values[i*strikes.length + j] = europeanMatrix[i][j].getValue(
															new TwoFactorBatesModelCF(
																	Math.exp(bestParameters[0]),
																	Math.exp(bestParameters[1]),
																	Math.exp(bestParameters[2]),
																	Math.exp(bestParameters[3]),
																	Math.exp(bestParameters[4]),
																	Math.exp(bestParameters[5]),
																	bestParameters[6],
																	bestParameters[7],
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
													);
				} catch (CalculationException e) {
					e.printStackTrace();
				}
				squaredError +=100*(values[i*strikes.length + j]-targetValues[i*strikes.length + j])*100*(values[i*strikes.length + j]-targetValues[i*strikes.length + j]);
				System.out.print(df.format(values[i*strikes.length + j]*100) + "       \t");
			}
			System.out.print("\n  \t");
			for(int j = 0; j < strikes.length; j++){
				System.out.print("(" + df.format(100*Math.abs((values[i*strikes.length + j]-targetValues[i*strikes.length + j]))) + ")     \t");
			}
			System.out.print("\n\n");
		}
		System.out.println("Root mean squared error: " + df.format(Math.sqrt((squaredError/(maturities.length*strikes.length)))));

	}
	
	public static void printOneFactorBates(double[] bestParameters, SABRdata sheetdata) throws IOException{
		System.out.println( "\n\tThe best fit paramters for the One-factor Bates model are: \n\n" +
				"\talpha = \t\t" + Math.exp(bestParameters[0]) + 
				"\n\tbeta = \t\t\t" + Math.exp(bestParameters[1]) + 
				"\n\tsigma = \t\t" + Math.exp(bestParameters[2]) + 
				"\n\trho = \t\t\t" + bestParameters[3] + 
				"\n\tlambdaZero = \t\t" + Math.exp(bestParameters[4]) + 
				"\n\tlambdaOne = \t\t" + Math.exp(bestParameters[5]) + 
				"\n\tk = \t\t\t" + bestParameters[6] + 
				"\n\tdelta = \t\t" + bestParameters[7] +
				"\n\tvolatility = \t\t" + Math.exp(bestParameters[8]) 
						);


		System.out.println("\nSwaption prices calculated with the calibrated parameters and the error in brackets:\n");
		System.out.println("Maturty" + "\t\t\t\t\t\t\t\t Strikes\n" + "\t-1%      \t-0.5%    \t-0.25%    \t0%      \t0.25%      \t0.5%      \t1%      \t1.5%      \t2%\n");
		double squaredError = 0;
		for(int i = 0; i < maturities.length; i++){
			System.out.print(maturities[i] + "  \t");
			for(int j = 0; j < strikes.length; j++){
				try {
					values[i*strikes.length + j] = europeanMatrix[i][j].getValue(
														new TwoFactorBatesModelCF(
																Math.exp(bestParameters[0]),
																Math.exp(bestParameters[1]),
																Math.exp(bestParameters[2]),
																bestParameters[3],
																Math.exp(bestParameters[4]),
																Math.exp(bestParameters[5]),
																bestParameters[6],
																bestParameters[7],
																Math.exp(bestParameters[8]),//[i],
																sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																0,
																sheetdata.getSwapAnnuity(maturities[i], tenor))
													);
				} catch (CalculationException e) {
					e.printStackTrace();
				}
				squaredError +=100*(values[i*strikes.length + j]-targetValues[i*strikes.length + j])*100*(values[i*strikes.length + j]-targetValues[i*strikes.length + j]);
				System.out.print(df.format(values[i*strikes.length + j]*100) + "       \t");
			}
			System.out.print("\n  \t");
			for(int j = 0; j < strikes.length; j++){
				System.out.print("(" + df.format(100*Math.abs((values[i*strikes.length + j]-targetValues[i*strikes.length + j]))) + ")     \t");
			}
			System.out.print("\n\n");
		}
		System.out.println("Root mean squared error: " + df.format(Math.sqrt((squaredError/(maturities.length*strikes.length)))));

	}
	
	public static void printOneFactorHeston(double[] bestParameters, SABRdata sheetdata) throws IOException{
		System.out.println( "\n\tThe best fit paramters for the One-factor Heston model are: \n\n" +
				"\talpha = \t\t" + Math.exp(bestParameters[0]) + 
				"\n\tbeta = \t\t\t" + Math.exp(bestParameters[1]) + 
				"\n\tsigma = \t\t" + Math.exp(bestParameters[2]) + 
				"\n\trho = \t\t\t" + bestParameters[3] + 
				"\n\tvolatility = \t\t" + Math.exp(bestParameters[4]) 
						);
						


		System.out.println("\nSwaption prices calculated with the calibrated parameters and the error in brackets:\n");
		System.out.println("Maturty" + "\t\t\t\t\t\t\t\t Strikes\n" + "\t-1%      \t-0.5%    \t-0.25%    \t0%      \t0.25%      \t0.5%      \t1%      \t1.5%      \t2%\n");

		double squaredError = 0;
		for(int i = 0; i < maturities.length; i++){
			System.out.print(maturities[i] + "  \t");
			for(int j = 0; j < strikes.length; j++){
				try {
					values[i*strikes.length + j] = europeanMatrix[i][j].getValue(
														new TwoFactorBatesModelCF(
																Math.exp(bestParameters[0]),
																Math.exp(bestParameters[1]),
																Math.exp(bestParameters[2]),
																bestParameters[3],
																lambdaZero,
																lambdaOne,
																k,
																delta,
																Math.exp(bestParameters[4]),
																sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																0,
																sheetdata.getSwapAnnuity(maturities[i], tenor))
													);
				} catch (CalculationException e) {
					e.printStackTrace();
				}
				squaredError +=100*(values[i*strikes.length + j]-targetValues[i*strikes.length + j])*100*(values[i*strikes.length + j]-targetValues[i*strikes.length + j]);
				System.out.print(df.format(values[i*strikes.length + j]*100) + "       \t");
			}
			System.out.print("\n  \t");
			for(int j = 0; j < strikes.length; j++){
				System.out.print("(" + df.format(100*Math.abs((values[i*strikes.length + j]-targetValues[i*strikes.length + j]))) + ")     \t");
			}
			System.out.print("\n\n");
		}
		System.out.println("Root mean squared error: " + df.format(Math.sqrt((squaredError/(maturities.length*strikes.length)))));

	}
	
	public static double printSABR(double[] bestParameters, SABRdata sheetdata) throws IOException{
		System.out.println( "\n\tThe best fit paramters for the SABR model are: \n\n" +
				"\talpha = \t\t" + bestParameters[0] + 
				"\n\trho = \t\t\t" + bestParameters[1] + 
				"\n\tv = \t\t\t" + Math.exp(bestParameters[2])
						);
		
		System.out.println("\nSwaption prices calculated with the calibrated parameters and the error in brackets:\n");
		System.out.println("Maturty" + "\t\t\t\t\t\t\t\t Strikes\n" + "\t-1%      \t-0.5%    \t-0.25%    \t0%      \t0.25%      \t0.5%      \t1%      \t1.5%      \t2%\n");
	
		
		return printSABR(bestParameters, sheetdata, maturities, targetValues);
	}
	
	public static double printSABR(double[] bestParameters, SABRdata sheetdata, int[] maturities, double[] targetValues) throws IOException{
		
		double [] prices = new double[maturities.length*strikes.length];
		double squaredError = 0;
		for(int i = 0; i < maturities.length; i++){
			System.out.print(maturities[i] + "  \t");
			for(int j = 0; j < strikes.length; j++){
				values[i*strikes.length + j] = SABRVolatility.getVolatility(
													bestParameters[0],
													sabrBeta,
													bestParameters[1],
													Math.exp(bestParameters[2]),
													sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
													strikes[j]+shift,
													maturities[i]);
 				prices[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.blackScholesGeneralizedOptionValue(sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,values[i*strikes.length + j],maturities[i], strikes[j]+shift,sheetdata.getSwapAnnuity(maturities[i], tenor));
				squaredError +=100*(prices[i*strikes.length + j]-targetValues[i*strikes.length + j])*100*(prices[i*strikes.length + j]-targetValues[i*strikes.length + j]);
				System.out.print(df.format(prices[i*strikes.length + j]*100) + "       \t");
			}
			System.out.print("\n  \t");
			for(int j = 0; j < strikes.length; j++){
				System.out.print("(" + df.format(100*Math.abs((prices[i*strikes.length + j]-targetValues[i*strikes.length + j]))) + ")     \t");
			}
			System.out.print("\n\n");
		}
		return squaredError;

	}
	
	public static void maturityWiseSABR(SABRdata sheetdata) throws IOException{
		
		
		
	}
	
	
	
}
