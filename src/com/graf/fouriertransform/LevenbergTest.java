package com.graf.fouriertransform;

import java.io.IOException;

import net.finmath.fouriermethod.models.ProcessCharacteristicFunctionInterface;
import net.finmath.fouriermethod.products.AbstractProductFourierTransform;
import net.finmath.optimizer.LevenbergMarquardt;
import net.finmath.optimizer.SolverException;

public class LevenbergTest {

	public static void main(String[] args) throws IOException {
        
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
		
        //set known constants
        
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
			k,
			delta,
			volatilityOne,
			volatilityTwo
        };
        
        
		SABRdata sheetdata = new SABRdata();
		
		int tenor = 10;
		double shift;
		if(tenor == 2){shift = 2.65;}else if(tenor == 5){shift = 1.6;}else if(tenor == 10){shift = 1.5;}else{shift = 0;};
		double[][] volatilities = sheetdata.getSmileData(tenor);
		int[] maturities = {1, 2, 5, 10, 15, 20};
//		int[] maturities = {5};
		double[] strikes = {-1,-0.5, -0.25, 0, 0.25, 0.5, 1, 1.5, 2};
		AbstractProductFourierTransform[][] europeanMatrix = new AbstractProductFourierTransform[maturities.length][strikes.length];
		
		double [] targetValues = new double[maturities.length*strikes.length];
//		double [] bachelierValues = new double[maturities.length*strikes.length];
		
		double [] weights = new double[maturities.length*strikes.length];
		for(int w = 0; w < weights.length; w++){weights[w] = 1;};

		for(int i = 0; i < maturities.length; i++){
			for(int j = 0; j < strikes.length; j++){
				targetValues[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.blackScholesGeneralizedOptionValue(sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,volatilities[i][j],maturities[i], strikes[j]+shift,sheetdata.getSwapAnnuity(maturities[i], tenor));
//				bachelierValues[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.bachelierOptionValue(forwardcurve[index][2]+1.5,volatilities[i][j],maturities[i], strikes[j]+1.5,sheetdata.getSwapAnnuity(maturities[i], tenor));
				europeanMatrix[i][j] = new EuropeanOption(maturities[i], strikes[j]+shift, sheetdata.getSwapAnnuity(maturities[i], tenor));
				System.out.print(targetValues[i*strikes.length + j] + "\t");
			}
			System.out.println();
		}
        
        
        
        
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
																Math.exp(parameters[8]),
			 													Math.exp(parameters[9]),
			 													Math.exp(parameters[10]),
																parameters[11],
																parameters[12],
																parameters[13],
																parameters[14],
//																volatilityOne,
//																volatilityTwo,
																sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																-Math.log(sheetdata.getSwapAnnuity(maturities[i], tenor)/tenor)/maturities[i])
															)
							;
							
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
//						System.out.println((values[i*strikes.length + j]-targetValues[i*strikes.length + j])*(values[i*strikes.length + j]-targetValues[i*strikes.length + j]));
					}
		 		}
			}
		};
//		optimizer.setLambda(0.001);
	  	optimizer.setInitialParameters(initialParameters);
	  	optimizer.setWeights(weights);
	  	optimizer.setMaxIteration(5);
	  	
	  	optimizer.setTargetValues(targetValues);
	  
	  	try {
			optimizer.run();
		} catch (SolverException e) {
			e.printStackTrace();
		}
	  
	  	double[] bestParameters = optimizer.getBestFitParameters();
	  	

	  
	  	
	  	System.out.println("The solver for problem 1 required " + optimizer.getIterations() + " iterations. The best fit parameters are:");
		System.out.println( "\talphaOne = \t\t" + Math.exp(bestParameters[0]) + 
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
							"\n\tk = \t\t\t" + bestParameters[11] + 
							"\n\tdelta = \t\t" + bestParameters[12] +
							"\n\tvolatilityOne = \t" + bestParameters[13] + 
							"\n\tvolatilityTwo = \t" + bestParameters[14]);
		
		double [] values = new double[maturities.length*strikes.length];
		double mse = 0;
		for(int i = 0; i < maturities.length; i++){
			for(int j = 0; j < strikes.length; j++){
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
 													bestParameters[11],
 													bestParameters[12],
 													bestParameters[13],
 													bestParameters[14],
// 													volatilityOne,
// 													volatilityTwo,
 													sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
 													-Math.log(sheetdata.getSwapAnnuity(maturities[i], tenor)/tenor)/maturities[i])
 												)
 				;
				mse +=(values[i*strikes.length + j]-targetValues[i*strikes.length + j])*(values[i*strikes.length + j]-targetValues[i*strikes.length + j]);
				System.out.print(values[i*strikes.length + j] + "\t");
			}
			System.out.println();
 		}
		System.out.println("RMSE =" + Math.sqrt((mse/(maturities.length*strikes.length))));
	}
}


