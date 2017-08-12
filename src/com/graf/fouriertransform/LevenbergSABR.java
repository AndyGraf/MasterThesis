package com.graf.fouriertransform;

import java.io.IOException;

import net.finmath.fouriermethod.models.ProcessCharacteristicFunctionInterface;
import net.finmath.fouriermethod.products.AbstractProductFourierTransform;
import net.finmath.optimizer.LevenbergMarquardt;
import net.finmath.optimizer.SolverException;

public class LevenbergSABR {

	public static void main(String[] args) throws IOException {
        
		SABRdata sheetdata = new SABRdata();
		
		//initial parameters  
        
        double alpha	    = 0.926;
        double rho      	= -0.52;
        double v    		= 0.73;
 
        
        
		
        //set known constants
        double beta	    	= 1;
        
        
 
        
		
		
		int tenor = 2;
		double shift;
		if(tenor == 2){shift = 2.65;}else if(tenor == 5){shift = 1.6;}else if(tenor == 10){shift = 1.5;}else{shift = 0;};
		double[][] volatilities = sheetdata.getSmileData(tenor);
		
//		double[] volatilityOne = new double[volatilities[3].length];
//		for(int i=0; i<volatilities[3].length;i++){
//			volatilityOne[i]   = volatilities[3][i]*volatilities[3][i];
//		}
		
        double[] initialParameters = new double[]{			
        		alpha,
        		rho,
        		Math.log(v)
            };
		
		int[] maturities = {1, 2, 5, 10, 15, 20};
//		int[] maturities = {5};
		double[] strikes = {-1,-0.5, -0.25, 0, 0.25, 0.5, 1, 1.5, 2};
		AbstractProductFourierTransform[][] europeanMatrix = new AbstractProductFourierTransform[maturities.length][strikes.length];
		
		double [] targetValues = new double[maturities.length*strikes.length];
		double [] targetPrices = new double[maturities.length*strikes.length];
		
		double [] weights = new double[maturities.length*strikes.length];
		for(int w = 0; w < weights.length; w++){weights[w] = 1;};

		for(int i = 0; i < maturities.length; i++){
			for(int j = 0; j < strikes.length; j++){
				targetValues[i*strikes.length + j] = volatilities[i][j];
				targetPrices[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.blackScholesGeneralizedOptionValue(sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,volatilities[i][j],maturities[i], strikes[j]+shift,sheetdata.getSwapAnnuity(maturities[i], tenor));
				System.out.print(targetPrices[i*strikes.length + j] + "\t");
			}
			System.out.println();
		}
        
        
        
        
	 	LevenbergMarquardt optimizer = new LevenbergMarquardt() {

		 	public void setValues(double[] parameters, double[] values) {
		 		for(int i = 0; i < maturities.length; i++){
					for(int j = 0; j < strikes.length; j++){
							try {
								values[i*strikes.length + j] = SABRVolatility.getVolatility(
																	parameters[0],
																	beta,
																	parameters[1],
																	Math.exp(parameters[2]),
																	sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																	strikes[j]+shift,
																	maturities[i]);
							} catch (IOException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
					}
		 		}
			}
		};
//		optimizer.setLambda(0.001);
	  	optimizer.setInitialParameters(initialParameters);
	  	optimizer.setWeights(weights);
	  	optimizer.setMaxIteration(100);
	  	
	  	optimizer.setTargetValues(targetValues);
	  
	  	try {
			optimizer.run();
		} catch (SolverException e) {
			e.printStackTrace();
		}
	  
	  	double[] bestParameters = optimizer.getBestFitParameters();
	  	

	  
	  	
	  	System.out.println("The solver for problem 1 required " + optimizer.getIterations() + " iterations. The best fit parameters are:");
		System.out.println( "\talpha = \t\t" + bestParameters[0] + 
							"\n\trho = \t\t\t" + bestParameters[1] + 
							"\n\tv = \t\t\t" + Math.exp(bestParameters[2])
									);
		
		double [] values = new double[maturities.length*strikes.length];
		double [] prices = new double[maturities.length*strikes.length];
		double mse = 0;
		for(int i = 0; i < maturities.length; i++){
			for(int j = 0; j < strikes.length; j++){
				values[i*strikes.length + j] = SABRVolatility.getVolatility(
													bestParameters[0],
													beta,
													bestParameters[1],
													Math.exp(bestParameters[2]),
													sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
													strikes[j]+shift,
													maturities[i]);
 				;
 				prices[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.blackScholesGeneralizedOptionValue(sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,values[i*strikes.length + j],maturities[i], strikes[j]+shift,sheetdata.getSwapAnnuity(maturities[i], tenor));
				mse +=(prices[i*strikes.length + j]-targetPrices[i*strikes.length + j])*(prices[i*strikes.length + j]-targetPrices[i*strikes.length + j]);
				System.out.print(prices[i*strikes.length + j] + "\t");
			}
			System.out.println();
 		}
		System.out.println("RMSE =" + Math.sqrt((mse/(maturities.length*strikes.length))));
		
	}
	
	
	
	
}


