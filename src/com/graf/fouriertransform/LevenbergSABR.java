package com.graf.fouriertransform;

import java.io.IOException;
import java.text.DecimalFormat;

import net.finmath.optimizer.LevenbergMarquardt;
import net.finmath.optimizer.SolverException;

public class LevenbergSABR {

	
		static DecimalFormat df = new DecimalFormat("#.###");
		
		
		
		 //initial parameters  
        
        static double alpha	    = 0.926;
        static double rho      	= -0.52;
        static double v    		= 0.73;
		
        //set known constants
        static double beta	    	= 1;
        
		
		static int tenor = 10;
		static double shift;
		
	public static void main(String[] args) throws IOException {
			
			
		int[] maturities = {1, 2, 5, 10, 15, 20};
		double[] strikes = {-0.01,-0.005, -0.0025, 0, 0.0025, 0.005, 0.01, 0.015, 0.02};
		
		
		double mse = calculate(maturities, strikes, true);
		
		double maturityWiseMSE = 0;
		for(int a = 0; a<maturities.length;a++){
			int[] maturity = {maturities[a]};
			maturityWiseMSE += calculate(maturity, strikes, false);
		}
		
		System.out.println("RMSE = " + Math.sqrt((mse/(maturities.length*strikes.length))));
		
		System.out.println("RMSE of maturity recalibrated SABR = " + Math.sqrt((maturityWiseMSE/(maturities.length*strikes.length))));
	}
		
	
	
	
	
	
	
	public static double calculate(int[] maturities, double[]strikes, boolean print) throws IOException{
		
		
		
		SABRdata sheetdata = new SABRdata();
		if(tenor == 2){shift = 0.0265;}else if(tenor == 5){shift = 0.016;}else if(tenor == 10){shift = 0.015;}else{shift = 0;};
		double[][] volatilities = sheetdata.getSmileData(tenor);
		
		double[] initialParameters = new double[]{			
        		alpha,
        		rho,
        		Math.log(v)
            };
		
		
		double [] targetValues = new double[maturities.length*strikes.length];
		double [] targetPrices = new double[maturities.length*strikes.length];
		
		double [] weights = new double[maturities.length*strikes.length];

		for(int w = 0; w < weights.length; w++){weights[w] = 1;};
		if (print) System.out.println("Maturty");
		for(int i = 0; i < maturities.length; i++){
			if(print) System.out.print(maturities[i] + "  \t");
			for(int j = 0; j < strikes.length; j++){
				targetValues[i*strikes.length + j] = volatilities[i][j];
				targetPrices[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.blackScholesGeneralizedOptionValue(sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,volatilities[i][j],maturities[i], strikes[j]+shift,sheetdata.getSwapAnnuity(maturities[i], tenor));
				if(print) System.out.print(df.format(targetPrices[i*strikes.length + j]*100) + "    \t");
			}
			if(print) System.out.println();
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
	  	

	  	if(print){
	  	
	  	System.out.println("The solver for problem 1 required " + optimizer.getIterations() + " iterations. The best fit parameters are:");
		System.out.println( "\talpha = \t\t" + bestParameters[0] + 
							"\n\trho = \t\t\t" + bestParameters[1] + 
							"\n\tv = \t\t\t" + Math.exp(bestParameters[2])
									);
		System.out.println("\n Factor 1:" + df.format(bestParameters[0]) + "&\t" +
				df.format(bestParameters[1]) + "&\t" +
				df.format(bestParameters[2]) + "&&&&&&&&&"			
				);
	  	}
		
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
				if(print)System.out.print("&" + df.format(prices[i*strikes.length + j]*100) + "\t");
			}
			if(print)System.out.print("\\\\ \n");
 		}
		
		return mse;
	}
		
		
	
	
	
	
	
}


