package com.graf.fouriertransform;

import java.io.IOException;
import java.text.DecimalFormat;

import net.finmath.fouriermethod.models.ProcessCharacteristicFunctionInterface;
import net.finmath.fouriermethod.products.AbstractProductFourierTransform;
import net.finmath.optimizer.LevenbergMarquardt;
import net.finmath.optimizer.SolverException;

public class LevenbergTwoHeston {

	public static void main(String[] args) throws IOException {
		DecimalFormat df = new DecimalFormat("#.###");
        
		SABRdata sheetdata = new SABRdata();
		
		//initial parameters  
        double alphaOne     = 0.028;
        double alphaTwo		= 0.130;
        double betaOne      = 0.0;
        double betaTwo      = 5.58;
//        double sigmaOne     = 0.582;
//        double sigmaTwo     = 0.346;
//        double rhoOne       = -0.848;
//        double rhoTwo       = -0.402;
        double sigmaOne     = 1.039;
        double sigmaTwo     = 0.667;
        double rhoOne       = -0.775;
        double rhoTwo       = -0.382;
        double lambdaZero   = 0;
        double lambdaOne    = 0;
        double lambdaTwo    = 0;
 
        double k       		= 0;
        double delta   		= 0;
        
        

		
        //set known constants
        
        double volatilityOne   = 0.00963;
        double volatilityTwo   = 0.01352;      
        
		
		
		int tenor = 10;
		double shift;
		if(tenor == 2){shift = 0.0265;}else if(tenor == 5){shift = 0.016;}else if(tenor == 10){shift = 0.015;}else{shift = 0;};
//		if(tenor == 2){shift = 2.65;}else if(tenor == 5){shift = 1.6;}else if(tenor == 10){shift = 1.5;}else{shift = 0;};
		double[][] volatilities = sheetdata.getSmileData(tenor);
		
		
        double[] initialParameters = new double[]{			
            	Math.log(alphaOne),
            	Math.log(alphaTwo),
    			Math.log(betaOne),
    			Math.log(betaTwo),
    			Math.log(sigmaOne),
    			Math.log(sigmaTwo),
    			rhoOne,
    			rhoTwo,
//    			volatilityOne,
//    			volatilityTwo
            };
		
		int[] maturities = {1, 2, 5, 10, 15, 20};
//		int[] maturities = {5};
		double[] strikes = {-0.01,-0.005, -0.0025, 0, 0.0025, 0.005, 0.01, 0.015, 0.02};
		AbstractProductFourierTransform[][] europeanMatrix = new AbstractProductFourierTransform[maturities.length][strikes.length];
		
		double [] targetValues = new double[maturities.length*strikes.length];
//		double [] bachelierValues = new double[maturities.length*strikes.length];
		
		double [] weights = new double[maturities.length*strikes.length];
		for(int w = 0; w < weights.length; w++){weights[w] = 1;};

		for(int w = 0; w < weights.length; w++){weights[w] = 1;};
		System.out.println("Maturty");
		for(int i = 0; i < maturities.length; i++){
			System.out.print(maturities[i] + "\t");
			for(int j = 0; j < strikes.length; j++){
				targetValues[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.blackScholesGeneralizedOptionValue(sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,volatilities[i][j],maturities[i], strikes[j]+shift,sheetdata.getSwapAnnuity(maturities[i], tenor));
//				bachelierValues[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.bachelierOptionValue(forwardcurve[index][2]+1.5,volatilities[i][j],maturities[i], strikes[j]+1.5,sheetdata.getSwapAnnuity(maturities[i], tenor));
				europeanMatrix[i][j] = new EuropeanOption(maturities[i], strikes[j]+shift, sheetdata.getSwapAnnuity(maturities[i], tenor));
				System.out.print(df.format(targetValues[i*strikes.length + j]*100) + "    \t");
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
																lambdaZero,
																lambdaOne,
																lambdaTwo,
																k,
																delta,
//																parameters[13],
//																parameters[14],
																volatilityOne,//[i],
																volatilityTwo,//[i],
																sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
																0)
															)
							;
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
		System.out.println( "\talphaOne = \t\t" + Math.exp(bestParameters[0]) + 
							"\n\talphaTwo = \t\t" + Math.exp(bestParameters[1]) + 
							"\n\tbetaOne = \t\t" + Math.exp(bestParameters[2]) + 
							"\n\tbetaTwo = \t\t" + Math.exp(bestParameters[3]) + 
							"\n\tsigmaOne = \t\t" + Math.exp(bestParameters[4]) + 
							"\n\tsigmaTwo = \t\t" + Math.exp(bestParameters[5]) + 
							"\n\trhoOne = \t\t" + bestParameters[6] + 
							"\n\trhoTwo = \t\t" + bestParameters[7] 
//							"\n\tvolatilityOne = \t" + bestParameters[13] + 
//							"\n\tvolatilityTwo = \t" + bestParameters[14]
									);
		
		System.out.println("\n Factor 1:" + df.format(Math.exp(bestParameters[0])) + "&\t" +
				df.format(Math.exp(bestParameters[2])) + "&\t" +
				df.format(Math.exp(bestParameters[4])) + "&\t" +
				df.format(bestParameters[6]) + "&\t"
				);
		System.out.println("\n Factor 2:" + df.format(Math.exp(bestParameters[1])) + "&\t" +
				df.format(Math.exp(bestParameters[3])) + "&\t" +
				df.format(Math.exp(bestParameters[5])) + "&\t" +
				df.format(bestParameters[7]) + "& &\t" + "& & &\n" 
				);
		
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
													lambdaZero,
													lambdaOne,
													lambdaTwo,
													k,
													delta,
// 													bestParameters[13],
// 													bestParameters[14],
 													volatilityOne,//[i],
 													volatilityTwo,//[i],
 													sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
 													0)
 												)
 				;
				mse +=(values[i*strikes.length + j]-targetValues[i*strikes.length + j])*(values[i*strikes.length + j]-targetValues[i*strikes.length + j]);
				System.out.print("&" + df.format(values[i*strikes.length + j]*100) + "\t");
			}
			System.out.print("\\\\ \n");
 		}
		System.out.println("RMSE =" + Math.sqrt((mse/(maturities.length*strikes.length))));
		
	}
	
	
	
	
}


