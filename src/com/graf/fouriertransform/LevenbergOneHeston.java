package com.graf.fouriertransform;

import java.io.IOException;
import java.text.DecimalFormat;

import net.finmath.fouriermethod.models.ProcessCharacteristicFunctionInterface;
import net.finmath.fouriermethod.products.AbstractProductFourierTransform;
import net.finmath.optimizer.LevenbergMarquardt;
import net.finmath.optimizer.SolverException;

public class LevenbergOneHeston {

	public static void main(String[] args) throws IOException {
		DecimalFormat df = new DecimalFormat("#.###");
        
		SABRdata sheetdata = new SABRdata();
		
		//initial parameters  
        
        double alphaOne     = 0.1;
        double betaOne      = 1.49;
        double sigmaOne     = 0.742;
        double rhoOne       = -0.571;
        double lambdaZero   = 0;
        double lambdaOne 	= 0;
 
        double k       		= 0;
        double delta   		= 0;
        
        
		
        //set known constants
        
        
        
        double volatilityOne   = 0.00963;    
        
		
		
		int tenor = 10;
		double shift;
		if(tenor == 2){shift = 0.0265;}else if(tenor == 5){shift = 0.016;}else if(tenor == 10){shift = 0.015;}else{shift = 0;};
		double[][] volatilities = sheetdata.getSmileData(tenor);
		
//		double[] volatilityOne = new double[volatilities[3].length];
//		for(int i=0; i<volatilities[3].length;i++){
//			volatilityOne[i]   = volatilities[3][i]*volatilities[3][i];
//		}
		
        double[] initialParameters = new double[]{			
            	Math.log(alphaOne),
    			Math.log(betaOne),
    			Math.log(sigmaOne),
    			rhoOne,
//    			volatilityOne
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
			System.out.print(maturities[i] + "  \t");
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
																parameters[3],
																lambdaZero,
																lambdaOne,
																k,
																delta,
																volatilityOne,
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
	  	optimizer.setMaxIteration(10);
	  	
	  	optimizer.setTargetValues(targetValues);
	  
	  	try {
			optimizer.run();
		} catch (SolverException e) {
			e.printStackTrace();
		}
	  
	  	double[] bestParameters = optimizer.getBestFitParameters();
	  	

	  
	  	
	  	System.out.println("The solver for problem 1 required " + optimizer.getIterations() + " iterations. The best fit parameters are:");
		System.out.println( "\talphaOne = \t\t" + Math.exp(bestParameters[0]) + 
							"\n\tbetaOne = \t\t" + Math.exp(bestParameters[1]) + 
							"\n\tsigmaOne = \t\t" + Math.exp(bestParameters[2]) + 
							"\n\trhoOne = \t\t" + bestParameters[3]
									);
		
		System.out.println("\n Factor 1:" + df.format(Math.exp(bestParameters[0])) + "&\t" +
				df.format(Math.exp(bestParameters[1])) + "&\t" +
				df.format(Math.exp(bestParameters[2])) + "&\t" +
				df.format(bestParameters[3]) + "&&&&&&&\t"				
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
 													bestParameters[3],
 													lambdaZero,
 													lambdaOne,
 													k,
 													delta,
 													volatilityOne,//[i],
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


