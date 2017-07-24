package com.graf.fouriertransform;

import java.io.IOException;

import net.finmath.fouriermethod.models.ProcessCharacteristicFunctionInterface;
import net.finmath.fouriermethod.products.AbstractProductFourierTransform;
import net.finmath.optimizer.LevenbergMarquardt;
import net.finmath.optimizer.SolverException;

public class LevenbergOneBates {

	public static void main(String[] args) throws IOException {
        
		SABRdata sheetdata = new SABRdata();
		
		//initial parameters  
//        double alphaOne     = 0.01;
//        double betaOne      = 0.91;
//        double sigmaOne     = 0.582;
//        double rhoOne       = -0.848;
//        double lambdaZero   = 0.1143;
// 
//        double k       		= -0.057;
//        double delta   		= 0.102;
        
        double alphaOne     = 0.049;
        double betaOne      = 2.45;
        double sigmaOne     = 0.378;
        double rhoOne       = -0.545;
        double lambdaZero   = 0.1;
 
        double k       		= -0.085;
        double delta   		= 0.113;
        
        
		
        //set known constants
        
        
        
        double volatilityOne   = 0.00963;    
        
		
		
		int tenor = 10;
		double shift;
		if(tenor == 2){shift = 2.65;}else if(tenor == 5){shift = 1.6;}else if(tenor == 10){shift = 1.5;}else{shift = 0;};
		double[][] volatilities = sheetdata.getSmileData(tenor);
		
//		double[] volatilityOne = new double[volatilities[3].length];
//		double[] volatilityTwo = new double[volatilities[3].length];
//		for(int i=0; i<volatilities[3].length;i++){
//			volatilityOne[i]   = volatilities[3][i]*volatilities[3][i];
//			volatilityTwo[i]   = volatilities[3][i]*volatilities[3][i];
//		}
		
        double[] initialParameters = new double[]{			
            	Math.log(alphaOne),
    			Math.log(betaOne),
    			Math.log(sigmaOne),
    			rhoOne,
    			Math.log(lambdaZero),
    			k,
    			delta
//    			volatilityOne
            };
		
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
																parameters[3],
																Math.exp(parameters[4]),
																parameters[5],
																parameters[6],
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
							"\n\tbetaOne = \t\t" + Math.exp(bestParameters[1]) + 
							"\n\tsigmaOne = \t\t" + Math.exp(bestParameters[2]) + 
							"\n\trhoOne = \t\t" + bestParameters[3] + 
							"\n\tlambdaZero = \t\t" + Math.exp(bestParameters[4]) + 
							"\n\tk = \t\t\t" + bestParameters[5] + 
							"\n\tdelta = \t\t" + bestParameters[6] 
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
 													Math.exp(bestParameters[4]),
 													bestParameters[5],
 													bestParameters[6],
 													volatilityOne,//[i],
 													sheetdata.getForwardSwapRate(maturities[i], tenor)+shift,
 													0)
 												)
 				;
				mse +=(values[i*strikes.length + j]-targetValues[i*strikes.length + j])*(values[i*strikes.length + j]-targetValues[i*strikes.length + j]);
				System.out.print(values[i*strikes.length + j] + "\t");
			}
			System.out.println();
 		}
		System.out.println("RMSE =" + Math.sqrt((mse/(maturities.length*strikes.length))));
		
		System.out.println("\n Same parameters but different tenor: targets\n");
		int tenor2 = 5;
		double shift2 = 0;
		if(tenor2 == 2){shift2 = 2.65;}else if(tenor2 == 5){shift2 = 1.6;}else if(tenor2 == 10){shift2 = 1.5;}else{shift2 = 0;};
		volatilities = sheetdata.getSmileData(tenor2);
//		for(int i=0; i<volatilities[3].length;i++){
//			volatilityOne[i]   = volatilities[3][i]*volatilities[3][i];
//			volatilityTwo[i]   = volatilities[3][i]*volatilities[3][i];
//		}
//		for(int i = 0; i < maturities.length; i++){
//			for(int j = 0; j < strikes.length; j++){
//				targetValues[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.blackScholesGeneralizedOptionValue(sheetdata.getForwardSwapRate(maturities[i], tenor2)+shift2,volatilities[i][j],maturities[i], strikes[j]+shift2,sheetdata.getSwapAnnuity(maturities[i], tenor2));
////				bachelierValues[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.bachelierOptionValue(forwardcurve[index][2]+1.5,volatilities[i][j],maturities[i], strikes[j]+1.5,sheetdata.getSwapAnnuity(maturities[i], tenor));
//				europeanMatrix[i][j] = new EuropeanOption(maturities[i], strikes[j]+shift2, sheetdata.getSwapAnnuity(maturities[i], tenor2));
//				System.out.print(targetValues[i*strikes.length + j] + "\t");
//			}
//			System.out.println();
//		}
//		System.out.println("\n result:\n");
//		for(int i = 0; i < maturities.length; i++){
//			for(int j = 0; j < strikes.length; j++){
//				values[i*strikes.length + j] = europeanMatrix[i][j].getValue(
// 													new TwoFactorBatesModelCF(
// 													Math.exp(bestParameters[0]),
// 													Math.exp(bestParameters[1]),
// 													Math.exp(bestParameters[2]),
// 													Math.exp(bestParameters[3]),
// 													Math.exp(bestParameters[4]),
// 													Math.exp(bestParameters[5]),
// 													bestParameters[6],
// 													bestParameters[7],
// 													Math.exp(bestParameters[8]),
// 													Math.exp(bestParameters[9]),
// 													Math.exp(bestParameters[10]),
// 													bestParameters[11],
// 													bestParameters[12],
//// 													bestParameters[13],
//// 													bestParameters[14],
// 													volatilityOne,//[i],
// 													volatilityTwo,//[i],
// 													sheetdata.getForwardSwapRate(maturities[i], tenor2)+shift2,
// 													0)
// 												)
// 				;
//				mse +=(values[i*strikes.length + j]-targetValues[i*strikes.length + j])*(values[i*strikes.length + j]-targetValues[i*strikes.length + j]);
//				System.out.print(values[i*strikes.length + j] + "\t");
//			}
//			System.out.println();
// 		}
//		System.out.println("RMSE =" + Math.sqrt((mse/(maturities.length*strikes.length))));
//	
	}
	
	
	
	
}


