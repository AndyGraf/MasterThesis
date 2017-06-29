package com.graf.fouriertransform;

import java.io.IOException;

import com.graf.sabrcall.CallOption;

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
        double lambdaOne    = 1.56;
        double lambdaTwo    = 0.28;
 
        double k       		= -0.057;
        double delta   		= 0.102;
        
        //set known constants
        
        double volatilityOne   = 0.00963;
        double volatilityTwo   = 0.01352;      
        
//        double initialPrice = 100;
//        double riskFreeRate = 0.05;
//        double maturity 	= 1;
//        double[] strikes		= new double[] {90, 95};
//	  	double[] targetValues = new double[] {50, 45};
//        
        double[] initialParameters = new double[]{			
        	alphaOne,
			alphaTwo,
			betaOne,
			betaTwo,
			sigmaOne,
			sigmaTwo,
			rhoOne,
			rhoTwo,
			lambdaZero,
			lambdaOne,
			lambdaTwo,
			k,
			delta};
        
        
		SABRdata sheetdata = new SABRdata();
		
		int tenor = 5;
		double[][] forwardcurve = sheetdata.getForwardCurveData();
		double[][] discount = sheetdata.getDiscountCurveData();
		double[][] volatilities = sheetdata.getSmileData(tenor);
		int[] maturities = {1, 2, 5, 10, 15, 20};
		double[] strikes = {-1, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5, 2};
		AbstractProductFourierTransform[][] europeanMatrix = new AbstractProductFourierTransform[maturities.length][strikes.length];
		
		double [] targetValues = new double[maturities.length*strikes.length];
		double [] bachelierValues = new double[maturities.length*strikes.length];
		
		double [] weights = new double[maturities.length*strikes.length];
		for(int w = 0; w < weights.length; w++){weights[w] = 1;};
		int index;
		

		for(int i = 0; i < maturities.length; i++){
			for(int j = 0; j < strikes.length; j++){
				index = maturities[i]+tenor;
				targetValues[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.blackScholesGeneralizedOptionValue(forwardcurve[index][2]+1.5,volatilities[i][j],maturities[i], strikes[j]+1.5,sheetdata.getSwapAnnuity(maturities[i], tenor));
				bachelierValues[i*strikes.length + j] = net.finmath.functions.AnalyticFormulas.bachelierOptionValue(forwardcurve[index][2]+1.5,volatilities[i][j],maturities[i], strikes[j]+1.5,sheetdata.getSwapAnnuity(maturities[i], tenor));
				europeanMatrix[i][j] = new EuropeanOption(maturities[i], strikes[j]+1.5, -Math.log(sheetdata.getSwapAnnuity(maturities[i], tenor))/maturities[i]);
				System.out.print(targetValues[i*strikes.length + j] + "\t");
			}
			System.out.println();
		}
        
        
        
        
	 	LevenbergMarquardt optimizer = new LevenbergMarquardt() {

		 	public void setValues(double[] parameters, double[] values) {
		 		for(int i = 0; i < maturities.length; i++){
					for(int j = 0; j < strikes.length; j++){
						int index = maturities[i]+tenor;
						try {
							values[i*strikes.length + j] = europeanMatrix[i][j].getValue(
																new TwoFactorBatesModelCF(
																parameters[0],
																parameters[1],
																parameters[2],
																parameters[3],
																parameters[4],
																parameters[5],
																parameters[6],
																parameters[7],
																parameters[8],
																parameters[9],
																parameters[10],
																parameters[11],
																parameters[12],
																volatilityOne,
																volatilityTwo,
																forwardcurve[index][2]+1.5,
																-Math.log(sheetdata.getSwapAnnuity(maturities[i], tenor))/maturities[i])
															)
							;
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						//System.out.println(values[i*strikes.length + j] + "\t");
					}
		 		}
			}
		};

	  	optimizer.setInitialParameters(initialParameters);
	  	optimizer.setWeights(weights);
	  	optimizer.setMaxIteration(100);
	  	
	  	optimizer.setTargetValues(targetValues);
	  
	  	try {
			optimizer.run();
		} catch (SolverException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	  
	  	double[] bestParameters = optimizer.getBestFitParameters();
	  	
	  	System.out.println("The solver for problem 1 required " + optimizer.getIterations() + " iterations. The best fit parameters are:");
		System.out.println( "\talphaOne: \t\t" + bestParameters[0] + 
							"\n\talphaTwo: \t\t" + bestParameters[1] + 
							"\n\tbetaOne: \t\t" + bestParameters[2] + 
							"\n\tbetaTwo: \t\t" + bestParameters[3] + 
							"\n\tsigmaOne: \t\t" + bestParameters[4] + 
							"\n\tsigmaTwo: \t\t" + bestParameters[5] + 
							"\n\trhoOne: \t\t" + bestParameters[6] + 
							"\n\trhoTwo: \t\t" + bestParameters[7] + 
							"\n\tlambdaZero: \t\t" + bestParameters[8] + 
							"\n\tlambdaOne: \t\t" + bestParameters[9] + 
							"\n\tlambdaTwo: \t\t" + bestParameters[10] + 
							"\n\tk: \t\t\t" + bestParameters[11] + 
							"\n\tdelta: \t\t\t" + bestParameters[12]);

//		System.out.println( "\tdouble alphaOne \t\t=" + bestParameters[0] +";"+ 
//							"\n\tdouble alphaTwo \t\t=" + bestParameters[1] +";"+ 
//							"\n\tdouble betaOne \t\t=" + bestParameters[2] +";"+ 
//							"\n\tdouble betaTwo \t\t=" + bestParameters[3] +";"+ 
//							"\n\tdouble sigmaOne \t\t=" + bestParameters[4] +";"+ 
//							"\n\tdouble sigmaTwo \t\t=" + bestParameters[5] +";"+ 
//							"\n\tdouble rhoOne \t\t=" + bestParameters[6] +";"+ 
//							"\n\tdouble rhoTwo \t\t=" + bestParameters[7] +";"+ 
//							"\n\tdouble lambdaZero \t\t=" + bestParameters[8] +";"+ 
//							"\n\tdouble lambdaOne \t\t=" + bestParameters[9] +";"+ 
//							"\n\tdouble lambdaTwo \t\t=" + bestParameters[10] +";"+ 
//							"\n\tdouble k \t\t\t=" + bestParameters[11] +";"+ 
//							"\n\tdouble delta \t\t\t=" + bestParameters[12] +";");
		
//		ProcessCharacteristicFunctionInterface bates = new TwoFactorBatesModelCF(
//					bestParameters[0],
//					bestParameters[1],
//					bestParameters[2],
//					bestParameters[3],
//					bestParameters[4],
//					bestParameters[5],
//					bestParameters[6],
//					bestParameters[7],
//					bestParameters[8],
//					bestParameters[9],
//					bestParameters[10],
//					bestParameters[11],
//					bestParameters[12],
//					volatilityOne,
//					volatilityTwo,
//					initialPrice,
//					riskFreeRate);
//        for(int i = 0; i < targetValues.length; i++ ){
//        	System.out.println("\nTargetvalue: " + targetValues[i] + "\t result: " + europeanMatrix[i].getValue(bates));
//        }
		double [] values = new double[maturities.length*strikes.length];
		double mse = 0;
		for(int i = 0; i < maturities.length; i++){
			for(int j = 0; j < strikes.length; j++){
				int index2 = maturities[i]+tenor;
				values[i*strikes.length + j] = europeanMatrix[i][j].getValue(
 													new TwoFactorBatesModelCF(
 													bestParameters[0],
 													bestParameters[1],
 													bestParameters[2],
 													bestParameters[3],
 													bestParameters[4],
 													bestParameters[5],
 													bestParameters[6],
 													bestParameters[7],
 													bestParameters[8],
 													bestParameters[9],
 													bestParameters[10],
 													bestParameters[11],
 													bestParameters[12],
 													volatilityOne,
 													volatilityTwo,
 													forwardcurve[index2][2]+1.5,
 													-Math.log(sheetdata.getSwapAnnuity(maturities[i], tenor))/maturities[i])
 												)
 				;
				mse +=(values[i*strikes.length + j]-targetValues[i*strikes.length + j])*(values[i*strikes.length + j]-targetValues[i*strikes.length + j]);
				System.out.print(values[i*strikes.length + j] + "\t");
			}
			System.out.println();
 		}
		System.out.println(Math.sqrt(mse));
	}
}


