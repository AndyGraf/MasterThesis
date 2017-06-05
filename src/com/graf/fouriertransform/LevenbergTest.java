package com.graf.fouriertransform;

import net.finmath.fouriermethod.models.ProcessCharacteristicFunctionInterface;
import net.finmath.fouriermethod.products.AbstractProductFourierTransform;
import net.finmath.optimizer.LevenbergMarquardt;
import net.finmath.optimizer.SolverException;

public class LevenbergTest {

	public static void main(String[] args) {
        
		//initial parameters  
        double alphaOne     = 0.1;
        double alphaTwo		= 0;
        double betaOne      = 0.5;
        double betaTwo      = 0.2;
        double sigmaOne     = 0.3;
        double sigmaTwo     = 0.4;
        double rhoOne       = -0.8;
        double rhoTwo       = -0.6;
        double lambdaZero   = 0.5;
        double lambdaOne    = 1;
        double lambdaTwo    = 2;
 
        double k       		= -0.5;
        double delta   		= Math.sqrt(0.3);
        
//        abc = BlackScholesFormula();
//        price = BlackScholesFormula(forwardSwaprate + 1.5%, strike + 1.5%, impliedvol, annuity, maturity);
        //set known constants
        
        double volatilityOne   = 0.2;
        double volatilityTwo   = 0.28;        
        
        double initialPrice = 100;
        double riskFreeRate = 0.05;
        double maturity 	= 1;
        double[] strikes		= new double[] {90, 95};
	  	double[] targetValues = new double[] {50, 45};
        
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
        
        AbstractProductFourierTransform[] europeanVector = new AbstractProductFourierTransform[10];

        for(int i = 0; i < strikes.length; i++ ){
        	europeanVector[i] = new EuropeanOption(maturity, strikes[i], riskFreeRate);
        }
        
        
	 	LevenbergMarquardt optimizer = new LevenbergMarquardt() {

		 	public void setValues(double[] parameters, double[] values) {
		 		for(int i = 0; i < values.length; i++ ){
		 			values[i] = europeanVector[i].getValue(
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
		 									initialPrice,
		 									riskFreeRate)
		 							)
		 						;
		 		}
			}
		};

	  	optimizer.setInitialParameters(initialParameters);
	  	optimizer.setWeights(new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
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
		
		ProcessCharacteristicFunctionInterface bates = new TwoFactorBatesModelCF(
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
					initialPrice,
					riskFreeRate);
        for(int i = 0; i < targetValues.length; i++ ){
        	System.out.println("\nTargetvalue: " + targetValues[i] + "\t result: " + europeanVector[i].getValue(bates));
        }
	}
}


