package com.graf.fouriertransform;

import net.finmath.fouriermethod.models.BlackScholesModel;
import net.finmath.fouriermethod.models.ProcessCharacteristicFunctionInterface;
import net.finmath.fouriermethod.products.AbstractProductFourierTransform;
import net.finmath.fouriermethod.products.EuropeanOption;


public class FourierTestRun {

    public static void main(String[] args) {
    	

//        double[] alpha     = {0.1, 0};
//        double[] beta    = {0.5, 0.00001};
//        double[] sigma    = {0.3, 0.00001};
//        double[] rho    = {-0.8, 0};
//        double[] lambda    = {0.5, 0, 0};
    	
        double theta        = 0.2;
        double kappa        = 0.5;
        double alpha        = theta*kappa;
        double beta         = kappa;
        double sigma        = 0.3;
        double rho          = -0.8;
        double lambda       = 0.5;
 
        double delta        = Math.sqrt(0.3);
        double k           = -0.5;
        double volatility   = 0.2;
        double initialPrice = 100;
        double riskFreeRate = 0.05;
        
        double maturity 	= 1;
        double strike		= 100;
        
//        double theta		= 0;
//    	double kappa		= 0;
//    	double alpha		= 0;
//    	double beta   	 	= 0.00000001;
//  	   	double sigma		= 0.000000001;
//  	   	double rho			= 0;
//  	   	double lambda		= 0;
//
//        double delta		= 0;
//        double mu			= 0;
//        double k			= 0;

        ProcessCharacteristicFunctionInterface bs = new BlackScholesModel(initialPrice, 0, Math.sqrt(volatility));
        ProcessCharacteristicFunctionInterface bates = new TwoFactorBatesModelCF(alpha, beta, sigma, rho, lambda, k, delta, volatility, initialPrice, riskFreeRate);
        AbstractProductFourierTransform european = new EuropeanOption(maturity, strike);

        System.out.println(Math.exp(-riskFreeRate*maturity)*european.getValue(bates));


    }

}