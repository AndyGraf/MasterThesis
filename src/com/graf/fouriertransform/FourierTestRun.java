package com.graf.fouriertransform;

import net.finmath.fouriermethod.models.BlackScholesModel;
import net.finmath.fouriermethod.models.ProcessCharacteristicFunctionInterface;
import net.finmath.fouriermethod.products.AbstractProductFourierTransform;


public class FourierTestRun {

    public static void main(String[] args) {
    	
    	
        double theta        = 14.229/248.49;
        double kappa        = 2.49;
        double alpha        = theta*kappa;
        double beta         = kappa;
        double sigma        = 46.63;
        double rho          = -0.454;
        double lambda       = 0;
 
        double delta        = 0;
        double k            = 0;
        double volatility   = 0.00963; 
        double initialPrice = 100;
        double riskFreeRate = 0.05;
        
        double maturity 	= 20;
        double strike		= 100;
        

//        ProcessCharacteristicFunctionInterface bs = new BlackScholesModel(initialPrice, 0, Math.sqrt(volatility));
        ProcessCharacteristicFunctionInterface bates = new TwoFactorBatesModelCF(alpha, beta, sigma, rho, lambda, 0, k, delta, volatility, initialPrice, riskFreeRate);
        ProcessCharacteristicFunctionInterface bates2 = new TwoFactorBatesModelCF(alpha, beta, sigma, rho, lambda, 0, k, delta, volatility, initialPrice, riskFreeRate);
        AbstractProductFourierTransform european = new EuropeanOption(maturity, strike, Math.exp(-riskFreeRate*maturity));

        System.out.println(european.getValue(bates));
        System.out.println(european.getValue(bates2));

    }

}