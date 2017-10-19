package com.graf.fouriertransform;

import net.finmath.exception.CalculationException;
import net.finmath.fouriermethod.models.BlackScholesModel;
import net.finmath.fouriermethod.models.ProcessCharacteristicFunctionInterface;
import net.finmath.fouriermethod.products.AbstractProductFourierTransform;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.assetderivativevaluation.AssetModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.assetderivativevaluation.HestonModel;
import net.finmath.montecarlo.assetderivativevaluation.HestonModelTest;
import net.finmath.montecarlo.assetderivativevaluation.MonteCarloAssetModel;
import net.finmath.montecarlo.assetderivativevaluation.HestonModel.Scheme;
import net.finmath.montecarlo.assetderivativevaluation.products.EuropeanOption;
import net.finmath.montecarlo.model.AbstractModel;
import net.finmath.montecarlo.process.AbstractProcess;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;


public class FourierTestRun {
	


    public static void main(String[] args) {
    	
    	
    	double	initialValue   = 1;
    	double	riskFreeRate   = 0.05;
    	double	volatility     = 0.30;

    	double theta = volatility*volatility;
    	double kappa = 1.0;
    	double xi = 0.00;
    	double rho = 0.0;

    	Scheme scheme = Scheme.FULL_TRUNCATION;

    	// Process discretization properties
    	int		numberOfPaths		= 1000;
    	int		numberOfTimeSteps	= 100;
    	double	deltaT				= 0.02;

    	int		seed				= 3815;

    	// Product properties
    	double	optionMaturity = 2.0;
    	double	optionStrike = 1.10;
    	
    	// parameter values for the fourier implementation
    	
        double alpha        = theta*kappa;
        double beta         = kappa;
        double sigma        = xi;
        
        double lambda       = 0;
        double delta        = 0;
        double k            = 0;
        
		TimeDiscretizationInterface timeDiscretization = new TimeDiscretization(0.0 /* initial */, numberOfTimeSteps, deltaT);

		BrownianMotionInterface brownianMotion = new BrownianMotion(timeDiscretization, 2 /* numberOfFactors */, numberOfPaths, seed);

        
        

//        ProcessCharacteristicFunctionInterface bs = new BlackScholesModel(initialPrice, 0, Math.sqrt(volatility));
        ProcessCharacteristicFunctionInterface bates = new TwoFactorBatesModelCF(alpha, beta, sigma, rho, lambda, 0, k, delta, volatility*volatility, initialValue, riskFreeRate, Math.exp(-riskFreeRate*optionMaturity));
        net.finmath.fouriermethod.products.EuropeanOption european = new net.finmath.fouriermethod.products.EuropeanOption(optionMaturity, optionStrike);

        double fourierValue = 0;
		try {
			fourierValue = european.getValue(bates);
		} catch (CalculationException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
        System.out.println("Analytic value using the Fourier method: \n" + fourierValue);

        
        
		AssetModelMonteCarloSimulationInterface monteCarloHestonModel;
		{
			// Create a model
			AbstractModel model = new HestonModel(initialValue, riskFreeRate, volatility, theta, kappa, xi, rho, scheme);

			// Create a corresponding MC process
			AbstractProcess process = new ProcessEulerScheme(brownianMotion);

			// Using the process (Euler scheme), create an MC simulation of a Black-Scholes model
			monteCarloHestonModel = new MonteCarloAssetModel(model, process);
		}
       


			EuropeanOption europeanOption = new EuropeanOption(optionMaturity, optionStrike);
			double value = 0;
			try {
				value = europeanOption.getValue(monteCarloHestonModel);
			} catch (CalculationException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.out.println("\nValue using the monte carlo method: \n" + value);
		
			System.out.println("\nError: \n" + Math.abs(fourierValue - value));
		

    }

}