package com.graf.fouriertransform;

import java.io.IOException;
import java.text.DecimalFormat;

import net.finmath.exception.CalculationException;
import net.finmath.fouriermethod.products.AbstractProductFourierTransform;
import net.finmath.fouriermethod.products.EuropeanOption;

// copy the printed values into the corresponding column Two-factor Bates of "vola smile sample plots.ods" 

public class VolatilitySmilePlot {

	public static void main(String[] args) throws IOException {
		DecimalFormat error = new DecimalFormat("#.#########");
		
		SABRdata sheetdata = new SABRdata();
		
        // parameters to be replaced by the calibrated parameters
		
        double alphaOne = 		0.01057829966771772;
        double alphaTwo = 		0.009128787405722911;
        double betaOne = 		0.0666006022640223;
        double betaTwo = 		1.6235837879480637;
        double sigmaOne = 		0.1523803477571126;
        double sigmaTwo = 		0.23316502414541135;
        double rhoOne = 		-0.9460013373386712;
        double rhoTwo = 		0.2771401117338455;
        double lambdaZero = 	0;
        double lambdaOne = 		0;
        double lambdaTwo = 		0;
        double k = 				0;
        double delta = 			0;
        double volatilityOne = 	0.007803940918576663;
        double volatilityTwo = 	0.16370502442716942;
        

        // maturity and tenor
        
        double maturities 	= 1;
        double tenor 		= 2;
        
        
        
        double[] strikes = new double[61];
        for(int j = 0; j < strikes.length; j++){
        	strikes[j] = 0.05*j-1;
        }
        
        AbstractProductFourierTransform[] europeanMatrix = new AbstractProductFourierTransform[strikes.length];
        
        double shift;
		if(tenor == 2){shift = 2.65;}else if(tenor == 5){shift = 1.6;}else if(tenor == 10){shift = 1.5;}else{shift = 0;};


		for(int j = 0; j < strikes.length; j++){
			europeanMatrix[j] = new EuropeanOption(maturities, strikes[j]+shift);
		}
		
        
        double implVol;
        double modelPrice;
        
        for(int j = 0; j < strikes.length; j++){
			try {
				modelPrice = 						europeanMatrix[j].getValue(
						new TwoFactorBatesModelCF(
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
						    	delta,
						    	volatilityOne,
						    	volatilityTwo,
								sheetdata.getForwardSwapRate(maturities, tenor)+shift,
								0,
								sheetdata.getSwapAnnuity(maturities, tenor))
							);
				implVol = net.finmath.functions.AnalyticFormulas.bachelierOptionImpliedVolatility(
						sheetdata.getForwardSwapRate(maturities, tenor)+shift,
						maturities,
						strikes[j]+shift,
						sheetdata.getSwapAnnuity(maturities, tenor),
						modelPrice
						);
				System.out.println(implVol 
						// optional check to see if the implied volatility is solved correctly
						
//						+ "\t\t\t" + error.format((modelPrice - net.finmath.functions.AnalyticFormulas.bachelierOptionValue(sheetdata.getForwardSwapRate(maturities, tenor)+shift,implVol, maturities,
//						strikes[j]+shift, sheetdata.getSwapAnnuity(maturities, tenor))))
						)
				;
			} catch (IOException e) {
				e.printStackTrace();
			} catch (CalculationException e) {
				e.printStackTrace();
			}
		}
        
        }

}
