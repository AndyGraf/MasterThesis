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
		
        double alphaOne = 		1.058168026131326E-4;
        double alphaTwo = 		0.006867109677130183;
        double betaOne = 		0.5021933073806344;
        double betaTwo = 		0.0011499094181474838;
        double sigmaOne = 		0.15514725394931037;
        double sigmaTwo = 		0.09489363562326633;
        double rhoOne = 		-0.3729490471652218;
        double rhoTwo = 		-0.9752770154360861;
        double lambdaZero = 	0.0046139184162780214;
        double lambdaOne = 		24.894455318935723;
        double lambdaTwo = 		0.17221778191218987;
        double k = 				0.23251509950953383;
        double delta = 			0.00237890974999215;
        double volatilityOne = 	0.04951194402156467;
        double volatilityTwo = 	4.10626281336595E-4;

        
        // maturity and tenor
        
        double maturities 	= 1;
        double tenor 		= 10;
        
        
        
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
