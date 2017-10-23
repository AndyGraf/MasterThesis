package com.graf.fouriertransform;

import java.io.IOException;
import net.finmath.functions.AnalyticFormulas;

//copy the printed values into the corresponding column SABR or Maturity-wise SABR of "vola smile sample plots.ods" 

public class VolatilitySmilePlotSABR {

	public static void main(String[] args) throws IOException {
		
		SABRdata sheetdata = new SABRdata();
		
		// parameters to be replaced by the calibrated parameters
		
        double alpha	    = 0.3222856232938358;
        double beta 		= 0.7492778681596671;
        double rho      	= -0.9999973982223452;
        double v   			= 0.03846594129745918;

        // maturity and tenor
        
        double maturities 	= 1;
        double tenor 		= 10;
        
        
        
        double[] strikes = new double[61];
        for(int j = 0; j < strikes.length; j++){
        	strikes[j] = 0.05*j-1;
        }
        
        
        double shift;
		if(tenor == 2){shift = 2.65;}else if(tenor == 5){shift = 1.6;}else if(tenor == 10){shift = 1.5;}else{shift = 0;};

    
				for(int j = 0; j < strikes.length; j++){
						try {
							System.out.println(AnalyticFormulas.sabrBerestyckiNormalVolatilityApproximation(
																alpha,
																beta,
																rho,
																v,
																0, //displacement
																sheetdata.getForwardSwapRate(maturities, tenor)+shift,
																strikes[j]+shift,
																maturities)
									);
						} catch (IOException e) {
							e.printStackTrace();
						}
				}
				
        
	}

}
