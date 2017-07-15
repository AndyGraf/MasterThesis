package com.graf.fouriertransform;

import java.io.File;
import java.io.IOException;

import org.jopendocument.dom.spreadsheet.Sheet;
import org.jopendocument.dom.spreadsheet.SpreadSheet;

import net.finmath.marketdata.model.curves.Curve;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.DiscountCurveInterface;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterface;
import net.finmath.marketdata.products.Swap;
import net.finmath.marketdata.products.SwapAnnuity;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretization.ShortPeriodLocation;
import net.finmath.time.TimeDiscretizationInterface;

public class SABRdata {
	
	private final Sheet volatilities;
	private final Sheet smile;
	DiscountCurveInterface discountCurve = null;
	ForwardCurveInterface forwardCurve = null;
	private double periodLength = 0.5;
	private final int rows = 107;
	private final int columns = 5;
	private final int rowsForward = 97;
	
	public SABRdata() throws IOException{
		
		File file = new File("src/Data.ods");
		volatilities = SpreadSheet.createFromFile(file).getSheet(0);
		smile = SpreadSheet.createFromFile(file).getSheet(1);
		


		double[][] discountCurveData = new double[rows][columns];
		for(int i=0; i< rows;i++){
			for(int j=0;j<columns; j++){
				discountCurveData[i][j] = Double.parseDouble(volatilities.getValueAt(j+8, i+27).toString());
			}
		}

		double[] times = new double[discountCurveData.length];
		double[] givenDiscountFactors = new double[discountCurveData.length];
		for(int i=0; i<times.length; i++) {
			times[i] = discountCurveData[i][3]*360/365;
			givenDiscountFactors[i] = discountCurveData[i][2];
		}
		discountCurve = DiscountCurve.createDiscountCurveFromDiscountFactors("discountCurve", times, givenDiscountFactors, Curve.InterpolationMethod.LINEAR, Curve.ExtrapolationMethod.LINEAR, Curve.InterpolationEntity.LOG_OF_VALUE);


		double[][] forwardCurveData = new double[rowsForward][columns];
		for(int i=0; i< rowsForward;i++){
			for(int j=0;j<columns; j++){
				forwardCurveData[i][j] = Double.parseDouble(volatilities.getValueAt(j+1, i+27).toString());

			}
		}

		double[] timesForward = new double[forwardCurveData.length];
		double[] givenForwardDiscountFactors = new double[forwardCurveData.length];
		for(int i=0; i<timesForward.length; i++) {
			timesForward[i] = discountCurveData[i][3]*360/365;
			givenForwardDiscountFactors[i] = forwardCurveData[i][2];
		}
		double paymentOffset = 0.5;
		forwardCurve = ForwardCurve.createForwardCurveFromDiscountFactors("forwardCurve", timesForward, givenForwardDiscountFactors, paymentOffset);
	}
	

	public double[][] getVolatilities() throws IOException {
		int rows = 10;
		int columns = 14;
		double[][] data = new double[rows][columns];
		for(int i=0; i< rows;i++){
			for(int j=0;j<columns; j++){
				data[i][j] = Double.parseDouble(volatilities.getValueAt(j+1, i+1).toString());
			}
		}
		
		return data;
	}

	
	public double getSwapAnnuity(double optionMaturity, double swapTenor) throws IOException {

		return SwapAnnuity.getSwapAnnuity(new TimeDiscretization(optionMaturity, optionMaturity+swapTenor, periodLength, ShortPeriodLocation.SHORT_PERIOD_AT_START), discountCurve);
		
	}
	
	public double getForwardSwapRate(double optionMaturity, double swapTenor) throws IOException {
		
		TimeDiscretizationInterface fixTenor = new TimeDiscretization(optionMaturity, optionMaturity+swapTenor, 1.0, ShortPeriodLocation.SHORT_PERIOD_AT_START);
		TimeDiscretizationInterface floatTenor = new TimeDiscretization(optionMaturity, optionMaturity+swapTenor, periodLength, ShortPeriodLocation.SHORT_PERIOD_AT_START);

		return Swap.getForwardSwapRate(fixTenor, floatTenor, forwardCurve, discountCurve);
		
	}

	
	public double[][] getSmileData(int tenor) throws IOException {
		int rows = 6;
		int columns = 9;
		double[][] data = new double[rows][columns];
		if(tenor == 2){
			for(int i=0; i< rows;i++){
				for(int j=0;j<columns; j++){
					data[i][j] = Double.parseDouble(smile.getValueAt(j+1, i+2).toString());
				}
			}
		}
		if(tenor == 5){
			for(int i=0; i< rows;i++){
				for(int j=0;j<columns; j++){
					data[i][j] = Double.parseDouble(smile.getValueAt(j+1, i+11).toString());
				}
			}
		}
		if(tenor == 10){
			for(int i=0; i< rows;i++){
				for(int j=0;j<columns; j++){
					data[i][j] = Double.parseDouble(smile.getValueAt(j+1, i+20).toString());
				}
			}
		}
		return data;
	}

}
