package com.graf.fouriertransform;

import java.io.File;
import java.io.IOException;

import org.jopendocument.dom.spreadsheet.Sheet;
import org.jopendocument.dom.spreadsheet.SpreadSheet;

public class SABRdata {
	
	private final Sheet volatilities;
	private final Sheet smile;
	private final Sheet curves;
	
	public SABRdata() throws IOException{
		
		File file = new File("src/Data.ods");
		volatilities = SpreadSheet.createFromFile(file).getSheet(0);
		smile = SpreadSheet.createFromFile(file).getSheet(1);
		curves = SpreadSheet.createFromFile(file).getSheet(2);
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
	
	public double[][] getForwardCurveData() throws IOException {
		int rows = 51;
		int columns = 5;
		double[][] data = new double[rows][columns];
		for(int i=0; i< rows;i++){
			for(int j=0;j<columns; j++){
				data[i][j] = Double.parseDouble(curves.getValueAt(j+1, i+3).toString());
			}
		}
		
		return data;
	}
	
	public double[][] getDiscountCurveData() throws IOException {
		int rows = 51;
		int columns = 5;
		double[][] data = new double[rows][columns];
		for(int i=0; i< rows;i++){
			for(int j=0;j<columns; j++){
				data[i][j] = Double.parseDouble(curves.getValueAt(j+8, i+3).toString());
			}
		}
		
		return data;
	}
	
	public double getSwapAnnuity(double optionMaturity, double swapTenor) throws IOException {
		int rows = 107;
		int columns = 5;
		double[][] discountCurveData = new double[rows][columns];
		for(int i=0; i< rows;i++){
			for(int j=0;j<columns; j++){
				discountCurveData[i][j] = Double.parseDouble(volatilities.getValueAt(j+8, i+27).toString());
			}
		}
		
		int startingIndex = 0;
		while( discountCurveData[startingIndex][3]< optionMaturity){startingIndex++;}
		
		double swapAnnuity=0;
		int index=0;
		while( discountCurveData[startingIndex+index][3]< optionMaturity + swapTenor){
			swapAnnuity+=discountCurveData[startingIndex+index][2]*(discountCurveData[startingIndex+index+1][3]-discountCurveData[startingIndex+index][3]);
			index++;
		}
		
		
		return swapAnnuity/(discountCurveData[startingIndex+index][3]-discountCurveData[startingIndex][3]);
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
