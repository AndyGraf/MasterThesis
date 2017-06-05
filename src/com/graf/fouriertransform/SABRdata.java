package com.graf.fouriertransform;

import java.io.File;
import java.io.IOException;

import org.jopendocument.dom.spreadsheet.Sheet;
import org.jopendocument.dom.spreadsheet.SpreadSheet;

public class SABRdata {
	
	private final Sheet volatilitiesAndCurves;
	private final Sheet smile;
	
	public SABRdata() throws IOException{
		
		File file = new File("src/Data.ods");
		volatilitiesAndCurves = SpreadSheet.createFromFile(file).getSheet(0);
		smile = SpreadSheet.createFromFile(file).getSheet(1);
		
	}
	

	public double[][] getVolatilities() throws IOException {
		int rows = 10;
		int columns = 14;
		double[][] data = new double[rows][columns];
		for(int i=0; i< rows;i++){
			for(int j=0;j<columns; j++){
				data[i][j] = Double.parseDouble(volatilitiesAndCurves.getValueAt(j+1, i+1).toString());
			}
		}
		
		return data;
	}

}
