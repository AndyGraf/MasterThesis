package com.graf.fouriertransform;

import java.io.File;
import java.io.IOException;

import org.jopendocument.dom.spreadsheet.Sheet;
import org.jopendocument.dom.spreadsheet.SpreadSheet;

public class dataTest {

	public static void main(String[] args) throws IOException {
		SABRdata sheetdata = new SABRdata();
		double[][] data = new double[10][14];
		data = sheetdata.getVolatilities();
		
		int rows = 10;
		int columns = 14;
		for(int i=0; i< rows;i++){
			for(int j=0;j<columns; j++){
				System.out.print(data[i][j] + "\t");
			}
			System.out.println();
		}
		
		
		
	}

}
