package vrp.Draw;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

public class ExcelUtils {
    public static void fillAvgToExcel(int max_iteration, double avg[], int rowIndex) throws IOException {
        String excelFilePath = "C:\\Users\\HOANG\\Desktop\\anh Nam\\vrp\\result_02122022_chay1000lan.xlsx";
        File file = new File(excelFilePath);
        FileInputStream inputStream = new FileInputStream(file);
        Workbook workbook = new XSSFWorkbook(inputStream);
        Sheet sheet = workbook.getSheet("avg_30customer_chiadeuxe");

        Row row = sheet.getRow(rowIndex);

        for (int i=0; i<max_iteration; i++){
            Cell avg_cell = row.getCell(1+i);
            avg_cell.setCellValue(avg[i]);
        }

        inputStream.close();
        FileOutputStream out = new FileOutputStream(file);
        workbook.write(out);
        out.close();
        System.out.println("Write to excel done!");
    }

    public static void fillForDrawFunctionToExcel(double [] Fmin, int startRowFmin, int maxIter) throws IOException {
        String excelFilePath = "C:\\Users\\HOANG\\Desktop\\anh Nam\\vrp\\result_26112022.xlsx";
        File file = new File(excelFilePath);
        FileInputStream inputStream = new FileInputStream(file);
        Workbook workbook = new XSSFWorkbook(inputStream);
        Sheet sheet = workbook.getSheet("avg_25customer");

        //F min
        int rowIndex = startRowFmin;
        Row row = sheet.getRow(rowIndex);
        for (int i=0; i<maxIter; i++){
            Cell cell = row.getCell(1+i);
            cell.setCellValue(Fmin[i]);
        }

        inputStream.close();
        FileOutputStream out = new FileOutputStream(file);
        workbook.write(out);
        out.close();
        System.out.println("Write to excel done!");
    }

    public static void fillX1X2ToExcelForDraw(double [] X_1, double [] X_2, int maxIter, int startRow) throws IOException {
        String excelFilePath = "C:\\Users\\HOANG\\Desktop\\anh Nam\\vrp\\x1x2.xlsx";
        File file = new File(excelFilePath);
        FileInputStream inputStream = new FileInputStream(file);
        Workbook workbook = new XSSFWorkbook(inputStream);
        Sheet sheet = workbook.getSheet("Sheet1");

        //X1
        int rowIndex = startRow;
        Row row = sheet.getRow(rowIndex);
        for (int i=0; i<maxIter; i++){
            Cell cell = row.getCell(1+i);
            cell.setCellValue(X_1[i]);
        }

        //X2
        rowIndex = startRow+1;
        row = sheet.getRow(rowIndex);
        for (int i=0; i<maxIter; i++){
            Cell cell = row.getCell(1+i);
            cell.setCellValue(X_2[i]);
        }

        inputStream.close();
        FileOutputStream out = new FileOutputStream(file);
        workbook.write(out);
        out.close();
        System.out.println("Write to excel done!");
    }
}
