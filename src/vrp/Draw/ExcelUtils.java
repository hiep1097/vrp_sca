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
        String excelFilePath = "C:\\Users\\HOANG\\Desktop\\anh Nam\\vrp\\result.xlsx";
        File file = new File(excelFilePath);
        FileInputStream inputStream = new FileInputStream(file);
        Workbook workbook = new XSSFWorkbook(inputStream);
        Sheet sheet = workbook.getSheet("Sheet2");

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
        String excelFilePath = "C:\\Users\\HOANG\\Desktop\\anh Nam\\vrp\\result.xlsx";
        File file = new File(excelFilePath);
        FileInputStream inputStream = new FileInputStream(file);
        Workbook workbook = new XSSFWorkbook(inputStream);
        Sheet sheet = workbook.getSheet("Sheet1");

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
}
