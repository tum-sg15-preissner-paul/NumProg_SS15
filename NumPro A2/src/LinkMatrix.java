import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;


public class LinkMatrix {
	public int[][] L;
	public String[] urls;

	// private static final String lineDelim = ";";
	private static final String delim = " ";

	public void read(String filename) throws FileNotFoundException, IOException {
		String curRow[];

		BufferedReader file = new BufferedReader(new FileReader(filename));

		try {
			int n = Integer.valueOf(file.readLine());

			L = new int[n][n];
			for (int i = 0; i < n; i++) {
				curRow = file.readLine().split(delim);
				for (int j = 0; j < n; j++)
					L[i][j] = Integer.valueOf(curRow[j]);
			}

			urls = new String[n];
			for (int i = 0; i < n; i++)
				urls[i] = file.readLine();
		} catch (IOException e) {
			throw e;
		} finally {
			file.close();
		}
	}

	public void write(String filename) throws FileNotFoundException,
			IOException {
		String temp;

		BufferedWriter file = new BufferedWriter(new FileWriter(filename));

		try {
			file.write(String.valueOf(L.length) + "\r\n");

			for (int i = 0; i < L.length; i++) {
				temp = "";
				for (int j = 0; j < L.length; j++)
					temp += String.valueOf(L[i][j]) + delim;
				file.write(temp + "\r\n");
			}

			for (int i = 0; i < urls.length; i++)
				file.write(urls[i] + "\r\n");
		} catch (FileNotFoundException e) {
			throw e;
		} catch (IOException e) {
			throw e;
		} finally {
			file.close();
		}
	}
}