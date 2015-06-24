
public class MatrixPlotter {

	
	public static  void plot(double[][] a, int size)
	{
		System.out.println("Plotting Matrix: ");
		for(int i =0; i<size; i++)
		{
			System.out.print("{");
			for(int j =0; j<size; j++)
			{
				System.out.print(""+a[i][j]+", ");
			}
			System.out.print("}");
			System.out.println("");
		}
	}
	
}
