

public class Gauss {

	/**
	 * Diese Methode soll die Loesung x des LGS R*x=b durch
	 * Rueckwaertssubstitution ermitteln.
	 * PARAMETER: 
	 * R: Eine obere Dreiecksmatrix der Groesse n x n 
	 * b: Ein Vektor der Laenge n
	 */
	public static double[] backSubst(double[][] R, double[] b) {
		/*formula (basically): Xi = [Bi - Sum(j=1, to n, of Aij*Xj)] / Aii*/
		int n = b.length;
		double[] x = new double[n];
		for(int i = n-1; i > -1; i--) {
			x[i] = b[i];
			System.out.print("do X"+i+" = ( B"+i+" ");
			for(int j = i+1; j < n; j++) {
				x[i] -= R[i][j] * x[j];
				System.out.print("- A"+i+j+"*X"+j+" ");
			}
			x[i] /= R[i][i];
			System.out.print(")/ A"+i+i+"\n");
		}
		return x;
	}

	/**
	 * Diese Methode soll die Loesung x des LGS A*x=b durch Gauss-Elimination mit
	 * Spaltenpivotisierung ermitteln. A und b sollen dabei nicht veraendert werden. 
	 * PARAMETER: A:
	 * Eine regulaere Matrix der Groesse n x n 
	 * b: Ein Vektor der Laenge n
	 */
	public static double[] solve(double[][] A, double[] b) {
		//declare vector/matrix length and matrix needed for solution, etc
		int n = b.length;
		int j = 0;
		double alpha = 0.0;
		double[][] L = new double[n][n];
		
		for(int k = 0; k < n-1; k++) {
			alpha = Math.abs(A[k][k]); j = k;
			for(int s = k+1; s < n; s++) {
				if(Math.abs(A[s][k]) > alpha) {
					alpha = Math.abs(A[s][k]); j=s;
				}
			}
			
			//# Pivotelement ist a(j,k) und Pivotzeile ist j
			for(int i = k; i < n; i++) {
				alpha = A[k][i]; A[k][i] = A[j][i]; A[j][i] = alpha;
			}
			alpha = b[j]; b[j] = b[k]; b[k] = alpha;
			
			//# Eliminationsschritt
			for(int s = k+1; s < n; s++) {
				L[s][k] = A[s][k] / A[k][k];
				b[s] = b[s] - L[s][k]*b[k];
				for(int i = k; i < n; i++) {
					A[s][i] = A[s][i] - L[s][k]*A[k][i];
				}
			}
		}
		
		return backSubst(A, b);
	}

	/**
	 * Diese Methode soll eine Loesung p!=0 des LGS A*p=0 ermitteln. A ist dabei
	 * eine nicht invertierbare Matrix. A soll dabei nicht veraendert werden.
	 * 
	 * Gehen Sie dazu folgendermassen vor (vgl.Aufgabenblatt): 
	 * -Fuehren Sie zunaechst den Gauss-Algorithmus mit Spaltenpivotisierung 
	 *  solange durch, bis in einem Schritt alle moeglichen Pivotelemente
	 *  numerisch gleich 0 sind (d.h. <1E-10) 
	 * -Betrachten Sie die bis jetzt entstandene obere Dreiecksmatrix T und
	 *  loesen Sie Tx = -v durch Rueckwaertssubstitution 
	 * -Geben Sie den Vektor (x,1,0,...,0) zurueck
	 * 
	 * Sollte A doch invertierbar sein, kann immer ein Pivot-Element gefunden werden(>=1E-10).
	 * In diesem Fall soll der 0-Vektor zurueckgegeben werden. 
	 * PARAMETER: 
	 * A: Eine singulaere Matrix der Groesse n x n 
	 */
	public static double[] solveSing(double[][] A) {
		//declare matrix/vector length, vectors/matrix needed for solving, etc
		int n = A.length;
		int j = 0;
		double alpha = 0.0;
		double[] p = new double[n];
		double[][] L = new double[n][n];
		double pseudoZero = 1E-10;
		int Tn = 0;
		
		/*Modified Gauss Elimination*/
		System.out.print("Mod. Gauss Elim. ...\n");
		for(int k = 0; k < n-1; k++) {
			
			System.out.print("\nGauss Matrix at step " + k + "\n");
			for(int x = 0; x < n; x++) {
				System.out.print(x + "{ ");
				for(int y = 0; y < n; y++) {
					System.out.print(A[x][y] + ", ");
				}
				System.out.print("}\n");
			}
			
			alpha = Math.abs(A[k][k]); j = k;
			for(int s = k+1; s < n; s++) {
				if(Math.abs(A[s][k]) > alpha) {
					alpha = Math.abs(A[s][k]); j=s;
				}
			}
			
			//check for whether the only "found" pivot element is effectively 0, if so, exit the loop
			if(alpha < pseudoZero) { System.out.print("breaking gauss elim\n"); pseudoZero = 1.0; break; }
			
			//# Pivotelement ist a(j,k) und Pivotzeile ist j
			for(int i = k; i < n; i++) {
				alpha = A[k][i]; A[k][i] = A[j][i]; A[j][i] = alpha;
			}

			//# Eliminationsschritt
			for(int s = k+1; s < n; s++) {
				L[s][k] = A[s][k] / A[k][k];
				for(int i = k; i < n; i++) {
					A[s][i] = A[s][i] - L[s][k]*A[k][i];
				}
			}
		}

		System.out.print("\nGauss Matrix after elim\n");
		for(int x = 0; x < n; x++) {
			System.out.print(x + "{ ");
			for(int y = 0; y < n; y++) {
				System.out.print(A[x][y] + ", ");
			}
			System.out.print("}\n");
		}
		
		//Matrix is invertible, i.e. no pivot element <1E-10 was found, 
		//so we never set pseudoZero to 1.0, so we return a vector with all zeroes
		if(pseudoZero < 1.0) { return new double[n]; }
		//otherwise we set pseudoZero correctly again and move on
		else { pseudoZero = 1E-10; }
		
		System.out.print("Determining T matrix size; n = " + n);
		//find out where T ends
		for(int i = 0; i < n; i++) {
			//since the matrix should be perfectly diagonal by now, 
			//our T matrix should end at the first element Aii that is effectively zero
			if(A[i][i] > pseudoZero) { 
				Tn = i; 
			} else {
				break;
			}
		}
		System.out.print("; Tn = " + Tn + "\n");
		
		//declare the T and v matrix/vector for doing Tx = -v
		double[][] T = new double[Tn][Tn];
		double[] v = new double[Tn];
		
		System.out.print("\nCreating T matrix\n");
		//copy values from A into T and v
		for(int x = 0; x < Tn; x++) {
			System.out.print(x + "{ ");
			for(int y = 0; y < Tn; y++) {
				T[x][y] = A[x][y];
				System.out.print(T[x][y] + ", ");
			}
			v[x] = -1 * A[x][Tn];
			System.out.print(" } v" +x+ "{" + v[x] + "}\n");
		}
		
		System.out.print("\nDo Tx = -v:\n");
		//calculate solution x for Tx = -v (write it straight into v again since we wont need that anymore afterwards
		v = backSubst(T, v);
		//copy values of said x (vector v) into p, then put 1 after it, as described in task description
		for(int i = 0; i < Tn; i++) {
			p[i] = v[i];
		}
		p[Tn] = 1;
		
		//just visual output of p solution vector
		System.out.print("p:{ ");
		for(int i = 0; i < n; i++) {
			System.out.print(p[i] + ", ");
		} System.out.print("}\n");
		
		return p;
	}

	/**
	 * Diese Methode berechnet das Matrix-Vektor-Produkt A*x mit A einer nxm
	 * Matrix und x einem Vektor der Laenge m. Sie eignet sich zum Testen der
	 * Gauss-Loesung
	 */
	public static double[] matrixVectorMult(double[][] A, double[] x) {
		int n = A.length;
		int m = x.length;

		double[] y = new double[n];

		for (int i = 0; i < n; i++) {
			y[i] = 0;
			for (int j = 0; j < m; j++) {
				y[i] += A[i][j] * x[j];
			}
		}

		return y;
	}
}
